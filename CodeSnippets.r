# This is a code snippet showing the code that runs the stage 1 model as described in the methods section of the 
# paper titled "Community Factors and excess mortality in the first wave of the COVID-19 pandemic?"

# The model uses the INLA library

# Data supplied to stage 1
# M & F analysed separately
# 135820 rows (6791 MSOAs * 5 years (comparion period 2015 - 2019) * 4 age groups)
#Data$count -  count of deaths for this MSOA,year,Age cat
#Data$MSOAOrder - index of MSOA
#Data$YEAR
#Data$AGECAT
#Data$value - population for this MSOA,year,Age cat
# IG is the inla.graph object containing the list of adjacencies of the MSOAs for England (islands are deemed to connect to the nearest MSOA on the mainland)


# INLA code
library(INLA)
#Set up the formula to use
formula = count~ 1 +
  f(MSOAOrder, model = 'bym', graph = IG, scale.model = TRUE,
    hyper=list(prec.unstruct=list(param=c(1,0.1)),
               prec.spatial=list(param=c(1,0.1)))) +
  f(YEAR, model = "iid") +
  relevel(as.factor(AGECAT), ref = 4)  #use the highest age cat (80+) as the reference (most deaths)

# run the stage 1 model
result=inla(formula, family = 'poisson', E=value, data = Data,  
             control.compute=list(dic=TRUE, config=TRUE), quantiles=c(0.025,0.5,0.975),verbose = TRUE)

# Save the model output
save(result, file = "results/AllCauseStage1M_INLA.Rdata")


#Extract the 100 stage 1 samples
draws <- inla.posterior.sample(50, result)
rate.draws <- lapply(draws, function(i)
  exp(i$latent[grep("Predictor",rownames(i$latent))]))
rate.drawsMed<-array(unlist(rate.draws), dim=c(dim(DataM)[1], 50)); dim(rate.drawsMed)
dM = as.data.frame(rate.drawsMed)
# Add to the data and save
Data= cbind(Data,dM)
rm(dM)
write.csv(DataM, file="results/AllCauseStage1M_data.csv")


# calculate the lambda from the results of stage1
DataGp = Data %>% group_by(AGECAT, MSOA11)  
DataGpMedian = DataGp %>% summarise_at(vars(ObjectID, V1:V100), mean)  

#merge with the 2020 population and mortlaity data
Data2020 = merge(x=Pop2020, y=DataGpMedian, by.x = c("MSOA11", "AGECAT"), by.y = c("MSOA11", "AGECAT"))
write.csv(Data2020, file="results/AllCauseData202M.csv",row.names = FALSE)


#At this stage we add in all the covariate data (quintiles) that we use. defined at MSOA level (see details in the paper)

#Now create 100 datafiles (1 for each of the samples we extracted) for easier loading in stage 2
for (i in 1:50){
  sampleM = Data2020[,c(1:14,(i+14))]
  write.csv(sampleM, file=paste("results/samples/AllCauseMsample",i,".csv",sep = ""),row.names = FALSE)
}


#The stage 2 models are run 100 times (50 males, 50 females) using the samples extracted in stage 1.
# Used a command line option with .bat files to run them

# extract from a stage 2 .r file:
#read in the sample
data = read.csv(paste(dataFile,".csv",sep = ""), header = TRUE, stringsAsFactors = FALSE)

#Set the expected column using the index we want to run
data$Exp = data$value * data[, 15] # the fifteenth column the lambda when there are 7 covariates present

library(INLA)

#Set up the formula, this one uses all 7 covariates
formula = count~ f(MSOAOrder, model = 'bym', graph = IG, scale.model = TRUE,
                   hyper=list(prec.unstruct=list(param=c(1,0.1)),
                              prec.spatial=list(param=c(1,0.1))))  +
  relevel(as.factor(AGECAT), ref = 4) +
  as.factor(OCR_Q) +
  as.factor(NWHITE_Q) +
  as.factor(INCOME_Q) +
  as.factor(NO2_Q) +
  as.factor(PD_Q) +
  as.factor(PM25_Q) +
  as.factor(CH_pop_Q)
  
cat("Starting INLA \n")
result2020=inla(formula, family = 'poisson', E=Exp, data = data, 
                control.compute=list(dic=TRUE, config=TRUE), quantiles=c(0.025,0.5,0.975),verbose = FALSE)
cat("INLA complete, extracting  draws \n")

draws<-inla.posterior.sample(100, result2020)
rate.draws <- lapply(draws, function(i)
  exp(i$latent[grep("Predictor",rownames(i$latent))]))
rate.drawsMed<-array(unlist(rate.draws), dim=c(dim(data)[1], 100)); dim(rate.drawsMed)
dF = as.data.frame(rate.drawsMed)
# Add to the data and save
data= cbind(data,dF)
rm(dF)
#Write the output
cat("Saving results \n")
write.csv(data, file=paste(dataFile,"res.csv",sep = ""),row.names = FALSE)
cat("Finish, time: ", (Sys.time() - start_time), "\n")

#Once the batch is complete there will be 5000 samples (50 from stage 1 * 100 from stage 2) which are used to extract the various measure of excess mortality 
# used in the paper. Plus the rate ratios for the covariates.



