##############################################################################################################
# Master script for "Wearable sensors enable personalized predictions of clinical laboratory measurements"   #
# Jessilyn Dunn & Lukasz Kidzinski                                                                           #
# initiated: 11/01/2017                                                                                      #
# last modified: 04/03/2021                                                                                  #
# Output: Figures for paper                                                                                  #
##############################################################################################################

#### LIBRARY DEPENDENCIES:
library(ggplot2)
library(data.table)
library(psych)
library(Rmisc)
library(gridExtra)
library(dplyr) #plyr
library(MASS)
library(ggthemes)
library(reshape2)
library(randomForest)
library(lme4)
library(fasttime) #fastPOSIXct
library(lubridate)
library(tidyr)
library(zoo) #rollapply
library(grid) 
library(nlme) 
library(MuMIn) 
library(glmnet)

# A theme for plots
weartals_theme = theme_classic() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

if(!dir.exists("plots")) dir.create("plots")

# FUNCTIONS
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

findOutlier <- function(data, cutoff = 2) {
  ## Calculate the sd
  sds <- apply(data, 2, sd, na.rm = TRUE)
  ## Identify the cells with value greater than cutoff * sd (column wise)
  result <- mapply(function(d, s) {
    which(d > cutoff * s)
  }, data, sds)
  result
}

fastDate <- function(x, tz=NULL)
  as.Date(x, tz=tz)

# Path to the directory with data
dir = "../SECURE_data/"

###################
#### READ DATA ####
###################

#iPOP wearables/clinical combined data
# wear <- read.csv(paste0(dir, "Basis2016_Cleaned_NotNorm0824_WeekPrior.csv"),
#                   header=TRUE,sep=',',stringsAsFactors=FALSE) 

timespans <-c("AllData",
              "MonthPrior",
              "2WeekPrior",
              "WeekPrior",
              "5DayPrior",
              "3DayPrior",
              "DayPrior" )

wear <- read.csv(paste0(dir, 
                "Basis2016_Clean_Norm_", timespans[7], "_20180504.csv"),
                 header=TRUE,sep=',',stringsAsFactors=FALSE)
wearables.people <- unique(wear$iPOP_ID)

# iPOP vitals 
iPOPvitals <- read.csv(paste0(dir, "vitals.csv"),
  header=TRUE,sep=',',stringsAsFactors=FALSE)
# restrict to wearables only people
iPOPvitals <- iPOPvitals[iPOPvitals$HIMCID %in% wearables.people,]

#iPOP Labs
iPOPlabs <- read.csv(
  paste0(dir, "lab_results_20170717.csv"),
  header=TRUE,sep=',',stringsAsFactors=FALSE)
# restrict to wearables only people
iPOPlabs <- iPOPlabs[iPOPlabs$HIMC_ID %in% wearables.people,]

# SEHR (aka 30k) vitals
vitals <- fread(paste0(dir, "all_vitals.csv"),
                header=TRUE,sep=',',stringsAsFactors=FALSE)
vitals <- data.frame(vitals)
length(na.omit(vitals$PULSE)) # 552,145
length(na.omit(vitals$TEMPERATURE)) # 333,821

# SEHR labs
labs <- fread(paste0(dir, "all_labs.csv"),
              header=TRUE,sep=',',stringsAsFactors=FALSE)
# SEHR labs/vitals combined file
corDf <- read.csv(paste0(dir, "20180412_Cleaned_joined_30k_labs_ALLvitals.csv"),
                   header=TRUE,sep=',',stringsAsFactors=FALSE)
# SEHR Codes
icd <- fread(paste0(dir, "all_icds.csv"),
             header=TRUE,sep=',',stringsAsFactors=FALSE)
cpt <- fread(paste0(dir, "all_cpts.csv"),
             header=TRUE,sep=',',stringsAsFactors=FALSE)

# iPOP demographics
iPOPdemographics <- read.csv(paste0(dir, "SECURE_ClinWearDemo_SamplePop.csv"),
                  header=TRUE,sep=',',stringsAsFactors=FALSE)
iPOPicd <-fread(paste0(dir, "diagnoses_110718.csv"),
                header=TRUE,sep=',',stringsAsFactors=FALSE)
iPOPcpt <-fread(paste0(dir, "procedures_110718.csv"),
                header=TRUE,sep=',',stringsAsFactors=FALSE)

# SEHR demographics
thirtyKdemog <- read.csv(paste0(dir, "SECURE_20180412_thirtyKDemog.csv"),
         header=TRUE,sep=',',stringsAsFactors=FALSE)

# names of top clinical tests used for the downstream analyses
top.names<-c("LYM", "NEUT", "LYMAB", "NEUTAB", "IGM", "HSCRP", "ALKP", "ALT", "HDL", "MCV", "TBIL", "CHOLHDL", "GLOB", "AG", "CO2", "CA", "LDLHDL", "BUN", "NHDL", "NA.", "UALB", "MONOAB", "CHOL", "MONO", "RDW", "HCT", "TP", "TGL", "EOS", "LDL", "GLU", "AST", "PLT", "K", "EOSAB", "BASOAB", "MCH", "ALB", "HGB", "A1C", "CL", "RBC", "BASO", "MCHC") 



###################
### CLEAN DATA ####
###################

### clean iPOP Vitals ###
names(iPOPvitals)[which(names(iPOPvitals)=="HIMCID")] <- "iPOP_ID"
names(iPOPvitals)[which(names(iPOPvitals)=="RECORDED_TIME")] <- "Clin_Result_Date"
names(iPOPvitals)[which(names(iPOPvitals)=="X.Pulse.")] <- "Pulse"
names(iPOPvitals)[which(names(iPOPvitals)=="X.Temp.")] <- "Temp"
names(iPOPvitals)[which(names(iPOPvitals)=="X.BP.")] <- "BP"
names(iPOPvitals)[which(names(iPOPvitals)=="X.Bmi.")] <- "BMI"

for (i in 1:length(iPOPvitals$BP)){
  iPOPvitals$systolic[i] <- strsplit(as.character(iPOPvitals$BP),'/')[[i]][1]
  iPOPvitals$diastolic[i] <- strsplit(as.character(iPOPvitals$BP),'/')[[i]][2]
}

#Reformat dates
iPOPvitals$Clin_Result_Date <- format(
  as.Date(iPOPvitals$Clin_Result_Date, "%d-%b-%Y"), "%Y-%m-%d")

#Make correlation variables numeric
iPOPvitals[,c("Pulse","Temp","BP","BMI","systolic", "diastolic")] <- apply(
  iPOPvitals[,c("Pulse","Temp","BP","BMI","systolic", "diastolic")], 2,
  function(x) as.numeric(as.character(x)))

#### CLEAN iPOP LABS DATA ####

#Rename columns
names(iPOPlabs)[which(names(iPOPlabs)=="HIMC_ID")] <- "iPOP_ID"
names(iPOPlabs)[which(names(iPOPlabs)=="RESULT_TIME")] <- "Clin_Result_Date"

#Reformat dates
iPOPlabs$Clin_Result_Date <- format(
  as.Date(iPOPlabs$Clin_Result_Date, "%d-%b-%Y"), "%Y-%m-%d")

allClin <- c("ALKP", "LYM", "HSCRP", "ALT", "NEUT", "TBIL", "IGM", "TGL", "MCV", "MCH", "CO2", "LYMAB", "NEUTAB", "UALB", "CHOL", "MONOAB", "ALB", "NA.", "HDL", "PLT", "AG", "HGB", "EOS", "CL", "BUN", "GLOB", "CA", "CHOLHDL", "HCT", "BASOAB", "A1C", "GLU", "LDLHDL", "TP", "EOSAB", "K", "NHDL", "RBC", "MONO", "AST", "MCHC", "RDW", "BASO", "LDL")

for(i in 1:length(allClin)){ #this removes non-numeric characters
  cache <- iPOPlabs[,c(allClin[i])]
  cache <- gsub("[^0-9.]","",cache) #this keeps decimals
  iPOPlabs[,c(allClin[i])] <- cache
}

#Make correlation variables numeric
iPOPlabs[,c(allClin)] <- apply(
  iPOPlabs[,c(allClin)], 2,
  function(x) as.numeric(as.character(x)))

#subset by allClin
iPOPlabs <- iPOPlabs[names(iPOPlabs) %in% c("iPOP_ID","Clin_Result_Date",allClin)]

#Merge data
iPOPcorDf <- merge(iPOPlabs,
                   iPOPvitals[,c("iPOP_ID","Clin_Result_Date",
                                 "Pulse","Temp","BMI","systolic","diastolic")],
                   by=c("iPOP_ID","Clin_Result_Date"))

iPOPcorDf2 <- iPOPcorDf
### clean iPOPcorDf ###
iPOPcorDf[, -c(1,2)] <- apply(iPOPcorDf[, -c(1,2)], 2, remove_outliers)
outliers <- findOutlier(iPOPcorDf)
library(compare)
comparison <- compare(iPOPcorDf,iPOPcorDf2,allowAll=TRUE)
difference <-
  data.frame(lapply(1:ncol(iPOPcorDf),function(i)setdiff(iPOPcorDf[,i],comparison$tM[,i])))
colnames(difference) <- colnames(iPOPcorDf)
#difference #no data points are removed.

### clean corDf ### 
summary(corDf)
corDf2 <- corDf
corDf[, -c(1,2)] <- apply(corDf[, -c(1,2)], 2, remove_outliers) 

comparison <- compare(corDf,corDf2,allowAll=TRUE)
difference <-
  data.frame(lapply(1:ncol(corDf),function(i)setdiff(corDf[,i],comparison$tM[,i])))
colnames(difference) <- colnames(corDf)
difference
outliers <- findOutlier(corDf) # http://rpubs.com/lcollado/7904

## clean SEHR icd and cpt dates
icd$ICD_DATE <- format(
  as.Date(icd$ICD_DATE, "%d-%b-%Y"), "%Y-%m-%d")
icd$ICD_DATE <- as.Date(icd$ICD_DATE, "%Y-%m-%d")
cpt$CPT_DATE <- format(
  as.Date(cpt$CPT_DATE, "%d-%b-%Y"), "%Y-%m-%d")
cpt$CPT_DATE <- as.Date(cpt$CPT_DATE, "%Y-%m-%d")

### clean wear ### for the CCA in Fig 2E
# wear2<-wear
# wear[,-c(1:6)] <- apply(wear[,-c(1:6)], 2, function(x) as.numeric(as.character(x)))
# wear[,-c(1:6)] <- apply(wear[,-c(1:6)], 2, remove_outliers) 

## merge iPOP and demographics
iPOPcorDf.demo <- merge(iPOPcorDf, iPOPdemographics[1:4], by="iPOP_ID")

###########################
#     Main Text Numbers   #
###########################

iPOPdaysMonitored <- read.csv(paste0(dir,"slide2_C_participant_data_summary.csv"),
                  header=TRUE,sep=',',stringsAsFactors=FALSE)
# restrict to wearables only people
iPOPdaysMonitored <- iPOPdaysMonitored[iPOPdaysMonitored$iPOP_ID %in% wearables.people,]
mean(iPOPdaysMonitored$Days_monitored_by_clinic)/365 # 3.3 years of clinic monitoring
mean(iPOPdaysMonitored$Total_NumOfClinMeasures) # avg of 42 clinic visits
mean(iPOPdaysMonitored$Days_monitored_by_basis) # avg of 343 clinic visits
mean(iPOPdaysMonitored$Days_overlapping_between_basis_and_clinic_monitoring) # avg of 313 clinic visits

#number of iPOPpers with > 50 clinic visits
sum(table(iPOPcorDf.demo$iPOP_ID) > 50)

iPOP.table <- fread(paste0(dir, "ClinWearDemo_SamplePop.csv"),
header=TRUE,sep=",",stringsAsFactors = FALSE)
# restrict to wearables only people
iPOP.table <- iPOP.table[iPOP.table$iPOP_ID %in% wearables.people,]
iPOP.table$iPOP_ID <- as.factor(iPOP.table$iPOP_ID)
iPOP.table$Gender <- as.factor(iPOP.table$Gender)
iPOP.table$Ethn <- as.factor(iPOP.table$Ethn)
iPOP.table$AgeIn2016 <- as.numeric(iPOP.table$AgeIn2016)
table(iPOP.table[1:54,]$Ethn)
female <- iPOP.table[iPOP.table$Gender=="F"]
max(female$AgeIn2016);  min(female$AgeIn2016); mean(female$AgeIn2016) # max age of female particpants

male <- iPOP.table[iPOP.table$Gender=="M"]
male <- male[-max(dim(male)),]
max(male$AgeIn2016);  min(male$AgeIn2016); mean(male$AgeIn2016) # max age of male particpants


df <- fread(paste0(dir, "/BasisData_20161111_PostSummerAddOns_Cleaned_NotNormalized_20180427.csv"),
                   header=TRUE,sep=",",stringsAsFactors = FALSE,
                   select=c("Timestamp_Local","Heart_Rate","Skin_Temperature_F",
                            "Steps","iPOP_ID")) # 38009228 observations
df<-df[!is.na(df$iPOP_ID),] # only keep wearables data that has an associated iPOP ID; 
df[,"Heart_Rate"] <- apply(df[,"Heart_Rate"], 2, remove_outliers) # clean data based on HR
df <- df[which(!is.na(df[,"Heart_Rate"])),] # remove observations where HR = NA
# restrict to daytime values only (between 6am-10pm)
df$Timestamp_Local<-fastPOSIXct(df$Timestamp_Local) # takes a very long time

# nighttime is currently represented below
daytime.df <- with( df, df[ hour( Timestamp_Local ) >= 19 & hour( Timestamp_Local ) < 21 , ] ) # pull data only from specific time window; store hourly resting data for boxplots

vitals <- iPOPvitals
vitals$Clin_Result_Date<-as.Date(vitals$Clin_Result_Date, "%Y-%m-%d")
colnames(vitals)[1:5] <- c("iPOP_ID", "Date", "BP", "Pulse", "Temp")

##################
# Fig 1B - Danny #
##################
#Set window and step maximums:
windows=c(10, 60) # define time windows with no steps for resting threshold
steps = c(1, 50) #define number max steps

#To check day-prior, supply TRUE or FALSE here:
#If TRUE, merges clinical based on day-prior, and also adds "dayPrior" to output file names.
dayPrior = FALSE

#Iterate through array of window sizes
for (window in windows){
  #Iterate through array of maximum steps
  for (maxsteps in steps){
    restingDf <- c() 
    restingDf.all <- list() # keep all resting data for boxplots later
    indiv.means <- c()
    indiv.sd <- c()
    
    for(i in unique(df$iPOP_ID)){
      subDf <- df[which(df$iPOP_ID %in% i),] #pull data per individual
      restingST<-c()
      restingST <- rollapplyr(subDf$Steps, width=window, by=1, sum, partial=TRUE)
      restingST[1:window-1]<-"NA" # remove 1st x observations because we dont know what happened prior to putting the watch on
      restingST <- as.numeric(restingST)  # Expected warning: "In as.numeric(restingST) : NAs introduced by coercion"
      restingDf <- subDf[restingST<maxsteps & !is.na(restingST)] # in the previous time window of X min plus the current minute,there are < maxsteps steps 
      indiv.means[i] <- mean(restingDf$Heart_Rate, na.rm=TRUE) # mean RHR for all days/times for individual i
      indiv.sd[i] <- sd(restingDf$Heart_Rate, na.rm=TRUE) # RHR var for all days/times for individual i
      restingDf.all[[i]] <- restingDf # store all resting data for boxplots
    }
    
    #Combine all subject IDs into restingDf.all
    restingDf.all <- rbindlist(restingDf.all)
    restingDf.all$Date <- as.Date(restingDf.all$Timestamp_Local)
    restingDf.all <- restingDf.all[,c("iPOP_ID","Date","Heart_Rate","Skin_Temperature_F","Steps","Timestamp_Local")]
    names(restingDf.all) <- c("iPOP_ID","Date","restingHR","restingSkinTemp","restingSteps","DateTime")
    
    #Optional: use day-prior rather than day-of wearable data for comparison:
    if(dayPrior){
      restingDf.all$Date <- restingDf.all$Date + days(1)
    }
    # ^ All this does is advance each day by 1 if we set dayPrior to TRUE.
    #   This simple change allows us to reuse most of our code while testing.
    
    # We do want only the dates that overlap for now, so:
    #merge vitals (i.e., clinical HR measurments) with resting HR data from wearable BASIS
    #device
    restingDf.vitals <- merge(restingDf.all,vitals,by=c("iPOP_ID","Date"))
    restingDf.vitals$DateTime<-as.POSIXct(restingDf.vitals$DateTime)
    restingDf.vitals <- restingDf.vitals[order(restingDf.vitals$DateTime),] 
    cols <- c("restingHR","restingSkinTemp","Pulse") #subset columns to convert
    restingDf.vitals[,(cols) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = cols ]
    df.name <- paste0("restingDf.vitals.window=", window,",msteps=", maxsteps)
    assign(df.name, restingDf.vitals) # to store data frame for each resting window definition
    
    #Check progress:
    print(paste0("Completed -- window:", window, " maxsteps:", maxsteps))
  }
  
}

# Gather aggregate daily wearable HR BASIS and CLINIC HR 

#Suppress warnings (but return to show warnings at end of this chunk):
original_warn <- getOption("warn")
options(warn = -1)

#Initialize subDfClinDay:

subDfClinDay <-list()
Delta_hr_parameters <- data.frame(matrix
                                  (vector(), 0, 13,
                                    dimnames=list(c(), c("numObs", 
                                                         "init_hour",
                                                         "end_hour",
                                                         "maxsteps",
                                                         "window",
                                                         "Mean_wHR",
                                                         "Median_wHR",
                                                         "Stdev_wHR",
                                                         "Mean_Clinic_HR",
                                                         "Mean_wTemp",
                                                         "Median_wTemp",
                                                         "Stdev_wTemp",
                                                         "Mean_Clinic_Temp"))),
                                  stringsAsFactors=F)

row_ind = 0

#Establish 3 part nested-loop with varying windows, maxsteps, and iterating across 6am-10am 2 hr time windows:
windows = c(10, 60)
steps = c(1, 50)
time = 'All'
hours = c(0:22)

for(window in windows){
  
  for(maxsteps in steps){
    
    restingDf.vitals <- eval(as.name(paste0("restingDf.vitals.window=", window,",msteps=", maxsteps)))
    
    #Set up sequence of starting time (hour) during day to iterate over for moving 2 hr time window
    #Note: we achieve 2 hr window with use of "init_hour" and (init_hour + 2) as window...
    
    par(mfrow=c(2,4))
    for (init_hour in hours){ # vary the hour of the day from 6am to 10am as start of window
      
      row_ind = row_ind + 1
      
      options(datatable.optimize=1)
      subDfClinDay[[row_ind]] <- with( restingDf.vitals, restingDf.vitals[ hour( DateTime ) >= init_hour & hour( DateTime ) < (init_hour+2) , ] ) # pull data only from specific time window; store hourly resting data for boxplots
      rhr.daily.means <- subDfClinDay[[row_ind]][, lapply(.SD, mean, na.rm= TRUE), by=c("iPOP_ID","Date")]
      restingDf.compare <- cbind(rhr.daily.means$restingHR, rhr.daily.means$Pulse) 
      restingDf.compare <- data.frame(restingDf.compare)
      colnames(restingDf.compare) <- c("Resting_wHR", "cHR")
      numObs <-dim(rhr.daily.means)[1]
      rhr.daily.means.id <- cbind(rhr.daily.means$iPOP_ID, rhr.daily.means$restingHR, rhr.daily.means$Pulse)
      rhr.daily.means.id <- as.data.frame(rhr.daily.means.id)
      rhr.daily.means.id[,1 ]<- as.factor(as.character(rhr.daily.means.id[,1 ])); 
      rhr.daily.means.id[,2 ]<- as.numeric(as.character(rhr.daily.means.id[,2 ])); 
      rhr.daily.means.id[,3 ]<- as.numeric(as.character(rhr.daily.means.id[,3 ])); 
      colnames(rhr.daily.means.id) <- c("iPOP_ID", "restingHR", "Pulse")
      
      title=(paste0((init_hour), " - ",(init_hour+2), " am\n Num Obs = ", numObs))
      sig.test<-t.test(restingDf.vitals$restingHR, 
                       restingDf.vitals$Pulse, 
                       var.equal=FALSE, 
                       paired=TRUE) # test wHR v cHR for significance
      
      #Calculate aggregate metrics per day at that given time window:
      
      mean_HR_BASIS = mean(rhr.daily.means[[3]], na.rm = TRUE)
      median_HR_BASIS = median(rhr.daily.means[[3]], na.rm = TRUE)
      stdev_HR_BASIS = sd(rhr.daily.means[[3]], na.rm = TRUE)
      
      mean_Temp_BASIS = mean(rhr.daily.means[[4]], na.rm = TRUE)
      median_Temp_BASIS = median(rhr.daily.means[[4]], na.rm = TRUE)
      stdev_Temp_BASIS = sd(rhr.daily.means[[4]], na.rm = TRUE)
      
      mean_HR_clinic = mean(rhr.daily.means[[8]], na.rm = TRUE)
      mean_Temp_clinic = mean(rhr.daily.means[[9]], na.rm = TRUE)
      
      #Assemble these metric statistics in a dataframe:
      
      Delta_hr_parameters[row_ind, "Mean_wHR"] <- mean_HR_BASIS
      Delta_hr_parameters[row_ind, "Median_wHR"] <- median_HR_BASIS
      Delta_hr_parameters[row_ind, "Stdev_wHR"] <- stdev_HR_BASIS
      Delta_hr_parameters[row_ind, "Mean_Clinic_HR"] <- mean_HR_clinic
      Delta_hr_parameters[row_ind, "Mean_wTemp"] <- mean_Temp_BASIS
      Delta_hr_parameters[row_ind, "Median_wTemp"] <- median_Temp_BASIS
      Delta_hr_parameters[row_ind, "Stdev_wTemp"] <- stdev_Temp_BASIS
      Delta_hr_parameters[row_ind, "Mean_Clinic_Temp"] <- mean_Temp_clinic
      Delta_hr_parameters[row_ind, "numObs"] <- numObs
      Delta_hr_parameters[row_ind, "init_hour"] <- init_hour
      Delta_hr_parameters[row_ind, "end_hour"] <- init_hour + 2
      Delta_hr_parameters[row_ind, "maxsteps"] <- maxsteps
      Delta_hr_parameters[row_ind, "window"] <- window
    }
  }
}

#Export file to csv and plot in Tableau:
write.csv(Delta_hr_parameters, file = "BASIS_mean_stdev_median_24_hr_wHR_and_Temp.csv")

############
#  Fig 1C  #
############
windows=c(60, 120, 180, 240) # define time windows with no steps for resting threshold (10,60,120, etc)

dayPrior = FALSE
idx=0
HR.personal.sd	<- c()
wRHR.mean	<- c()
wRHR.sd	<- c()
wRHR.num.obs	<- c()
wRTemp.mean	<- c()
wRTemp.sd	<- c()
wRTemp.num.obs	<- c()
Temp.personal.sd	<- c()
for (window in windows){
  idx=idx+1
  restingDf <- c() 
  restingDf.all <- list() # keep all resting data for boxplots later
  maxsteps <- 1 #define max steps for resting threshold
  indiv.means <- c()
  indiv.sd <- c()
  
  for(i in unique(daytime.df$iPOP_ID)){
    subDf <- daytime.df[which(daytime.df$iPOP_ID %in% i),] #pull data per individual
    if (dim(subDf)[1] > window) {
    print(i)
    restingST<-c()
    restingST <- rollapplyr(subDf$Steps, width=window, by=1, sum, partial=TRUE)
    restingST[1:window-1]<-"NA" # remove 1st x observations because we dont know what happened prior to putting the watch on
    restingST <- as.numeric(restingST)  # Expected warning: "In as.numeric(restingST) : NAs introduced by coercion"
    restingDf <- subDf[restingST<maxsteps & !is.na(restingST)] # in the previous time window of X min plus the current minute,there are < maxsteps steps 
    indiv.means[i] <- mean(restingDf$Heart_Rate, na.rm=TRUE) # mean RHR for all days/times for individual i
    indiv.sd[i] <- sd(restingDf$Heart_Rate, na.rm=TRUE) # RHR var for all days/times for individual i
    restingDf.all[[i]] <- restingDf # store all resting data for boxplots
    }
  }
  
  restingDf.all <- rbindlist(restingDf.all)
  restingDf.all$Date <- as.Date(restingDf.all$Timestamp_Local)
  restingDf.all <- restingDf.all[,c("iPOP_ID","Date","Heart_Rate","Skin_Temperature_F","Steps","Timestamp_Local")]
  names(restingDf.all) <- c("iPOP_ID","Date","restingHR","restingSkinTemp","restingSteps","DateTime")
  means.by.id<-aggregate(restingDf.all, list(restingDf.all$iPOP_ID), na.omit(mean)) 
  sd.by.id<-aggregate(restingDf.all, list(restingDf.all$iPOP_ID), sd) 
  HR.personal.sd[idx] <- mean(na.omit(sd.by.id$restingHR)) # personal / intra-individual SD
  wRHR.mean[idx] <- mean(na.omit(restingDf.all$restingHR)) # wRHR; mean 65.23 +/- 10.41, n=1,198,040 measurements
  wRHR.sd[idx] <-sqrt(var((na.omit(restingDf.all$restingHR))))  # sd wRHR
  wRHR.num.obs[idx] <- length(na.omit(restingDf.all$restingHR))
  wRTemp.mean[idx] <- mean(na.omit(restingDf.all$restingSkinTemp)) # wRTemp; 90.89 +/- 3.09, n=1,181,648 measurements
  wRTemp.sd[idx] <-sd(na.omit(restingDf.all$restingSkinTemp)) # sd wRTemp
  wRTemp.num.obs[idx] <- length(na.omit(restingDf.all$restingSkinTemp))
  Temp.personal.sd[idx] <- mean(na.omit(sd.by.id$restingSkinTemp))
  
  #Optional: use day-prior rather than day-of wearable data for comparison:
  if(dayPrior){
    restingDf.all$Date <- restingDf.all$Date + days(1)
  }

  restingDf.vitals <- merge(restingDf.all,vitals,by=c("iPOP_ID","Date"))
  restingDf.vitals$DateTime<-as.POSIXct(restingDf.vitals$DateTime)
  restingDf.vitals <- restingDf.vitals[order(restingDf.vitals$DateTime),] 
  cols <- c("restingHR","restingSkinTemp","Pulse","Temp") #subset columns to convert
  restingDf.vitals[,(cols) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = cols]
  df.name <- paste0("restingDf.vitals.", window)
  assign(df.name, restingDf.vitals) # to store data frame for each resting window definition
}
  
cols <- c("red","blue","green","purple")
ggplot()+
  geom_line(aes(x=windows, y=wRHR.sd), color="red") +
  geom_point(aes(x=windows, y=wRHR.sd), color="red") +
  geom_line(aes(x=windows, y=HR.personal.sd), color="blue") +
  geom_point(aes(x=windows, y=HR.personal.sd), color="blue") +
  #geom_line(aes(x=windows, y=rep(cHR.sd, length(windows))), color="coral") +
  #geom_point(aes(x=windows, y=rep(cHR.sd, length(windows))), color="coral") +
  #geom_line(aes(x=windows, y=rep(cHR.individual.sd, length(windows))), color="skyblue") +
  #geom_point(aes(x=windows, y=rep(cHR.individual.sd, length(windows))), color="skyblue") +
  xlab("Resting Time Window (seconds)") +
  ylab("HR SD") # or wRTemp

for(window in windows){
  restingDf.vitals <- eval(as.name(paste0("restingDf.vitals.",window)))

  # resting values to go in the text of the manuscript
  mean(restingDf.vitals$restingHR)
  sd(restingDf.vitals$restingHR)
  mean(na.omit(restingDf.vitals$restingSkinTemp))
  sd(na.omit(restingDf.vitals$restingSkinTemp))

  # wRHR from *specific day* of clinic visit only; aggregate wRHR into daily values corresponding to the clinic date
  # compare all resting watch data with all vitals data - this was fixed to aggregate wRHR into daily values corresponding to the clinic date
  options(datatable.optimize=1) #need this here; otherwsie won't handle character class
  rhr.daily.means <- restingDf.vitals[, lapply(.SD, mean, na.rm=TRUE), by=c("iPOP_ID","Date")]
  options(datatable.optimize=0)
  numObs <- dim(rhr.daily.means)[1]
  numPeople <- length(unique(rhr.daily.means$iPOP_ID))
  restingDf.compare <- cbind(rhr.daily.means$restingHR, rhr.daily.means$Pulse, rhr.daily.means$restingSkinTemp, rhr.daily.means$Temp)
  colnames(restingDf.compare) <- c("Resting wHR", "cHR", "Resting wSkinTemp", "cTemp")
  rhr.daily.means.id <- cbind(rhr.daily.means$iPOP_ID, rhr.daily.means$restingHR, rhr.daily.means$Pulse, rhr.daily.means$restingSkinTemp, rhr.daily.means$Temp)
  rhr.daily.means.id <- as.data.frame(rhr.daily.means.id)
  rhr.daily.means.id[,1 ]<- as.factor(as.character(rhr.daily.means.id[,1 ])); 
  rhr.daily.means.id[,2 ]<- as.numeric(as.character(rhr.daily.means.id[,2 ])); 
  rhr.daily.means.id[,3 ]<- as.numeric(as.character(rhr.daily.means.id[,3 ])); 
  rhr.daily.means.id[,4 ]<- as.numeric(as.character(rhr.daily.means.id[,4 ]));
  rhr.daily.means.id[,5 ]<- as.numeric(as.character(rhr.daily.means.id[,5 ])); 
  colnames(rhr.daily.means.id) <- c("iPOP_ID", "restingHR", "Pulse", "restingSkinTemp", "Temp")
  
  means<-aggregate(rhr.daily.means.id, list(rhr.daily.means.id$iPOP_ID), mean) # check this compare to indiv means
  sd<-aggregate(rhr.daily.means.id, list(rhr.daily.means.id$iPOP_ID), sd)
  delta.daily.mean.RHR.pulse<-means$Pulse - means$restingHR
  delta.daily.sd.RHR.pulse<-sd$Pulse - sd$restingHR 
  hist(delta.daily.mean.RHR.pulse, col="darkred", main = paste0("Delta Mean cHR - Mean wRHR)"))
  hist(delta.daily.sd.RHR.pulse, col="darkred", main = paste0("Delta StdDev cHR - StdDev wRHR)"))
  
  rhr.daily.means.id$idx <- as.numeric(rhr.daily.means.id$iPOP_ID)
}

###########################
#     Main Text Numbers   #
###########################
wRHR.mean
wRHR.sd
wRHR.num.obs
HR.personal.sd
wRTemp.mean
wRTemp.sd
wRTemp.num.obs
Temp.personal.sd
  
  
##########
# Fig 1D #
##########
  ## Figure 1D, left
  hist(iPOPdaysMonitored$Days_monitored_by_clinic, col="grey", breaks=10,
       xlab = "Time Monitored by Clinic (Days)", main = NULL, font.lab=2,lwd=2,font=2)
  rug(iPOPdaysMonitored$Days_monitored_by_clinic, ticksize = 0.1, lwd = 1)
  
  ## Figure 1D, right
  hist(iPOPdaysMonitored$Total_NumOfClinMeasures, col="grey", breaks=10,
       xlab = "Number of Clinic Visits / Person", main = NULL, font.lab=2,lwd=2,font=2)
  rug(iPOPdaysMonitored$Total_NumOfClinMeasures, ticksize = 0.1, lwd = 1)

############
# Fig 1E   #
############

## Fig 1E, top left
par(mfrow = c(2,2), mai = c(0.7, 0.7, 0.7, 0.7))
hist(iPOPvitals$Pulse, col="tomato3", , border="tomato4", breaks=1000,
     xlab = "cHR", xlim=c(50,110),
     main = NULL, font.lab=2,lwd=2,font=2)
abline(v=mean(iPOPvitals$Pulse[!is.na(iPOPvitals$Pulse)]), col="blue", lwd=2)

## Fig 1E, top right
hist(iPOPvitals$Temp, col="turquoise3", border="turquoise4", breaks=100,
     xlab = "cTemp", xlim=c(65,105),
     main = NULL, font.lab=2,lwd=2,font=2)
abline(v=mean(iPOPvitals$Temp[!is.na(iPOPvitals$Temp)]), col="blue", lwd=2)

#Fig 1E, bottom left
hist(restingDf.all$restingHR, col="tomato3", border="tomato4", breaks=100,
     xlab = "wRHR", xlim=c(50,100),
     main = NULL, font.lab=2,lwd=2,font=2)
scale_y_continuous()
abline(v=mean(restingDf.all$restingHR[!is.na(restingDf.all$restingHR)]), col="blue", lwd=2)

#Fig 1E, bottom right
hist(restingDf.all$restingSkinTemp, col="turquoise3", border="turquoise4", breaks=100,
     xlab = "wRTemp", xlim=c(65,105),
     main = NULL, font.lab=2,lwd=2,font=2)
abline(v=mean(restingDf.all$restingSkinTemp[!is.na(restingDf.all$restingSkinTemp)]), col="blue", lwd=2)

###########################
#     Main Text Numbers   #
###########################

length(iPOPvitals$Pulse[!is.na(iPOPvitals$Pulse)]) # number of cHR measurements in iPOP cohort
mean(iPOPvitals$Pulse[!is.na(iPOPvitals$Pulse)]) # mean cHR; 71.54 +/- 9.92, n=1644
cHR.sd <- sqrt(var(iPOPvitals$Pulse[!is.na(iPOPvitals$Pulse)])) # stdev of cHR
means<-aggregate(iPOPvitals$Pulse,list(iPOPvitals$iPOP_ID), mean) # check this compare to indiv means
sd<-aggregate(iPOPvitals$Pulse, list(iPOPvitals$iPOP_ID), sd)
# personal sd:
cHR.individual.sd <- mean(na.omit(sd$x)) # intra-individual SD 6.913

length(iPOPvitals$Temp[!is.na(iPOPvitals$Temp)]) # number of cTemp measurements in iPOP cohort
mean(iPOPvitals$Temp[!is.na(iPOPvitals$Temp)]) # mean cTemp; 97.84 +/- 0.38, n=1136
cTemp.sd <- sqrt(var(iPOPvitals$Temp[!is.na(iPOPvitals$Temp)])) # stdev of cTemp
means<-aggregate(iPOPvitals$Temp,list(iPOPvitals$iPOP_ID), mean) # check this compare to indiv means
sd<-aggregate(iPOPvitals$Temp, list(iPOPvitals$iPOP_ID), sd)
cTemp.individual.sd <- mean(na.omit(sd$x)) # intra-individual SD 0.2536

#characterize the iPOP data set
length(na.omit(iPOPvitals$Temp)) + length(na.omit(iPOPvitals$Pulse)) # total number of clinical vital signs measured
describe(iPOPlabs[names(iPOPlabs) %in% allClin]) # summary of clinical labs data
length(unique(wear$iPOP_ID)) # num people in iPOP wearables dataset

#############################
#    Suppl. Table 1A and B  #
#############################
# make table for vitals (Suppl. Table 1A)
describe(iPOPvitals)
# make table for labs (Suppl. Table 1B)
describe(iPOPlabs)


#################
#  Fig 3A - LM  #
#################

# creates ranked list of clinical laboratory tests by the %var explained in simple LM; LOO cross validation at the subject level 

# choose for during troubleshooting
use.Troubleshoot.mode = TRUE
#choose whether Demographics in models (supply TRUE or FALSE)
use.Demog <- FALSE
#choose whether iPOP_ID variable is used (supply TRUE or FALSE)
use.iPOP <- FALSE

if (use.Troubleshoot.mode){
  #   top.names<-c("MONOAB", "HGB", "HCT", "RBC") # "RBC", "PLT") # for testing model on small subset
  # top.names<-c("ALT")# "TGL", "BASOAB", "EOSAB") # for testing model on small subset
   top.names <- c("HGB", "TGL")#, "BASOAB", "EOSAB")
  #top.names <- c("PLT")
}

####
# CODE FOR SIMPLE LM
#
sum.vectors.in.list <- function(x) {
  c(sum(na.omit(x)))
}

rsq.all = c()
pct.var.all = c()
iPOPcorDf.demo <- merge(iPOPcorDf, iPOPdemographics[1:4], by="iPOP_ID")
iPOPcorDf.demo$Gender <- as.factor(iPOPcorDf.demo$Gender)
iPOPcorDf.demo$Ethn <- as.factor(iPOPcorDf.demo$Ethn)

#change gender and ethnicity to dummy variables
gender <- data.frame(model.matrix( ~ Gender - 1, data=iPOPcorDf.demo))
ethn <- data.frame(model.matrix( ~ Ethn - 1, data=iPOPcorDf.demo))

#remove the least populated gender and ethnicity (NCOL-1)
cache <- names(gender)[which(sapply(gender,sum)==max(sapply(gender,sum)))]
gender <- data.frame(cache=gender[which(sapply(gender,sum)==max(sapply(gender,sum)))])
ethn <- ethn[,-which(sapply(ethn,sum)==min(sapply(ethn,sum)))]
#cbind new gender and ethnicity variables to "wear"
iPOPcorDf.demo <- cbind(iPOPcorDf.demo,gender,ethn)

#store names as vitals.variables
#vitals.variables <- c("Pulse", "Temp", "AgeIn2016", names(gender), names(ethn)) # "BMI", "systolic", "diastolic", 
if(use.Demog){
  vitals.variables <- c("Pulse", "Temp", "AgeIn2016", "Gender", "Ethn") # "BMI", "systolic", "diastolic", 
} else if(!use.Demog) {
  vitals.variables <- c("Pulse", "Temp") #
}

patients = unique(iPOPcorDf$iPOP_ID)

val.true <- rep(list(NA),length(top.names)) # list of vectors to store true values; each vector is for 1 clinical lab
val.pred <- rep(list(NA),length(top.names)) # list of vectors to store trained model predicted values; each vector is for 1 clinical lab
val.null.pred <- rep(list(NA),length(top.names)) # list of vectors to store null model predicted values; each vector is for 1 clinical lab
num.true <- rep(list(NA),length(top.names)) # number of observations per individual per clinic test

for (k in 1:length(patients)){
  train <- patients[patients != patients[k]]
  test <- patients[patients == patients[k]]
  cat("Patient",patients[k],"\n") # LOO
  dat.train.unsorted<-iPOPcorDf.demo[ iPOPcorDf.demo$iPOP_ID %in% train, ] # subset input data by training set
  dat.train <- dat.train.unsorted[order(dat.train.unsorted$iPOP_ID),] #order by iPOP_ID in order to supply correct nfolds arg to glmnet
  dat.test<-iPOPcorDf.demo[ iPOPcorDf.demo$iPOP_ID %in% test, ] # subset input data by testing set
  
  for (l in 1:length(top.names)){
    cat("Test",top.names[l],"\n")
    # create training set
    x.train<-dat.train[,colnames(dat.train) %in% c("iPOP_ID",top.names[l], vitals.variables)] # subset input data by lab: only take current lab test of interest
    x.train<- na.omit(x.train) # skip nas and nans ## TODO: the way this script is written, you will lose a lot of data because you take the number of lab visits down to the test with the minimum number of visits. However, if you do na.omit after the next line, you have to change your matrix to accept dynamic number of row entries. 
    x.train.ids<-x.train$iPOP_ID
    x.train<-x.train[,-1]
    predictors <- as.data.frame(x.train[,colnames(x.train) %in% c(vitals.variables)]) # later add in demographics
    outcome <- as.matrix(x.train[,colnames(x.train) %in% top.names[l]]) # matrix of outcome for model building # tried adding as.numeric after as.matrix() but that introduced new issues
    # create test set
    x.test<-dat.test[,colnames(dat.test) %in% c(top.names[l], vitals.variables)] # subset input data by lab: only take current lab test of interest
    x.test<- na.omit(x.test) # skip nas and nans ## TODO: SEE ABOVE na.omit FOR ISSUE WITH THIS
    val.true[[l]] = c(val.true[[l]], x.test[,top.names[l]]) # true values of left out person
    num.true[[l]]<-c(num.true[[l]],length(x.test[,top.names[l]])) # number of test observations for left out person
    fml = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(vitals.variables,collapse=" + "))
    fml.null = paste(top.names[l]," ~ 1")
    bivar.lm.model = lm(as.formula(fml), data = x.train) # build the model
    val.pred[[l]] = c(val.pred[[l]], predict(bivar.lm.model, newdata = x.test)) # predict on trained model
    bivar.null.lm.model<-lm(as.formula(fml.null), data = x.train) # create null model for hypothesis testing and for calculating RSS0
    val.null.pred[[l]] = c(val.null.pred[[l]], predict(bivar.null.lm.model, newdata = x.test)) # predict on null model
    # t<- anova(bivar.null.lm.model, bivar.lm.model) # to get p-values for model
    # p.value[[l]] <- as.numeric(t[2,][["Pr(>F)"]])  # to get p-values for model
  }
}
num.test.obs <- lapply(num.true, sum.vectors.in.list)

rsq.vitals = c()
p.val.rsq.vitals = c()
rssm.vitals = c()
rss0.vitals = c()
pct.var.explained = c()
num.Records.check <- c()
for (j in 1:length(top.names)){
  rsq.vitals = c(rsq.vitals, cor(val.pred[[j]], val.true[[j]], use = "complete.obs"))
  p.val.rsq.vitals = c(p.val.rsq.vitals, cor.test(val.pred[[j]], val.true[[j]], use = "complete.obs")$p.value)
  rssm.vitals = sum(na.omit((val.true[[j]] - val.pred[[j]])^2))
  rss0.vitals = sum(na.omit((val.true[[j]] - val.null.pred[[j]])^2))
  pct.var.explained = c(pct.var.explained, (1 - ( rssm.vitals / rss0.vitals )))
  num.Records.check <- c(num.Records.check, (length(val.pred[[j]])-1)) # same as num.test.obs
}
names(rsq.vitals) = top.names
names(pct.var.explained) = top.names
sqrt.pct.var <- sqrt(pct.var.explained)
as.matrix(sort(sqrt(pct.var.explained)))


##############################
# Fig 3A: CODE FOR LASSO, RF # 
##############################

#clean wear data frame
wear[,8:length(names(wear))] <- apply(
  wear[,8:length(names(wear))], 2,
  function(x) as.numeric(as.character(x)))
wear$Gender <- as.factor(wear$Gender)
wear$Ethn <- as.factor(wear$Ethn)
wear.variables <- unlist(read.table("FinalLasso_153WearableFactors.csv", stringsAsFactors = FALSE)) # the table of model features we want to work with

#change gender and ethnicity to dummy variables
gender <- data.frame(model.matrix( ~ Gender - 1, data=wear))
ethn <- data.frame(model.matrix( ~ Ethn - 1, data=wear))

#remove the least populated gender and ethnicity (NCOL-1)
cache <- names(gender)[which(sapply(gender,sum)==max(sapply(gender,sum)))]
gender <- data.frame(cache=gender[which(sapply(gender,sum)==max(sapply(gender,sum)))])
ethn <- ethn[,-which(sapply(ethn,sum)==min(sapply(ethn,sum)))]

#store names as demo.variables
demo.variables <- c("AgeIn2016", names(gender), names(ethn))

#cbind new gender and ethnicity variables to "wear"
wear <- cbind(wear,gender,ethn)

val.true <- rep(list(NA),length(top.names)) # list of vectors to store true values; each vector is for 1 clinical lab
null.val.pred <- rep(list(NA),length(top.names))  # list of vectors to store nullmodel-predicted values; each vector is for 1 clinical lab
lasso.val.pred.lambda.manual <- rep(list(NA),length(top.names)) # list of vectors to store lasso-trainedmodel-predicted values; each vector is for 1 clinical lab
lasso.val.pred.lambda.min <- rep(list(NA),length(top.names)) # list of vectors to store lasso-trainedmodel-predicted values; each vector is for 1 clinical lab
lasso.val.pred.lambda.1se <- rep(list(NA),length(top.names)) # list of vectors to store lasso-trainedmodel-predicted values; each vector is for 1 clinical lab
rf.val.pred <- rep(list(NA),length(top.names))  # list of vectors to store rf-trainedmodel-predicted values; each vector is for 1 clinical lab

num.Records <- c()
lasso.features.lambda.manual <- data.frame("test"=character(),"cv.run"=character(),
                                           "left.out.person"=character(),"lasso.feature"=character(),
                                           "lasso.coef.value"=character())
lasso.features.lambda.min <- data.frame("test"=character(),"cv.run"=character(),
                                        "left.out.person"=character(),"lasso.feature"=character(),
                                        "lasso.coef.value"=character())
lasso.features.lambda.1se <- data.frame("test"=character(),"cv.run"=character(),
                                        "left.out.person"=character(),"lasso.feature"=character(),
                                        "lasso.coef.value"=character())
rf.features <- data.frame("test"=character(),"cv.run"=character(),
                          "left.out.person"=character(),"rf.feature"=character(),
                          "rf.coef.value"=character())

if(use.Demog){
    demo.variables <- c("AgeIn2016", names(gender), names(ethn))
  } else if(!use.Demog) {
      demo.variables <- c()
}

for (k in 1:length(patients)){
  train <- patients[patients != patients[k]]
  test <- patients[patients == patients[k]]
  cat("Patient",patients[k],"\n") # LOO
  #set up iPOP dummy variable if use.iPOP=TRUE
  if(use.iPOP){
    cache <- data.frame(model.matrix( ~ iPOP_ID - 1, data=wear))
    cache <- cache[,-which(names(cache)==paste0("iPOP_ID",gsub("-",".",test)))]
    wear <- cbind(wear,cache)
    demo.variables <- c(demo.variables,names(cache)) #append iPOPs to demo
  } 
  dat.train.unsorted <- wear[ wear$iPOP_ID %in% train, ] # subset input data by training set
  dat.train <- dat.train.unsorted[order(dat.train.unsorted$iPOP_ID),] #order by iPOP_ID in order to supply correct nfolds arg to glmnet
  dat.test<-wear[ wear$iPOP_ID %in% test, ] # subset input data by testing set
  for (l in 1:length(top.names)){
    cat("Test",top.names[l],"\n")
    x.train<-dat.train[,colnames(dat.train) %in% c("iPOP_ID", top.names[l], wear.variables, demo.variables)] # subset input data by lab: only take current lab test of interest
    x.train<-na.omit(x.train) # skip nas and nans ## TODO: the way this script is written, you will lose a lot of data because you take the number of lab visits down to the test with the minimum number of visits. However, if you do na.omit after the next line, you have to change your matrix to accept dynamic number of row entries. Not sure how to do this yet, so for now just reducing the data amount by a lot. 
    x.train.ids<-x.train$iPOP_ID
    x.train<-x.train[,-1] 
    # if(!NROW(x.train)){ #if x.train is empty
    #   print(paste0("The x.train data was empty for ",patients[k],"'s ",top.names[l]," test."))
    # } else {
    #   print(paste0("The x.train data for ",patients[k],"'s ",top.names[l]," test had ",NROW(x.train)," observations."))
    # }
    predictors <- as.data.frame(x.train[,colnames(x.train) %in% c(wear.variables, demo.variables)]) # later add in demographics
    outcome <- as.matrix(x.train[,colnames(x.train) %in% top.names[l]]) # matrix of outcome for model building # tried adding as.numeric after as.matrix() but that introduced new issues
    
    # create test set
    x.test<-dat.test[,colnames(dat.test) %in% c(top.names[l], wear.variables, demo.variables)] # subset input data by lab: only take current lab test of interest
    x.test<- na.omit(x.test) # skip nas and nans ## TODO: SEE ABOVE na.omit FOR ISSUE WITH THIS
    # if(!NROW(x.test)){ #if x.test is empty
    #   print(paste0("The x.test data was empty for ",patients[k],"'s ",top.names[l]," test."))
    # } else {
    #   print(paste0("The x.test data for ",patients[k],"'s ",top.names[l]," test had ",NROW(x.test)," observations."))
    # }
    val.true[[l]] = c(val.true[[l]], x.test[,top.names[l]]) # true values of left out person
    
    num.Records <- rbind(num.Records, c(IPOP_ID=patients[k], test=top.names[l], TrainingObs=length(outcome), TestObs=length(x.test[,top.names[l]])))
    
    rf.variables.to.use = c(wear.variables, demo.variables) # rf variables (use all)
    
    #decide on number for nfolds from number of obs per subject
    frq <- as.vector(table(x.train.ids))
    
    #optional argument for leave-one-out CV method
    n <- length(frq)
    
    #optional argument to specify folds for CV
    folds <- rep(1:length(frq),frq[1:length(frq)])
    
    ## run lasso for variable selection
    # n <- as.numeric(length(outcome)) #optional argument for leave-one-out CV method for nfold
    x_train <- model.matrix( ~ .-1, as.data.frame(predictors))
    glm.res = cv.glmnet(x=x_train,y=outcome,
                        standardize.response=FALSE,
                        family="gaussian",
                        nfolds=n,
                        foldid=folds,
                        nlambda=100)
 
#    lasso.model.lambda.min = lm(as.formula(lasso.fml.lambda.min), data = x.train) # , weights = labs.wear$weight)   
    
    #store all lasso variable coefs (lambda specific: manual, min, and 1se)
    factors.lambda.manual = glm.res$glmnet.fit$beta[,25] # TODO: this is an arbitrary rule for now
    lasso.variables.lambda.manual = factors.lambda.manual#[abs(factors.lambda.manual)!=0]
    factors.lambda.min <- glm.res$glmnet.fit$beta[,which(glm.res$glmnet.fit$lambda==glm.res$lambda.min)]
    lasso.variables.lambda.min = factors.lambda.min#[abs(factors.lambda.min)!=0]
    factors.lambda.1se <- glm.res$glmnet.fit$beta[,which(glm.res$glmnet.fit$lambda==glm.res$lambda.1se)]
    lasso.variables.lambda.1se = factors.lambda.1se#[abs(factors.lambda.1se)!=0]
    
    ## pull out features from lasso models ##
    
    tmp <- data.frame("test"=rep(top.names[l],length(lasso.variables.lambda.manual)),
                      "cv.run"=rep(k,length(lasso.variables.lambda.manual)),
                      "left.out.person"=rep(patients[k],length(lasso.variables.lambda.manual)),
                      "lasso.feature"=names(lasso.variables.lambda.manual),
                      "lasso.coef.value"=as.numeric(lasso.variables.lambda.manual))
    lasso.features.lambda.manual <- rbind(lasso.features.lambda.manual,tmp)
    
    tmp <- data.frame("test"=rep(top.names[l],length(lasso.variables.lambda.min)),
                      "cv.run"=rep(k,length(lasso.variables.lambda.min)),
                      "left.out.person"=rep(patients[k],length(lasso.variables.lambda.min)),
                      "lasso.feature"=names(lasso.variables.lambda.min),
                      "lasso.coef.value"=as.numeric(lasso.variables.lambda.min))
    lasso.features.lambda.min <- rbind(lasso.features.lambda.min,tmp)
    
    tmp <- data.frame("test"=rep(top.names[l],length(lasso.variables.lambda.1se)),
                      "cv.run"=rep(k,length(lasso.variables.lambda.1se)),
                      "left.out.person"=rep(patients[k],length(lasso.variables.lambda.1se)),
                      "lasso.feature"=names(lasso.variables.lambda.1se),
                      "lasso.coef.value"=as.numeric(lasso.variables.lambda.1se))
    lasso.features.lambda.1se <- rbind(lasso.features.lambda.1se,tmp)
    
    #store lasso variable names based on coef threshold (lambda specific: manual, min, and 1se)
    lasso.variables.to.use.lambda.manual = names(factors.lambda.manual[abs(factors.lambda.manual)>1e-10]) # TODO: this is an arbitrary rule for now
    lasso.variables.to.use.lambda.min = names(factors.lambda.min[abs(factors.lambda.min)>1e-10])
    lasso.variables.to.use.lambda.1se = names(factors.lambda.1se[abs(factors.lambda.1se)>1e-10])
    
    # build null, lasso, and rf models
    set.seed(1)
    null.fml = paste(top.names[l]," ~ 1")
    null.model = lm(as.formula(null.fml), data = x.train) # create null model for hypothesis testing and for calculating RSS0
    null.val.pred[[l]] = c(null.val.pred[[l]], predict(null.model, newdata = x.test)) # predict on null model
    
    lasso.fml.lambda.manual = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(lasso.variables.to.use.lambda.manual,collapse=" + "))
    # lasso.model.lambda.manual = lm(as.formula(lasso.fml.lambda.manual), data = x.train) # , weights = labs.wear$weight) 
    # lasso.val.pred.lambda.manual[[l]] = c(lasso.val.pred.lambda.manual[[l]], predict(lasso.model.lambda.manual, newdata = x.test)) # predict on trained model
    # check that the formula is valid (i.e. not empty)
    if(lasso.fml.lambda.manual!=paste0("cbind( ",top.names[l]," ) ~ ")){
      lasso.model.lambda.manual = lm(as.formula(lasso.fml.lambda.manual), data = x.train) # , weights = labs.wear$weight)
      lasso.val.pred.lambda.manual[[l]] = c(lasso.val.pred.lambda.manual[[l]], predict(lasso.model.lambda.manual, newdata = x.test)) # predict on trained model
    } else {
      lasso.model.lambda.manual = null.model # if lasso sets all coeffs = 0, then use the null model
      lasso.val.pred.lambda.manual[[l]] = c(lasso.val.pred.lambda.manual[[l]], predict(lasso.model.lambda.manual, newdata = x.test)) # predict on trained model
      # cache <-  length(val.true[[l]])-length(lasso.val.pred.lambda.manual[[l]])
      # lasso.val.pred.lambda.manual[[l]] = c(lasso.val.pred.lambda.manual[[l]], rep(NA,cache)) # fill with NA(s) if invalid model was supplied
    }
    
    lasso.fml.lambda.min = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(lasso.variables.to.use.lambda.min,collapse=" + "))
    #check that the formula is valid (i.e. not empty)
    if(lasso.fml.lambda.min!=paste0("cbind( ",top.names[l]," ) ~ ")){
      lasso.model.lambda.min = lm(as.formula(lasso.fml.lambda.min), data = x.train) # , weights = labs.wear$weight) 
      lasso.val.pred.lambda.min[[l]] = c(lasso.val.pred.lambda.min[[l]], predict(lasso.model.lambda.min, newdata = x.test)) # predict on trained model
    } else {
      lasso.model.lambda.min = null.model # if lasso sets all coeffs = 0, then use the null model
      lasso.val.pred.lambda.min[[l]] = c(lasso.val.pred.lambda.min[[l]], predict(lasso.model.lambda.min, newdata = x.test)) # predict on trained model
      # cache <-  length(val.true[[l]])-length(lasso.val.pred.lambda.min[[l]])
      # lasso.val.pred.lambda.min[[l]] = c(lasso.val.pred.lambda.min[[l]], rep(NA,cache)) # fill with NA(s) if invalid model was supplied
    }
    
    lasso.fml.lambda.1se = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(lasso.variables.to.use.lambda.1se,collapse=" + "))
    #check that the formula is valid (i.e. not empty)
    if(lasso.fml.lambda.1se!=paste0("cbind( ",top.names[l]," ) ~ ")){
      lasso.model.lambda.1se = lm(as.formula(lasso.fml.lambda.1se), data = x.train) # , weights = labs.wear$weight) # 
      lasso.val.pred.lambda.1se[[l]] = c(lasso.val.pred.lambda.1se[[l]], predict(lasso.model.lambda.1se, newdata = x.test)) # predict on trained model
    } else {
      lasso.model.lambda.1se = null.model # if lasso sets all coeffs = 0, then use the null model
      lasso.val.pred.lambda.1se[[l]] = c(lasso.val.pred.lambda.1se[[l]], predict(lasso.model.lambda.1se, newdata = x.test)) # predict on trained model
      # cache <-  length(val.true[[l]])-length(lasso.val.pred.lambda.1se[[l]])
      # lasso.val.pred.lambda.1se[[l]] = c(lasso.val.pred.lambda.1se[[l]], rep(NA,cache)) # fill with NA(s) if invalid model was supplied
    }
    
    rf.fml = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(rf.variables.to.use,collapse=" + "))
    #check that the formula is valid (i.e. not empty)
    if(rf.fml!=paste0("cbind( ",top.names[l]," ) ~ ")){
      rf.model = randomForest(as.formula(rf.fml), data = x.train)  #weights = labs.wear$weight)
      rf.val.pred[[l]] = c(rf.val.pred[[l]], predict(rf.model, newdata = x.test)) # predict on left out person
    } else {
      # rf.model = null.model # if lasso sets all coeffs = 0, then use the null model
      # rf.val.pred[[l]] = c(rf.val.pred[[l]], predict(rf.model, newdata = x.test)) # predict on left out person
    }

    ## pull out features from rf models ##
    rf.features.list <- as.matrix(importance(rf.model)[order(importance(rf.model), decreasing=TRUE),])
    

    tmp <- data.frame("test"=rep(top.names[l],length(rf.features.list)),
                      "cv.run"=rep(k,length(rf.features.list)),
                      "left.out.person"=rep(patients[k],length(rf.features.list)),
                      "rf.feature"=unlist(dimnames(rf.features.list)),
                      "rf.coef.value"=as.numeric(rf.features.list))
    rf.features <- rbind(rf.features,tmp)
    # t<- anova(bivar.null.lm.model, bivar.lm.model) # to get p-values for model
    # p.value[[l]] <- as.numeric(t[2,][["Pr(>F)"]])  # to get p-values for model
  }
}

## calculate correlation coefficients and pct var explained by the models (lambda manual)
rsq.lasso.lambda.manual = c()
num.cor.pairs.lasso.lambda.manual = c()
p.val.rsq.lasso.manual = c()
rssm.lasso.lambda.manual = c()
rss0.lasso.lambda.manual = c()
lasso.pct.var.explained.lambda.manual = c()
lasso.num.Records.lambda.manual <- c()

rsq.lasso.lambda.min = c()
num.cor.pairs.lasso.lambda.min = c()
p.val.rsq.lasso.min = c()
rssm.lasso.lambda.min = c()
rss0.lasso.lambda.min = c()
lasso.pct.var.explained.lambda.min = c()
lasso.num.Records.lambda.min <- c()

rsq.lasso.lambda.1se = c()
num.cor.pairs.lasso.lambda.1se = c()
p.val.rsq.lasso.1se = c()
rssm.lasso.lambda.1se = c()
rss0.lasso.lambda.1se = c()
lasso.pct.var.explained.lambda.1se = c()
lasso.num.Records.lambda.1se <- c()

rsq.rf = c()
num.cor.pairs.rf = c()
p.val.rsq.rf = c()
rssm.rf = c()
rss0.rf = c()
rf.pct.var.explained = c()
rf.num.Records <- c()

#note: rf only runs once below

for (j in 1:length(top.names)){
  #lasso (lambda.manual)
  if(length(which(!is.na(lasso.val.pred.lambda.manual[[j]])))>=3){
    num.cor.pairs.lasso.lambda.manual <- c(num.cor.pairs.lasso.lambda.manual, length(which(!is.na(lasso.val.pred.lambda.manual[[j]]))))
    #insert step here to check if sample size for correlation test meets minimum threshold
    rsq.lasso.lambda.manual = c(rsq.lasso.lambda.manual, cor(lasso.val.pred.lambda.manual[[j]], val.true[[j]], use = "complete.obs"))
    p.val.rsq.lasso.manual <- c(p.val.rsq.lasso.manual, cor.test(lasso.val.pred.lambda.manual[[j]], val.true[[j]], use = "complete.obs")$p.value)
    #insert step here to check if correlation was significant?
    rssm.lasso.lambda.manual = sum(na.omit((val.true[[j]] - lasso.val.pred.lambda.manual[[j]])^2))
    rss0.lasso.lambda.manual = sum(na.omit((val.true[[j]] - null.val.pred[[j]])^2))
    lasso.pct.var.explained.lambda.manual = c(lasso.pct.var.explained.lambda.manual, (1 - ( rssm.lasso.lambda.manual / rss0.lasso.lambda.manual )))
    lasso.num.Records.lambda.manual <- c(num.Records, length(lasso.val.pred.lambda.manual[[j]]))
  } else {
    p.val.rsq.lasso.manual <- c(p.val.rsq.lasso.manual, NA)
    num.cor.pairs.lasso.lambda.manual <- c(num.cor.pairs.lasso.lambda.manual, 0)
    rsq.lasso.lambda.manual = c(rsq.lasso.lambda.manual, NA)
    lasso.pct.var.explained.lambda.manual = c(lasso.pct.var.explained.lambda.manual, NA)
    lasso.num.Records.lambda.manual <- c(num.Records, NA)
  }
  #lasso (lambda.min)
  if(length(which(!is.na(lasso.val.pred.lambda.min[[j]])))>=3){
    num.cor.pairs.lasso.lambda.min <- c(num.cor.pairs.lasso.lambda.min, length(which(!is.na(lasso.val.pred.lambda.min[[j]]))))
    #insert step here to check if sample size for correlation test meets minimum threshold
    rsq.lasso.lambda.min = c(rsq.lasso.lambda.min, cor(lasso.val.pred.lambda.min[[j]], val.true[[j]], use = "complete.obs"))
    p.val.rsq.lasso.min <- c(p.val.rsq.lasso.min, cor.test(lasso.val.pred.lambda.min[[j]], val.true[[j]], use = "complete.obs")$p.value)
    #insert step here to check if correlation was significant?
    rssm.lasso.lambda.min = sum(na.omit((val.true[[j]] - lasso.val.pred.lambda.min[[j]])^2))
    rss0.lasso.lambda.min = sum(na.omit((val.true[[j]] - null.val.pred[[j]])^2))
    lasso.pct.var.explained.lambda.min = c(lasso.pct.var.explained.lambda.min, (1 - ( rssm.lasso.lambda.min / rss0.lasso.lambda.min )))
    lasso.num.Records.lambda.min <- c(num.Records, length(lasso.val.pred.lambda.min[[j]]))    
  } else {
    p.val.rsq.lasso.min <- c(p.val.rsq.lasso.min, NA)
    num.cor.pairs.lasso.lambda.min <- c(num.cor.pairs.lasso.lambda.min, 0)
    rsq.lasso.lambda.min = c(rsq.lasso.lambda.min, NA)
    lasso.pct.var.explained.lambda.min = c(lasso.pct.var.explained.lambda.min, NA)
    lasso.num.Records.lambda.min <- c(num.Records, NA)   
  }
  #lasso (lambda.1se)
  if(length(which(!is.na(lasso.val.pred.lambda.1se[[j]])))>=3){
    num.cor.pairs.lasso.lambda.1se <- c(num.cor.pairs.lasso.lambda.1se, length(which(!is.na(lasso.val.pred.lambda.min[[j]]))))
    #insert step here to check if sample size for correlation test meets minimum threshold
    rsq.lasso.lambda.1se = c(rsq.lasso.lambda.1se, cor(lasso.val.pred.lambda.1se[[j]], val.true[[j]], use = "complete.obs"))
    p.val.rsq.lasso.1se <- c(p.val.rsq.lasso.1se, cor.test(lasso.val.pred.lambda.1se[[j]], val.true[[j]], use = "complete.obs")$p.value)
    #insert step here to check if correlation was significant?
    rssm.lasso.lambda.1se = sum(na.omit((val.true[[j]] - lasso.val.pred.lambda.1se[[j]])^2))
    rss0.lasso.lambda.1se = sum(na.omit((val.true[[j]] - null.val.pred[[j]])^2))
    lasso.pct.var.explained.lambda.1se = c(lasso.pct.var.explained.lambda.1se, (1 - ( rssm.lasso.lambda.1se / rss0.lasso.lambda.1se )))
    lasso.num.Records.lambda.1se <- c(num.Records, length(lasso.val.pred.lambda.1se[[j]]))    
  } else {
    p.val.rsq.lasso.1se <- c(p.val.rsq.lasso.1se, NA)
    num.cor.pairs.lasso.lambda.1se <- c(num.cor.pairs.lasso.lambda.1se, 0)
    rsq.lasso.lambda.1se = c(rsq.lasso.lambda.1se, NA)
    lasso.pct.var.explained.lambda.1se = c(lasso.pct.var.explained.lambda.1se, NA)
    lasso.num.Records.lambda.1se <- c(num.Records, NA)   
  }
  #rf
  if(length(which(!is.na(rf.val.pred[[j]])))>=3){
    num.cor.pairs.rf <- c(num.cor.pairs.rf, length(which(!is.na(rf.val.pred[[j]]))))
    #insert step here to check if sample size for correlation test meets minimum threshold
    rsq.rf = c(rsq.rf, cor(rf.val.pred[[j]], val.true[[j]], use = "complete.obs"))
    p.val.rsq.rf <- c(p.val.rsq.rf, cor.test(rf.val.pred[[j]], val.true[[j]], use = "complete.obs")$p.value)
    #insert step here to check if correlation was significant?
    rssm.rf = sum(na.omit((val.true[[j]] - rf.val.pred[[j]])^2))
    rss0.rf = sum(na.omit((val.true[[j]] - null.val.pred[[j]])^2))
    rf.pct.var.explained = c(rf.pct.var.explained, (1 - ( rssm.rf / rss0.rf )))
    rf.num.Records <- c(num.Records, length(rf.val.pred[[j]]))
  } else {
    p.val.rsq.rf <- c(p.val.rsq.rf, NA)
    num.cor.pairs.rf <- c(num.cor.pairs.rf, 0)
    rsq.rf = c(rsq.rf, NA)
    rf.pct.var.explained = c(rf.pct.var.explained, NA)
    rf.num.Records <- c(num.Records, NA)
  }
}
names(rsq.lasso.lambda.manual) = top.names
names(lasso.pct.var.explained.lambda.manual) = top.names
lasso.sqrt.pct.var.lambda.manual <- sqrt(lasso.pct.var.explained.lambda.manual)
#^ Error about NaNs may occur from trying to sqrt a negative number

names(rsq.lasso.lambda.min) = top.names
names(lasso.pct.var.explained.lambda.min) = top.names
lasso.sqrt.pct.var.lambda.min <- sqrt(lasso.pct.var.explained.lambda.min)
#^ Error about NaNs may occur from trying to sqrt a negative number

names(rsq.lasso.lambda.1se) = top.names
names(lasso.pct.var.explained.lambda.1se) = top.names
lasso.sqrt.pct.var.lambda.1se <- sqrt(lasso.pct.var.explained.lambda.1se)
#^ Error about NaNs may occur from trying to sqrt a negative number

names(rsq.rf) = top.names
names(rf.pct.var.explained) = top.names
rf.sqrt.pct.var <- sqrt(rf.pct.var.explained)
#^ Error about NaNs may occur from trying to sqrt a negative number

fig.2c.df <- cbind(rownames(as.data.frame(sqrt.pct.var)), 
                            as.data.frame(sqrt.pct.var), 
                            as.data.frame(lasso.sqrt.pct.var.lambda.manual),
                            as.data.frame(lasso.sqrt.pct.var.lambda.min),
                            as.data.frame(lasso.sqrt.pct.var.lambda.1se),
                            as.data.frame(rf.sqrt.pct.var), 
                            as.data.frame(num.cor.pairs.lasso.lambda.manual),
                            as.data.frame(num.cor.pairs.lasso.lambda.min),
                            as.data.frame(num.cor.pairs.lasso.lambda.1se),
                            as.data.frame(num.cor.pairs.rf),
                   as.data.frame(p.val.rsq.lasso.manual),
                   as.data.frame(p.val.rsq.lasso.min),
                   as.data.frame(p.val.rsq.lasso.1se),
                   as.data.frame(p.val.rsq.rf),
                            row.names=NULL)


colnames(fig.2c.df)<-c("test", "vitals", "lasso.manual", "lasso.min", "lasso.1se", "rf", 
                       "num.obs.lasso.manual","num.obs.lasso.min","num.obs.lasso.1se","num.obs.rf",
                       "p.val.lasso.manual","p.val.lasso.min","p.val.lasso.1se","p.val.rf")
#fig.2c.df$test = factor(fig.2c.df$test, levels = as.factor(names(sqrt.pct.var)[order(-sqrt.pct.var)]))
fig.2c.corr.coefs <- cbind(rownames(as.data.frame(rsq.vitals)), 
                           as.data.frame(rsq.vitals), 
                           as.data.frame(rsq.lasso.lambda.manual), 
                           as.data.frame(rsq.lasso.lambda.min), 
                           as.data.frame(rsq.lasso.lambda.1se), 
                           as.data.frame(rsq.rf), 
                           row.names=NULL)
colnames(fig.2c.corr.coefs)<-c("test", "vitals", "lasso.manual", "lasso.min", "lasso.1se", "rf")

#fig.2c.corr.coefs$test = factor(fig.2c.corr.coefs$test, levels = as.factor(names(rsq.vitals)[order(-rsq.vitals)]))
# fig.2c.corr.coefs[fig.2c.corr.coefs<0]=0 # clamp to zero
# ^This line was throwing an error; revised version below:
fig.2c.corr.coefs[,c("vitals","lasso.manual","lasso.min","lasso.1se","rf")] <- sapply(fig.2c.corr.coefs[,c("vitals","lasso.manual","lasso.min","lasso.1se","rf")],function(x) ifelse(x<0,0,x)) # clamp to zero
#fig.2c.corr.coefs$test = factor(fig.2c.corr.coefs$test, levels = as.factor(fig.2c.corr.coefs$test[order(-fig.2c.corr.coefs$lasso.min)]))

#choose lambda to plot ("lasso.manual", "lasso.min", or "lasso.1se")

##UNCOMMENT IF RUNNING LOCALLY
lambda.choice <- "lasso.min"

fig.2c.plot <- melt(fig.2c.corr.coefs,id.vars="test")
fig.2c.plot[,3][is.na(fig.2c.plot[,3])] <- 0 #replace % var explained of NaN w/ 0
fig.2c.plot$test = factor(fig.2c.plot$test, levels = as.factor(fig.2c.plot$test[order(-fig.2c.plot$value)]))
#^ Ran out of time, but I can simplify this later, which will probably rid the error.
fig.2c <- fig.2c.plot
#fig.2c <- fig.2c.plot[order(-fig.2c.plot[,3]),] # reorder by LM Vitals
## DONE UNCOMMENT

num.Records <- as.data.frame(num.Records)
# num.Records.2 <- transform(num.Records, TrainingObs = as.numeric(TrainingObs),
#                          TestObs = as.numeric(TestObs))
# Plot the % var explained

##UNCOMMENT IF RUNNING LOCALLY
ggplot(fig.2c[fig.2c$variable %in% c("vitals",lambda.choice,"rf"),], aes(x=test, y=value, color = variable)) + geom_point(size = 5, aes(shape=variable, color=variable)) +
  weartals_theme +
  ylim(0,1) +
  scale_shape_discrete(breaks=c("vitals", lambda.choice, "rf"),
                       labels=c("LM vitals", "LASSO", "RF")) +
  scale_color_discrete(breaks=c("vitals", lambda.choice, "rf"),
                       labels=c("LM vitals", "LASSO", "RF")) +
  labs(x = "Lab tests",y = expression(paste("Sqrt of % Variance Explained")))
## DONE UNCOMMENT

## calculate correlation coefficients and pct var explained by the models (lambda min)

# remove coefficients associated with participants that had zero x.test data
num.Records <- data.frame(num.Records)
cache <- num.Records[num.Records$TestObs=="0",]
cache <- cache[,c(2,1)] #list of clin tests and iPOPs that had zero x.test data
names(cache) <- c("test","left.out.person")
cache$delete <- 1
tmp <- merge(lasso.features.lambda.manual,cache,all=TRUE)
tmp <- tmp[which(is.na(tmp$delete)),] #remove rows selected for deletion
lasso.features.lambda.manual <- tmp[,-which(names(tmp)=="delete")] #remove column called "delete"
tmp <- merge(lasso.features.lambda.min,cache,all=TRUE)
tmp <- tmp[which(is.na(tmp$delete)),] #remove rows selected for deletion
lasso.features.lambda.min <- tmp[,-which(names(tmp)=="delete")] #remove column called "delete"
tmp <- merge(lasso.features.lambda.1se,cache,all=TRUE)
tmp <- tmp[which(is.na(tmp$delete)),] #remove rows selected for deletion
lasso.features.lambda.1se <- tmp[,-which(names(tmp)=="delete")] #remove column called "delete"

# store the results
write.csv(fig.2c.df, "../SECURE_data/20180608/20180608_pct_var_Dayprior_noDemog_ThreeLambdas.csv",row.names=FALSE)
write.csv(fig.2c.corr.coefs, "../SECURE_data/20180608/20180608_corr_coefs_Dayprior_noDemog_ThreeLambdas.csv",row.names=FALSE)
write.table(num.Records, "../SECURE_data/20180608/20180608_Dayprior_noDemog_num_Records.csv",row.names=FALSE,col.names=FALSE, sep=",")
write.table(num.Records.check, "../SECURE_data/20180608/20180608_Dayprior_noDemog_num_Records_check.csv",row.names=FALSE,col.names=FALSE, sep=",")
write.table(lasso.features.lambda.manual, "../SECURE_data/20180608/20180608_Dayprior_noDemog_LassoFeaturesLambdaManual.csv",row.names=FALSE,col.names=FALSE, sep=",")
write.table(lasso.features.lambda.1se, "../SECURE_data/20180608/20180608_Dayprior_noDemog_LassoFeaturesLambda1se.csv",row.names=FALSE,col.names=FALSE, sep=",")
write.table(lasso.features.lambda.min, "../SECURE_data/20180608/20180608_Dayprior_noDemog_LassoFeaturesLambdaMin.csv",row.names=FALSE,col.names=FALSE, sep=",")
write.table(rf.features, "../SECURE_data/20180608/20180608_Dayprior_noDemog_RF_Features.csv",row.names=FALSE,col.names=FALSE, sep=",")


## categories to color the x-axis
## UALAB data formatting problem - removed from diabetes group
clinical.groups = list()
clinical.groups[["Electrolytes"]] =c("CA","K","CL","CO2","NA.","AG")
clinical.groups[["Diabetes"]] =c("A1C","ALB","GLU","CR","ALCRU") #"UALB",
clinical.groups[["Cardiovascular.Disease"]]=c("CHOL","LDLHDL","HDL","CHOLHDL","NHDL","TGL","LDL")
clinical.groups[["Liver Function"]]=c("ALKP","BUN","ALT","TBIL","AST")
clinical.groups[["Inflammation"]]=c("BASO","LYM","LYMAB","MONO","MONOAB","NEUT","NEUTAB","IGM","EOS","EOSAB","BASOAB","WBC","HSCRP")
clinical.groups[["Blood"]] = c("PLT","GLOB","TP","HGB","HCT","RDW","MCH","MCV","RBC","MCHC")

## RF only
withDemog <- read.csv("20180621_pct_var_DayPrior_ThreeLambdas.csv",
                      header=TRUE,sep=',',stringsAsFactors=FALSE)[,c("test","rf" )]
noDemog <- read.csv("20180621_pct_var_DayPrior_noDemog_ThreeLambdas.csv",
                    header=TRUE,sep=',',stringsAsFactors=FALSE)[,c("test","rf" )]
noDemog[is.na(noDemog)] <- 0; withDemog[is.na(withDemog)] <- 0
withDemog$demog <-"withDemog"
noDemog$demog <- "noDemog"
df <- melt(rbind(withDemog, noDemog), id.vars=c("test", "demog", "shapes"))
withDemog$shapes <- rep("25", length(withDemog$demog))
noDemog$shapes <- rep("24", length(withDemog$demog))
df <- melt(rbind(withDemog, noDemog), id.vars=c("test", "demog", "shapes"))
df$value <- as.numeric(df$value)
df <- df[df$value>0,] # only show those models that actually do well.

#colors for the x-axis labels by clinical.groups
library(RColorBrewer)
myColors <- c("peachpuff2", "lightskyblue3", "firebrick", "mediumpurple4", "sandybrown", "palegreen4", "salmon")
df$col <- rep("NA", length(df$test))
for (i in 1:length(clinical.groups)){
  df$col[df$test %in% clinical.groups[[i]]] <- as.character(myColors[i])
}
df$col[df$col == "NA"] <- "black"
df$test = factor(df$test, levels = as.factor(df$test[order(-df$value)]))

ggplot(df) + 
  geom_point(aes(x=test, y=value, color=variable, shape=as.character(shapes)), size = 4) +
  weartals_theme +
  scale_color_discrete(breaks=c("vitals", "lasso.min", "rf"),
                       labels=c("LM vitals (cVS)", "LASSO (wVS)", "RF (wVS)"),
                       name="Model") +
  scale_shape_discrete(breaks=c("24", "25"),
                       labels=c("Without Demographics", "With Demographics"),
                       name=NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = df$col))

## LM, Lasso, RF
withDemog <- read.csv("20180622_pct_var_DayPrior_ThreeLambdas.csv",
                      header=TRUE,sep=',',stringsAsFactors=FALSE)[,c("test","vitals","lasso.min","rf" )]
noDemog <- read.csv("20180622_pct_var_DayPrior_noDemog_ThreeLambdas.csv",
                    header=TRUE,sep=',',stringsAsFactors=FALSE)[,c("test","vitals","lasso.min","rf" )]
noDemog[is.na(noDemog)] <- 0; withDemog[is.na(withDemog)] <- 0
withDemog$demog <-"withDemog"
noDemog$demog <- "noDemog"
withDemog$shapes <- rep("25", length(withDemog$demog))
noDemog$shapes <- rep("24", length(withDemog$demog))
df <- melt(rbind(withDemog, noDemog), id.vars=c("test", "demog", "shapes"))
df$test <- factor(df$test, levels = test_levels)
df$value <- as.numeric(df$value)
df$test = factor(df$test, levels = as.factor(df$test[order(-df$value)]))

# figure with upside down & rightside up triangles
ggplot(df) + 
  geom_point(aes(x=test, y=value, color=variable, shape=as.integer(as.character(shapes))), size = 3) +
  weartals_theme +
  scale_color_discrete(breaks=c("vitals", "lasso.min", "rf"),
                       labels=c("LM vitals (cVS)", "LASSO (wVS)", "RF (wVS)"),
                       name="Model") +
  scale_shape_identity()
                       #  breaks=as.integer(c("24", "25")),
                       #  labels=c("Without Demographics", "With Demographics"),
                       #  guide="legend"
                       # )

# figure with filled triangles and circles
ggplot(df) + 
  geom_point(aes(x=test, y=value, color=variable, shape=as.character(shapes)), size = 4) +
  weartals_theme +
  scale_color_discrete(breaks=c("vitals", "lasso.min", "rf"),
                       labels=c("LM vitals (cVS)", "LASSO (wVS)", "RF (wVS)"),
                       name="Model") +
  scale_shape_discrete(breaks=c("24", "25"),
                       labels=c("Without Demographics", "With Demographics"),
                       name=NULL) 

###
# Calculate difference between top models with and without demographics
model.diff <-max()
demog.diff <- withDemog[,2:4] - noDemog[,2:4]
demog.diff[order(demog.diff$vitals, decreasing=TRUE),]

###
# Plot the top Lasso and RF top features

lasso.features.lambda.manual <- read.csv("../SECURE_data/20180530/20180530_Dayprior_LassoFeaturesLambdaManual.csv",
                                         header=FALSE,sep=',',stringsAsFactors=FALSE)
colnames(lasso.features.lambda.manual) <- c("test", "cv.run", "left.out.person", "lasso.feature",
                                            "lasso.coef.value")
lasso.feature.summaries <- as.data.frame(summarise(group_by(lasso.features.lambda.manual, test, lasso.feature),
          mean=mean(lasso.coef.value), sd=sd(lasso.coef.value)))

# pick top 10 features from each test
top.features <- as.data.frame(lasso.feature.summaries %>% 
  group_by_(~ test) %>% 
  top_n(n = 10, wt = abs(mean)))

# figure below is probably better served with a table given the dramatically different scales for the different tests
ggplot(top.features, aes(x=lasso.feature, y=mean, color = test)) + 
  geom_point(size = 3, aes(color=test)) +
  weartals_theme +
  scale_color_discrete(breaks=unique(top.features$test),
                       labels=unique(top.features$test)) +
  labs(x = "Model Features",y = expression(paste("Coefficient")))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.5)

############
#  Fig 3D  #
############

# library("devtools)
# install_url("https://cran.r-project.org/src/contrib/Archive/impute/impute_1.26.0.tar.gz")
# install.packages("PMA)

library("PMA")
library("Hmisc")

## UALAB data formatting problem - removed from diabetes group
clinical.groups = list()
clinical.groups[["Electrolytes"]] =c("CA","K","CL","CO2","NA.","AG")
clinical.groups[["Diabetes"]] =c("A1C","ALB","GLU","CR","ALCRU") #"UALB",
clinical.groups[["Cardiovascular Disease"]]=c("CHOL","LDLHDL","HDL","CHOLHDL","NHDL","TGL","LDL")
clinical.groups[["Cardiometabolic Disease"]]=c("A1C","ALB","GLU","UALB","CR","ALCRU","CHOL"," LDLHDL","HDL","CHOLHDL","NHDL","TGL","LDL")
clinical.groups[["Liver Function"]]=c("ALKP","BUN","ALT","TBIL","AST")
clinical.groups[["Inflammation"]]=c("BASO","LYM","LYMAB","MONO","MONOAB","NEUT","NEUTAB","IGM","EOS","EOSAB","BASOAB","WBC","HSCRP")
clinical.groups[["Hematologic"]] = c("PLT","GLOB","TP","HGB","HCT","RDW","MCH","MCV","RBC","MCHC")

cca.corr.coefs <- c()
patients <- unique(wear$iPOP_ID)

wear$AG = as.numeric(wear$AG)
wear$AST = as.numeric(wear$AST)
wear$TBIL = as.numeric(wear$TBIL)
wear$HSCRP = as.numeric(wear$HSCRP)
wear$IGM = as.numeric(wear$IGM)
wear$LDL = as.numeric(wear$LDL)
wear$TGL = as.numeric(wear$TGL)

best.weigths = list(
  Hematologic = c(0.7,0.7), # ~ 0.5
  Inflammation = c(0.1,0.9), # ~ 0.37
  Electrolytes = c(0.5,0.5), # ~ 0.12
  Diabetes = c(0.1,0.7), # ~ 0.21
  Cardiovascular.Disease = c(0.1,0.7), # ~ 0.24
  'Liver Function' = c(0.3,0.7) # ~ 0.26
)

set.seed(123)
results = c()
for (nm in names(clinical.groups)){
  print(nm)
  # Remove rows with NAs
  data.clin = wear[,which(colnames(wear) %in% c("iPOP_ID", clinical.groups[[nm]]))]
  data.wear = wear[,which(colnames(wear) %in% c(wear.variables))]
  # data.wear$AgeIn2016 = wear$AgeIn2016
  # data.wear$male = wear$Gender == "M"
  # data.wear$ethnA = wear$Ethn == "A"
  # data.wear$ethnB = wear$Ethn == "B"
  # data.wear$ethnC = wear$Ethn == "C"
  
  d <- cbind(data.clin, data.wear)
  d <- na.omit(d)
  iPOP.idx <- d[,1]
  d<-d[-1]

  # remove correlated columns
  tmp <- cor(d)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  d <- d[,!apply(tmp,2,function(x) any(x > 0.99999999999))] 

  d = scale(d,center = TRUE, scale = TRUE) # center data
  indexX = c()
  indexY = c()
  
  # leave one person out CV
  for (i in 1:length(patients)){
    train <- d[!iPOP.idx %in% patients[i],,drop=FALSE]
    test <- d[iPOP.idx %in% patients[i],,drop=FALSE]
    if (nrow(test) != 1){ #maybe make this > 0?

  # scheme = c(1:3 / 4.0)
  # pensx = rep(scheme, each=3)
  # pensz = rep(scheme, times=3)
      # pensx = best.weigths[[nm]][1]
      # pensz = best.weigths[[nm]][2]

      # pensx = c(0.1,0.4,0.7,0.1,0.5,0.1,0.1,0.1)
      # pensz = c(0.1,0.4,0.7,0.9,0.1,0.3,0.5,0.9)
      pensx = c(0.7,0.1,0.7,0.5)
      pensz = c(0.7,0.7,0.1,0.5)
      
  # build the CCA model
  model.cc.cv = CCA.permute(train[,(ncol(data.clin)):(ncol(train))],  
                             train[,1:(ncol(data.clin)-1)],
                            trace = FALSE,
                            standardize = FALSE,
                            penaltyxs = pensx,
                            penaltyzs = pensz)
  model.cc = CCA(train[,(ncol(data.clin)):(ncol(train))],  
                        train[,1:(ncol(data.clin)-1)],trace = FALSE,K=1,standardize = FALSE,
                        penaltyx = model.cc.cv$bestpenaltyx,
                        penaltyz = model.cc.cv$bestpenaltyz)
      
  #plug in test data using coefficients from CCA model and compare right and left sides
  indexX = c(indexX, as.matrix(test[,(ncol(data.clin)):(ncol(test))]) %*% model.cc$u)
  indexY = c(indexY, as.matrix(test[,1:(ncol(data.clin)-1)]) %*% model.cc$v)
  
  #cca.corr <- cor(indexX, indexY)
  #print(model.cc$cor[1])
  #cca.corr.coefs <- rbind(cca.corr.coefs, c(nm, model.cc$cor[1], patients[i]))
  #cca.corr.coefs <- rbind(cca.corr.coefs, c(nm, cca.corr, patients[i]))
    }
  }
  
  cca.corr <- abs(cor(indexX, indexY))
  print(ggplot(data.frame(indexX = indexX, indexY = indexY),aes(indexX,indexY)) +
    weartals_theme +
    geom_point() +
    ggtitle(nm) +
    stat_summary(fun.data=mean_cl_normal) +
    geom_smooth(method='lm',formula=y~x))
  print(cca.corr)
  results = c(results, cca.corr)
}
df.res = data.frame(name = names(clinical.groups), cor = results)

# library(dplyr)
# data <- (cca.corr.coefs %>%
#                        group_by(nm) %>% 
#                        summarise_at(vars("cca.corr"), funs(mean,sd)))
df.res = df.res[order(-df.res$cor),]
df.res$name = factor(as.character(df.res$name), levels = as.character(df.res$name))

p=ggplot(df.res, aes(x=name, y=cor)) +
  theme(legend.title = element_blank()) +
  geom_point(size=3, shape =1) +
  weartals_theme + 
  ylim(0,0.5) +
  labs(x = NULL, y ="Correlation Coefficient")
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.5)
ggsave(paste0("plots/Figure3D.png"),p,width=5,height=4)


###################
#  Suppl. Table 3 #
###################

rf.features <-read.table("20180621_DayPrior_noDemog_RF_Features.csv",
                  header=FALSE,sep=',',stringsAsFactors=FALSE)
colnames(rf.features) <- c("test", "cv.run", "iPOP_ID", "feature", "coefficient")
rf.feature.summaries <- as.data.frame(summarise(group_by(rf.features, test, feature),
                                                   mean=mean(coefficient), sd=sd(coefficient)))
top.models <- c("HCT", "HGB", "RBC", "MONOAB", "GLU", "UALB",  "CL", "A1C")
top.rf.features <- c()
for (i in top.models) {
  rf.feature.subset<-rf.feature.summaries[rf.feature.summaries$test %in% i,] 
  rf.feature.sorted <- rf.feature.subset[order(rf.feature.subset$mean, decreasing=TRUE),]
  top.rf.features <- rbind(top.rf.features, head(rf.feature.sorted[1:10,]))
}

rf.feature.summaries[rf.feature.summaries$test %in% top.models,]


#########################################################
#     Main Text Numbers  - like Fig1D but for SEHR Data #
#########################################################

hist(table(corDf$ANON_ID), col="red", breaks=200, xlab = "Number of Clinic Visits / Person",
     main = NULL, font.lab=2,lwd=2,font=2,lty="blank",
     xlim = c(0,350)) # dist. of clinic visits in SEHR cohort
hist(table(corDf$ANON_ID)[table(corDf$ANON_ID)>50], col="red", breaks=200, xlab = "Number of Clinic Visits / Person",
     main = NULL, font.lab=2,lwd=2,font=2,lty="blank",
     xlim = c(50,350)) # dist. of clinic visits in SEHR cohort
describe(as.matrix(table(corDf$ANON_ID))) # mean & median number visits in SEHR cohort
mean(na.omit(corDf$Pulse)) # mean pulse 77.51
sd(na.omit(corDf$Pulse)) # sd pulse 14.12
mean(na.omit(corDf$Temp)) # mean temp 97.96
sd(na.omit(corDf$Temp)) # sd temp 0.50
length(na.omit(corDf$Pulse)) # 86,515
length(na.omit(corDf$Temp)) # 75,187

# duration of time monitored in SEHR dataset:
maxDate <-as.Date(as.matrix(tapply(corDf$Clin_Result_Date, corDf$ANON_ID, max)))
minDate <- as.Date(tapply(corDf$Clin_Result_Date, corDf$ANON_ID, min))
duration <- as.numeric(maxDate-minDate)
withDuration <- cbind(as.data.frame(table(corDf$ANON_ID)), duration)
describe(duration) # mean & median number of days of monitoring in SEHR cohort
hist(withDuration$duration, col="red", breaks=200, xlab = "Time Monitored by Clinic (Days)",
     main = NULL, font.lab=2,lwd=2,font=2, lty="blank")
hist(withDuration$duration[withDuration$Freq > 50], col="red", breaks=200, xlab = "Time Monitored by Clinic (Days)",
     main = NULL, font.lab=2,lwd=2,font=2, lty="blank")

#characterize the SEHR data set
length(unique(corDf$ANON_ID)) # num people in SEHR dataset where both labs and vitals exist
length(unique(labs$ANON_ID)) # num people in SEHR dataset
length(na.omit(labs$Clin_Result_Date)) # num lab tests (in the 50 labs we explored) in SEHR dataset
as.matrix(table(labs$LAB_NAME)) # number of each clinical lab
length(na.omit(vitals$Temp)) + length(na.omit(vitals$Pulse)) # total number of clinical vital signs measured
#304 people have more than 50 observations per person
length(table(corDf$ANON_ID)[table(corDf$ANON_ID)>50])


##############################
#    Suppl. Table 2 and 3    #
##############################
# create ranked list of clinical laboratory tests by the correlation coefficients between observed and predicted values; 
# predicted values from simple bivariate models of (lab test ~ pulse + temp) using SEHR dataset
# Do 10-fold cross validation at the subject level (e.g. each test set contains 1/10 of the people in the SEHR dataset)
# RUN SEHR CORRELATIONS between labs and vitals

names(corDf)[names(corDf) %in% "GLU_SerPlas"] <-"GLU"  # fix names to be same between iPOP and SEHR datasets ; number of NAs for each GLU: GLU_nonFasting (113472), GLU_wholeBld (111726), GLU_SerPlas (30949), GLU_byMeter (NA = 101012), GLU_fasting (110303)
names(corDf)[names(corDf)  %in% "LDL_Calc"] <-"LDL"  # fix names to be same between iPOP and SEHR datasets ; corDf$LDL_Calc range = wear$LDL range
options("scipen"=100, "digits"=4)
models=c(" ~ Pulse", # univariate with pulse only
         " ~ Temp",   # univariate with temp only
         " ~ Pulse + Temp", # bivariate with pulse + temp
         " ~ Pulse + I(Pulse^2)",
         " ~ Temp + I(Temp^2)", " ~ Pulse + I(Pulse^2) + Temp + I(Temp^2)" )
cv.runs <- 50
models.corr.coefs <- c()
models.pct.var <- c()

for (i in 1:cv.runs){ #50 fold cross validation (10% test set; 90% training set)
  print(i)
  ANON_ID = corDf$ANON_ID # Remember the list of subjects
  corDf.tmp = corDf[,-c(1,2)]  #remove ANON_ID and Clin_Result_Date
  corDf.tmp <- subset(corDf.tmp, select=-c(ALCRU, CR)) # all values for ALCRU tests are NA, only 20 values for CR are not NA
  nms = names(subset(corDf.tmp, select=-c(Pulse, Temp)))

  # Do cross-validation
  subjects = unique(ANON_ID)
  n = length(subjects) # total num of observations
  test = sample(n)[1:floor(n*0.1)] # 10% of subjects are held for testing
  test.subj = subjects[test]
  test.mask = ANON_ID %in% test.subj

    for (nm in top.names){ # for each of the 50 clinical lab tests
      print(nm)
      tmp=0
      corDf2 = data.frame(labtest = corDf.tmp[[nm]], Pulse = corDf.tmp$Pulse, Temp = corDf.tmp$Temp) # prepare data for LM
      #df <- cbind(corDf2[[i]], corDf2[,c("Pulse", "Temp")])
      corDf2 <- na.omit(corDf2)
      test.data <- na.omit(corDf2[test.mask,])
      train.data <-na.omit(corDf2[!test.mask,])
        for (k in 1:length(models)){
          model<-lm(as.formula(paste0("labtest",models[k])),data=train.data)
          m <- summary(model) # quadratic univariate with pulse or temp only
          model.null <- lm(as.formula(paste0("labtest"," ~ 1")),data=train.data)
          # r[tmp,tmp2]<-m$adj.r.squared # matrix of r-squared values for each left-one-out model
          # p[tmp,tmp2]<-1-pf(m$fstatistic[1],m$fstatistic[2],m$fstatistic[3]) # matrix of p-squared values for each left-one-out model
          numTrainObs<-length(train.data$Pulse) # train: the number of each clinical lab test that has corresponding vital signs
          numTestObs<-length(test.data[,1]) #  test: the number of each clinical lab test that has corresponding vital signs
          pred=predict(model, newdata=test.data)# prediction on test data set
          pred.null=predict(model.null, newdata=test.data)# prediction on test data set
          #rsq.pred = 1 - (mean( pred - test.data[,1])**2 ) / var( (test.data[,1]) ) # test r.sq
          if (length(pred)<1){next}
          r.pred = cor(pred, test.data[,1], use = "complete.obs") # test r.sq
          rssm <- sum((test.data[,1] - pred)^2)
          rss0 <- sum((test.data[,1]- pred.null)^2)
          sqrt.pct.var <- sqrt(1- (rssm/rss0))
          name.rsq <- paste("model.mean.rsq", k, sep = ".")
          models.corr.coefs <- rbind(models.corr.coefs,
                                     c(model = name.rsq, cv.step = i, test = nm, corr.coef = r.pred, sqrt.pct.var = sqrt.pct.var, numTestObs = numTestObs, numTrainObs = numTrainObs))
    }
  }
}

# models.corr.coefs <-read.csv("20180503_pct_var_30k_noDemog.csv")
corr.coefs <- as.data.frame(models.corr.coefs)
corr.coefs$cv.step <- as.numeric(as.character(corr.coefs$cv.step))
corr.coefs$corr.coef <- as.numeric(as.character(corr.coefs$corr.coef))
corr.coefs$sqrt.pct.var <- as.numeric(as.character(corr.coefs$sqrt.pct.var))

library(dplyr)
model.corr.coefs <- (corr.coefs %>%
  group_by(test, model) %>% 
  summarise_at(vars("sqrt.pct.var"), funs(mean,sd)))
model.corr.coefs$model <- mapvalues(model.corr.coefs$model, from = c("model.mean.rsq.1", "model.mean.rsq.2", "model.mean.rsq.3", "model.mean.rsq.4", "model.mean.rsq.5", "model.mean.rsq.6"), to = c("~ Pulse", "~ Temp", "~ Pulse + Temp", " ~ Pulse + P^2", " ~ Temp + T^2", " ~ Pulse + P^2 + Temp + T^2"))
model.corr.coefs <- na.omit(model.corr.coefs)
model.corr.coefs$test  = factor(model.corr.coefs$test, levels=pull(model.corr.coefs[order(-model.corr.coefs$mean),][,1]))
model.corr.coefs$mean <- pmax(model.corr.coefs$mean, 0)

# plot how rsq changes with the different models, and add in error bars from sd.plot
ggplot(model.corr.coefs, aes(x=test, y=mean, group=model, col=as.factor(model.corr.coefs$model))) +
  theme(legend.title = element_blank())+
  geom_point() +
  #guides(fill=guide_legend(title="Model")) +
  xlab("Clinical Laboratory Test") + 
  #ylab(expression(atop("Cross-Validated", paste( "Cor Coef (+/- SD)")))) +
  ylab(expression(atop("Cross-Validated", paste( "Sqrt of % Variance Explained (+/- SD)")))) +
  theme(axis.title=element_text(face="bold",size="12"),axis.text=element_text(size=12,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text.y = element_text(hjust = 1)) +
  ylim(0,0.5) +
  scale_fill_discrete(name="Model")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.5)
  #, position=position_dodge(.7))

write.table(models.corr.coefs, "20180504_pct_var_30k_noDemog.csv",row.names=FALSE,col.names=TRUE, sep=",")

#################################
#   Fig 5A  - Population Models #
#################################
library(caret)
library(plyr)
#function to pull out lm model p-values
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
names(corDf)[names(corDf) %in% "GLU_SerPlas"] <-"GLU"  # fix names to be same between iPOP and SEHR datasets ; number of NAs for each GLU: GLU_nonFasting (113472), GLU_wholeBld (111726), GLU_SerPlas (30949), GLU_byMeter (NA = 101012), GLU_fasting (110303)
names(corDf)[names(corDf)  %in% "LDL_Calc"] <-"LDL"  # fix names to be same between iPOP and SEHR datasets ; corDf$LDL_Calc range = wear$LDL range
corDf.demog <- merge(thirtyKdemog, corDf, by="ANON_ID")
corDf.demog$Gender <- as.factor(corDf.demog$Gender)
corDf.demog$Ethn <- as.factor(corDf.demog$Ethn)

# make sure you have sufficient # of tests
summaries <- summary(corDf.demog) 
to.remove <-c() 
for (i in 6:dim(summaries)[2]){
  if ( #as.numeric(unlist(strsplit(summaries[7,][i], ":"))[2]) != "NA" &
    as.numeric(unlist(strsplit(summaries[7,][i], ":"))[2]) > (dim(corDf.demog)[1] - .005*dim(corDf.demog)[1])){ #remove anything that is missing X% of our total # of observations 
    print(i)
    to.remove <- c(to.remove, names(summaries[7,][i]))
  }
}
to.remove <- gsub("\\s", "", to.remove)
corDf.demog <- corDf.demog[ , -which(names(corDf.demog) %in% c(to.remove))]

cv.runs <- 50
folds <- createFolds(factor(corDf.demog$Ethn), k = cv.runs, list = FALSE) # break data into (cv.runs) folds with equal proportion of ethnicities in each fold - if it becomes unbalanced sometimes one ethnicity will appear in training and not in test and it breaks the pipeline

#check that your folds work the way you expect
corDf.demog$cv.folds <- folds; # ddply(corDf.demog, 'cv.folds', summarise, prop=sum(Ethn=="White")) # check that this equals table(corDf.tmp$Ethn) / cv.runs

options("scipen"=100, "digits"=4)
models=c(" ~ Pulse", # univariate with pulse only
         " ~ Temp", # univariate with temp only
         " ~ Systolic", # univariate with sys only
         " ~ Diastolic",   # univariate with dias only
         " ~ Respiration", # univariate with resp only
         " ~ Pulse + Temp + Systolic + Diastolic + Respiration", # this is the total possible info we can gain from vitals
         " ~ Age + Gender + Ethn", # this is the total possible info we can gain from demog
         " ~ Pulse + Temp + Systolic + Diastolic + Respiration + Age + Gender + Ethn") # this is the total possible info we can gain from vitals and demog

         # " ~ Systolic + Diastolic + Respiration + Age + Gender + Ethn", # trivariate - this is the info we are losing by not having a wearable that measures these things
         # " ~ Pulse + I(Pulse^2) + Temp + I(Temp^2)
         # + Systolic + I(Systolic^2) + Diastolic + I(Diastolic^2) + Respiration + I(Respiration^2) + Age + Gender + Ethn" ) # this is the total possible info we can gain from vitals

# for random forest:
models=c(" ~ Pulse + Temp + Systolic + Diastolic + Respiration", # this is the total possible info we can gain from vitals
         " ~ Pulse + Temp + Systolic + Diastolic + Respiration + Age + Gender + Ethn") # this is the total possible info we can gain from vitals and demog

models.corr.coefs <- c()
cv.runs = 1
nms
for (i in 1:cv.runs){ #50 fold cross validation (10% test set; 90% training set)
  print(i)
  corDf.tmp = corDf.demog[corDf.demog$cv.folds==i,]  #remove ANON_ID and Clin_Result_Date & demographics
  # ANON_ID = corDf.tmp$ANON_ID # Remember the list of subjects
  corDf.tmp = corDf.tmp[,-c(1,5)]  #remove ANON_ID and Clin_Result_Date & demographics
  nms = names(subset(corDf.tmp, select=-c(Pulse, Temp, Systolic, Diastolic, Respiration, Age, Gender, Ethn)))
  
  # Do stratified cross-validation per subject
  # split into training and test
  #folds <- createFolds(factor(corDf.tmp$Ethn), k = 10, list = FALSE) # break data into (10% training, 90% test) folds with equal proportion of ethnicities in each fold - if it becomes unbalanced sometimes one ethnicity will appear in training and not in test and it breaks the pipeline
  # for random forest
  folds <- createFolds(factor(corDf.tmp$Ethn), k = 3000, list = FALSE) # ~50 entries per fold; folds with equal proportion of ethnicities in each fold - if it becomes unbalanced sometimes one ethnicity will appear in training and not in test and it breaks the pipeline
  corDf.tmp$test.train <- folds
  
  # subjects = unique(ANON_ID)
  # n = length(subjects) # total num of observations
  # test = sample(n)[1:floor(n*0.1)] # 10% of subjects are held for testing
  # test.subj = subjects[test]
  # test.mask = ANON_ID %in% test.subj
  
  for (nm in nms){ # for each of the 50 clinical lab tests
    print(nm)
    #tmp=0
    corDf2 = data.frame(labtest = corDf.tmp[[nm]], Pulse = corDf.tmp$Pulse, Temp = corDf.tmp$Temp,
                        Systolic = corDf.tmp$Systolic, Diastolic = corDf.tmp$Diastolic, Respiration = corDf.tmp$Respiration,
                        Age = corDf.tmp$Age, Gender = corDf.tmp$Gender, Ethn = corDf.tmp$Ethn, test.train = corDf.tmp$test.train
                         ) # prepare data for LM
    #df <- cbind(corDf2[[i]], corDf2[,c("Pulse", "Temp")])
    corDf2 <- na.omit(corDf2)
    # train.data <- corDf2[corDf2$test.train<9,] 
    # test.data <-corDf2[corDf2$test.train==10,] # training set is ~10% of total set, but not exactly because it is balancing by ethnicities
    # 
    # for random forest
    train.data <- corDf2[corDf2$test.train<180,] # 60% training, 40% test data
    test.data <-corDf2[corDf2$test.train>=180 & corDf2$test.train<300,] # training set is ~10% of total set, but not exactly because it is balancing by ethnicities
    #train.data<-sample(train.data, 450, replace=FALSE) # comment out if not doing RF
    
    t<-as.data.frame(table(train.data$Ethn)) # if there is an ethnicity that has zero entries in the training data
    test.data <- test.data[!(test.data$Ethn %in% as.character(t[t[,2]<1,][,1])),] #remove that ethnicity from the test data
    for (k in 1:length(models)){
      # model<-lm(as.formula(paste0("labtest",models[k])),data=train.data)
      # p.val <- lmp(model)
      # m <- summary(model) # quadratic univariate with pulse or temp only
      # model.null <- lm(as.formula(paste0("labtest"," ~ 1")),data=train.data)
      
      #for random forest
      model<-randomForest(as.formula(paste0("labtest",models[k])),data=train.data)
      model.null <- lm(as.formula(paste0("labtest"," ~ 1")),data=train.data)
      
      # r[tmp,tmp2]<-m$adj.r.squared # matrix of r-squared values for each left-one-out model
      # p[tmp,tmp2]<-1-pf(m$fstatistic[1],m$fstatistic[2],m$fstatistic[3]) # matrix of p-squared values for each left-one-out model
      numTrainObs<-length(train.data$Pulse) # train: the number of each clinical lab test that has corresponding vital signs
      numTestObs<-length(test.data[,1]) # test: the number of each clinical lab test that has corresponding vital signs
      pred=predict(model, newdata=test.data)# prediction on test data set
      pred.null=predict(model.null, newdata=test.data)# prediction on test data set
      #rsq.pred = 1 - (mean( pred - test.data[,1])**2 ) / var( (test.data[,1]) ) # test r.sq
      if (length(pred)<1){next}
      r.pred = cor(pred, test.data[,1], use = "complete.obs") # test r.sq
      rssm <- sum((test.data[,1] - pred)^2)
      rss0 <- sum((test.data[,1]- pred.null)^2)
      sqrt.pct.var <- sqrt(1- (rssm/rss0))
      if ((1- (rssm/rss0)) <0) {sqrt.pct.var <- 0}
      name.rsq <- paste("model.mean.rsq", k, sep = ".")
      # for lm, add in pval = p.val below
      models.corr.coefs <- rbind(models.corr.coefs,
                                 c(model = name.rsq, cv.step = i, test = nm, corr.coef = r.pred, sqrt.pct.var = sqrt.pct.var, numTestObs = numTestObs, numTrainObs = numTrainObs))
    }
  }
}

write.table(models.corr.coefs, "20180804_cVS_noID_RF_train_test_60_40.csv",row.names=FALSE,col.names=TRUE, sep=",")

corr.coefs <- as.data.frame(models.corr.coefs)
corr.coefs$cv.step <- as.numeric(as.character(corr.coefs$cv.step))
corr.coefs$corr.coef <- as.numeric(as.character(corr.coefs$corr.coef))
corr.coefs$sqrt.pct.var <- as.numeric(as.character(corr.coefs$sqrt.pct.var))
corr.coefs$numTestObs <- as.numeric(as.character(corr.coefs$numTestObs))
corr.coefs$numTrainObs <- as.numeric(as.character(corr.coefs$numTrainObs))
mean.numTrainObs <- (corr.coefs %>%
                       group_by(test, model) %>% 
                       summarise_at(vars("numTrainObs"), funs(mean)))
mean.numTestObs <- (corr.coefs %>%
                       group_by(test, model) %>% 
                       summarise_at(vars("numTestObs"), funs(mean)))
meanNumObs <- merge(mean.numTrainObs, mean.numTestObs, by = c("test", "model"))
write.table(meanNumObs, "20180804_NumObs_summarized_cVS_noID_RF_train_test_60_40.csv",row.names=FALSE,col.names=TRUE, sep=",")

library(dplyr)
## plot corr coefs or sqrt pct var explained (need to change code to plot second one)
model.corr.coefs <- (corr.coefs %>%
                       group_by(test, model) %>% 
                       summarise_at(vars("sqrt.pct.var"), funs(mean,sd))) # change to summarise_at(vars("sqrt.pct.var") or "corr.coef"
model.corr.coefs$model <- mapvalues(model.corr.coefs$model, from = c("model.mean.rsq.1", "model.mean.rsq.2"), 
                                    #from = c("model.mean.rsq.1", "model.mean.rsq.2", "model.mean.rsq.3", "model.mean.rsq.4", "model.mean.rsq.5", "model.mean.rsq.6", "model.mean.rsq.7", "model.mean.rsq.8"), 
                                    to = c("~ All Vitals", "~ All Vitals + Demographics"))
                                    # "~ Pulse","~ Temp","~ Systolic", "~ Diastolic", "~ Respiration", "~ All Vitals", "~ Demographics", "~ All Vitals + Demographics"))
                                    # to = c("~ Systolic", "~ Diastolic", "~ Respiration", "~ Systolic + Diastolic + Respiration", "~ Pulse + P^2 + Temp + T^2 + Systolic + S^2) + Diastolic + D^2 + Respiration + R^2"))
model.corr.coefs <- na.omit(model.corr.coefs)
model.corr.coefs$test  = factor(model.corr.coefs$test, levels=pull(model.corr.coefs[order(-model.corr.coefs$mean),][,1]))
model.corr.coefs$mean <- pmax(model.corr.coefs$mean, 0)
model.corr.coefs.top.names <- model.corr.coefs[model.corr.coefs$test %in% top.names,]

write.table(model.corr.coefs.top.names, "20180804_summarized_cVS_noID_RF_train_test_60_40.csv",row.names=FALSE,col.names=TRUE, sep=",")

# plot how rsq changes with the different models, and add in error bars from sd.plot
ggplot(model.corr.coefs.top.names, aes(x=test, y=mean, group=model, col=as.factor(model.corr.coefs.top.names$model))) +
  theme(legend.title = element_blank())+
  geom_point(position=position_dodge(.05)) +
  #guides(fill=guide_legend(title="Model")) +
  xlab("Clinical Laboratory Test") + 
  #ylab(expression(atop("Cross-Validated", paste( "Cor Coef (+/- SD)")))) + # if corr coefs
  ylab(expression(atop("Cross-Validated", paste( "Sqrt of % Variance Explained (+/- SD)")))) + # if sqrt pct var
  theme(axis.title=element_text(face="bold",size="12"),axis.text=element_text(size=12,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text.y = element_text(hjust = 1)) +
  ylim(0,0.7) +
  scale_fill_discrete(name="Model")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.5, position=position_dodge(.05))


# Delta Corr Coeff between Demographics Only and Vitals + Demographics:
# Which clinical lab models benefit the most by the addition of vital signs?
delta.corr.coef <- model.corr.coefs.top.names[model.corr.coefs.top.names$model %in% "~ All Vitals + Demographics",][,3] - model.corr.coefs.top.names[model.corr.coefs.top.names$model %in% "~ Demographics",][,3]
delta.corr.coef<-cbind(as.data.frame(model.corr.coefs.top.names[model.corr.coefs.top.names$model %in% "~ All Vitals + Demographics",][,1]), delta.corr.coef)
delta.corr.coef$test  = factor(delta.corr.coef$test, levels=delta.corr.coef[order(-delta.corr.coef$mean),][,1])
ggplot(delta.corr.coef, aes(x=test, y=mean))+
  geom_point() +
  theme(legend.title = element_blank())+
  xlab("Clinical Laboratory Test") + 
  #ylab(expression(atop("Delta Cor Coef ", paste( "(~ Vitals + Demog ) - (~ Demog )")))) + # if corr coefs
  ylab(expression(atop("Delta Sqrt % Var Explained ", paste( "(~ Vitals + Demog ) - (~ Demog )")))) + 
  theme(axis.title=element_text(face="bold",size="12"),axis.text=element_text(size=12,face="bold"), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.text.y = element_text(hjust = 1))
  #ylim(0,0.5)
  
## get mean RPVE for each model:

allModelResults <- read.csv("20180505_model_compare_30k_withDemog.csv",
         header =TRUE, stringsAsFactors = FALSE)
models <- unique(allModelResults$model)

# mean RVPE for the univariate models 
for (i in models){
print(i)
sub <- allModelResults[allModelResults$model %in% i,]
print(c("mean" , mean(sub$sqrt.pct.var)))
print(c("sd" , sd(sub$sqrt.pct.var)))
print(c("numTestObs", mean(sub$numTestObs)))
print(c("numTrainingObs", mean(sub$numTrainObs)))
}

# increase in RVPE from demog-only to cVS + demog for top models:
delta.corr.coef[order(delta.corr.coef$mean, decreasing = TRUE),]

###############################
# Fig 5 - Personalized Models #
###############################

library("lme4")

# Extra analysis for the text
dist = aggregate(corDf$Clin_Result_Date, list(ANON_ID = corDf$ANON_ID), length)
sum(dist$x > 50)

getTestStats = function(test){
  hct = na.omit(corDf[,c("ANON_ID",test)])
  hct.mn = aggregate(hct[[test]], list(ANON_ID = hct$ANON_ID), mean)
  hct.var = aggregate(hct[[test]], list(ANON_ID = hct$ANON_ID), var)
  hct.cnt = aggregate(hct[[test]], list(ANON_ID = hct$ANON_ID), length)

  bs = var(hct.mn$x) # between-subject
  ws = mean(hct.var$x[hct.cnt$x > 50]) # within-subject
  list(bs=bs,ws=ws,tv=var(hct[,2]),cnt=hct.cnt$x)
}
res = getTestStats("HCT")
res$bs
res$ws
res$tv
sum(res$cnt >= 50)

getDiastolicSlope = function(test){
  hct = na.omit(corDf[,c("ANON_ID","Diastolic",test)])
  hct.mn = aggregate(list(HCT = hct[[test]], diastolic=hct$Diastolic), list(ANON_ID = hct$ANON_ID), mean)

  model = lm(paste0(test," ~ Diastolic"),data=hct)
  model
}
res = getDiastolicSlope("HCT")
summary(res)

# Get ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Compute stats of an LM model of a certain test given the coefs
# The following script cross-validates by taking one observation from each patient with at least 4 observations. 

getLMresults = function(corDf4A, test.name, threshold, identifier, cap, model_coefs, threshold_hi = 1e7, mixed=FALSE, type="LM"){
  # TODO: that's an ugly way to remove NAs but correct
  corDf.tmp = corDf4A[!is.na(corDf4A[,test.name]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Temp"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Pulse"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Systolic"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Diastolic"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Respiration"]),]

  ids = sort(table(corDf.tmp[[identifier]]))
  atleastafew = names(ids[(ids >= threshold) & (ids <= threshold_hi)]) # select patients with at least 10 tests
  atleastafew = atleastafew[1:min(cap,length(atleastafew))] # troubles with training bigger models
  corDf.tmp = corDf.tmp[corDf.tmp[[identifier]] %in% atleastafew,]
  testids = c()
  
  # compute correlation if we have at least 10 patients 
  if (length(atleastafew) > 10){
    # Select the last observation from each patient who qualified (at least n0 obs)
    for (id in atleastafew){
      lastvisit = tail(which(corDf.tmp[[identifier]] == id),1)
      testids = c(testids, lastvisit )
    }
    
    # Build individual models and check correlation predicted vs true
    # model = lm(paste(test.name,"~",identifier), data=corDf.tmp[-testids,],na.action = na.omit)
    
    frm = paste(test.name,"~",model_coefs)
    print(frm)
    
    
    if (mixed)
      model = lmer(frm, data=corDf.tmp[-testids,],na.action = na.omit, verbose = 1, control = lmerControl(
        boundary.tol = 1e-2
      ))
    else if (type=="RF"){
      vars = c(test.name, attr(terms(formula(frm)),"term.labels"))
      corDf.tmp.rf = corDf.tmp[,vars]
      corDf.tmp.rf[[identifier]] = as.factor(corDf.tmp.rf[[identifier]])
      
      # one-hot encoding
      onehot = model.matrix(~ ANON_ID - 1, corDf.tmp.rf)
      colnames(onehot) = paste0("ID",1:ncol(onehot))
      corDf.tmp.rf = cbind(corDf.tmp.rf, onehot)
      corDf.tmp.rf[[identifier]] = NULL

      model = randomForest(as.formula(paste(test.name,"~ .")), data=corDf.tmp.rf[-testids,])
    }
    else
      model = lm(frm, data=corDf.tmp[-testids,],na.action = na.omit)
    
    if (type=="RF")
      preds = predict(model, newdata = corDf.tmp.rf[testids,])
    else
      preds = predict(model, newdata = corDf.tmp[testids,])
    
#    print(summary(model))
    
    # Correlation
    # res.meanpred = c(res.meanpred, cor(corDf.tmp[testids,test.name],preds))
    
    # Sqrd root of variance explained
    var.exp = sum( (corDf.tmp[testids,test.name] - preds)**2)
    var.null = sum( (corDf.tmp[testids,test.name] - mean(corDf.tmp[-testids,test.name]))**2)
    if (test.name=="LYM")
      plot(corDf.tmp[testids,test.name], preds)
    if (var.exp/var.null > 1)
      res = 0
    else 
      res = cor(preds, corDf.tmp[testids,test.name])#sqrt(1 - var.exp/var.null)
  }
  else {
    res = NA
    testids = NULL
  }
  print(res)
  list(res=res,n=length(testids))
}

# Loop over all tests and all models for the personal medels comparsion. Plot results
generate4A = function(dataset, threshold = 4, cap = 10, threshold_hi = 1e7, ntest = NULL){
  set.seed(0)
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf4A = iPOPcorDf
    test.names = allClin
  }
  else{
    identifier = "ANON_ID"
    corDf4A = corDf
    test.names = intersect(allClin,colnames(corDf4A))
  }
  
  res.values = c()
  res.n = c()
  res.models = c()
  res.tests = c()
  if (is.null(ntest))
    ntest = length(test.names)
  
  vits = c("Pulse","Temp","Systolic","Diastolic","Respiration")
  
  for (test.name in test.names[1:ntest]){
    models = list()
    
    # population vitals model
    models[["vitals"]] = paste(vits, collapse = " + ")
    
    # population vitals + personal intercept model
    models[["vitals + personal mean"]] = paste0(identifier," + ",paste(vits, collapse = " + "))
    
    # population vitals + personal intercept model
    models[["vitals + personal mean and slope"]] = paste(models[["vitals + personal mean"]]," + ",
    #    "(Pulse + Temp | ANON_ID)") 
      "(Temp + Pulse + Systolic + Diastolic + Respiration | ANON_ID)")
    #  "Temp * ANON_ID + Pulse * ANON_ID + Systolic * ANON_ID + Diastolic * ANON_ID + Respiration * ANON_ID")
    
    # population vitals + personal intercept model
    models[["personal mean"]] = paste0(identifier)
    
    res.tests  = c(res.tests, test.name)
    mdl = "vitals + personal mean"
    model = getLMresults(corDf4A, test.name, threshold, identifier, cap, models[[mdl]], threshold_hi = threshold_hi, mixed = FALSE, type = "RF")
    res.values = c(res.values, model$res)
    res.models = c(res.models, paste(mdl,"(RF)"))
    res.n = c(res.n, model$n)
    print("RF done")

    for (mdl in names(models)){
      res.tests  = c(res.tests, test.name)
      mixed = "vitals + personal mean and slope" == mdl
      model = getLMresults(corDf4A, test.name, threshold, identifier, cap, models[[mdl]], threshold_hi = threshold_hi, mixed = mixed)
      res.values = c(res.values, model$res)
      res.models = c(res.models, mdl)
      res.n = c(res.n, model$n)
    }
  }
  res = data.frame(test = res.tests, model = res.models, value = res.values, n = res.n)
  print(res)
  res = res[order(-res$value),]
  res = res[order(res$model),] # order by the population vitals model
  res$test = factor(as.character(res$test), levels = unique(as.character(res$test)))
  res = na.omit(res)
  pp = ggplot(res, aes(test, value, group = model, color = model)) +
    ylab(expression(sqrt("Variance explained"))) +
    xlab("Lab test") +
    geom_jitter(size = 2, height = 0.0, width=0.2) + 
    weartals_theme + 
    theme(text = element_text(size=14))
  ggsave(paste0("plots/Figure-4A-",dataset,".png"), plot = pp, width = 12, height = 3)
  print(pp)
  write.table(res, file=paste0("data/Figure-4A-",dataset,".csv"))
  res
}
# Increase cap in the final version to get better accuracy. Lower caps are for speeding up
#  cap - cut of number of patients for building the individual models (the higher the better population slope estimates)
#  threshold - minimum number of visits for being included in the model (the higher the more accurate personal models)
res = generate4A("30k",threshold = 50, cap = 500)

res = res[order(as.numeric(row.names(res))),]
res$model = as.character(res$model)

## Reproduce Figure 4A without rerunning the code (read from a file)
if ("from.file" %in% ls()){
res = read.table(file=paste0("plots/Figure-4A-30k.csv"))
res = res[order(-res$value),]
res = res[order(res$model),] # order by the population vitals model
res$test = factor(as.character(res$test), levels = unique(as.character(res$test)))
pp = ggplot(res, aes(test, value, group = model, color = model)) +
  ylab(expression(sqrt("Variance explained"))) +
  xlab("Lab test") +
  geom_jitter(size = 2, height = 0.0, width=0.2) + 
  weartals_theme + 
#  geom_point(size = 1) + 
  theme(text = element_text(size=14))
print(pp)
ggsave(paste0("plots/Figure-4A-30k.png"), plot = pp, width = 12, height = 3)
}

#res[1:30 * 5 - 4,]$model = paste(res[1:30 * 5 - 4,]$model,"(RF)")

# res = read.table("plots/Figure-4A-30k.csv")
res

#we select people with the largest number of observations
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Temp"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Pulse"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Systolic"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Diastolic"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Respiration"]),]
  
  toppat = -sort(-table(corDf.tmp[[identifier]]))
  toppat = toppat[toppat > 50]
  toppat = names(toppat)
  
  vitals = c("Pulse","Temp","Systolic","Diastolic","Respiration")
  dd = corDf.tmp[corDf.tmp[[identifier]] %in% toppat ,c(identifier,vitals,clin)]
  
  # Use lm to estimate the population model
  frm = paste0(clin," ~ ",paste(vitals, collapse=" + "))
  model = lm(frm,corDf.tmp)
  
  # Compute R for individual models
  dd$accuracy = dd[[identifier]]
  
  res.pat = c()
  res.err = c()
  for (pat in toppat){
    cat(paste("Building a model for",pat,"\n"))
    frm = paste0(clin," ~ ",paste(vitals, collapse=" + "))
    corDf.ind = corDf.tmp[corDf.tmp[[identifier]] == pat,]
    
    model = lm(frm, data = corDf.ind)
    res.err = c(res.err, summary(model)$adj.r)
    res.pat = c(res.pat, pat)
  }
  
  data.frame(patient = res.pat, error = res.err)
}
res = identifyNiceCase("HCT","SEHR")
res[order(-res$error),]

####################################################
## Figure 4D and 4F: Individual models Lab ~ Vital #
####################################################

generate4DF = function(clin,vit,dataset = "SEHR"){
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
  }
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Temp"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Pulse"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Systolic"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Diastolic"]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,"Respiration"]),]
  
  # Here we can select a few ANON_ID
  if (dataset == "iPOP")
    toppat = c("1636-70-1005","1636-70-1008","1636-70-1014","1636-69-001")
  else
    toppat = c("PD-4419","D-6050","N-8819","N-5362","D-3086")
  
  # Alternatively, we select people with the largest number of observations
  # toppat = names(sort(-table(corDf.tmp[[identifier]]))[1:5])
  
  dd = corDf.tmp[corDf.tmp[[identifier]] %in% toppat ,c(identifier,vit,clin)]
  
  # Use loess to estimate the population model
  vitals = c("Pulse","Temp","Systolic","Diastolic","Respiration")

  # Use lm to estimate the population model
  frm = paste0(clin," ~ ",vit)
  ww = lm(frm, corDf.tmp[sample(nrow(corDf.tmp))[1:10000],])
  grid = seq(min(corDf.tmp[[vit]], na.rm = T),max(corDf.tmp[[vit]],na.rm = T),length.out = 100)
  df = data.frame(vit = grid)
  df[[vit]] = grid
  ff = approxfun(grid, predict(ww,df ))
  
  # Compute R for individual models
  dd$accuracy = dd[[identifier]]
  
  stats = c()
  for (pat in toppat){
    cat(paste("Building a model for",pat,"\n"))
    frm = paste0(clin," ~ ",paste(vitals, collapse=" + ")) #," + ",vit,"^2")
    corDf.ind = corDf.tmp[corDf.tmp[[identifier]] == pat,]
    corDf.ind = corDf.ind[!is.na(corDf.ind[,vit]),]
    
    # Var explained (CV)
    errors = c()
    for (i in 1:20){
      testids = sample(nrow(corDf.ind))[1:floor(nrow(corDf.ind)*0.2)]
      model = lm(frm, data = corDf.ind[-testids,])
      preds = predict(model, newdata = corDf.ind[testids,])
      var.exp = sum( (corDf.ind[testids,clin] - preds) ** 2)
      var.null = sum( (corDf.ind[testids,clin] - mean(corDf.ind[-testids,clin])) ** 2)
      err = max(1 - var.exp / var.null,0)
      errors = c(errors, sqrt(err))
    }
    err = mean(errors)

    # Correlation
    # model.sum = summary(model)
    # err = sqrt(model.sum$r.squared)
    
    dd[dd[[identifier]] == pat,]$accuracy = paste0(pat," (r=",sprintf("%.1f", err),")")
    
    slope = lm(paste0(clin," ~ ",vit),corDf.ind)$coefficients[vit]
    stats = rbind(stats, c(err,slope,nrow(corDf.ind)))
  }
  rownames(stats) = toppat
  colnames(stats) = c("sqvarexp",vit,"visits")
  
  # Plot individual models with accuracies
  ggplot(dd, aes_string(vit, clin, group = "accuracy", colour = "accuracy")) + 
    weartals_theme + theme(text = element_text(size=25)) +
    #  geom_point(size=0) + 
    theme(legend.position="none") +
    geom_smooth(method="lm", formula = y ~ x, size=1, fill=NA) +
    stat_function(fun = ff, size=0.7, color="black", linetype="dashed")
  ggsave(paste0("plots/Figure-4D-",dataset,".png"),width = 8, height = 6)
  write.table(dd, file=paste0("data/Figure-4D-",dataset,".csv"))
  
  ## TRUE VS PREDICTED
  preds = predict(model, newdata = corDf.ind)

  dd = data.frame(true = corDf.ind$HCT, predicted = preds)

  pp = ggplot(dd, aes(x=true, y=predicted)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_abline(intercept = 0, size=1) + 
    geom_point(size=3, colour = gg_color_hue(5)[3])

  cor(pp$data$true, pp$data$predicted)
  ggsave(paste0("plots/Figure-4F-true-predicted-",corDf.ind$ANON_ID[1],".png"),width = 8, height = 6)
  write.table(dd, file=paste0("data/Figure-4F-",corDf.ind$ANON_ID[1],".csv"))
  print(pp)
  
  stats
}
stats = generate4DF("HCT","Diastolic","SEHR")

generate4E = function(clin,vit,dataset = "SEHR"){
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
    tocmp = c("1636-70-1005","1636-70-1008","1636-70-1014","1636-69-001")
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
    tocmp = c("N-5362","D-3086")
  }
  dd = corDf.tmp[,c(clin,vit,identifier)]

  # Use lm to estimate the population model
  frm = paste0(clin," ~ ",vit)
  ww = lm(frm, corDf.tmp[sample(nrow(corDf.tmp))[1:10000],])
  grid = seq(min(corDf.tmp[[vit]], na.rm = T),max(corDf.tmp[[vit]],na.rm = T),length.out = 100)
  df = data.frame(vit = grid)
  df[[vit]] = grid
  ff = approxfun(grid, predict(ww,df ))
  # frm = paste0(clin," ~ ",vit)
  # ww = loess(frm, corDf.tmp[sample(nrow(corDf.tmp))[1:10000],])
  # grid = seq(min(corDf.tmp[[vit]], na.rm = T),max(corDf.tmp[[vit]],na.rm = T),length.out = 100)
  # ff = approxfun(grid, predict(ww,grid))

  cols = gg_color_hue(5)
    ggplot(dd[dd[[identifier]] %in% tocmp,], aes_string(vit, clin, color = identifier)) + 
      weartals_theme + theme(text = element_text(size=25)) +
      scale_fill_manual(values = cols[c(1,3)]) +
      scale_color_manual(values = cols[c(1,3)]) +
      geom_point(size=2, aes_string(color = identifier)) + 
      geom_smooth(method="lm", formula = y ~ x, size=1) +
      theme(legend.position="none") +
      stat_function(fun = ff, size=0.7, color="black", linetype="dashed")
    ggsave(paste0("plots/Figure-4E-",dataset,".png"), width = 8, height = 6)
    write.table(dd[dd[[identifier]] %in% tocmp,],
                file=paste0("data/Figure-4E.csv"))
}
generate4E("HCT","Diastolic","SEHR")

library("grid")

generate4C = function(dataset){
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf
  }
  
  png(paste0("plots/Figure-4C-",dataset,".png"),height=300,width=1200,units="px")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, 5)))
  
  annotate.pat = c("D-3086","D-6050","N-5362","N-8819","PD-4419")
  cols = gg_color_hue(length(annotate.pat))
  
  ## Distribution of slopes in individual models
  # !! Only patients with at least min_visits = 10
  min_visits = 10
  top.vars = c("HCT")
  vits = c("Temp","Pulse","Systolic","Diastolic","Respiration")

  for (i in 1:length(top.vars)){
    for (j in 1:length(vits)){
      vit = vits[j]
      clin = top.vars[i]
      #2*(i - 1) + j + 1
      #matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      patients.tmp = sort(table(corDf.tmp[!is.na(corDf.tmp[[clin]]),][[identifier]]))
      corDf.tmp = corDf.tmp[corDf.tmp[[identifier]] %in% names(patients.tmp[patients.tmp > min_visits]),]
      if (nrow(corDf.tmp) > 2){
        frm = paste0(clin," ~ ",vit," + (",vit,"|",identifier,")")
        print(frm)
        mm = lmer(frm, data = corDf.tmp)
        cf = coef(mm)
        
        
        qq = qplot(cf[[identifier]][vit], geom="histogram")  + weartals_theme + xlab(paste0(clin," ~ ",vit)) +
          ylab("count") + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
        write.table(cf[[identifier]],
                    file=paste0("data/Figure-4C-",vit,".csv"))
        
        for (k in 1:length(annotate.pat)){
          val = cf[[identifier]][annotate.pat[k],vit]
          print(val)
          
          qq = qq + geom_vline(xintercept=val,color=cols[k],size=1.5)
        }
        print(qq, vp = viewport(layout.pos.row = i, layout.pos.col = j))
        
      }
    }
  }
  dev.off()
  
  cf
}
res = generate4C("SEHR")

###############
#  Figure 4C  #
###############
# Visits vs R^2
generate5A = function(clin,dataset = "SEHR",min_visits=10,cap=200){
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
  }
  
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Pulse"]]),]
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Temp"]]),]
  
  # Here we select people with the largest number of observations
  toppat = table(corDf.tmp[[identifier]])
  toppat = toppat[toppat > min_visits]
  toppat = names(toppat)
  toppat = toppat[1:min(cap,length(toppat))]
  
  dd = corDf.tmp[corDf.tmp[[identifier]] %in% toppat ,c(identifier,"Temp","Pulse",clin,"Clin_Result_Date")]
  
  # Compute R for individual models
  res = c()
  for (pat in toppat){
    frm = paste0(clin," ~ Pulse + Pulse^2 + Temp + Temp^2")
    
    # Crossval
    allpreds = c()
    for (i in 1:10){
      ind.data = dd[dd[[identifier]] == pat,]
      nsmpl = nrow(ind.data)
      train.idx = sample(nsmpl)[1:floor(0.8 * nsmpl)]
      model = lm(frm, data = ind.data[train.idx,])
      preds = predict(model,ind.data[-train.idx,])
      
      allpreds = rbind(allpreds, cbind(ind.data[-train.idx,clin], preds))
    }
    
    err = sqrt(1 - var(allpreds[,1] - allpreds[,2])/var(allpreds[,1] - mean(allpreds[,1])))
    vis = sum(dd[[identifier]] == pat)
    span = max(as.Date(ind.data$Clin_Result_Date)) - min(as.Date(ind.data$Clin_Result_Date))
    res = rbind(res, c(vis, span, err))
  }
  dres = data.frame(span = res[,2]/365, r = res[,3])
  ggplot(dres, aes(span, r)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_point(size=2) + 
    geom_smooth(method="lm", formula = y ~ x, size=1)
  ggsave(paste0("plots/Figure-5A-",clin,"-",dataset,".png"),width = 9,height = 6,units = "in")
}
generate5A("CHOL","SEHR",cap=100,min_visits = 10)

#########################################
#   Figure 4A & Extended Data Figure 2A #
#########################################

#Run after running individual time course from 2C (reading in various files for wear by timespans) that produced pct_var, corr_coeffs, & num_Records files   #
weartals_theme = theme_bw() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
# read in each of the corr_coeffs from the different time windows
# save as pdf 4x12.5"
fig.2c.df <-read.csv("with_restingbugfix_demog_pctdev/20180507_AllData_pct_var.csv",
                     header=TRUE,sep=',',stringsAsFactors=FALSE)
fig.2c.plot <- melt(fig.2c.df)
fig.2c.plot[,3][is.nan(fig.2c.plot[,3])] <- 0 #replace % var explained of NaN w/ 0
fig.2c$test <- fig.2c.plot[order(-fig.2c.plot[,3]),] # reorder by LM Vitals
fig.2c$test = factor(fig.2c$test, levels = order(-fig.2c.plot[,3]))

# Plot the % var explained
ggplot(fig.2c, aes(x=test, y=value, color = variable)) + geom_point(size = 5, aes(shape=variable, color=variable)) +
  weartals_theme +
  ylim(0,1) +
  scale_shape_discrete(breaks=c("vitals", "lasso", "rf"),
                       labels=c("LM vitals", "LASSO", "RF")) +
  scale_color_discrete(breaks=c("vitals", "lasso", "rf"),
                       labels=c("LM vitals", "LASSO", "RF")) +
  labs(x = "Lab tests",y = expression(paste("Sqrt of % Variance Explained"))) 

################################
#  Fig 4B and Suppl. Figure 2B #
################################
getEvents = function(dres, codes)
{
  dres$date = as.Date(dres$date)
  dres$change = c(0, dres$slope[2:length(dres$slope)] - dres$slope[2:length(dres$slope) - 1])
  dres$change = abs(dres$change)
  dres$time_diff = c(0, dres$date[2:length(dres$slope)] - dres$date[2:length(dres$slope) - 1])
  newpat = c(TRUE,dres$identifier[2:length(dres$slope)] != dres$identifier[2:length(dres$slope) - 1])
  dres = dres[!newpat,]
  dres = dres[order(-dres$change),]
  
  res.events = c()
  
  for(i in 1:10){
    events = codes[as.character(dres$identifier[i]) == codes$ANON_ID,]
    events$dist = abs(events$date - dres$date[i])
    events = events[order(events$dist),]
    res.events = rbind(res.events, events[1:10,])
  }
  res.events
}


# Temporal evolution of the estimate of the mean
generate5C = function(clin,vit,dataset = "SEHR",window=50,filter = NULL,col = 1)
  {
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
  }
  
  corDf.tmp = corDf.tmp[!is.na(corDf.tmp[[vit]]),]
  if (!is.null(filter)){
    rows = corDf.tmp[[identifier]] %in% filter
    corDf.tmp = corDf.tmp[rows,]
  }

  # Here we select people with the largest number of observations
  toppat = table(corDf.tmp[[identifier]])
  toppat = names(sort(-toppat))[1:length(filter)]

  dd = corDf.tmp[corDf.tmp[[identifier]] %in% toppat, c(identifier,"Clin_Result_Date",vit,clin)]
  
  dates = c()
  slopes = c()
  rsquared = c()
  ids = c()
  vtrue = c()
  vpred = c()
  vitals = c("Pulse","Temp","Systolic","Diastolic","Respiration")
  
  for (pat in toppat){
    dd.pat = dd[dd[[identifier]] == pat,]
    for (i in (window+1):(nrow(dd.pat)-1) ){
      model = lm(formula(paste(clin,"~",paste(vitals,sep=" + "))), data = dd.pat[(i-window):i, ])
      slopes = c(slopes, model$coefficients[1])
      rsquared = c(rsquared, sqrt(summary(model)$r.squared))
      dates = c(dates, dd.pat$Clin_Result_Date[i])
      vtrue = c(vtrue, dd.pat[i+1,clin])
      vpred = c(vpred, predict(model, newdata = dd.pat[i+1,]) )
      ids = c(ids, pat)
    }
  }
  
  dres = data.frame(date = as.Date(as.POSIXct(dates)), slope = slopes, identifier = ids, rsquared = rsquared, true = vtrue,
                    predicted = vpred)
  dres$time = as.numeric(difftime(dres$date,min(dres$date),units="days")) / 365.0
  
  dres.pred = gather(dres, values, measurement, true:predicted, factor_key=TRUE)
  
  plt_pred = ggplot(dres.pred, aes(time, measurement, group = values, color = values)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_line(size=1.3) +
    geom_point(size=2) 
  print(plt_pred)
  
  plt_slope = ggplot(dres, aes(time, slope, group = identifier, color = identifier)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    geom_line(size=1.3) +
    geom_point(size=2) 
  write.table(dres, file=paste0("data/Figure-5C-",clin,"-",vit,"-",pat,".csv"),sep = ',')
  
  cols = gg_color_hue(3)
  
  
  plt_rsq = ggplot(dres, aes(time, rsquared, group = identifier, color = identifier)) + 
    weartals_theme + theme(text = element_text(size=25)) +
    ylab(expression(sqrt("Variance explained")))
  
  if (dataset != "iPOP"){
    plt_rsq = plt_rsq + geom_point(size=4, colour=cols[col]) +
      geom_line(size=2, colour=cols[col]) 
  }
  else{
    plt_rsq = plt_rsq + geom_point(size=4) +
      geom_line(size=1.3) + ylim(c(0,0.6)) + theme(legend.position = "none")
  }
  
  events = NULL
  if (dataset != "iPOP")
    events = getEvents(dres, codes)
  
  if (dataset == "iPOP"){
    ggsave(paste0("plots/Figure-5C-",clin,'-',vit,"-",window,"-",dataset,"-slopes.png"), 
          plot = plt_slope, width = 6, height = 6,units = "in")
    ggsave(paste0("plots/Figure-5C-",clin,'-',vit,"-",window,"-",dataset,"-rsqured.png"),
           plot=plt_rsq,width = 6, height = 6,units = "in")
    ggsave(paste0("plots/Figure-5C-",clin,'-',vit,"-",window,"-",dataset,"-predictions.png"),
           plot=plt_pred,width = 6, height = 6,units = "in")
  }
  list(dres = dres, plt_slope = plt_slope, plt_rqs = plt_rsq, plt_pred = plt_pred, events = events)
}

visits = aggregate(iPOPcorDf$iPOP_ID,by=list(iPOPcorDf$iPOP_ID),length)
visits = visits[order(-visits$x),]
pats = visits[1:1,1]
dres = generate5C("HCT","Pulse","iPOP",
                  window = 30,filter=pats,1)

## Load codes (initial)
codes = read.csv("../SECURE_data/SECURE/initial_MI.csv",stringsAsFactors=FALSE,header = FALSE)
colnames(codes) = c("ANON_ID","ICD_DATE","ICD9","ICD10","DX_NAME")
codes$date = as.Date(as.POSIXct(codes$ICD_DATE,format="%d-%b-%Y")) + 1 # POSIX no time mapped by Date to the previous day so add 1

# Find Patients with enough visits
getTop10Patients = function(){
  pats = unique(codes$ANON_ID)
  corDf.tmp = corDf[corDf$ANON_ID %in% pats,]
  counts = aggregate(corDf.tmp$ANON_ID, by=list(corDf.tmp$ANON_ID), FUN=length)
  counts = counts[order(-counts$x),]
  
  npats = 10
  pats = counts[counts[,2] > 50,1][1:10] # top 5 with at least 50 visits
  pats
}
#pats = getTop10Patients()

generate5Cevents = function(pats,col,dataset="SEHR"){
  dres = generate5C("HCT","Pulse",dataset,
                    window = 30,filter=pats,col)
  
  if(dataset == "SEHR"){
    codes_pats = codes[codes$ANON_ID %in% pats,]
    ids.tmp = c("!",codes_pats$ANON_ID)
    first = c(ids.tmp[-1] != ids.tmp[-length(ids.tmp)])
    last = c(first[-1],TRUE)
    
    # correct for D-148 (last event instead of first)
    if (pats == "D-148"){
      first[] = FALSE
      first[6] = TRUE
    }
    codes_pats$time = as.numeric(codes_pats$date - min(dres$dres$date))/365
    codes_pats = codes_pats[first,]
  }
  else{
    times = c(7,43,52,64)
    first = rep(FALSE, length(dres$dres$time))
    first[times] = TRUE
    codes_pats = list(time = as.numeric(dres$dres$time[times]), ICD10 = NULL)
    print(dres$dres$date[times])
  }
  filename = paste0("plots/Figure-5C-rsqured-",pats,".png")
  
  plt_cur = dres$plt_rqs + theme(legend.position="none")
  for (evid in 1:sum(first)){
    plt_cur = plt_cur + geom_vline(xintercept = codes_pats$time[evid],color="red",size=1) +
    geom_text(aes_q(x=codes_pats$time[evid], label=paste0("\n",codes_pats$ICD10[evid]),
                    y=max(dres$dres$rsquared) - sd(dres$dres$rsquared)/2),
                    colour="black", angle=90, text=element_text())
  }
  
  ## Other events (TODO: very manual for now...)
  print(pats)
  if (pats == "D-148"){
    eid = 28
    plt_cur = plt_cur + geom_vline(xintercept = dres$dres$time[eid],color="blue",size=1)
    print(dres$dres$date[eid])
  }
  if (pats == "D-145"){
    for (eid in c(12, 17)){
      plt_cur = plt_cur + geom_vline(xintercept = dres$dres$time[eid],color="blue",size=1)
      print(dres$dres$date[eid])
    }
  }
  if (pats == "PD-6145"){
    for (eid in c(25)){
      plt_cur = plt_cur + geom_vline(xintercept = dres$dres$time[eid],color="blue",size=1)
      print(dres$dres$date[eid])
    }
  }
  
  print(plt_cur)
  ggsave(paste0(filename),
         plot=plt_cur,width = 6,height = 6,units = "in")
  filename = paste0("plots/Figure-5C-pred-",pats,".png")
  ggsave(paste0(filename),
         plot=dres$plt_pred,width = 6,height = 6,units = "in")
  
  dres
}

pats = c("D-145","D-148","PD-6145") #, "N-3683")
for (i in 1:length(pats))
  generate5Cevents(pats[i],i)

generate5D = function(pat = "1636-69-001",lab.test = "HCT"){
## 1636-69-001


nobs = 25
neval = 3
nstart = 6
shift = 2

## How much data we need for the estimates
corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[lab.test]]),]
corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Pulse"]]),]
corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["Temp"]]),]
corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["systolic"]]),]
corDf.tmp = corDf.tmp[!is.na(corDf.tmp[["diastolic"]]),]
tbl = table(corDf.tmp$iPOP_ID)
patids = names(tbl[tbl >= nobs + neval])
wall = which(corDf.tmp$iPOP_ID == patids[1])
half2 = wall[ceiling(length(wall)/2):length(wall)]
corDf.tmp[half2,]$iPOP_ID = "1636-69-001-late"
patids = c("1636-69-001-late", patids)
patids = patids[1:2]

# Here we select people with the largest number of observations
vts = c("Pulse","Temp","systolic","diastolic") #,"Respiration")

res = matrix(1, length(patids), nobs)
for (i in nstart:nobs){
  for (j in 1:length(patids)){
    pat = patids[j]
    d = corDf.tmp[corDf.tmp$iPOP_ID %in% pat, c(lab.test,vts)]
    ntest = nrow(d) - shift
    
    train = shift + (ntest - neval + 1 - i):(ntest-neval)
    test = shift + (ntest - neval + 1):ntest
    
    model = lm(paste0(lab.test," ~ ."),data=d[train,])
    
    sm = summary(model)
    pp = predict(model, newdata=d[ntest,])
    err = mean((pp - d[(ntest - neval + 1):ntest,lab.test])**2) / mean((mean(d[train,lab.test]) - d[test,lab.test])**2)
    res[j,i] = err
  }
}
mse = colMeans(res[,1:nobs])
mse[mse > 1] = 1
df = data.frame(observations = 1:nobs, RPVE = sqrt(1 - mse))
#plot(df,ylab="RPVE of HCT prediction",xlab="observation used for the model")

plt = ggplot(df, aes(observations, RPVE)) + 
  geom_point(size=3) +
  weartals_theme + theme(text = element_text(size=25)) +
  theme(legend.position="none") +
  scale_x_continuous(breaks = pretty(df$observations, n = 13)) 
print(plt)
ggsave(paste0("plots/Figure-5D-",pat,"-",lab.test,".png"),width = 6, height = 6)
write.table(df, file=paste0("data/Figure-5D-",pat,"-",lab.test,".csv"))
df
}
res = generate5D("1636-69-001", "HCT")

###############
#   Figure 5 #
###############
# run after reading in and cleaning data and running Figure 2D section to get top.names

weartals_theme = theme_bw() + theme(text = element_text(size=18), panel.border = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

corr.coefs <-read.table("../SECURE_data/20180403_ranked_models_ipop_lm.csv",row.names=1, sep=",")
top.names<-rownames(corr.coefs) # names of lab tests from either the SEHR or the iPOP simple bivariate models
top.names<-top.names[top.names %in% names(wear)] # only keep the lab names that are also present in the iPOP data

## Univariate Mixed-effect: True vs predicted 
# !! Only patients with at least min_visits = 20

min_visits = 20
mm.corr.coefs <- c()
lr.corr.coefs <- c()
id.corr.coefs <- c()
clin.idx <- c()
#for (i in 1:length(top.names)){
for (i in 1:4){
  clin = top.names[i]
  patients = sort(table(corDf[!is.na(corDf[[clin]]),]$ANON_ID))
  labs.vitals.tmp = corDf[corDf$ANON_ID %in% names(patients[patients > min_visits]),]
  labs.vitals.tmp$ANON_ID = factor(labs.vitals.tmp$ANON_ID)
  
  nn = nrow(labs.vitals.tmp)
  smp = sample(nn)
  test = smp[(1+floor(nn*0.9)):nn]
  train = smp[1:floor(nn*0.9)]
  
  frm = paste0(clin," ~ Pulse + Temp + (Pulse + Temp|ANON_ID)")
  print(frm)
  if (nrow(labs.vitals.tmp[train,]) && length(unique(labs.vitals.tmp[train,]$ANON_ID)) > 1){
    clin.idx <-c(clin.idx, clin)
    mm = lmer(frm, data = labs.vitals.tmp[train,])
    cf = coef(mm)
    vit = "Pulse"
    #qq = qplot(cf$ANON_ID[vit], geom="histogram")  + weartals_theme + xlab(paste0(top8[i]," ~ ",vit)) + ylab("count")
    #print(qq) 
    #, vp = viewport(layout.pos.row = matchidx$row,
    #                         layout.pos.col = matchidx$col))
    tt = labs.vitals.tmp[test,clin]
  
    # Evaluate LR model
    frm = paste0(clin," ~ Pulse + Temp")
    m0 = lm(frm, labs.vitals.tmp[train,])
    pp = predict(m0, newdata = labs.vitals.tmp[test,])
    #plot(pp, tt)
    lr.corr.coefs <- c(lr.corr.coefs, cor(pp,tt,use = "na.or.complete")) # corr coef of LR model
    
    # Evaluate MM model
    pp = predict(mm, newdata = labs.vitals.tmp[test,])
    #plot(pp, tt)
    mm.corr.coefs<-c(mm.corr.coefs, cor(pp,tt,use = "na.or.complete")) # corr coef of MM model
    
    # Evaluate LR model with ID
    frm = paste0(clin," ~ ANON_ID")
    m0 = lm(frm, labs.vitals.tmp[train,])
    pp = predict(m0, newdata = labs.vitals.tmp[test,])
    #plot(pp, tt)
    id.corr.coefs <- c(id.corr.coefs, cor(pp,tt,use = "na.or.complete")) # corr coef of LR model only with patient ID
  }
}
indiv.corr.coefs <- cbind(lr.corr.coefs, mm.corr.coefs, id.corr.coefs)
rownames(indiv.corr.coefs) <- clin.idx
#write.table(indiv.corr.coefs, "../SECURE_data/20180329_indiv_SEHR_corr_coeffs.csv",row.names=TRUE,col.names=TRUE, sep=",")
data <-read.table("../SECURE_data/20180330/20180329_indiv_SEHR_corr_coeffs.csv",
                  header=TRUE,sep=',',stringsAsFactors=FALSE)
d <- melt(data, id.vars="X")
ggplot(d, aes(x=X, y=value, col=variable, shape=variable))+
  geom_point(cex=2.5) + 
  weartals_theme

########################
#    Suppl. Table 4A   #
########################

allClin <- c("A1C","AG","ALB","ALKP","ALT","AST","BASO",
             "BASOAB","BUN","CA","CHOL","CHOLHDL","CL","CO2",
             "CR","EOS","EOSAB","ESR", "GLOB",
             "GLU_fasting","HCT","HDL",
             "HGB","HSCRP","IGM","K","LDL_Direct","LDLHDL","LYM","LYMAB",
             "MCH","MCHC","MCV","MONO","MONOAB","NA.","NEUT",
             "NEUTAB","NHDL","PLT","RBC","RDW","TBIL","TGL","TP","WBC")# RUN SEHR CORRELATIONS

# for each lab run a multiple regression:
r<-c()
p<-c()
fstat <-c()
degfree <- c()
tot= 0 
for (i in allClin){
  call <-paste0("corDf$",i)
  df <- cbind(corDf[[i]], corDf[,c("Pulse", "Temp")])
  df <- na.omit(df)
  #tot = tot + length(df$Pulse)
  print(c(i , length(df$Pulse))) # the number of each clinical lab test that has corresponding vital signs
  m <- summary(lm(df[,1] ~ df$Pulse + df$Temp)) # bivariate with pulse + temp
  #m <- summary(lm(df[,1] ~ df$Pulse)) # univariate with pulse or temp only
  #m <- summary(lm(df[,1] ~ df$Temp + I(df$Temp^2))) # quadratic univariate with pulse or temp only
  r[i]<-m$adj.r.squared 
  fstat[i]<-m$fstatistic
  p[i]<-1-pf(m$fstatistic[1],m$fstatistic[2],m$fstatistic[3])
  degfree[i]<-m$df
}

options("scipen"=100, "digits"=4)
str(data.frame(as.list(r)))
str(data.frame(as.list(p))); p<-sort(p) 
str(data.frame(as.list(fstat)))
tot # total number of labs that have clin vitals measures corresponding to it

# num lab tests in iPOP dataset
tot <- 0; for (i in 7:56){
  tmp <- length(as.matrix(na.omit(wear[i]))); tot <- tot + tmp}; tot 

# num vital signs in iPOP dataset
tot <- 0; for (i in 7:56){tmp <- length(as.matrix(na.omit(wear[i]))); tot <- tot + tmp}; tot

###################
# Suppl. Table 1B #
###################
# get mean +/- SE for each of the 50 clin values
thirtyk.summary <- describe(corDf)
ipop.corDf.summary <- describe(iPOPcorDf)
ipop.demog.summary <- describe(iPOPdemographics$AgeIn2016)
table(iPOPdemographics)
write.csv(ipop.corDf.summary, "supplTable2.csv")

###################
# Suppl. Table 5B #
###################
write.csv(thirtyk.summary, "supplTable1.csv")

#### END ####

