source("load-data.R")
df <- fread(paste0(dir, "/BasisData_20161111_PostSummerAddOns_Cleaned_NotNormalized_20180427.csv"),
            header=TRUE,sep=",",stringsAsFactors = FALSE,
            select=c("Timestamp_Local","Heart_Rate","Skin_Temperature_F",
                     "Steps","iPOP_ID")) # 38009228 observations
df<-df[!is.na(df$iPOP_ID),] # only keep wearables data that has an associated iPOP ID; 
df[,"Heart_Rate"] <- apply(df[,"Heart_Rate"], 2, remove_outliers) # clean data based on HR (TODO: later also clean on Skin Temp, Steps)
df <- df[which(!is.na(df[,"Heart_Rate"])),] # remove observations where HR = NA
# restrict to daytime values only (between 6am-10pm)
df$Timestamp_Local<-fastPOSIXct(df$Timestamp_Local) # takes a very long time

# keep only resting times: can use the same time as the clinic visits or the resting definition from the PlosBio paper for consistency with the wVS dataframes used throughout the paper
daytime.df <- with( df, df[ hour( Timestamp_Local ) >= 7 & hour( Timestamp_Local ) < 9 , ] ) # pull data only from specific time window; store hourly resting data for boxplots

vitals <- iPOPvitals
vitals$Clin_Result_Date<-as.Date(vitals$Clin_Result_Date, "%Y-%m-%d")
colnames(vitals)[1:5] <- c("iPOP_ID", "Date", "BP", "Pulse", "Temp")

windows=c(60, 10)
#, 30, 60) # define time windows with no steps for resting threshold (10,60,120, etc)
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
fig1b.list <- list()
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
  
  groupColumns = c("iPOP_ID", "Date")
  dataColumns = c("restingHR", "restingSkinTemp","restingSteps")
  dailyRestingMeans = ddply(restingDf.all, groupColumns, function(x) colMeans(x[dataColumns])) # calculate daily means
  
  # if wearables data is within X days of the clinic date, keep it.
  X = c(30, 14, 7)
  wRHR<-c()
  wTemp<-c()
  mean.cHR<-c()
  mean.wRHR<-c()
  mean.cTemp<-c()
  mean.wTemp<-c()
  sd.cHR<-c()
  sd.wRHR<-c()
  num.Obs<-c()
  index=0
  for (i in X){
    for (d in 1:dim(vitals)[1]){ # mean wRHR for the clinic date is the mean wRHR for X days prior to the clinic visit
      wRHR[d]<- mean(dailyRestingMeans[dailyRestingMeans$iPOP_ID %in% vitals$iPOP_ID[d] & dailyRestingMeans$Date < vitals$Date[d] & dailyRestingMeans$Date > vitals$Date[d]-i,]$restingHR)
      wTemp[d]<- mean(dailyRestingMeans[dailyRestingMeans$iPOP_ID %in% vitals$iPOP_ID[d] & dailyRestingMeans$Date < vitals$Date[d] & dailyRestingMeans$Date > vitals$Date[d]-i,]$restingSkinTemp)
    }
    index=index+1
    wRHR.v.cHR <- data.frame(vitals$iPOP_ID, vitals$Date, vitals$Pulse, vitals$Temp, wRHR, wTemp) #use this for the new histograms
    wRHR.v.cHR.full<-wRHR.v.cHR[complete.cases(wRHR.v.cHR), ]
    mean.cHR[index] <- mean(wRHR.v.cHR.full$vitals.Pulse) #MEAN cHR for PAPER
    mean.wRHR[index] <- mean(wRHR.v.cHR.full$wRHR)       #MEAN wRHR for PAPER
    sd.cHR[index] <- sd(wRHR.v.cHR.full$vitals.Pulse)  #SD cHR for PAPER
    sd.wRHR[index] <- sd(wRHR.v.cHR.full$wRHR)         #SD wRHR for PAPER
    num.Obs[index]<-dim(wRHR.v.cHR.full)[1]
    
  }
  fig1b <-data.frame(rep(window,length(X)), X, mean.cHR, mean.wRHR, sd.cHR, sd.wRHR, num.Obs)
  fig1b.list[[idx]]<-fig1b
}
fig1b.df <- do.call("rbind",fig1b.list)

par(mfrow = c(2,2), mai = c(0.7, 0.7, 0.7, 0.7))
options(scipen=10)
hist(wRHR.v.cHR.full$vitals.Pulse, col="tomato3", border="tomato4", breaks=50,
     xlab = "cHR", xlim=c(50,120),
     main = NULL, font.lab=2,lwd=2,font=2)
scale_y_continuous()
hist(wRHR.v.cHR.full$vitals.Temp, col="turquoise4", border="turquoise3", breaks=46,
     xlab = "cTemp", xlim=c(65,105),
     main = NULL, font.lab=2,lwd=2,font=2)
hist(wRHR.v.cHR.full$wRHR, col="tomato3", border="tomato4", breaks=50,
     xlab = "wRHR", xlim=c(50,120),
     main = NULL, font.lab=2,lwd=2,font=2)
scale_y_continuous()
hist(wRHR.v.cHR.full$wTemp, col="turquoise4", border="turquoise3", breaks=46,
     xlab = "wRTemp", xlim=c(65,105),
     main = NULL, font.lab=2,lwd=2,font=2)
dev.off()
