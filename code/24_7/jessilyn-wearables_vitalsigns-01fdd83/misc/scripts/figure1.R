####################
#### Figure 1  #####
####################
iPOPdaysMonitored <- read.csv("/Users/jessilyn/Desktop/framework_paper/Figure1/Slide 2/slide2_C_participant_data_summary.csv",
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

###############
#  Figure 1B  #
############### 
# expanded analysis in 20180605_Fig1D_analysis.R in the whitecoat github directory
df <- fread(paste0(dir, "/BasisData_20161111_PostSummerAddOns_Cleaned_NotNormalized_20180427.csv"),
            header=TRUE,sep=",",stringsAsFactors = FALSE,
            select=c("Timestamp_Local","Heart_Rate","Skin_Temperature_F",
                     "Steps","iPOP_ID")) # 38009228 observations
df<-df[!is.na(df$iPOP_ID),] # only keep wearables data that has an associated iPOP ID; 
df[,"Heart_Rate"] <- apply(df[,"Heart_Rate"], 2, remove_outliers) # clean data based on HR (TODO: later also clean on Skin Temp, Steps)
df <- df[which(!is.na(df[,"Heart_Rate"])),] # remove observations where HR = NA
# restrict to daytime values only (between 6am-10pm)
#df$Timestamp_Local<-as.POSIXct(df$Timestamp_Local) # takes forever
df$Timestamp_Local<-fastPOSIXct(df$Timestamp_Local) # takes a very long time
#daytime.df <- with( df, df[ hour( Timestamp_Local ) >= 6 & hour( Timestamp_Local ) < 22 , ] ) # pull data only from specific time window; store hourly resting data for boxplots
#daytime.df <- with( df, df[ hour( Timestamp_Local ) >= 20 & hour( Timestamp_Local ) < 24 , ] ) # pull data only from specific time window; store hourly resting data for boxplots
#daytime.df <- with( df, df[ hour( Timestamp_Local ) >= 24 | hour( Timestamp_Local ) >= 0 & hour( Timestamp_Local ) < 1, ] ) # pull data only from specific time window; store hourly resting data for boxplots

#below is actually nighttime for the test mike asked for
daytime.df <- with( df, df[ hour( Timestamp_Local ) >= 19 & hour( Timestamp_Local ) < 21 , ] ) # pull data only from specific time window; store hourly resting data for boxplots

vitals <- iPOPvitals
vitals$Clin_Result_Date<-as.Date(vitals$Clin_Result_Date, "%Y-%m-%d")
colnames(vitals)[1:5] <- c("iPOP_ID", "Date", "BP", "Pulse", "Temp")
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


ggplot()+
  geom_line(aes(x=windows, y=wRHR.mean), color="red") +
  geom_point(aes(x=windows, y=wRHR.mean), color="red") +
  # geom_line(aes(x=windows, y=wRTemp.mean), color="blue") +
  # geom_point(aes(x=windows, y=wRTemp.mean), color="blue") +
  xlab("Resting Time Window (seconds)") +
  ylab("Mean wRHR") # or wRTemp

cols <- c("red","blue","green","purple")
ggplot()+
  geom_line(aes(x=windows, y=wRHR.sd), color="red") +
  geom_point(aes(x=windows, y=wRHR.sd), color="red") +
  geom_line(aes(x=windows, y=HR.personal.sd), color="blue") +
  geom_point(aes(x=windows, y=HR.personal.sd), color="blue") +
  geom_line(aes(x=windows, y=rep(cHR.sd, length(windows))), color="coral") +
  geom_point(aes(x=windows, y=rep(cHR.sd, length(windows))), color="coral") +
  geom_line(aes(x=windows, y=rep(cHR.individual.sd, length(windows))), color="skyblue") +
  geom_point(aes(x=windows, y=rep(cHR.individual.sd, length(windows))), color="skyblue") +
  xlab("Resting Time Window (seconds)") +
  ylab("HR SD") # or wRTemp

wRHR.mean
wRHR.sd
wRHR.num.obs
HR.personal.sd
wRTemp.mean
wRTemp.sd
wRTemp.num.obs
Temp.personal.sd

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
  
  # Scatter plot HR
  p1 <- ggplot(rhr.daily.means.id,
               aes(x=restingHR,y=Pulse,col=as.factor(rhr.daily.means.id$idx))) +
    #col=as.factor(substr(iPOP_ID,9,12)))) +
    geom_point() +
    #labs(title=paste0("Clinical Pulse vs. Wearable RHR Mean (",window,"min resting)"),x="Resting Heart Rate",y="Pulse") +
    labs(title=NULL,x="wRHR",y="cHR") +
    annotate("segment",x=-Inf,xend=Inf,y=-Inf,yend=Inf,
             lwd=1, color="blue", alpha=.25) +
    #guides(col=guide_legend("Subject ID")) +
    xlim(40, 100) +
    ylim(40, 100)+
    theme(plot.title=element_text(face="bold",colour="black",size=16),
          axis.title.x=element_text(face="bold",colour="black",size=16),
          axis.text.x=element_text(face="bold",colour="black",size=16,angle=55,vjust=0.9,hjust=1),
          axis.title.y=element_text(face="bold",colour="black",size=16),
          axis.text.y=element_text(face="bold",colour="black",size=16),
          axis.ticks.length = unit(.2,"cm"),
          #legend.position="none",
          #legend.title=element_text(face="bold", colour="black", size=16),
          #legend.text=element_text(face="bold", colour="black", size=16),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #panel.background=element_rect(fill="grey94"))
  
  ## scatter plot skin temp
  p2 <- ggplot(rhr.daily.means.id,
               aes(x=restingSkinTemp,y=Temp,col=as.factor(rhr.daily.means.id$idx))) +
    #col=as.factor(substr(iPOP_ID,9,12)))) +
    geom_point() +
    #labs(title=paste0("Clinical Temp vs. Wearable Temp (",window,"min resting)"), x="Resting Skin Temp",y="Core Temperature") +
    labs(title=NULL, x="wRTemp",y="cTemp") +
    annotate("segment",x=-Inf,xend=Inf,y=-Inf,yend=Inf,
             lwd=1, color="blue", alpha=.25) +
    guides(col=guide_legend("Subject ID")) +
    xlim(92, 100) +
    ylim(92, 100)+
    theme(plot.title=element_text(face="bold",colour="black",size=16),
          axis.title.x=element_text(face="bold",colour="black",size=16),
          axis.text.x=element_text(face="bold",colour="black",size=16,angle=55,vjust=0.9,hjust=1),
          axis.title.y=element_text(face="bold",colour="black",size=16),
          axis.text.y=element_text(face="bold",colour="black",size=16),
          axis.ticks.length = unit(.2,"cm"),
          #legend.position="none",
          legend.title=element_text(face="bold", colour="black", size=16),
          legend.text=element_text(face="bold", colour="black", size=16),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #panel.background=element_rect(fill="grey94"))
  
  grid.arrange(p1, p2, nrow = 2)
  # save as a 7x4.5 pdf
}

##########
# Fig 1C #
##########
hist(iPOPdaysMonitored$Days_monitored_by_clinic, col="grey", breaks=10,
     xlab = "Time Monitored by Clinic (Days)", main = NULL, font.lab=2,lwd=2,font=2)
hist(iPOPdaysMonitored$Days_monitored_by_clinic, col="red", breaks=10,
     xlab = "Time Monitored by Clinic (Days)", main = NULL, font.lab=2,lwd=2,font=2, lty="blank")
hist(iPOPdaysMonitored$Total_NumOfClinMeasures, col="grey", breaks=10,
     xlab = "Number of Clinic Visits / Person", main = NULL, font.lab=2,lwd=2,font=2)
hist(iPOPdaysMonitored$Days_monitored_by_basis, col="grey", breaks=20,
     xlab = "Time Monitored by Watch (Days)", main = NULL, font.lab=2,lwd=2,font=2)
mean(iPOPdaysMonitored$Total_NumOfClinMeasures)
mean(iPOPdaysMonitored$Days_monitored_by_clinic)

###############
# Fig 1D  top #
###############

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

pdf(file = paste0(dir, "../Figure1/Figure1D_hists.pdf"))
par(mfrow = c(2,2), mai = c(0.7, 0.7, 0.7, 0.7))
hist(iPOPvitals$Pulse, col="tomato3", , border="tomato4", breaks=50,
     xlab = "cHR", xlim=c(50,200),
     main = NULL, font.lab=2,lwd=2,font=2)
hist(iPOPvitals$Temp, col="turquoise3", border="turquoise4", breaks=10,
     xlab = "cTemp", xlim=c(65,105),
     main = NULL, font.lab=2,lwd=2,font=2)
#####################
#  Figure 1D Bottom #
#####################

# resting HR and ST from Fig 1B data - need to run code for Fig 1B first
options(scipen=10)
hist(restingDf.all$restingHR, col="tomato3", border="tomato4", breaks=50,
     xlab = "wRHR", xlim=c(50,200),
     main = NULL, font.lab=2,lwd=2,font=2)
scale_y_continuous()
hist(restingDf.all$restingSkinTemp, col="turquoise3", border="turquoise4", breaks=46,
     xlab = "wRTemp", xlim=c(65,105),
     main = NULL, font.lab=2,lwd=2,font=2)
dev.off()

dfFigOne <- fread(paste0(paste0(dir, "BasisData_20161111_PostSummerAddOns_Cleaned_NotNormalized_20180427.csv")),
                  header=TRUE,sep=",",stringsAsFactors = FALSE)

dfFigOne$Heart_Rate <- remove_outliers(dfFigOne$Heart_Rate) # clean data based on HR (TODO: later also clean on Skin Temp, Steps)
length(dfFigOne$Heart_Rate[!is.na(dfFigOne$Heart_Rate)]) # number of wHR measurements in iPOP cohort
mean(dfFigOne$Heart_Rate[!is.na(dfFigOne$Heart_Rate)]) # mean wHR mean 74.31 +/- 15.17, n=25,341,508
sqrt(var(dfFigOne$Heart_Rate[!is.na(dfFigOne$Heart_Rate)])) # stdev of wHR

# hist(dfFigOne$Heart_Rate, col="darkred", breaks=100,
#      xlab = "wHR",
#      main = NULL, font.lab=2,lwd=2,font=2,
#      xlim=c(0,200))

dfFigOne$Skin_Temperature_F <- remove_outliers(dfFigOne$Skin_Temperature_F) # clean data based on HR (TODO: later also clean on Skin Temp, Steps)
length(dfFigOne$Skin_Temperature_F[!is.na(dfFigOne$Skin_Temperature_F)]) # number of wTemp measurements in iPOP cohort
mean(dfFigOne$Skin_Temperature_F[!is.na(dfFigOne$Skin_Temperature_F)]) # mean wTemp mean 88.57 +/- 3.74, n=27,136,802
sqrt(var(dfFigOne$Skin_Temperature_F[!is.na(dfFigOne$Skin_Temperature_F)])) # stdev of wTemp 
# hist(dfFigOne$Skin_Temperature_F, col="darkgrey", breaks=100,
#      xlab = "wTemp", xlim=c(65,105),
#      main = NULL, font.lab=2,lwd=2,font=2)

#characterize the iPOP data set
length(na.omit(iPOPvitals$Temp)) + length(na.omit(iPOPvitals$Pulse)) # total number of clinical vital signs measured
describe(iPOPlabs[names(iPOPlabs) %in% allClin]) # summary of clinical labs data
length(unique(wear$iPOP_ID)) # num people in iPOP wearables dataset
