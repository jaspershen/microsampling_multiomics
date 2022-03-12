##############
#  Figure 3A #
##############

hist(table(corDf$ANON_ID), col="red", breaks=200, xlab = "Number of Clinic Visits / Person",
     main = NULL, font.lab=2,lwd=2,font=2,lty="blank",
     xlim = c(0,350)) # dist. of clinic visits in 30k cohort
hist(table(corDf$ANON_ID)[table(corDf$ANON_ID)>50], col="red", breaks=200, xlab = "Number of Clinic Visits / Person",
     main = NULL, font.lab=2,lwd=2,font=2,lty="blank",
     xlim = c(50,350)) # dist. of clinic visits in 30k cohort
describe(as.matrix(table(corDf$ANON_ID))) # mean & median number visits in 30k cohort
mean(na.omit(corDf$Pulse)) # mean pulse 77.51
sd(na.omit(corDf$Pulse)) # sd pulse 14.12
mean(na.omit(corDf$Temp)) # mean temp 97.96
sd(na.omit(corDf$Temp)) # sd temp 0.50
length(na.omit(corDf$Pulse)) # 86,515
length(na.omit(corDf$Temp)) # 75,187

# duration of time monitored in 30K dataset:
maxDate <-as.Date(as.matrix(tapply(corDf$Clin_Result_Date, corDf$ANON_ID, max)))
minDate <- as.Date(tapply(corDf$Clin_Result_Date, corDf$ANON_ID, min))
duration <- as.numeric(maxDate-minDate)
withDuration <- cbind(as.data.frame(table(corDf$ANON_ID)), duration)
describe(duration) # mean & median number of days of monitoring in 30k cohort
hist(withDuration$duration, col="red", breaks=200, xlab = "Time Monitored by Clinic (Days)",
     main = NULL, font.lab=2,lwd=2,font=2, lty="blank")
hist(withDuration$duration[withDuration$Freq > 50], col="red", breaks=200, xlab = "Time Monitored by Clinic (Days)",
     main = NULL, font.lab=2,lwd=2,font=2, lty="blank")

#characterize the 30k data set
length(unique(corDf$ANON_ID)) # num people in 30k dataset where both labs and vitals exist
length(unique(labs$ANON_ID)) # num people in 30k dataset
length(na.omit(labs$Clin_Result_Date)) # num lab tests (in the 50 labs we explored) in 30k dataset
as.matrix(table(labs$LAB_NAME)) # number of each clinical lab
length(na.omit(vitals$Temp)) + length(na.omit(vitals$Pulse)) # total number of clinical vital signs measured
#304 people have more than 50 observations per person
length(table(corDf$ANON_ID)[table(corDf$ANON_ID)>50])

##############
#  Figure 3B #
##############
pdf(file = paste0(dir, "Figure3B_hists.pdf"))
par(mfrow = c(2,2))
hist(corDf$Pulse, col="tomato3", border="tomato4", breaks=50,
     xlab = "cHR", xlim=c(50,200),
     main = NULL, font.lab=2,lwd=2,font=2)
hist(corDf$Temp, col="turquoise3", border="turquoise4", breaks=5,
     xlab = "cTemp", xlim=c(65,105),
     main = NULL, font.lab=2,lwd=2,font=2)
dev.off()

##############
#  Figure 3XX #
##############
#30 K Univariate Correlation Fit Plots by Lukasz/Jessie
vitalVars <- which(names(corDf) %in% c("Pulse","Temp"))
allClin <- c("A1C","AG","ALB","ALKP","ALT","AST","BASO",
             "BASOAB","BUN","CA","CHOL","CHOLHDL","CL","CO2",
             "CR","EOS","EOSAB","ESR", "GLOB","GLU_byMeter",
             "GLU_fasting","GLU_nonFasting","GLU_SerPlas",
             "GLU_wholeBld","HCT","HDL",
             "HGB","HSCRP","IGM","K","LDL_Calc", "LDL_Direct","LDLHDL","LYM","LYMAB",
             "MCH","MCHC","MCV","MONO","MONOAB","NA.","NEUT",
             "NEUTAB","NHDL","PLT","PROCALCITONIN", "RBC","RDW","TBIL","TGL","TP","TroponinI","WBC")
clinVars <- which(names(corDf) %in% allClin)
#clin subset of the top 10 most predictive models from bivariate analysis:
# clinTopTen <- c("GLU_fasting","CR","HSCRP", "NEUTAB","NEUT","LYM", "RDW","ALB","AG", "PLT","PROCALCITONIN", "ESR")
# clinTopTen <- c("NA." , "NEUT", "HSCRP", "RBC", "LDLHDL", "ALB", "NHDL", "HGB", "GLU_fasting", "CL", "LYM")
clin.WBCs<- c("NEUT", "LYM", "BASO","MONO","EOS",
              "NEUTAB", "LYMAB", "BASOAB","MONOAB","EOSAB",) 
summary.pulse<-list()
summary.Temp<-list()
r.squared <-c()
plots <- list()
idx=0
for (j in clin.WBCs){
  idx=idx+1
  ## 30k All data scatterplots for Fig 3D and 3E
  pname <- paste0("Plot-",i)
  p<-ggplot(corDf, aes_string(y = corDf[[j]], x = corDf$Pulse)) +
    #ggplot(corDf, aes(y = corDf[[j]], x = corDf$Temp)) +
    stat_density_2d(aes(fill = ..level..), geom = 'polygon') +
    #scale_fill_viridis_c(name = "density") +
    geom_point(shape = '.') +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1.5, col="darkred") +
    theme(axis.title=element_text(face="bold",size="14"),axis.text=element_text(size=16,face="bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ylab(paste0(c(j ," Bin")))+
    xlab("cHR")
  ggsave(paste0(pname,".png"),p)
  plots[[idx]] = p   
}
p <- grid.arrange(grobs=plots,ncol=5)

clin.WBCs<- c("NEUT", "LYM", "BASO","MONO","EOS",
              "NEUTAB", "LYMAB", "BASOAB","MONOAB","EOSAB") 
#clin.WBCs<- c("HSCRP")
summary.pulse<-list()
summary.Temp<-list()
r.squared <-c()


## All points with curves
plots <- list()
idx=0
for (j in clin.WBCs){
  idx=idx+1
  ## 30k All data scatterplots for Fig 3D and 3E
  #pname <- paste0("Plot-",i)
  p<-ggplot(
    corDf, aes_string(y = corDf[[j]], x = corDf$Pulse)) +
    #corDf, aes_string(y = corDf[[j]], x = corDf$Temp)) +
    stat_density_2d(aes(fill = ..level..), geom = 'polygon') +
    #scale_fill_viridis_c(name = "density") +
    geom_point(shape = '.') +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1.5, col="darkred") +
    theme(axis.title=element_text(face="bold",size="14"),axis.text=element_text(size=16,face="bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ylab(paste0(c(j ," Bin")))+
    xlab("cHR")
  #ggsave(paste0(pname,".png"),p)
  plots[[idx]] = p   
}

p <- grid.arrange(grobs=plots,ncol=5)

## Just the binned curves
plots1 <- list()
plots2 <- list()
idx=0
for (j in clin.WBCs){
  idx=idx+1
  corDf$bin2<-ntile(corDf[[j]], 40)
  # for Temp - point plot of bin values (looks like line of points)
  # corDf2 <- summarySE(corDf, measurevar="Temp", groupvars="bin2", na.rm=TRUE)
  # p1<- ggplot(corDf2, aes(x=bin2, y=Temp)) +
  #   geom_point(stat="identity", fill="darkblue") +
  #   geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), width=.4) +
  #   xlab(paste(c(j, "bins", sep=" ")))+
  #   scale_y_continuous(limits = c(97,99))+
  #   theme(text = element_text(size=9),
  #         axis.text.x = element_text(angle = 60, hjust = 1))
  # # For Pulse
  corDf2 <- summarySE(corDf, measurevar="Pulse", groupvars="bin2", na.rm=TRUE)
  # barplot of bin values (looks like the line of points but with bars instead)
  # p1 <- ggplot(corDf2, aes(x=bin2, y=Pulse)) +
  #   geom_bar(stat="identity", fill="darkred") +
  #   geom_errorbar(aes(ymin=Pulse-se, ymax=Pulse+se), width=.2) +
  #   xlab(paste(c(j, "bins", sep=" "))) +
  #   theme(text = element_text(size=9),
  #         axis.text.x = element_text(angle = 60, hjust = 1))
  # print(paste0(j, ": number of data points in bin = ", sum(corDf$bin2 %in% "2")))
  # quadratic models of bins
  model <-lm(corDf2$Pulse  ~ corDf2$bin2 + I((corDf2$bin2)^2))
  p1 <- ggplot(corDf2, aes(y = bin2, x = Pulse)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1.5, col="darkred") +
    theme(axis.title=element_text(face="bold",size="11"), axis.text.x = element_text(angle = 60, hjust = 1), axis.text=element_text(size=11,face="bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    #geom_point(col="black") +
    ylab(paste0(c(j ," Bin")))
  # summary.pulse <- summary(lm(corDf2$Pulse ~ corDf2$bin2 + I(corDf2$bin2^2)))
  # r.squared[j] <- summary.pulse$adj.r.squared
  corDf2 <- summarySE(corDf, measurevar="Temp", groupvars="bin2", na.rm=TRUE)
  p2<-ggplot(corDf2, aes(y = bin2, x = Temp)) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1.5, col="darkblue") +
    #geom_point(col="black") +
    #ylim(c(96,98.5))+
    theme(axis.title=element_text(face="bold",size="11"), axis.text.x = element_text(angle = 60, hjust = 1), axis.text=element_text(size=11,face="bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ylab(paste0(c(j ," Bin")))
  # # summary.Temp <- summary(lm(corDf2$Temp ~ corDf2$bin2 + I(corDf2$bin2^2)))
  # # r.squared[j] <- summary.Temp$adj.r.squared
  plots1[[idx]] = p1
  plots2[[idx]] = p2
}
p <- grid.arrange(grobs=plots1,ncol=5)
p2 <- grid.arrange(grobs=plots2,ncol=5)
as.matrix(r.squared)


#####################
#  Figure 3F and 3G #
#####################

#for (i in clinTopTen){
i <- "NEUT"
pulse.diff <- c()
temp.diff <- c()
pulse.fourth.quartile <- c()
pulse.num.fourth.quartile <- c()
pulse.third.quartile <- c()
pulse.num.third.quartile <- c()
pulse.second.quartile <- c()
pulse.num.second.quartile <- c()
pulse.first.quartile <- c()
pulse.num.first.quartile <- c()
temp.fourth.quartile <- c()
temp.num.fourth.quartile <- c()
temp.third.quartile <- c()
temp.num.third.quartile <- c()
temp.second.quartile <- c()
temp.num.second.quartile <- c()
temp.first.quartile <- c()
temp.num.first.quartile <- c()
idx=0
ptm <- proc.time()
for (j in unique(corDf$ANON_ID)){
  idx=idx+1
  #create personalized quartiles for each person/measurement type; this step takes a very very long time
  person <- corDf[corDf$ANON_ID == j,]
  if (sum(!is.na(person[,i])) >= 4 & sum(!is.na(person$Pulse)) >= 4){
    print(paste0(idx, " : ", j))
    person$bins2 <- ntile(person[,i], 4)
    #get pulse values when the lab measurement for that person is in their lowest or highest quartile
    pulse.fourth.quartile[j] <- mean(person$Pulse[person$bins2 >= 4])
    pulse.num.fourth.quartile[j] <-length(person$Pulse[person$bins2 >= 4 ])
    pulse.third.quartile[j] <- mean(person$Pulse[person$bins2 >= 3  &  person$bins2 < 4 ])
    pulse.num.third.quartile[j] <- length(person$Pulse[person$bins2 >= 3  &  person$bins2 < 4 ])
    pulse.second.quartile[j] <- mean(person$Pulse[person$bins2 >= 2  &  person$bins2 < 3 ])
    pulse.num.second.quartile[j] <- length(person$Pulse[person$bins2 >= 2  &  person$bins2 < 3 ])
    pulse.first.quartile[j] <- mean(person$Pulse[person$bins2 <= 1 ])
    pulse.num.first.quartile[j] <-length(person$Pulse[person$bins2 <= 1 ])
    # make a way to save this for each i
    
    #get temp values when the lab measurement for that person is in their lowest or highest quantile
    temp.fourth.quartile[j] <- mean(person$Temp[person$bins2 >= 4 ])
    temp.num.fourth.quartile[j] <-length(person$Temp[person$bins2 >= 4 ])
    temp.third.quartile[j] <- mean(person$Temp[person$bins2 >= 3  &  person$bins2 < 4 ])
    temp.num.third.quartile[j] <- length(person$Temp[person$bins2 >= 3  &  person$bins2 < 4 ])
    temp.second.quartile[j] <- mean(person$Temp[person$bins2 >= 2  &  person$bins2 < 3 ])
    temp.num.second.quartile[j] <- length(person$Temp[person$bins2 >= 2  &  person$bins2 < 3 ])
    temp.first.quartile[j] <- mean(person$Temp[person$bins2 <= 1 ])
    temp.num.first.quartile[j] <-length(person$Temp[person$bins2 <= 1 ])
    # make a way to save this for each i
    
  }
}
proc.time() - ptm

personalQuartiles<-cbind(as.matrix(pulse.first.quartile), as.matrix(pulse.second.quartile), as.matrix(pulse.third.quartile),as.matrix(pulse.fourth.quartile),
                         as.matrix(temp.first.quartile), as.matrix(temp.second.quartile), as.matrix(temp.third.quartile), as.matrix(temp.fourth.quartile))
personalQuartiles <- personalQuartiles[(!is.na(personalQuartiles[,2]) & !is.na(personalQuartiles[,5])),] #remove NAs
#names(personalQuartiles) <- c("pulse.first.quartile", "pulse.second.quartile", "pulse.third.quartile", "pulse.fourth.quartile", "temp.first.quartile", "temp.second.quartile", "temp.third.quartile","temp.fourth.quartile")
pulse.diff<-as.matrix(personalQuartiles[,4] - personalQuartiles[,1])
hist(pulse.diff, breaks=100, col="darkred", main=paste0("Mean Pulse (1st - 4th quartile of ",i, " values)"))
boxplot(pulse.diff, col="darkred", outline=FALSE, main=paste0("Mean Pulse (1st - 4th quartile of ",i, " values)"))
temp.diff<-as.matrix(personalQuartiles[,8] - personalQuartiles[,5])
hist(temp.diff, breaks=100, col="darkblue", main=paste0("Mean Temp (1st - 4th quartile of ",i, " values)"))
boxplot(temp.diff, col="darkblue", outline=FALSE, main=paste0("Mean Temp (1st - 4th quartile of ",i, " values)"))


# data in play: temp.diff.neut and pulse.diff.neut and temp.diff.lym and pulse.diff.lym
hist(temp.diff.neut, breaks=100, main="Temperature Difference Between Personalized  4th and 1st Quartile of Neutrophil Levels", xlab="Temperature Difference", ylab="Number of Individuals", border="black", col="darkblue")
hist(pulse.diff.neut, breaks=100, main="Pulse Difference Between Personalized 4th and 1st Quartile of Neutrophil Levels", xlab="Pulse Difference", ylab="Number of Individuals", border="black", col="darkred")
hist(temp.diff.lym, breaks=100, main="Temperature Difference Between Personalized 4th and 1st Quartile of Lymphocyte Levels", xlab="Temperature Difference", ylab="Number of Individuals", border="black", col="darkblue")
hist(pulse.diff.lym, breaks=100, main="Pulse Difference Between Personalized 4th and 1st Quartile of Lymphocyte Levels", xlab="Pulse Difference", ylab="Number of Individuals", border="black", col="darkred")

write.csv(temp.diff.lym, "~/Desktop/tempdifflym.csv")
write.csv(pulse.diff.lym, "~/Desktop/pulsedifflym.csv")

length(temp.diff.neut[!is.na(temp.diff.neut)])
#  }
#}