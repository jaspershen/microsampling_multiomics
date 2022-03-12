####################
# Suppl. Figure 1  #
####################
iPOPtopTen <- c( "PLT",
                 "ALKP",
                 "BUN",
                 "MONOAB",
                 "HSCRP",
                 "GLU",
                 "GLOB",
                 "A1C",
                 "WBC",
                 "IGM" )

pList <- list(); j=0  
for (i in iPOPtopTen){
  j=j+1
  call <-paste0("iPOPcorDf$",i)
  df <- cbind(iPOPcorDf[[i]], iPOPcorDf[,c("Pulse", "Temp")])
  df <- na.omit(df)
  #pList[[j]] <- 
  print(ggplot(df, aes(x = df$Pulse, y = df[,1])) +
          geom_point(col="black", pch=19, cex=0.5) +
          stat_smooth(method = "lm", formula = y  ~ x + I(x^2), size = 1.5, col="darkred") +
          theme(axis.title=element_text(face="bold",size="14"),axis.text=element_text(size=16,face="bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          xlab("cHR") + ylab(i)) }
#grid.arrange(pList[[1]],pList[[2]],pList[[3]],pList[[4]],pList[[5]],pList[[6]], ncol=2,top="Main Title")
#####
# iPOP binning plot figures
#####
summary.pulse<-list()
summary.Temp<-list()
r.squared <-c()
for (j in clinTopTen){
  iPOPcorDf$bin2<-ntile(iPOPcorDf[[j]], 40)
  # for Temp
  # corDf2 <- summarySE(corDf, measurevar="Temp", groupvars="bin2", na.rm=TRUE)
  # print(ggplot(corDf2, aes(x=bin2, y=Temp)) +
  # geom_point(stat="identity", fill="darkblue") +
  #   geom_errorbar(aes(ymin=Temp-se, ymax=Temp+se), width=.4) +
  #   xlab(paste(c(j, "bins", sep=" ")))+
  #   scale_y_continuous(limits = c(97,99)) 
  # + theme(text = element_text(size=9),
  #         axis.text.x = element_text(angle = 60, hjust = 1)))
  # For Pulse
  iPOPcorDf2 <- summarySE(iPOPcorDf, measurevar="Pulse", groupvars="bin2", na.rm=TRUE)
  # print(ggplot(corDf2, aes(x=bin2, y=Pulse)) +
  #  geom_bar(stat="identity", fill="darkred") +
  # geom_errorbar(aes(ymin=Pulse-se, ymax=Pulse+se), width=.2) +
  # xlab(paste(c(j, "bins", sep=" ")))
  # + theme(text = element_text(size=9),
  #         axis.text.x = element_text(angle = 60, hjust = 1)))
  print(paste0(j, ": number of data points in bin = ", sum(iPOPcorDf$bin2 %in% "2")))
  print(ggplot(iPOPcorDf2, aes(x = bin2, y = Pulse)) +
          stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1.5, col="darkred") +
          theme(axis.title=element_text(face="bold",size="14"),axis.text=element_text(size=16,face="bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          geom_point(col="black") +
          xlab(paste0(c(j ," Bin"))))
  summary.pulse <- summary(lm(iPOPcorDf2$Pulse ~ iPOPcorDf2$bin2 + I(iPOPcorDf2$bin2^2)))
  r.squared[j] <- summary.pulse$adj.r.squared
  # print(ggplot(corDf2, aes(x = bin2, y = Temp)) +
  #         stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1.5, col="darkblue") +
  #         geom_point(col="black") +
  #         #ylim(c(96,98.5))+
  #         theme(axis.title=element_text(face="bold",size="14"),axis.text=element_text(size=16,face="bold"), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #         xlab(paste0(c(j ," Bin"))))
  # summary.Temp <- summary(lm(corDf2$Temp ~ corDf2$bin2 + I(corDf2$bin2^2)))
  # r.squared[j] <- summary.Temp$adj.r.squared
  
}
as.matrix(r.squared)



##################
# Suppl. Fig 2A  #
##################
plot(corDf$Pulse ~ corDf$NEUT, pch='.',
     xlab="Neutrophils", ylab="Pulse", font.lab=2,lwd=2,font=2)
abline(lm(corDf$Pulse ~ corDf$NEUT+ (corDf$NEUT)^2), col="blue",lwd=4)

#############################
#    Suppl. Table 2 and 3   #
#############################
# allClin <- c("A1C","AG","ALB","ALKP","ALT","AST","BASO",
#              "BASOAB","BUN","CA","CHOL","CHOLHDL","CL","CO2",
#              "CR","EOS","EOSAB","ESR", "GLOB","GLU_byMeter",
#              "GLU_fasting","GLU_nonFasting","GLU_SerPlas",
#              "GLU_wholeBld","HCT","HDL",
#              "HGB","HSCRP","IGM","K","LDL_Calc", "LDL_Direct","LDLHDL","LYM","LYMAB",
#              "MCH","MCHC","MCV","MONO","MONOAB","NA.","NEUT",
#              "NEUTAB","NHDL","PLT","PROCALCITONIN", "RBC","RDW","TBIL","TGL","TP","TroponinI","WBC")# RUN 30K CORRELATIONS 
allClin <- c("A1C","AG","ALB","ALKP","ALT","AST","BASO",
             "BASOAB","BUN","CA","CHOL","CHOLHDL","CL","CO2",
             "CR","EOS","EOSAB","ESR", "GLOB",
             "GLU_fasting","HCT","HDL",
             "HGB","HSCRP","IGM","K","LDL_Direct","LDLHDL","LYM","LYMAB",
             "MCH","MCHC","MCV","MONO","MONOAB","NA.","NEUT",
             "NEUTAB","NHDL","PLT","RBC","RDW","TBIL","TGL","TP","WBC")# RUN 30K CORRELATIONS

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


##############
#  Figure 2B #
##############
# make boxplots for the top best correlated clinic values with vital signs in the 30k cohort

#topHits <- c("GLU_fasting","HSCRP","ESR", "NEUT","RDW","LYM", "ALB", "PROCALCITONIN")
#topHits <- c("HSCRP","NEUTAB","NEUT","LYM", "PLT", "TroponinI", "ESR", "PROCALCITONIN")

for (i in clinTopTen){
  #quartile
  #below<-corDf$Temp[corDf[i] < summary(corDf[i])[5]]
  #above<-corDf$Temp[corDf[i] > summary(corDf[i])[5]]
  #decile
  below<- corDf$Temp[sapply(ntile(corDf[i], 40) <= 1, isTRUE)]
  #fifth<- corDf$Temp[ntile(corDf[i], 40) <= 6 & ntile(corDf[i], 40) >= 4 & !is.na(ntile(corDf[i], 40))]
  #tenth<- corDf$Temp[ntile(corDf[i], 40) <= 11 & ntile(corDf[i], 40) >= 9 & !is.na(ntile(corDf[i], 40))]
  above<- corDf$Temp[sapply(ntile(corDf[i], 40) >= 40, isTRUE)]
  print(paste0(i, ": number in upper is ", length(below), " and number of lower is ", length(above)))
  pval<-round(unlist(t.test(below, above)[3]), digits=3)
  #boxplot(below,above, outline=FALSE, main=paste(i,", P=", pval, sep=""), names=c("Below 3rd Quartile","Above 3rd Quartile"), ylab="Temp")
  #boxplot(below,fifth, tenth, above, outline=FALSE, main=paste(i,", P=", pval, sep=""), names=c("Below 1st ventile","fifth ventile", "tenth ventile","Above 10th ventile"), ylab="Temp")
  boxplot(below, above, outline=FALSE, main=paste(i,", P=", pval, sep=""), names=c("Lowest","Highest"), ylab="Temp", col="lavenderblush4")
  #color for skin temp is lavenderblush4
}

# check diabetes ones:
for (i in topHits){
  normal<- corDf$Temp[corDf$GLU_fasting < 100]
  prediabetes<- corDf$Temp[corDf$GLU_fasting > 100 & corDf$GLU_fasting < 110]
  diabetes<- corDf$Temp[ntile(corDf[i], 40) <= 11 & ntile(corDf[i], 40) >= 9 & !is.na(ntile(corDf[i], 40))]
  print(paste0(i, ": number in upper is ", length(below), " and number of lower is ", length(above)))
  pval<-round(unlist(t.test(below, above)[3]), digits=3)
  #boxplot(below,above, outline=FALSE, main=paste(i,", P=", pval, sep=""), names=c("Below 3rd Quartile","Above 3rd Quartile"), ylab="Temp")
  #boxplot(below,fifth, tenth, above, outline=FALSE, main=paste(i,", P=", pval, sep=""), names=c("Below 1st ventile","fifth ventile", "tenth ventile","Above 10th ventile"), ylab="Temp")
  boxplot(below, above, outline=FALSE, main=paste(i,", P=", pval, sep=""), names=c("Below 1st ventile","Above 10th ventile"), ylab="Temp")
}



############
# Figure 4 #
############
library(lme4)
topNeut <- names(head(sort(pulse.num.top.quartile, decreasing = TRUE), 20)) #the 20 people that had the most measurements of NEUT (between 67 and 171 measurements)
p1 <- ggplot(data = corDf[corDf$ANON_ID %in% topNeut,], aes(x = NEUT, y = Temp, colour = ANON_ID)) +       
  #geom_point() + 
  #geom_smooth(method='lm',formula=y~x, se = FALSE)
  geom_smooth(method='loess',formula=y~x, se = FALSE)
#  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 2), se = FALSE)
model <- lmList(NEUT ~ Temp | ANON_ID, data=corDf)

na.omit(corDf) %>% 
  group_by(ANON_ID) %>% 
  do({
    mod = lm(NEUT ~ Temp, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })



