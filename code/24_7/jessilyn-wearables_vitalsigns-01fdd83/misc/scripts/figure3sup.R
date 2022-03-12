##########################################
#    Figure 3XX; Suppl. Table 2 and 3        #
##########################################
# create ranked list of clinical laboratory tests by the correlation coefficients between observed and predicted values; checked by Jessie on 2017-12-20
# predicted values from simple bivariate models of (lab test ~ pulse + temp) using 30k dataset
# Do 10-fold cross validation at the subject level (e.g. each test set contains 1/10 of the people in the 28k dataset)
# RUN 30K CORRELATIONS between labs and vitals

names(corDf)[names(corDf) %in% "GLU_SerPlas"] <-"GLU"  # fix names to be same between iPOP and 30K datasets ; number of NAs for each GLU: GLU_nonFasting (113472), GLU_wholeBld (111726), GLU_SerPlas (30949), GLU_byMeter (NA = 101012), GLU_fasting (110303)
names(corDf)[names(corDf)  %in% "LDL_Calc"] <-"LDL"  # fix names to be same between iPOP and 30K datasets ; corDf$LDL_Calc range = wear$LDL range
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

# models.corr.coefs <-read.csv("/Users/jessilyn/Desktop/framework_paper/Figure3/Fig3C/20180503_pct_var_30k_noDemog.csv")
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

write.table(models.corr.coefs, "/Users/jessilyn/Desktop/framework_paper/Figure3/Fig3C/20180504_pct_var_30k_noDemog.csv",row.names=FALSE,col.names=TRUE, sep=",")
#write.table(models.corr.coefs, "../SECURE_data/20180503_pct_var_30k_noDemog.csv",row.names=FALSE,col.names=FALSE, sep=",")

#############################################################
#    Figure 3C + Demographics + BloodPressure + Respiration #
#############################################################
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
names(corDf)[names(corDf) %in% "GLU_SerPlas"] <-"GLU"  # fix names to be same between iPOP and 30K datasets ; number of NAs for each GLU: GLU_nonFasting (113472), GLU_wholeBld (111726), GLU_SerPlas (30949), GLU_byMeter (NA = 101012), GLU_fasting (110303)
names(corDf)[names(corDf)  %in% "LDL_Calc"] <-"LDL"  # fix names to be same between iPOP and 30K datasets ; corDf$LDL_Calc range = wear$LDL range
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
folds <- createFolds(factor(corDf.demog$Ethn), k = cv.runs, list = FALSE) # break data into (cv.runs) folds with equal proportion of ethnicities in each fold - if it becomes unbalanced sometimes one ethnicity will appear in training and note in test and it breaks the pipeline

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
  #folds <- createFolds(factor(corDf.tmp$Ethn), k = 10, list = FALSE) # break data into (10% training, 90% test) folds with equal proportion of ethnicities in each fold - if it becomes unbalanced sometimes one ethnicity will appear in training and note in test and it breaks the pipeline
  # for random forest
  folds <- createFolds(factor(corDf.tmp$Ethn), k = 3000, list = FALSE) # ~50 entries per fold; folds with equal proportion of ethnicities in each fold - if it becomes unbalanced sometimes one ethnicity will appear in training and note in test and it breaks the pipeline
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

#20180506_model_compare_30k_withDemog.csv
# 20180801_30k_RF_noIDnoDemog.csv
write.table(models.corr.coefs, "/Users/jessilyn/Desktop/framework_paper/Figure3/Fig3C/20180804_cVS_noID_RF_train_test_60_40.csv",row.names=FALSE,col.names=TRUE, sep=",")

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
write.table(meanNumObs, "/Users/jessilyn/Desktop/framework_paper/Figure3/Fig3C/20180804_NumObs_summarized_cVS_noID_RF_train_test_60_40.csv",row.names=FALSE,col.names=TRUE, sep=",")

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

write.table(model.corr.coefs.top.names, "/Users/jessilyn/Desktop/framework_paper/Figure3/Fig3C/20180804_summarized_cVS_noID_RF_train_test_60_40.csv",row.names=FALSE,col.names=TRUE, sep=",")

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

allModelResults <- read.csv("/Users/jessilyn/Desktop/framework_paper/Figure3/Fig3C/20180505_model_compare_30k_withDemog.csv",
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
