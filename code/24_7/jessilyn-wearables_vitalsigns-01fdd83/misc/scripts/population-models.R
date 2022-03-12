###############
#  Figure 2D  #
###############

# creates ranked list of clinical laboratory tests by the %var explained in simple LM; LOO cross validation at the subject level 
# use.Troubleshoot.mode - choose for during troubleshooting
# use.Demog - choose whether Demographics in models (supply TRUE or FALSE)
# use.iPOP - choose whether iPOP_ID variable is used (supply TRUE or FALSE)
cVS.lm = function(data, use.Troubleshoot.mode = TRUE, use.Demog = FALSE, use.iPOP = FALSE, npatients = NULL, model = "lm"){
  if (use.Troubleshoot.mode){
    top.names <- c("HGB", "TGL")
  }
  
  ####
  # CODE FOR SIMPLE LM
  #
  sum.vectors.in.list <- function(x) {
    c(sum(na.omit(x)))
  }
  
  rsq.all = c()
  pct.var.all = c()
  iPOPcorDf.demo <- merge(data, iPOPdemographics[1:4], by="iPOP_ID")
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
    vitals.variables <- c("Pulse", "systolic", "diastolic", "Temp", "AgeIn2016", "Gender", "Ethn") # "BMI", "systolic", "diastolic", 
  } else if(!use.Demog) {
    vitals.variables <- c("Pulse", "Temp", "systolic", "diastolic") #
  }
  if (use.iPOP){
    vitals.variables = c(vitals.variables,"iPOP_ID")
  }
  
  patients = unique(data$iPOP_ID)
  
  res = list()
  res$val.true <- rep(list(NA),length(top.names)) # list of vectors to store true values; each vector is for 1 clinical lab
  res$val.pred <- rep(list(NA),length(top.names)) # list of vectors to store trained model predicted values; each vector is for 1 clinical lab
  res$val.null.pred <- rep(list(NA),length(top.names)) # list of vectors to store null model predicted values; each vector is for 1 clinical lab
  # p.value<-list()  # TODO: decide if this is a relevant parameter to collect - I think not...?
  res$num.true <- rep(list(NA),length(top.names)) # number of observations per individual per clinic test
  
  if (is.null(npatients))
    npatients = length(patients)
  
  #patients.test = c("1636-69-001","1636-70-1005","1636-70-1008","1636-70-1014")
  patients.test = patients[1:npatients]
  
  for (pat in patients.test[1:npatients]){
    cat("Patient",pat,"\n") # LOO
    for (l in 1:length(top.names)){
      
      train <- patients[patients != pat]
      test <- patients[patients == pat]
      train.ids = iPOPcorDf.demo$iPOP_ID %in% train
      test.ids = iPOPcorDf.demo$iPOP_ID %in% test
      
      if (use.iPOP){
        test.ids = iPOPcorDf.demo$iPOP_ID %in% test
        test.idx = which(test.ids)
        
        nfrac = length(test.idx)
        if (nfrac < 6)
          next
        nfrac = floor(length(test.idx)*0.2)
        
        test.idx = sample(test.idx)[1:nfrac]
        test.ids = (1:length(test.ids) %in% test.idx)
        train.ids = !test.ids
      }
      
      iPOPcorDf.demo.tmp = iPOPcorDf.demo
      if (use.iPOP){
        iPOPcorDf.demo.tmp$iPOP_ID = 1*(iPOPcorDf.demo.tmp$iPOP_ID == pat)
      }
      
      dat.train.unsorted = iPOPcorDf.demo.tmp[train.ids, ] # subset input data by training set
      dat.train = dat.train.unsorted[order(dat.train.unsorted$iPOP_ID),] #order by iPOP_ID in order to supply correct nfolds arg to glmnet
      dat.test = iPOPcorDf.demo.tmp[test.ids, ] # subset input data by testing set
      
      cat("Test",top.names[l],"\n")
      # create training set
      x.train <- dat.train[,colnames(dat.train) %in% c("iPOP_ID",top.names[l], vitals.variables)] # subset input data by lab: only take current lab test of interest
      x.train <- na.omit(x.train) # skip nas and nans ## TODO: the way this script is written, you will lose a lot of data because you take the number of lab visits down to the test with the minimum number of visits. However, if you do na.omit after the next line, you have to change your matrix to accept dynamic number of row entries. Not sure how to do this yet, so for now just reducing the data amount by a lot. 
      x.train.ids <- x.train$iPOP_ID
      
      if (!use.iPOP){
        x.train <- x.train[,-1]
      }
      
      predictors <- as.data.frame(x.train[,colnames(x.train) %in% c(vitals.variables)]) # later add in demographics
      outcome <- as.matrix(x.train[,colnames(x.train) %in% top.names[l]]) # matrix of outcome for model building # tried adding as.numeric after as.matrix() but that introduced new issues
      # create test set
      x.test<-dat.test[,colnames(dat.test) %in% c(top.names[l], vitals.variables)] # subset input data by lab: only take current lab test of interest
      x.test<- na.omit(x.test) # skip nas and nans ## TODO: SEE ABOVE na.omit FOR ISSUE WITH THIS
      res$val.true[[l]] = c(res$val.true[[l]], x.test[,top.names[l]]) # true values of left out person
      res$num.true[[l]]<-c(res$num.true[[l]],length(x.test[,top.names[l]])) # number of test observations for left out person
      fml = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(vitals.variables,collapse=" + "))
      fml.null = paste(top.names[l]," ~ 1")
      
      if (model=="lm"){
        bivar.lm.model = lm(as.formula(fml), data = x.train) # build the model
        res$val.pred[[l]] = c(res$val.pred[[l]], predict(bivar.lm.model, newdata = x.test)) # predict on trained model
      }
      if (model=="rf"){
        bivar.lm.model = lm(as.formula(fml), data = x.train) # build the model
        res$val.pred[[l]] = c(res$val.pred[[l]], predict(bivar.lm.model, newdata = x.test)) # predict on trained model
      }
      
      bivar.null.lm.model<-lm(as.formula(fml.null), data = x.train) # create null model for hypothesis testing and for calculating RSS0
      res$val.null.pred[[l]] = c(res$val.null.pred[[l]], predict(bivar.null.lm.model, newdata = x.test)) # predict on null model
      # t<- anova(bivar.null.lm.model, bivar.lm.model) # to get p-values for model
      # p.value[[l]] <- as.numeric(t[2,][["Pr(>F)"]])  # to get p-values for model
    }
  }
  num.test.obs <- lapply(res$num.true, sum.vectors.in.list)
  
  res$rsq.vitals = c()
  res$p.val.rsq.vitals = c()
  res$rssm.vitals = c()
  res$rss0.vitals = c()
  res$pct.var.explained = c()
  res$num.Records.check <- c()
  for (j in 1:length(top.names)){
    res$rsq.vitals = c(res$rsq.vitals, cor(res$val.pred[[j]], res$val.true[[j]], use = "complete.obs"))
    res$p.val.rsq.vitals = c(res$p.val.rsq.vitals, cor.test(res$val.pred[[j]], res$val.true[[j]], use = "complete.obs")$p.value)
    res$rssm.vitals = sum(na.omit((res$val.true[[j]] - res$val.pred[[j]])^2))
    res$rss0.vitals = sum(na.omit((res$val.true[[j]] - res$val.null.pred[[j]])^2))
    res$pct.var.explained = c(res$pct.var.explained, (1 - ( res$rssm.vitals / res$rss0.vitals )))
    res$num.Records.check <- c(res$num.Records.check, (length(res$val.pred[[j]])-1)) # same as num.test.obs
  }
  names(res$rsq.vitals) = top.names
  names(res$pct.var.explained) = top.names
  
  res$pct.var.explained[res$pct.var.explained<0]=0
  res$sqrt.pct.var <- sqrt(res$pct.var.explained)
  res
}

####
# CODE FOR LASSO, RF
####

wVS.models = function(data, use.Troubleshoot.mode = FALSE, use.Demog = FALSE, use.iPOP = FALSE, npatients = NULL){
  if (use.Troubleshoot.mode){
    top.names <- c("HGB", "TGL")
  }
  #clean wear data frame
  data[,8:length(names(data))] <- apply(
    data[,8:length(names(data))], 2,
    function(x) as.numeric(as.character(x)))
  data$Gender <- as.factor(data$Gender)
  data$Ethn <- as.factor(data$Ethn)
  wear.variables <- unlist(read.table(paste0(dir,"FinalLasso_153WearableFactors.csv"), stringsAsFactors = FALSE)) # the table of model features we want to work with
  wear.variables = wear.variables[1:20]
  #wear.variables = c("hr_mean","st_mean","gsr_mean","rhr_mean")
  
  if (use.iPOP){
    wear.variables = c(wear.variables,"iPOP_ID")
  }
  
  #change gender and ethnicity to dummy variables
  gender <- data.frame(model.matrix( ~ Gender - 1, data=data))
  ethn <- data.frame(model.matrix( ~ Ethn - 1, data=data))
  
  #remove the least populated gender and ethnicity (NCOL-1)
  cache <- names(gender)[which(sapply(gender,sum)==max(sapply(gender,sum)))]
  gender <- data.frame(cache=gender[which(sapply(gender,sum)==max(sapply(gender,sum)))])
  ethn <- ethn[,-which(sapply(ethn,sum)==min(sapply(ethn,sum)))]
  
  #store names as demo.variables
  demo.variables <- c("AgeIn2016", names(gender), names(ethn))
  
  #cbind new gender and ethnicity variables to "data"
  data <- cbind(data,gender,ethn)
  
  res = list()
  res$val.true <- rep(list(NA),length(top.names)) # list of vectors to store true values; each vector is for 1 clinical lab
  res$null.val.pred <- rep(list(NA),length(top.names))  # list of vectors to store nullmodel-predicted values; each vector is for 1 clinical lab
  res$lasso.val.pred.lambda.manual <- rep(list(NA),length(top.names)) # list of vectors to store lasso-trainedmodel-predicted values; each vector is for 1 clinical lab
  res$lasso.val.pred.lambda.min <- rep(list(NA),length(top.names)) # list of vectors to store lasso-trainedmodel-predicted values; each vector is for 1 clinical lab
  res$lasso.val.pred.lambda.1se <- rep(list(NA),length(top.names)) # list of vectors to store lasso-trainedmodel-predicted values; each vector is for 1 clinical lab
  res$rf.val.pred <- rep(list(NA),length(top.names))  # list of vectors to store rf-trainedmodel-predicted values; each vector is for 1 clinical lab
  
  res$num.Records <- c()
  res$lasso.features.lambda.manual <- data.frame("test"=character(),"cv.run"=character(),
                                                 "left.out.person"=character(),"lasso.feature"=character(),
                                                 "lasso.coef.value"=character())
  res$lasso.features.lambda.min <- data.frame("test"=character(),"cv.run"=character(),
                                              "left.out.person"=character(),"lasso.feature"=character(),
                                              "lasso.coef.value"=character())
  res$lasso.features.lambda.1se <- data.frame("test"=character(),"cv.run"=character(),
                                              "left.out.person"=character(),"lasso.feature"=character(),
                                              "lasso.coef.value"=character())
  res$rf.features <- data.frame("test"=character(),"cv.run"=character(),
                                "left.out.person"=character(),"rf.feature"=character(),
                                "rf.coef.value"=character())
  
  if(use.Demog){
    demo.variables <- c("AgeIn2016", names(gender), names(ethn))
  } else if(!use.Demog) {
    demo.variables <- c()
  }
  
  patients = unique(iPOPcorDf$iPOP_ID)
  if (is.null(npatients))
    npatients = length(patients)
  
  #patients.test = c("1636-69-001") #,"1636-70-1008")
  #patients.test = c("1636-69-001","1636-70-1005","1636-70-1008","1636-70-1014")
  patients.test = patients[1:npatients]
  
  k = 0
  for (pat in patients.test){
    k=k+1
    cat("Patient",pat,"\n") # LOO
    
    for (l in 1:length(top.names)){
      wear.tmp = data
      if (use.iPOP){
        wear.tmp = filter.nas(wear.tmp, c(top.names[l], wear.variables, demo.variables))
        patients = unique(wear.tmp$iPOP_ID)
      }
      
      train <- patients[patients != pat]
      test <- patients[patients == pat]
      train.ids = wear.tmp$iPOP_ID %in% train
      test.ids = wear.tmp$iPOP_ID %in% test
      
      if (use.iPOP){
        test.ids = iPOPcorDf.demo$iPOP_ID %in% test
        test.idx = which(test.ids)
        
        nfrac = length(test.idx)
        if (nfrac < 6)
          next
        nfrac = floor(length(test.idx)*0.2)
        
        test.idx = sample(test.idx)[1:nfrac]
        test.ids = (1:length(test.ids) %in% test.idx)
        train.ids = !test.ids
        
        wear.tmp$iPOP_ID = 1*(wear.tmp$iPOP_ID == pat)
      }
      
      
      
      dat.train.unsorted <- wear.tmp[ train.ids, ] # subset input data by training set
      dat.train <- dat.train.unsorted[order(dat.train.unsorted$iPOP_ID),] #order by iPOP_ID in order to supply correct nfolds arg to glmnet
      dat.test<-wear.tmp[ test.ids, ] # subset input data by testing set
      
      cat("Test",top.names[l],"\n")
      x.train<-dat.train[,colnames(dat.train) %in% c("iPOP_ID", top.names[l], wear.variables, demo.variables)] # subset input data by lab: only take current lab test of interest
      x.train<-na.omit(x.train) # skip nas and nans ## TODO: the way this script is written, you will lose a lot of data because you take the number of lab visits down to the test with the minimum number of visits. However, if you do na.omit after the next line, you have to change your matrix to accept dynamic number of row entries. Not sure how to do this yet, so for now just reducing the data amount by a lot. 
      x.train.ids<-x.train$iPOP_ID
      if (!use.iPOP){
        x.train <- x.train[,-1,drop=FALSE]
      }
      
      # if(!NROW(x.train)){ #if x.train is empty
      #   print(paste0("The x.train data was empty for ",pat,"'s ",top.names[l]," test."))
      # } else {
      #   print(paste0("The x.train data for ",pat,"'s ",top.names[l]," test had ",NROW(x.train)," observations."))
      # }
      predictors <- as.data.frame(x.train[,colnames(x.train) %in% c(wear.variables, demo.variables)]) # later add in demographics
      outcome <- as.matrix(x.train[,colnames(x.train) %in% top.names[l]]) # matrix of outcome for model building # tried adding as.numeric after as.matrix() but that introduced new issues
      
      # create test set
      x.test<-dat.test[,colnames(dat.test) %in% c(top.names[l], wear.variables, demo.variables)] # subset input data by lab: only take current lab test of interest
      x.test<- na.omit(x.test) # skip nas and nans ## TODO: SEE ABOVE na.omit FOR ISSUE WITH THIS
      # if(!NROW(x.test)){ #if x.test is empty
      #   print(paste0("The x.test data was empty for ",pat,"'s ",top.names[l]," test."))
      # } else {
      #   print(paste0("The x.test data for ",pat,"'s ",top.names[l]," test had ",NROW(x.test)," observations."))
      # }
      res$val.true[[l]] = c(res$val.true[[l]], x.test[,top.names[l]]) # true values of left out person
      
      res$num.Records <- rbind(res$num.Records, c(IPOP_ID=pat, test=top.names[l], TrainingObs=length(outcome), TestObs=length(x.test[,top.names[l]])))
      
      rf.variables.to.use = c(wear.variables, demo.variables) # rf variables (use all)
      
      #decide on number for nfolds from number of obs per subject
      frq <- as.vector(table(x.train.ids))
      
      #optional argument for leave-one-out CV method
      n <- max(3,length(frq))
      
      #optional argument to specify folds for CV
      folds <- rep(1:length(frq),frq[1:length(frq)])
      
      ## run lasso for variable selection
      # n <- as.numeric(length(outcome)) #optional argument for leave-one-out CV method for nfold
      x_train <- model.matrix( ~ .-1, as.data.frame(predictors))
      glm.res = cv.glmnet(x=x_train,y=outcome,
                          standardize.response=FALSE,
                          family="gaussian",
                          nfolds=n,
                          #                        foldid=folds,
                          nlambda=100)
      
      #    lasso.model.lambda.min = lm(as.formula(lasso.fml.lambda.min), data = x.train) # , weights = labs.wear$weight) # TODO: do we need to include weights?    
      
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
                        "left.out.person"=rep(pat,length(lasso.variables.lambda.manual)),
                        "lasso.feature"=names(lasso.variables.lambda.manual),
                        "lasso.coef.value"=as.numeric(lasso.variables.lambda.manual))
      res$lasso.features.lambda.manual <- rbind(res$lasso.features.lambda.manual,tmp)
      
      tmp <- data.frame("test"=rep(top.names[l],length(lasso.variables.lambda.min)),
                        "cv.run"=rep(k,length(lasso.variables.lambda.min)),
                        "left.out.person"=rep(pat,length(lasso.variables.lambda.min)),
                        "lasso.feature"=names(lasso.variables.lambda.min),
                        "lasso.coef.value"=as.numeric(lasso.variables.lambda.min))
      res$lasso.features.lambda.min <- rbind(res$lasso.features.lambda.min,tmp)
      
      tmp <- data.frame("test"=rep(top.names[l],length(lasso.variables.lambda.1se)),
                        "cv.run"=rep(k,length(lasso.variables.lambda.1se)),
                        "left.out.person"=rep(pat,length(lasso.variables.lambda.1se)),
                        "lasso.feature"=names(lasso.variables.lambda.1se),
                        "lasso.coef.value"=as.numeric(lasso.variables.lambda.1se))
      res$lasso.features.lambda.1se <- rbind(res$lasso.features.lambda.1se,tmp)
      
      #store lasso variable names based on coef threshold (lambda specific: manual, min, and 1se)
      lasso.variables.to.use.lambda.manual = names(factors.lambda.manual[abs(factors.lambda.manual)>1e-10]) # TODO: this is an arbitrary rule for now
      lasso.variables.to.use.lambda.min = names(factors.lambda.min[abs(factors.lambda.min)>1e-10])
      lasso.variables.to.use.lambda.1se = names(factors.lambda.1se[abs(factors.lambda.1se)>1e-10])
      
      # # Remove iPop
      # if (use.iPOP){
      #   lasso.variables.to.use.lambda.manual = c("iPOP_ID",lasso.variables.to.use.lambda.manual[!grepl("iPOP", lasso.variables.to.use.lambda.manual)])
      #   lasso.variables.to.use.lambda.min = c("iPOP_ID",lasso.variables.to.use.lambda.min[!grepl("iPOP", lasso.variables.to.use.lambda.min)])
      #   lasso.variables.to.use.lambda.1se = c("iPOP_ID",lasso.variables.to.use.lambda.1se[!grepl("iPOP", lasso.variables.to.use.lambda.1se)])
      # }
      if (use.iPOP){
        lasso.variables.to.use.lambda.manual = c("iPOP_ID",lasso.variables.to.use.lambda.manual)
        lasso.variables.to.use.lambda.min = c("iPOP_ID",lasso.variables.to.use.lambda.min)
        lasso.variables.to.use.lambda.1se = c("iPOP_ID",lasso.variables.to.use.lambda.1se)
      }
      
      # build null, lasso, and rf models
      set.seed(1)
      null.fml = paste(top.names[l]," ~ 1")
      null.model = lm(as.formula(null.fml), data = x.train) # create null model for hypothesis testing and for calculating RSS0
      res$null.val.pred[[l]] = c(res$null.val.pred[[l]], predict(null.model, newdata = x.test)) # predict on null model
      
      lasso.fml.lambda.manual = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(lasso.variables.to.use.lambda.manual,collapse=" + "))
      # lasso.model.lambda.manual = lm(as.formula(lasso.fml.lambda.manual), data = x.train) # , weights = labs.wear$weight) # TODO: do we need to include weights?
      # res$lasso.val.pred.lambda.manual[[l]] = c(res$lasso.val.pred.lambda.manual[[l]], predict(lasso.model.lambda.manual, newdata = x.test)) # predict on trained model
      # check that the formula is valid (i.e. not empty)
      if(lasso.fml.lambda.manual!=paste0("cbind( ",top.names[l]," ) ~ ")){
        lasso.model.lambda.manual = lm(as.formula(lasso.fml.lambda.manual), data = x.train) # , weights = labs.wear$weight) # TODO: do we need to include weights?
        res$lasso.val.pred.lambda.manual[[l]] = c(res$lasso.val.pred.lambda.manual[[l]], predict(lasso.model.lambda.manual, newdata = x.test)) # predict on trained model
      } else {
        lasso.model.lambda.manual = null.model # if lasso sets all coeffs = 0, then use the null model
        res$lasso.val.pred.lambda.manual[[l]] = c(res$lasso.val.pred.lambda.manual[[l]], predict(lasso.model.lambda.manual, newdata = x.test)) # predict on trained model
        # cache <-  length(res$val.true[[l]])-length(res$lasso.val.pred.lambda.manual[[l]])
        # res$lasso.val.pred.lambda.manual[[l]] = c(res$lasso.val.pred.lambda.manual[[l]], rep(NA,cache)) # fill with NA(s) if invalid model was supplied
      }
      
      lasso.fml.lambda.min = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(lasso.variables.to.use.lambda.min,collapse=" + "))
      #check that the formula is valid (i.e. not empty)
      if(lasso.fml.lambda.min!=paste0("cbind( ",top.names[l]," ) ~ ")){
        lasso.model.lambda.min = lm(as.formula(lasso.fml.lambda.min), data = x.train) # , weights = labs.wear$weight) # TODO: do we need to include weights?
        res$lasso.val.pred.lambda.min[[l]] = c(res$lasso.val.pred.lambda.min[[l]], predict(lasso.model.lambda.min, newdata = x.test)) # predict on trained model
      } else {
        lasso.model.lambda.min = null.model # if lasso sets all coeffs = 0, then use the null model
        res$lasso.val.pred.lambda.min[[l]] = c(res$lasso.val.pred.lambda.min[[l]], predict(lasso.model.lambda.min, newdata = x.test)) # predict on trained model
        # cache <-  length(res$val.true[[l]])-length(res$lasso.val.pred.lambda.min[[l]])
        # res$lasso.val.pred.lambda.min[[l]] = c(res$lasso.val.pred.lambda.min[[l]], rep(NA,cache)) # fill with NA(s) if invalid model was supplied
      }
      
      lasso.fml.lambda.1se = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(lasso.variables.to.use.lambda.1se,collapse=" + "))
      #check that the formula is valid (i.e. not empty)
      if(lasso.fml.lambda.1se!=paste0("cbind( ",top.names[l]," ) ~ ")){
        lasso.model.lambda.1se = lm(as.formula(lasso.fml.lambda.1se), data = x.train) # , weights = labs.wear$weight) # TODO: do we need to include weights?
        res$lasso.val.pred.lambda.1se[[l]] = c(res$lasso.val.pred.lambda.1se[[l]], predict(lasso.model.lambda.1se, newdata = x.test)) # predict on trained model
      } else {
        lasso.model.lambda.1se = null.model # if lasso sets all coeffs = 0, then use the null model
        res$lasso.val.pred.lambda.1se[[l]] = c(res$lasso.val.pred.lambda.1se[[l]], predict(lasso.model.lambda.1se, newdata = x.test)) # predict on trained model
        # cache <-  length(res$val.true[[l]])-length(res$lasso.val.pred.lambda.1se[[l]])
        # res$lasso.val.pred.lambda.1se[[l]] = c(res$lasso.val.pred.lambda.1se[[l]], rep(NA,cache)) # fill with NA(s) if invalid model was supplied
      }
      
      rf.fml = paste("cbind(",paste(top.names[l],collapse=" , "),") ~",paste(rf.variables.to.use,collapse=" + "))
      #check that the formula is valid (i.e. not empty)
      if(rf.fml!=paste0("cbind( ",top.names[l]," ) ~ ")){
        rf.model = randomForest(as.formula(rf.fml), data = x.train)  #weights = labs.wear$weight) # TODO: do we need to include weights?
        res$rf.val.pred[[l]] = c(res$rf.val.pred[[l]], predict(rf.model, newdata = x.test)) # predict on left out person
      } else {
        # rf.model = null.model # if lasso sets all coeffs = 0, then use the null model
        # res$rf.val.pred[[l]] = c(res$rf.val.pred[[l]], predict(rf.model, newdata = x.test)) # predict on left out person
      }
      
      ## pull out features from rf models ##
      rf.features.list <- as.matrix(importance(rf.model)[order(importance(rf.model), decreasing=TRUE),])
      
      
      tmp <- data.frame("test"=rep(top.names[l],length(rf.features.list)),
                        "cv.run"=rep(k,length(rf.features.list)),
                        "left.out.person"=rep(pat,length(rf.features.list)),
                        "rf.feature"=unlist(dimnames(rf.features.list)),
                        "rf.coef.value"=as.numeric(rf.features.list))
      res$rf.features <- rbind(res$rf.features,tmp)
      # t<- anova(bivar.null.lm.model, bivar.lm.model) # to get p-values for model
      # p.value[[l]] <- as.numeric(t[2,][["Pr(>F)"]])  # to get p-values for model
    }
  }
  res
}

# ## QUANTILES OF LASSO FOR EACH TEST
# # Percent var explained (not root because root might be undefined)
# PVE = function(predicted, true){
#   mse = sum((predicted - true)**2)
#   mse_null = sum((mean(true) - true)**2)
#   1 - mse / mse_null
# }
# # Sample some glm results 'under the null' (random outputs)
# glm_sample = function(x,y,reps = 100,nfolds=30){
#   samples = c()
#   for (i in 1:reps){
#     y_random = sample(y)
#     glm.res = cv.glmnet(x=x,y=y_random,
#                         standardize.response=FALSE,
#                         family="gaussian",
#                         nfolds=nfolds,
#                         nlambda=100)
#     lasso.variables.lambda.min <- glm.res$glmnet.fit$beta[,which(glm.res$glmnet.fit$lambda==glm.res$lambda.min)]
#     lasso.fml.lambda.min = paste("cbind(",paste(colnames(x)[1],collapse=" , "),") ~",paste(lasso.variables.to.use.lambda.min,collapse=" + "))
#     lasso.model.lambda.min = lm(as.formula(lasso.fml.lambda.min), data = data.frame(x)) # , weights = labs.wear$weight) # TODO: do we need to include weights?
#     
#     y_model_pred = predict(lasso.model.lambda.min, newdata = data.frame(x))
#     r_model = RPVE(y_model_pred, y_random)
#     samples = c(samples, r_model)
#   }
#   samples
# }
# # Estimate distribution of the glm 'under the null'
# glm_dist = function(test_name, data, variables, reps = 100){
#   dt = as.matrix(wear[,colnames(wear) %in% c(test_name, variables)])
#   dt = na.omit(dt)
#   samples = glm_sample(dt[,-1],dt[,1],reps=5)
#   ecdf(samples)
# }
# # Elementary two-sided test given the null distribution
# emp_2side_pval = function(dist, val){
#   p = dist(val)
#   min(1 - p, p)*2
# }
# 
# # Execute example
# emp_quant = glm_dist("GLU", wear, c(wear.variables, demo.variables), reps = 5)
# emp_2side_pval(emp_quant, -31.15)


combineResults = function(cVS.results, wVS.results, use.Troubleshoot.mode = FALSE){ 
  if (use.Troubleshoot.mode){
    top.names <- c("HGB", "TGL")
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
    if(length(which(!is.na(wVS.results$lasso.val.pred.lambda.manual[[j]])))>=3){
      num.cor.pairs.lasso.lambda.manual <- c(num.cor.pairs.lasso.lambda.manual, length(which(!is.na(wVS.results$lasso.val.pred.lambda.manual[[j]]))))
      #insert step here to check if sample size for correlation test meets minimum threshold
      rsq.lasso.lambda.manual = c(rsq.lasso.lambda.manual, cor(wVS.results$lasso.val.pred.lambda.manual[[j]], wVS.results$val.true[[j]], use = "complete.obs"))
      p.val.rsq.lasso.manual <- c(p.val.rsq.lasso.manual, cor.test(wVS.results$lasso.val.pred.lambda.manual[[j]], wVS.results$val.true[[j]], use = "complete.obs")$p.value)
      #insert step here to check if correlation was significant?
      rssm.lasso.lambda.manual = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$lasso.val.pred.lambda.manual[[j]])^2))
      rss0.lasso.lambda.manual = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$null.val.pred[[j]])^2))
      lasso.pct.var.explained.lambda.manual = c(lasso.pct.var.explained.lambda.manual, (1 - ( rssm.lasso.lambda.manual / rss0.lasso.lambda.manual )))
      lasso.num.Records.lambda.manual <- c(wVS.results$num.Records, length(wVS.results$lasso.val.pred.lambda.manual[[j]]))
    } else {
      p.val.rsq.lasso.manual <- c(p.val.rsq.lasso.manual, NA)
      num.cor.pairs.lasso.lambda.manual <- c(num.cor.pairs.lasso.lambda.manual, 0)
      rsq.lasso.lambda.manual = c(rsq.lasso.lambda.manual, NA)
      lasso.pct.var.explained.lambda.manual = c(lasso.pct.var.explained.lambda.manual, NA)
      lasso.num.Records.lambda.manual <- c(wVS.results$num.Records, NA)
    }
    #lasso (lambda.min)
    if(length(which(!is.na(wVS.results$lasso.val.pred.lambda.min[[j]])))>=3){
      num.cor.pairs.lasso.lambda.min <- c(num.cor.pairs.lasso.lambda.min, length(which(!is.na(wVS.results$lasso.val.pred.lambda.min[[j]]))))
      #insert step here to check if sample size for correlation test meets minimum threshold
      rsq.lasso.lambda.min = c(rsq.lasso.lambda.min, cor(wVS.results$lasso.val.pred.lambda.min[[j]], wVS.results$val.true[[j]], use = "complete.obs"))
      p.val.rsq.lasso.min <- c(p.val.rsq.lasso.min, cor.test(wVS.results$lasso.val.pred.lambda.min[[j]], wVS.results$val.true[[j]], use = "complete.obs")$p.value)
      #insert step here to check if correlation was significant?
      rssm.lasso.lambda.min = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$lasso.val.pred.lambda.min[[j]])^2))
      rss0.lasso.lambda.min = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$null.val.pred[[j]])^2))
      lasso.pct.var.explained.lambda.min = c(lasso.pct.var.explained.lambda.min, (1 - ( rssm.lasso.lambda.min / rss0.lasso.lambda.min )))
      lasso.num.Records.lambda.min <- c(wVS.results$num.Records, length(wVS.results$lasso.val.pred.lambda.min[[j]]))    
    } else {
      p.val.rsq.lasso.min <- c(p.val.rsq.lasso.min, NA)
      num.cor.pairs.lasso.lambda.min <- c(num.cor.pairs.lasso.lambda.min, 0)
      rsq.lasso.lambda.min = c(rsq.lasso.lambda.min, NA)
      lasso.pct.var.explained.lambda.min = c(lasso.pct.var.explained.lambda.min, NA)
      lasso.num.Records.lambda.min <- c(wVS.results$num.Records, NA)   
    }
    #lasso (lambda.1se)
    if(length(which(!is.na(wVS.results$lasso.val.pred.lambda.1se[[j]])))>=3){
      num.cor.pairs.lasso.lambda.1se <- c(num.cor.pairs.lasso.lambda.1se, length(which(!is.na(wVS.results$lasso.val.pred.lambda.min[[j]]))))
      #insert step here to check if sample size for correlation test meets minimum threshold
      rsq.lasso.lambda.1se = c(rsq.lasso.lambda.1se, cor(wVS.results$lasso.val.pred.lambda.1se[[j]], wVS.results$val.true[[j]], use = "complete.obs"))
      p.val.rsq.lasso.1se <- c(p.val.rsq.lasso.1se, cor.test(wVS.results$lasso.val.pred.lambda.1se[[j]], wVS.results$val.true[[j]], use = "complete.obs")$p.value)
      #insert step here to check if correlation was significant?
      rssm.lasso.lambda.1se = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$lasso.val.pred.lambda.1se[[j]])^2))
      rss0.lasso.lambda.1se = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$null.val.pred[[j]])^2))
      lasso.pct.var.explained.lambda.1se = c(lasso.pct.var.explained.lambda.1se, (1 - ( rssm.lasso.lambda.1se / rss0.lasso.lambda.1se )))
      lasso.num.Records.lambda.1se <- c(wVS.results$num.Records, length(wVS.results$lasso.val.pred.lambda.1se[[j]]))    
    } else {
      p.val.rsq.lasso.1se <- c(p.val.rsq.lasso.1se, NA)
      num.cor.pairs.lasso.lambda.1se <- c(num.cor.pairs.lasso.lambda.1se, 0)
      rsq.lasso.lambda.1se = c(rsq.lasso.lambda.1se, NA)
      lasso.pct.var.explained.lambda.1se = c(lasso.pct.var.explained.lambda.1se, NA)
      lasso.num.Records.lambda.1se <- c(wVS.results$num.Records, NA)   
    }
    #rf
    if(length(which(!is.na(wVS.results$rf.val.pred[[j]])))>=3){
      num.cor.pairs.rf <- c(num.cor.pairs.rf, length(which(!is.na(wVS.results$rf.val.pred[[j]]))))
      #insert step here to check if sample size for correlation test meets minimum threshold
      rsq.rf = c(rsq.rf, cor(wVS.results$rf.val.pred[[j]], wVS.results$val.true[[j]], use = "complete.obs"))
      p.val.rsq.rf <- c(p.val.rsq.rf, cor.test(wVS.results$rf.val.pred[[j]], wVS.results$val.true[[j]], use = "complete.obs")$p.value)
      #insert step here to check if correlation was significant?
      rssm.rf = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$rf.val.pred[[j]])^2))
      rss0.rf = sum(na.omit((wVS.results$val.true[[j]] - wVS.results$null.val.pred[[j]])^2))
      rf.pct.var.explained = c(rf.pct.var.explained, (1 - ( rssm.rf / rss0.rf )))
      rf.num.Records <- c(wVS.results$num.Records, length(wVS.results$rf.val.pred[[j]]))
    } else {
      p.val.rsq.rf <- c(p.val.rsq.rf, NA)
      num.cor.pairs.rf <- c(num.cor.pairs.rf, 0)
      rsq.rf = c(rsq.rf, NA)
      rf.pct.var.explained = c(rf.pct.var.explained, NA)
      rf.num.Records <- c(wVS.results$num.Records, NA)
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
  
  fig.2c.df <- cbind(rownames(as.data.frame(cVS.results$sqrt.pct.var)), 
                     as.data.frame(cVS.results$sqrt.pct.var), 
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
  fig.2c.corr.coefs <- cbind(rownames(as.data.frame(cVS.results$rsq.vitals)), 
                             as.data.frame(cVS.results$rsq.vitals), 
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
  list(fig.2c.corr.coefs = fig.2c.corr.coefs, fig.2c.df = fig.2c.df)
}

post.process.results = function(wVS.results)
{
  ##UNCOMMENT IF RUNNING LOCALLY
  lambda.choice <- "lasso.min"
  
  ## calculate correlation coefficients and pct var explained by the models (lambda min)
  
  # remove coefficients associated with participants that had zero x.test data
  wVS.results$num.Records <- data.frame(wVS.results$num.Records)
  cache <- wVS.results$num.Records[wVS.results$num.Records$TestObs=="0",]
  cache <- cache[,c(2,1)] #list of clin tests and iPOPs that had zero x.test data
  names(cache) <- c("test","left.out.person")
  if (is.null(nrow(cache$delete)))
    return()
  cache$delete <- 1
  tmp <- merge(wVS.results$lasso.features.lambda.manual,cache,all=TRUE)
  tmp <- tmp[which(is.na(tmp$delete)),] #remove rows selected for deletion
  wVS.results$lasso.features.lambda.manual <- tmp[,-which(names(tmp)=="delete")] #remove column called "delete"
  tmp <- merge(wVS.results$lasso.features.lambda.min,cache,all=TRUE)
  tmp <- tmp[which(is.na(tmp$delete)),] #remove rows selected for deletion
  wVS.results$lasso.features.lambda.min <- tmp[,-which(names(tmp)=="delete")] #remove column called "delete"
  tmp <- merge(wVS.results$lasso.features.lambda.1se,cache,all=TRUE)
  tmp <- tmp[which(is.na(tmp$delete)),] #remove rows selected for deletion
  wVS.results$lasso.features.lambda.1se <- tmp[,-which(names(tmp)=="delete")] #remove column called "delete"
  wVS.results
}

bootstrap.experiment = function(clin, wear, debug = FALSE, personalized = FALSE, demographics = FALSE, bootstrap = FALSE){
  if (bootstrap){
    clin.boot = bootstrap.dataset(clin)
    wear.boot = bootstrap.dataset(wear)
  }
  else {
    clin.boot = clin
    wear.boot = wear
  }
  
  npatients = NULL
  if (debug)
    npatients = 5
  
  cVS.results = cVS.lm(clin.boot, use.Troubleshoot.mode = debug, use.iPOP = personalized, npatients = npatients, use.Demog = demographics)
  wVS.results = wVS.models(wear.boot, use.Troubleshoot.mode = debug, use.iPOP = personalized, npatients = npatients, use.Demog = demographics)
  fig.tables = combineResults(cVS.results,wVS.results,use.Troubleshoot.mode = debug)
  #  wVS.results = post.process.results(wVS.results)
  
  list(cVS.results = cVS.results, wVS.results = wVS.results, fig.tables = fig.tables)
}

combine.experiments = function(experiments,pernsonal){
  res = list(fig.2c.corr.coefs = c(), fig.2c.df = c())
  for (i in 1:length(personal) ){
    experiment = experiments[[i]]
    coefs = data.frame(experiment$fig.tables$fig.2c.corr.coefs)
    personal.coefs = personal[[i]][,1:2]
    colnames(personal.coefs)[2] = "personal.mean"
    coefs = merge(coefs, personal.coefs, by="test")
    coefs$experiment = i
    res$fig.2c.corr.coefs = rbind(res$fig.2c.corr.coefs, coefs)
    
    df = data.frame(experiment$fig.tables$fig.2c.df)
    df$experiment = i
    res$fig.2c.df = rbind(res$fig.2c.df, df)
  }
  res
}

plot.comparison = function(fig.tables){
  lambda.choice <- "lasso.min"
  
  fig.2c.plot <- gather(fig.tables$fig.2c.corr.coefs, variable, value, -test)
  fig.2c.plot = fig.2c.plot[fig.2c.plot$variable!="experiment",]
  #  fig.2c.plot[,3][is.na(fig.2c.plot[,3])] <- 0 #replace % var explained of NaN w/ 0
  fig.2c.plot = fig.2c.plot %>% group_by(test,variable) %>% summarise(mean = mean(value), sd = sd(value))
  
  tmp = fig.2c.plot[fig.2c.plot$variable == "personal.mean",]
  lvls = as.character(tmp$test[order(-tmp$mean)])
  
  fig.2c.plot$test = factor(fig.2c.plot$test, levels = lvls)
  #^ Ran out of time, but I can simplify this later, which will probably rid the error.
  fig.2c <- fig.2c.plot
  #fig.2c <- fig.2c.plot[order(-fig.2c.plot[,3]),] # reorder by RF
  ## DONE UNCOMMENT
  
  # wVS.results$num.Records <- as.data.frame(wVS.results$num.Records)
  # num.Records.2 <- transform(num.Records, TrainingObs = as.numeric(TrainingObs),
  #                          TestObs = as.numeric(TestObs))
  # Plot the % var explained
  
  ##UNCOMMENT IF RUNNING LOCALLY
  ggplot(fig.2c[fig.2c$variable %in% c("rf.pers","personal.mean"),], aes(x=test, y=mean, color = variable)) +
    geom_errorbar(size = 0.8, aes(ymin=mean-sd, ymax=mean+sd), width=.8, position=position_dodge(width=0.7)) +
    geom_point(size = 3, position=position_dodge(width=0.7)) + #, aes(shape=variable)
    weartals_theme +
    ylim(0,1) +
    scale_color_manual(values=gg_color_hue(5)[c(3,1,2,5,4)]) +
    #    scale_shape_manual(breaks=c("vitals", lambda.choice, "rf","personal.mean"),
    #                         labels=c("LM", "LASSO", "RF","Pers. mean")) +
    #    scale_color_discrete(breaks=c("vitals", lambda.choice, "rf","personal.mean"),
    #                         labels=c("LM", "LASSO", "RF","Pers. mean")) +
    labs(x = "Lab tests",y = expression(paste("Sqrt of % Variance Explained")))
  ## DONE UNCOMMENT
}

personal.mean.model = function(){
  dset = wear
  idvar = "iPOP_ID"
  
  for (test in allClin){
    dset[,test] = as.numeric(gsub("<","",dset[,test]))
    dset[,test] = as.numeric(gsub(">","",dset[,test]))
  }
  res = data.frame(test=NULL, error=NULL)
  
  
  patients = unique(wear$iPOP_ID)
  for (patient in patients){
    test.ids = dset$iPOP_ID %in% patient
    test.idx = which(test.ids)
    
    nfrac = length(test.idx)
    if (nfrac < 6)
      next
    nfrac = floor(length(test.idx)*0.2)
    
    test.idx = sample(test.idx)[1:nfrac]
    test.ids = (1:length(test.ids) %in% test.idx)
    train.ids = !test.ids
    
    
    for (test in allClin){
      model.our = lm(paste(test,"~",idvar),na.omit(dset[train.ids,c(idvar,test)]))
      
      model.null = lm(paste(test,"~ 1"),na.omit(dset[train.ids,c(idvar,test)]))
      dset.test = na.omit(dset[test.ids,c(idvar,test)])
      preds.our = predict(model.our, newdata = dset.test)
      preds.null = predict(model.null, newdata = dset.test)
      r.squared = max(0, 1 - sum((dset.test[,2] - preds.our)**2) / sum((dset.test[,2] - preds.null)**2))
      
      res = rbind(res, data.frame(test=test, error=sqrt(r.squared)))
    }
  }
  res <- transform(res, test=reorder(test, -error) ) 
  
  # ggplot(res, aes(x=test, y=error)) +
  #   ylim(0,1) +
  #   weartals_theme +
  #   geom_point(size = 5)
  res = na.omit(res) %>% group_by(test) %>% summarise_all(funs(mean,sd))
  #  res <- transform(res, test=reorder(test, -mean) ) 
  res
}

plot_tree <- function(final_model, 
                      tree_num) {
  
  # get tree by index
  tree <- randomForest::getTree(final_model, 
                                k = tree_num, 
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(`split var`), NA, `split point`)) %>%
    mutate(`prediction` = ifelse(is.na(`split var`), `prediction`, NA))
  
  tree[["split point"]][is.na(tree[["split var"]])] = NA
  
  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))
  
  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  
  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- format(tree$prediction, digits = 4)
  V(graph)$leaf_label[is.na(tree$prediction)] = NA
  V(graph)$split <- format(tree[["split point"]],digits = 2)
  V(graph)$split[is.na(tree[["split point"]])] <- NA
  print(tree[["split point"]])
  print(V(graph)$split)
  
  # plot
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, vjust = -0.5, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 1.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label),  vjust=1.15, na.rm = TRUE,
                    colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  
  print(plot)
}

library(dplyr)
library(ggraph)
library(igraph)
demo_rf = function(){
  wear.variables <- unlist(read.table(paste0(dir,"FinalLasso_153WearableFactors.csv"), stringsAsFactors = FALSE)) # the table of model features we want to work with
  #wear.variables = c("hr_mean","st_mean","gsr_mean","rhr_mean")
  data = na.omit(wear[,c("HCT",wear.variables)])
  frm = paste0("HCT ~ ",paste(wear.variables,collapse=" + "))
  
  model = randomForest(as.formula(frm), data = data, importance=TRUE, proximity=TRUE, maxnodes=20)
  
  plot_tree(model,2)
}
plt = demo_rf()
ggsave("hct-rf.png",plt,width = 12,height = 5)

reps = 6
cores = 6
debug = TRUE

# reps = 1
# cores = 1
# debug = TRUE

personal = mclapply(1:reps, function(x) { personal.mean.model() },
                    mc.cores = cores)

experiments = mclapply(1:reps, function(x) { bootstrap.experiment(iPOPcorDf, wear, debug = debug, bootstrap = TRUE, demographics = FALSE, personalized = FALSE) },
                       mc.cores = cores)
experiments.nopers = mclapply(1:reps, function(x) { bootstrap.experiment(iPOPcorDf, wear, debug = debug, bootstrap = TRUE, demographics = FALSE, personalized = FALSE) },
                              mc.cores = cores)

# Save experiments
#save(experiments, experiments, file="data/experiments-ipop-pers-dem.Rda")
#save(experiments.nopers, experiments, file="data/experiments-ipop-nopers-dem.Rda")
save(personal,file="data/personal-ipop-100.Rda")

load("data/experiments-ipop-pers.Rda")
load("data/experiments-ipop-nopers.Rda")
load("data/personal-ipop-6.Rda")

fig.tables = combine.experiments(experiments.nopers,personal)
fig.tables.pers = combine.experiments(experiments,personal)

colnames(fig.tables.pers$fig.2c.corr.coefs)[6]="rf.pers"
colnames(fig.tables.pers$fig.2c.df)[c(6,10,14)] = c("rf.pers","num.obs.rf.pers","p.val.rf.pers")

fig.tables$fig.2c.corr.coefs = merge(fig.tables$fig.2c.corr.coefs, fig.tables.pers$fig.2c.corr.coefs[,c("test","rf.pers")])
fig.tables$fig.2c.df = merge(fig.tables$fig.2c.df, fig.tables.pers$fig.2c.df[,c("test","rf.pers","num.obs.rf.pers","p.val.rf.pers")])

#fig.tables = list(fig.2c.df = aggregate(. ~ test, data=fig.tables$fig.2c.df, function(x,na.rm=TRUE){ c(mean(x),sd(x))}, na.rm=TRUE),
#fig.2c.corr.coefs = aggregate(. ~ test, data=fig.tables$fig.2c.corr.coefs, function(x,na.rm=TRUE){ c(mean(x),sd(x))}, na.rm=TRUE))

res.30k = fig.tables$fig.2c.corr.coefs
res = fig.tables
save(res,file = "figure5c.Rda")

plt = plot.comparison(fig.tables)
plt

ggsave("iPOP-wear-model-comparison.png",width=20,height=5)
results = fig.tables$fig.2c.corr.coefs
results$diff = results$rf - results$personal.mean
results = results[,c("test","diff")]
results = results %>% group_by(test) %>% summarise_all(funs(mean,sd))
lvls = as.character(results$test[order(-results$mean)])
results$test = factor(results$test, levels = lvls)

ggplot(results, aes(x=test, y=mean)) +
  geom_errorbar(size = 0.8, aes(ymin=mean-sd, ymax=mean+sd), width=.3, position=position_dodge(width=0.5)) +
  weartals_theme +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(size = 5)

cols = gg_color_hue(3)
testId = 2
allClin[testId]
tvp = data.frame(predicted=experiments[[1]]$wVS.results$rf.val.pred[[testId]], true=experiments[[1]]$wVS.results$val.true[[testId]])
plt = ggplot(tvp, aes(x=true, y=predicted)) +
  weartals_theme +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(size = 3)
plt
ggsave("hct-true-vs-pred.png",plt,width = 6,height = 5)

# number of HCT observations per person
wear.tmp = wear
wear.tmp = cbind(wear.tmp[,c("iPOP_ID", "HCT", "highAct_sk_min")],wear.tmp[,60:100])
wear.tmp = wear.tmp[wear.tmp$iPOP_ID == "1636-69-001",]
wear.tmp = na.omit(wear.tmp)
dim(na.omit(wear.tmp))

