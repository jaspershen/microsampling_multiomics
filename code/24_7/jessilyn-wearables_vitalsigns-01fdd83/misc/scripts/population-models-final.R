run.on.patient = function(data, patient.id, test.id, personalized = FALSE, model = "RF", mode = "all", randomized = FALSE, min.visits = 10)
{
  # We assume that the data that went in contains all the variables to use in the model and iPOP_ID
  wear.tmp = filter.nas(data,names(data))
  wear.tmp = wear.tmp[order(wear.tmp$iPOP_ID),]
  
  if (randomized)
    wear.tmp[,test.id] = sample(wear.tmp[,test.id])
  
  # List all patients and remove patient.id for LOO cross-validation
  patients = unique(wear.tmp$iPOP_ID)
  train = patients[!(patients %in% patient.id)]
  test = patients[patients %in% patient.id]
  train.ids = wear.tmp$iPOP_ID %in% train
  test.ids = wear.tmp$iPOP_ID %in% test
  
  if (personalized){
    test.idx = which(test.ids)
    
    # Check if the patient has at least six observations    
    nfrac = length(test.idx)
    if (nfrac < min.visits){
      return(list(pred = c(),
                  true = c(),
                  null = c()))
    }
    
    # if yes, take 20% of their observation and keep them as a test set
    # the rest goes to the training set
    nfrac = floor(length(test.idx)*0.2)
    
    test.idx = sample(test.idx)[1:nfrac]
    test.ids = (1:length(test.ids) %in% test.idx)
    train.ids = !test.ids
    
    # We will use iPOP ID to build a personalized model for that subject
    wear.tmp$iPOP_ID = 1*(wear.tmp$iPOP_ID == patient.id)
  }
  
  # subsect
  x.train = wear.tmp[train.ids,] 
  
  # if we are not building personalized models, remove subject id
  if (!personalized){
    x.train = x.train[,-1,drop=FALSE]
  }
  
  # create test set
  x.test = wear.tmp[test.ids,] # subset input data by lab: only take current lab test of interest
  
  # List all variables going into the model
  vars = names(x.train)
  vars = vars[-which(vars %in% test.id)]
  
  # Choose variables using lasso
  if (mode == "LASSO"){
    glm.res = cv.glmnet(x=as.matrix(x.train[,vars]),y=as.matrix(x.train[,test.id]),
                        standardize.response=FALSE,
                        family="gaussian",
                        nfolds=10,
                        nlambda=100)
    beta.1se = which(glm.res$lambda < glm.res$lambda.1se)[1] # lambda 1se
    vars = rownames(glm.res$glmnet.fit$beta[abs(glm.res$glmnet.fit$beta[,beta.1se]) > 1e-10,]) # TODO: this is an arbitrary rule for now
  }
  
  if (length(vars) == 0){
    vars = "1"
    if (personalized)
      vars = "iPOP_ID"
  }

  # Build the RF model 
  rf.fml = paste(test.id, "~", paste(vars, collapse=" + "))
  if (model == "RF")
    mdl = randomForest(as.formula(rf.fml), data = x.train) 
  
  if (model == "LM")
    mdl = lm(as.formula(rf.fml), data = x.train)
  
  # Build the NULL model
  null.fml = paste(test.id, "~ 1")
  null.model = lm(as.formula(null.fml), data = x.train)
  
  # Return predictions and true values
  list(pred = predict(mdl, newdata = x.test),
       true = x.test[,test.id],
       null = predict(null.model, newdata = x.test))
}

population.loo = function(data, model = "RF", mode="all", personalized = FALSE,
                          debug = FALSE, vars = NULL, randomized = FALSE, num.patients = NULL,
                          kfolds = NULL){
  patients = unique(iPOPcorDf$iPOP_ID)
  if (!is.null(num.patients))
    patients = patients[1:num.patients]
  else
    num.patients = length(patients)

  if (is.null(kfolds))
    kfolds = length(patients)
  
  if (debug){
    top.names <- c("HCT")
    patients = patients[1:50]
  }
  demo.variables = c()
  
  k = 0

  fold_size = num.patients %/% kfolds
  
  res = list()
  for (test.id in top.names){
    cat("Test",test.id,"\n")
    res[[test.id]] = list()
    
#    for (patient.id in patients){
#      cat("Patient",patient.id,"\n") # LOO
#      vars.subset = c("iPOP_ID", test.id, vars, demo.variables)
#      res[[test.id]][[patient.id]] = run.on.patient(data[,vars.subset], patient.id, test.id, personalized = personalized, model = model, mode = mode, randomized = randomized)
#    }
    
    for (i in 1:kfolds){
      test.patients = patients[((i-1)*fold_size):min(num.patients,(i*fold_size-1))]
      cat("Fold",i,"\n") # LOO
      vars.subset = c("iPOP_ID", test.id, vars, demo.variables)
      res[[test.id]][[i]] = run.on.patient(data[,vars.subset], test.patients, test.id, personalized = personalized, model = model, mode = mode, randomized = randomized)
    }
  }
  res
}

wear.data.preprocess = function(wear){
  wear[,top.names] <- apply(
    wear[,top.names], 2,
    function(x) as.numeric(as.character(x)))
  wear
}

reduce = function(res.list){
  ret = c()
  for (res in res.list){
    ret = c(ret, res)
  }
  ret
} 

# Summary metrics
get.stats = function(res){
  stats = data.frame(model = c(), test = c(), rve = c(), ve = c(), pval = c())
  for (r in names(res)){
    for (res.test in names(res[[r]]) ){
      true = reduce(lapply(res[[r]][[res.test]], function(res){ res$true} ))
      null = reduce(lapply(res[[r]][[res.test]], function(res){ res$null} ))
      pred = reduce(lapply(res[[r]][[res.test]], function(res){ res$pred} ))
      
      df = data.frame(model = r,
#                      pval = cor.test(true,pred)$p.val,
#                      pval = t.test( (true-pred)**2, (true-null)**2 )$p.value,
                      test = res.test,
                      rve = sqrt(max(0,1 - sum((true - pred)**2) / sum((true - null)**2))),
                      ve = 1 - sum((true - pred)**2) / sum((true - null)**2)
      )
      stats = rbind(stats, df)
    }
  }
  stats
}

bootstrap.experiment.2d = function(clin, wear, debug = FALSE, bootstrap = FALSE, randomized = FALSE, kfolds = NULL){
  if (bootstrap){
    clin = bootstrap.dataset(clin, replace = TRUE)
    wear = bootstrap.dataset(wear, replace = TRUE)
  }
  
  res = list()
  
  vars.all = unlist(read.table(paste0(dir,"FinalLasso_153WearableFactors.csv"), stringsAsFactors = FALSE))
  
  res[["wear_nopers_lm_lasso"]] = population.loo(wear, debug = debug, personalized = FALSE, randomized = randomized, vars = vars.all, model = "LM", mode = "LASSO", kfolds = kfolds)
  res[["wear_nopers_rf"]] = population.loo(wear, debug = debug, personalized = FALSE, randomized = randomized, vars = vars.all, model = "RF", kfolds = kfolds)
  res[["clin_nopers_rf"]] = population.loo(clin, debug = debug, personalized = FALSE, randomized = randomized, vars = c("Pulse","Temp"), model = "RF", kfolds = kfolds)
  res[["clin_nopers_lm"]] = population.loo(clin, debug = debug, personalized = FALSE, randomized = randomized, vars = c("Pulse","Temp"), model = "LM", kfolds = kfolds)
  
  res
}

bootstrap.experiment.4.5a = function(clin, wear, debug = FALSE, bootstrap = FALSE, randomized = FALSE){
  if (bootstrap){
    clin = bootstrap.dataset(clin, replace = FALSE)
    wear = bootstrap.dataset(wear, replace = FALSE)
  }
  
  res = list()
  
  vars.all = unlist(read.table(paste0(dir,"FinalLasso_153WearableFactors.csv"), stringsAsFactors = FALSE))
  # vars.all = c("hr_mean",#"hr_sd",
  #              "st_mean",#"st_sd",
  #              "gsr_mean",#"gsr_sd",
  #              "sk_mean",#"sk_sd",
  #              "rhr_mean"#"rhr_sd"
  # )
  
  # Figure 4.5
  res[["wear_pers_lm_null"]] = population.loo(wear, debug = debug, personalized = TRUE, vars = numeric(0), model = "LM", num.patients = 1)
  res[["wear_pers_lm"]] = population.loo(wear, debug = debug, personalized = TRUE, vars = vars.all[1:20], num.patients = 1, mode = "LASSO")
  res[["wear_nopers_rf"]] = population.loo(wear, debug = debug, personalized = FALSE, vars = vars.all, model = "RF",num.patients = 1)
  
  res
}

### EXPERIMENTS
## NULL distributions (save summary in null.res)
#res = mclapply(1:100, function(i){bootstrap.experiment.2d(iPOPcorDf, wear.data.preprocess(wear), debug = FALSE, randomized = FALSE, bootstrap = FALSE)}, mc.cores = 100)
#res.randomized = mclapply(1:1000, function(i){bootstrap.experiment.2d(iPOPcorDf, wear.data.preprocess(wear), debug = FALSE, randomized = TRUE, bootstrap = FALSE)}, mc.cores = 100)

## Final runs for 2d (save summary in all.res)
# res = mclapply(1:100, function(i){bootstrap.experiment.2d(iPOPcorDf, wear.data.preprocess(wear), debug = FALSE, randomized = FALSE, bootstrap = TRUE)}, mc.cores = 50)
res = mclapply(1:100, function(i){
  bootstrap.experiment.2d(iPOPcorDf, wear.data.preprocess(wear), debug = FALSE, randomized = FALSE, bootstrap = TRUE, kfolds = 25)
}, mc.cores = 50)
res = mclapply(1:6, function(i){
  bootstrap.experiment.2d(iPOPcorDf, wear.data.preprocess(wear), debug = TRUE, randomized = FALSE, bootstrap = FALSE, kfolds = 25)
  }, mc.cores = 6)

## Experiments for 4.5
# res = mclapply(1:6, function(i){bootstrap.experiment.4.5a(iPOPcorDf, wear.data.preprocess(wear), debug = FALSE, bootstrap = TRUE)}, mc.cores = 6)

### Summary stats for the randomized trials (null hypothesis)
# null.res = data.frame()
# 
# for (r in res.randomized){
#   null.res = rbind(null.res, get.stats(r))
# }

### Summary stats for the actual experiment
all.res = data.frame()

for (r in res){
  all.res = rbind(all.res, get.stats(r))
}

df = data.frame(all.res) %>% group_by(model, test) %>%
  summarise(mean=mean(rve), mean.ve=mean(ve), sd=sd(rve))  %>%
  arrange(desc(model),desc(mean))
df$test = factor(df$test, levels = unique(df$test))

df$min = df$mean - df$sd
df$min[df$min < 0] = 0
df$max = df$mean + df$sd
p = ggplot(df, mapping= aes(x=test,y=mean,color=model)) + 
  geom_errorbar(aes(ymin=min, ymax=max), width=.2) +
  geom_point(size=3) +
  weartals_theme + 
  labs(x = NULL, y ="RPVE")
p

# Save results
save(res, file="experiment5b.Rda")
save(all.res, file="fig5b.Rda")

#source("scripts/extra-plotting/fig2d.R")
