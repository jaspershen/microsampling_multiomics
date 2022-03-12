############
# Figure 4 #
############
## Figure 4A: Check how individual means perform

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

# For the response letter
getTestStatsDemo = function(test, th = 10){
  hct = na.omit(iPOPcorDf[,c("iPOP_ID",test)])
  hct.mn = aggregate(hct[[test]], list(iPOP_ID = hct$iPOP_ID), mean)
  hct.var = aggregate(hct[[test]], list(iPOP_ID = hct$iPOP_ID), var)
  hct.cnt = aggregate(hct[[test]], list(iPOP_ID = hct$iPOP_ID), length)
  
  mn = mean(hct.var$x[hct.cnt$x >= th])
  wmn = sum(hct.var$x[hct.cnt$x >= th] * sqrt(hct.cnt$x[hct.cnt$x >= th]) / sum(sqrt(hct.cnt$x[hct.cnt$x >= th])))
  ssd = sd(hct.mn$x[hct.cnt$x >= th])
  list(diff=mn-wmn,ssd=ssd)
}
res.all = data.frame()
for (test in c("A1C")){
  res = data.frame(getTestStatsDemo(test))
  res$test = test
  res.all = rbind(res.all, res)
}

plot(hct.mn$x)
th=1

mn = mean(hct.mn$x[hct.cnt$x >= th])
wmn = sum(hct.mn$x[hct.cnt$x >= th] * hct.cnt$x[hct.cnt$x >= th] / sum(hct.cnt$x[hct.cnt$x >= th]))


data = na.omit(iPOPcorDf[,c("iPOP_ID",test)])
a <- data %>% 
  group_by(iPOP_ID) %>% 
  summarise(Response=mean(A1C),se=sd(A1C),cnt=length(A1C)) %>% ungroup %>%
  ggplot(aes(x=iPOP_ID,y=Response)) + 
  geom_point(show.legend = FALSE) +
  weartals_theme + 
  geom_errorbar(aes(ymin = Response - se, ymax = Response + se), width = 0.1) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("A1C") + 
  xlab("iPOP cohort participant") + 
  geom_hline(yintercept=mn,color="red") + 
  geom_hline(yintercept=wmn,color="blue")
a

getDiastolicSlope = function(test){
  hct = na.omit(corDf[,c("ANON_ID","Diastolic",test)])
  hct.mn = aggregate(list(HCT = hct[[test]], diastolic=hct$Diastolic), list(ANON_ID = hct$ANON_ID), mean)
  
  model = lm(paste0(test," ~ Diastolic"),data=hct)
  model
}
res = getDiastolicSlope("HCT")
summary(res)

# Compute stats of an LM model of a certain test given the coefs
# The following script cross-validates by taking one observation from each
# patient with at least 4 observations. The model simply 
getLMresults = function(corDf4A, test.name, threshold, identifier, cap, model_coefs, vits = c("Temp","Pulse","Systolic","Diastolic","Respiration"), model_vars=c(), threshold_hi = 1e7, mixed=FALSE, type="LM", min_patients = 10){
  # TODO: that's an ugly wayL to remove NAs but correct
  corDf.tmp.full = filter.nas(corDf4A, c(test.name,model_vars,vits))
  
  ids = sort(table(corDf.tmp.full[[identifier]]))
  atleastafew = names(ids[(ids >= threshold) & (ids <= threshold_hi)]) # select patients with at least 10 tests
  atleastafew = atleastafew[1:min(cap,length(atleastafew))] # troubles with training bigger models

    corDf.tmp = corDf.tmp.full[corDf.tmp.full[[identifier]] %in% atleastafew,]
  testids = c()

  print(dim(corDf.tmp))
  
  # compute correlation if we have at least 10 patients 
  if (length(atleastafew) > min_patients){
    # Select the last observation from each patient who qualified (at least n0 obs)
    for (id in atleastafew){
      lastvisit = sample(which(corDf.tmp[[identifier]] == id),ceiling(length(atleastafew)*0.2))
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
      onehot = model.matrix(as.formula(paste0("~ ",identifier," - 1")), corDf.tmp.rf)
      colnames(onehot) = paste0("ID",1:ncol(onehot))
      corDf.tmp.rf = cbind(corDf.tmp.rf, onehot)
      corDf.tmp.rf[[identifier]] = NULL
      
      model = randomForest(as.formula(paste(test.name,"~ .")), data=corDf.tmp.rf[-testids,])
    }
    else{
      corDf.tmp[[identifier]] = factor(corDf.tmp[[identifier]])
      model = lm(frm, data=corDf.tmp[-testids,])
    }
    
    if (type=="RF"){
      preds = predict(model, newdata = corDf.tmp.rf[testids,])
      print(preds)
    }
    else{
      preds = predict(model, newdata = corDf.tmp[testids,])
    }
    #    print(summary(model))
    
    # Correlation
    # res.meanpred = c(res.meanpred, cor(corDf.tmp[testids,test.name],preds))
    
    # Sqrd root of variance explained
    var.exp = sum( (corDf.tmp[testids,test.name] - preds)**2)
    var.null = sum( (corDf.tmp[testids,test.name] - mean(corDf.tmp.full[-testids,test.name]))**2)
    if (test.name=="LYM")
      plot(corDf.tmp[testids,test.name], preds)
    if (var.exp/var.null > 1)
      res = 0
    else 
      res = sqrt(1 - var.exp/var.null)
  }
  else {
    res = NA
    testids = NULL
  }
  print(res)
  list(res=res,n=length(testids))
}

# Loop over all tests and all models for the personal medels comparsion. Plot results
generate4A = function(dataset, dataset_type, threshold = 4, cap = 10, threshold_hi = 1e7, ntest = NULL){
  if (dataset_type == "iPOP"){
    identifier = "iPOP_ID"
    corDf4A = dataset
    test.names = allClin
    
    corDf4A = correct.vars(corDf4A)
    vits = c("Pulse","Temp","systolic","diastolic")
  }
  else{
    identifier = "ANON_ID"
    corDf4A = dataset
    test.names = intersect(allClin,colnames(corDf4A))
    vits = c("Pulse","Temp","Systolic","Diastolic","Respiration")
  }
  
  res.values = c()
  res.n = c()
  res.models = c()
  res.tests = c()
  if (is.null(ntest))
    ntest = length(test.names)

  for (test.name in test.names[1:ntest]){
    models = list()
    model_vars = c()

    if (dataset_type == "iPOP"){
      model_vars = names(wear)[57:563] #c("gsr_mean")
    }
      # population vitals model
      models[["vitals"]] = paste(vits, collapse = " + ")
      
      # population vitals + personal intercept model
      models[["vitals + personal mean"]] = paste0(identifier," + ",paste(vits, collapse = " + "))
      
      # population vitals + personal intercept model
      models[["vitals + personal mean and slope"]] = paste(models[["vitals + personal mean"]],
                                                           " + (",paste(vits,sep = " + ")," | ",identifier,")")
      
      # population vitals + personal intercept model
      models[["personal mean"]] = paste0(identifier)

      # RF MODELS (only vitals & personal)
    models_rf = list("vitals + personal mean","personal mean")
    for (mdl in models_rf){
      res.tests  = c(res.tests, test.name)
      model = getLMresults(corDf4A, test.name, threshold, identifier, cap, models[[mdl]],
                           model_vars = model_vars,threshold_hi = threshold_hi, mixed = FALSE, type = "RF",
                           vits = vits)
      res.values = c(res.values, model$res)
      res.models = c(res.models, paste(mdl,"(RF)"))
      res.n = c(res.n, model$n)
      print("RF done")
    }
    
    for (mdl in names(models)){
      res.tests  = c(res.tests, test.name)
      mixed = "vitals + personal mean and slope" == mdl
      model = getLMresults(corDf4A, test.name, threshold, identifier, cap, models[[mdl]],
                           model_vars = model_vars,threshold_hi = threshold_hi, mixed = mixed,
                           vits = vits)
      res.values = c(res.values, model$res)
      res.models = c(res.models, mdl)
      res.n = c(res.n, model$n)
    }
  }
  res = data.frame(test = res.tests, model = res.models, value = res.values, n = res.n)
  res = na.omit(res)
  res
}
# TODO: Increase cap in the final version to get better accuracy. Lower caps are for speeding up
#  cap - cut of number of patients for building the individual models (the higher the better population slope estimates)
#  threshold - minimum number of visits for being included in the model (the higher the more accurate personal models)
#res = generate4A("30k",threshold = 50, cap = 500)
#res = generate4A(wear, "iPOP",threshold = 5, cap = 50, ntest=2)
experiments5A = mclapply(1:4, function(x) { 
  # IPOP
  corDf.boot = bootstrap.dataset(iPOPcorDf)
  generate4A(corDf.boot, "iPOP",threshold = 4, cap = 100, ntest = 2)
}, mc.cores = 4)
#save(experiments5A,file="experiments5A.Rda")
load("experiments5A.Rda")

experiments5B = mclapply(1:4, function(x) { 
  # 30k
  corDf.boot = bootstrap.dataset(corDf)
  generate4A(corDf.boot, "30k",threshold = 50, cap = 50, ntest = 2)
}, mc.cores = 4)
#save(experiments5B,file="experiments5B.Rda")
load("experiments5B.Rda")

boxplot(iPOPcorDf$IGM ~ iPOPcorDf$iPOP_ID) 
boxplot(iPOPcorDf$CO2 ~ iPOPcorDf$iPOP_ID) 

experiments = experiments5B
res = data.frame()
for (i in 1:length(experiments)){
  experiment = experiments[[i]]
  experiment$experiment_id = i
  res = rbind(res, experiment)
}
res.ipop = res

toplot = group_by(res,test,model) %>% summarise(mean = mean(value), sd = sd(value))
toplot = toplot[order(-toplot$mean),]
toplot = toplot[order(toplot$model),] # order by the population vitals model
toplot$test = factor(as.character(toplot$test), levels = unique(as.character(toplot$test)))

pp = ggplot(toplot, aes(test, mean, group = model, color = model)) +
  ylab(expression(sqrt("Variance explained"))) +
  xlab("Lab test") +
  geom_point(size = 3, position=position_dodge(width=0.7)) +
  geom_errorbar(size = 0.8, aes(ymin=mean-sd, ymax=mean+sd), width=.8, position=position_dodge(width=0.7)) +
  scale_color_manual(values=gg_color_hue(6)[c(1,2,3,5,4,6)]) +
  weartals_theme + 
  theme(text = element_text(size=14))
print(pp)
ggsave(paste0("plots/Figure-4A.png"), plot = pp, width = 20, height = 4.5)


res.wide = spread(res, model, value)
res.wide$diff = res.wide$`vitals + personal mean and slope` - res.wide$`personal mean`
data.frame(group_by(res.wide, test) %>% summarise(p.val = t.test(diff,y=rep(0,100),"greater")$p.value,
                                                  improvement.mean = mean(diff),
                                                  improvement.sd = sd(diff),
                                                  times.better = mean(diff > 0)))

tt = t.test(res.wide[res.wide$test == "HCT","diff"])

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

# Find patients and tests with the following conditions:
#  * patients with > 20 observations
#  * personal models work great (high R^2)
#  * coefficients are far apart
identifyNiceCase = function(clin,dataset){
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
  }
  # Alternatively, we select people with the largest number of observations
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
res = identifyNiceCase("HCT","30k")
res[order(-res$error),]

## Figure 4D and 4F: Individual models Lab ~ Vital
generate4DF = function(clin,vit,dataset = "30k"){
  if (dataset == "iPOP"){
    identifier = "iPOP_ID"
    corDf.tmp = iPOPcorDf[!is.na(iPOPcorDf[[clin]]),]
  }
  else{
    identifier = "ANON_ID"
    corDf.tmp = corDf[!is.na(corDf[[clin]]),]
  }
  # Use loess to estimate the population model
#  vitals = c("Pulse","Temp","Systolic","Diastolic","Respiration")
  vitals = c("Pulse","Temp","Systolic","Diastolic","Respiration")
  
  
  for (vital in vitals)
    corDf.tmp = corDf.tmp[!is.na(corDf.tmp[,vital]),]

  # Here we can select a few ANNON_ID
  if (dataset == "iPOP")
    toppat = c("1636-70-1005","1636-70-1008","1636-70-1014","1636-69-001")
  else
    toppat = c("PD-4419","D-6050","N-8819","N-5362","D-3086")
  
  toppat = unique(corDf.tmp[[identifier]])
  pat.analyzed = c()
  
  # Alternatively, we select people with the largest number of observations
  # toppat = names(sort(-table(corDf.tmp[[identifier]]))[1:5])
  
  dd = corDf.tmp[corDf.tmp[[identifier]] %in% toppat ,c(identifier,vit,clin)]
  
  # Use lm to estimate the population model
  frm = paste0(clin," ~ ",vit)
  ww = lm(frm, corDf.tmp[sample(nrow(corDf.tmp)),])
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
    
    if (nrow(corDf.ind) < 60 || pat == "PD-9999")
      next
    pat.analyzed = c(pat.analyzed, pat)
    
    # Var explained (CV)
    errors = c()
    errors.persmean = c()
    for (i in 1:1000){
      #testids = (nrow(corDf.ind) - 10):nrow(corDf.ind)
      testids = sample(nrow(corDf.ind))[1:ceiling(nrow(corDf.ind)*0.2)]
      model = lm(frm, data = corDf.ind[-testids,])
      
      # personal model
      preds = predict(model, newdata = corDf.ind[testids,])
      
      var.null = sum( (corDf.ind[testids,clin] - mean(corDf.tmp[,clin])) ** 2)

      # compute var exp
      var.exp = sum( (corDf.ind[testids,clin] - preds) ** 2)
      err = max(1 - var.exp / var.null,0)
      
      errors = c(errors, sqrt(err))
      
      var.exp = sum( (corDf.ind[testids,clin] - mean(corDf.ind[testids,clin])) ** 2)
      err = max(1 - var.exp / var.null,0)
      
      errors.persmean = c(errors.persmean, sqrt(err))
    }
    err = mean(errors)
    err.pm = mean(errors.persmean)
    
    plot(corDf.ind[testids,clin])
    abline(h=mean(corDf.tmp[,clin]),col="red")
    abline(h=mean(corDf.ind[testids,clin]),col="blue")
    
    # Correlation
    # model.sum = summary(model)
    # err = sqrt(model.sum$r.squared)
    
    dd[dd[[identifier]] == pat,]$accuracy = paste0(pat," (r=",sprintf("%.1f", err),")")
    
    slope = lm(paste0(clin," ~ ",vit),corDf.ind)$coefficients[vit]
    stats = rbind(stats, c(err,err.pm,slope,nrow(corDf.ind)))
  }
  rownames(stats) = pat.analyzed
  colnames(stats) = c("sqvarexp","sqvarexp.pm",vit,"visits")
  
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
set.seed(1)
stats = generate4DF("HCT","Diastolic","30k")
print(stats)

generate4E = function(clin,vit,dataset = "30k"){
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
generate4E("HCT","Diastolic","30k")

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
#        print(qq, vp = viewport(layout.pos.row = i, layout.pos.col = j))
        
      }
    }
  }
  dev.off()
  
  cf
}
res = generate4C("30k")
