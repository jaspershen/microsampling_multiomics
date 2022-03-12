################################################
#  Figure 2E  - Canonical Correlation Analysis #
################################################

# should we make this regularized - what is the right proportion between observations and features to decide this (I think there we have ~2X the numver of obs as features)

# library("devtools)
# install_url("https://cran.r-project.org/src/contrib/Archive/impute/impute_1.26.0.tar.gz")
# install.packages("PMA)

library("PMA")
library("Hmisc")
wear.variables <- unlist(read.table(paste0(dir,"FinalLasso_153WearableFactors.csv"), stringsAsFactors = FALSE)) # the table of model features we want to work with


#TODO: UALAB coding causes problems! I removed it from diabetes
clinical.groups = list()
clinical.groups[["Electrolytes"]] =c("CA","K","CL","CO2","NA.","AG")
clinical.groups[["Diabetes"]] =c("A1C","ALB","GLU","CR","ALCRU") #"UALB",
clinical.groups[["Cardiovascular Disease"]]=c("CHOL","LDLHDL","HDL","CHOLHDL","NHDL","TGL","LDL")
clinical.groups[["Cardiometabolic Disease"]]=c("A1C","ALB","GLU","UALB","CR","ALCRU","CHOL"," LDLHDL","HDL","CHOLHDL","NHDL","TGL","LDL")
clinical.groups[["Liver Function"]]=c("ALKP","BUN","ALT","TBIL","AST")
clinical.groups[["Inflammation"]]=c("BASO","LYM","LYMAB","MONO","MONOAB","NEUT","NEUTAB","IGM","EOS","EOSAB","BASOAB","WBC","HSCRP")
clinical.groups[["Hematologic"]] = c("PLT","GLOB","TP","HGB","HCT","RDW","MCH","MCV","RBC","MCHC")
#clinical.groups[["Cardiometabolic.Disease"]]=c("A1C","ALB","GLU","UALB","CR","ALCRU","CHOL"," LDLHDL","HDL","CHOLHDL","NHDL","TGL","LDL")
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

get.cca = function(wear,name){
  nm = name
  # print(nm)
  # Remove rows with NAs
  data.clin = correct.vars(wear[,which(colnames(wear) %in% c("iPOP_ID", clinical.groups[[nm]]))])
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
  d <- d[,!apply(tmp,2,function(x) any(x > 0.99999999999))] # how does it choose which variable to get rid of? Does it matter which one bc they are linear combos of eachother?
  
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
      model.cc.cv = CCA.permute(train[,(ncol(data.clin)):(ncol(train))],  # TODO: cehck that this is choosing the correct columns for the right and left hand sides
                                train[,1:(ncol(data.clin)-1)],
                                trace = FALSE,
                                standardize = FALSE,
                                penaltyxs = pensx,
                                penaltyzs = pensz)
      model.cc = CCA(train[,(ncol(data.clin)):(ncol(train))],  # TODO: cehck that this is choosing the correct columns for the right and left hand sides
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
  list(indexX = indexX, indexY = indexY)
}

# Bootstrap standard deviation
cca.experiment = function(i){
  set.seed(i)
  results = c()
  for (nm in names(clinical.groups)){
    res = get.cca(bootstrap.dataset(wear), nm)
    indexX = res$indexX
    indexY = res$indexY
    cca.corr <- abs(cor(indexX, indexY))
    # print(ggplot(data.frame(indexX = indexX, indexY = indexY),aes(indexX,indexY)) +
    #         weartals_theme +
    #         geom_point() +
    #         ggtitle(nm) +
    #         stat_summary(fun.data=mean_cl_normal) +
    #         geom_smooth(method='lm',formula=y~x))
    # print(cca.corr)
    results = c(results, cca.corr)
  }
  results
}
res.cca = mclapply(1:12, cca.experiment, mc.cores = 6)
res.cca = t(simplify2array(res.cca))
colnames(res.cca) = names(clinical.groups)

# Get means from the experiment on full data
res = c()
for (nm in names(clinical.groups)){
  r = get.cca(wear, nm)
  res = c(res,abs(cor(r$indexX, r$indexY)))
}
res.cca.means = data.frame(test=names(clinical.groups),value=res)

save(res.cca, res.cca.means, file="res.cca.Rda")
load("res.cca.Rda")

source("scripts/extra-plotting/fig2f.R")
