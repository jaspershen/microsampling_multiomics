#########################
# Supplementary Table 1 #
#########################
# get mean +/- SE for each of the 50 clin values
library("psych")
thirtyk.summary <- describe(corDf)
ipop.corDf.summary <- describe(iPOPcorDf)
ipop.demog.summary <- describe(iPOPdemographics$AgeIn2016)
table(iPOPdemographics)
write.csv(thirtyk.summary, "/Users/jessilyn/Desktop/framework_paper/Suppl_Table_1/supplTable1.csv")
write.csv(ipop.corDf.summary, "/Users/jessilyn/Desktop/framework_paper/Suppl_Table_1/supplTable2.csv")

raw.wear <- read.csv(paste0(dir,
                            "BasisData_20161111_PostSummerAddOns_Cleaned_NotNormalized_20180427.csv"),
                     header=TRUE,sep=',',stringsAsFactors=FALSE)

describe(raw.wear$GSR)
describe(raw.wear$Heart_Rate)
describe(raw.wear$Steps)

# get p-values and rsquared of univariate correlations
vitalVars <- which(names(corDf) %in% c("Pulse","Temp"))
results <- corr.test(corDf[,3:56],
                     corDf[,c(vitalVars)],
                     method="pearson",adjust="fdr")
#### CREATE DATA FRAME(S) SUMMARIZING RESULTS ####
rCols <- paste0(dimnames(results$r)[[2]],"_R")
pCols <- paste0(dimnames(results$p)[[2]],"_P")
nCols <- paste0(dimnames(results$n)[[2]],"_N")
allCors_R <- data.frame(placeHolder=rep(NA,length(allClin)))
rownames(allCors_R) <- dimnames(results$r)[[1]]
for(i in 1:length(rCols)){
  allCors_R[,i] <- results$r[,i]
  names(allCors_R)[i] <- rCols[i]
}

allCors_P <- data.frame(placeHolder=rep(NA,length(allClin)))
rownames(allCors_P) <- dimnames(results$p)[[1]]
for(i in 1:length(pCols)){
  allCors_P[,i] <- results$p[,i]
  names(allCors_P)[i] <- pCols[i]
}

allCors_N <- data.frame(placeHolder=rep(NA,length(allClin)))
rownames(allCors_N) <- dimnames(results$n)[[1]]
for(i in 1:length(nCols)){
  allCors_N[,i] <- results$n[,i]
  names(allCors_N)[i] <- nCols[i]
}

allCors <- cbind(allCors_R,allCors_P,allCors_N)

newColOrder <- c("Pulse","Temp")

tmp <- allCors
for(i in 1:length(newColOrder)){
  if(i == 1){
    cache <- grep(newColOrder[i],names(allCors))
    tmp[,c(i:sum(0+3))] <- allCors[,c(cache)]
    names(tmp)[c(i:sum(0+3))] <- c(names(allCors)[c(cache)])
    count <- sum(0+3)
  } else {
    if(i != 1){
      cache <- grep(paste0("^",newColOrder[i]),names(allCors))
      tmp[,c(sum(count+1):sum(count+3))] <- allCors[,c(cache)]
      names(tmp)[c(sum(count+1):sum(count+3))] <- c(names(allCors)[c(cache)])
      count <- sum(count+3)
    }
  }
}

allCors <- tmp
sigCors_Ranked <- allCors

#class(sigCors_Ranked$Pulse_P)
sigCors_Ranked <- sigCors_Ranked[
  which(sigCors_Ranked$Pulse_P < 0.05),]

sigCors_Ranked <- sigCors_Ranked[
  order(abs(sigCors_Ranked$Pulse_R),decreasing = TRUE),]

#### SAVE DATA FRAME(S) SUMMARIZING RESULTS ####

write.csv(allCors,
          "~/Desktop/20170810_allCors_VitalsVsLabs.csv",
          row.names=TRUE)

write.csv(sigCors_Ranked,
          "~/Desktop/20170810_sigCorsRanked_VitalsVsLabs.csv",
          row.names=TRUE)
