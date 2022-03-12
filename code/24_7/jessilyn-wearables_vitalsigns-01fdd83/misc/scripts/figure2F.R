##############
#  Figure 2F #
##############

rf.features <-read.table(paste0(dir,"20180622/20180621_DayPrior_noDemog_RF_Features.csv"),
                         header=FALSE,sep=',',stringsAsFactors=FALSE)
colnames(rf.features) <- c("test", "cv.run", "iPOP_ID", "feature", "coefficient")
rf.feature.summaries <- as.data.frame(summarise(group_by(rf.features, test, feature),
                                                mean=mean(coefficient), sd=sd(coefficient)))
top.models <- c("HCT", "HGB", "RBC", "MONOAB", "GLU", "UALB",  "CL", "A1C")
top.rf.features <- c()
for (i in top.models) {
  rf.feature.subset<-rf.feature.summaries[rf.feature.summaries$test %in% i,] 
  rf.feature.sorted <- rf.feature.subset[order(rf.feature.subset$mean, decreasing=TRUE),]
  top.rf.features <- rbind(top.rf.features, head(rf.feature.sorted[1:10,]))
}

rf.feature.summaries[rf.feature.summaries$test %in% top.models,]
