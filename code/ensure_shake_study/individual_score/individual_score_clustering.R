##
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

###load data
load("data/ensure_shake_study/3_omics/individual_variation/cluster")
cluster <- 
data.frame(subject_id = names(cluster),
           cluster)

masstools::setwd_project()
setwd("data/ensure_shake_study/3_omics/individual_scores")

load("all_score")

dir.create("consensus_clustering")
setwd("consensus_clustering/")

###consensus clustering
library(CancerSubtypes)

result <-
  ExecuteCC(
    clusterNum = 2,
    d = as.matrix(all_score),
    maxK = 6,
    reps = 1000,
    pItem = 0.8,
    pFeature = 0.8,
    title = "k_means_consensus",
    clusterAlg = "km",
    distance = "euclidean",
    plot = "png",
    writeTable = TRUE
  )

save(result, file = "result")
load("result")

cluster
group <- result$group
cluster2 <-
  data.frame(subject_id = names(group), group = group)

final_result <- 
cluster %>% 
  dplyr::left_join(cluster2, by = "subject_id")

sum(final_result$cluster == final_result$group)/23


