##
no_function()

sxtTools::setwd_project()
setwd("data/")
library(tidyverse)
rm(list = ls())
load("shake_study/shake_omics.RData")

sxtTools::setwd_project()
setwd("data/shake_study/metabolome_data_analysis/data_preparation/")

##metabolomics data
data <- 
  omics %>% 
  dplyr::filter(Class == "Metabolite")

sample_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                subject_id = PID,
                TP = TP) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame()

variable_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                pathway = Pathway,
                subclass = Subclass,
                p.value) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame()

expression_data <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                mol_name = MolName,
                intensity = Intensity) %>%
  tidyr::pivot_wider(names_from = sample_id,
                     values_from = intensity) %>% 
  as.data.frame()

expression_data$mol_name[1:100] == variable_info$mol_name[1:100]

variable_info <- 
  variable_info %>% 
  dplyr::mutate(variable_id = paste("peak", 1:nrow(variable_info), sep = "_")) %>% 
  dplyr::select(variable_id, everything())

rownames(expression_data) <- variable_info$variable_id

expression_data <- 
  expression_data %>% 
  dplyr::select(-mol_name) %>% 
  as.data.frame()

dim(expression_data)
dim(variable_info)
dim(sample_info)

colnames(expression_data) == sample_info$sample_id


###information from metID
metid1 <- readr::read_csv("Annotated_metabolomics_whole_blood_Shake_Study.csv")

metid3 <- readr::read_csv("Ryan_Shake_MetID.csv")

metid1$Level <- metid3$Confidence_Level

expression_data$S1_T0 == metid1$S1_T0

metid1 <- 
  metid1 %>% 
  dplyr::select(-contains("_T"))

metid1 <-
metid1 %>% 
  dplyr::select(Compound:Retention.time..min., Metabolite_metID:Mode, Level) %>% 
  dplyr::select(variable_id = Compound, mz = m.z, rt = Retention.time..min., everything()) %>% 
  dplyr::mutate(rt = rt * 60)

colnames(metid1) <- stringr::str_replace(colnames(metid1), "\\_metID", "")

metid1 
variable_info <- 
  cbind(metid1, variable_info[,-1])

rownames(expression_data) <- variable_info$variable_id


##remove duplicated metabolites
dir.create("peaks")
dir.create("metabolites")

save(expression_data, file = "peaks/expression_data")
save(variable_info, file = "peaks/variable_info")
save(sample_info, file = "peaks/sample_info")

idx <- 
  which(duplicated(variable_info$Metabolite[!is.na(variable_info$Metabolite)]))

duplicated_name <- 
  variable_info$Metabolite[!is.na(variable_info$Metabolite)][idx] %>% 
  unique()

remove_name <- NULL

for(name in duplicated_name){
  cat(name, " ")
  temp_idx <- which(variable_info$Metabolite == name)
  remain_name <- 
  variable_info[temp_idx,] %>% 
    dplyr::filter(Confidence_Level == 1) %>% 
    dplyr::filter(Total_Score == max(Total_Score)) %>% 
    dplyr::pull(variable_id)
  if(length(remain_name) > 1){
    stop("error")
  }
  
  temp_remove_name <- setdiff(variable_info[temp_idx,]$variable_id, remain_name)
  remove_name <- c(remove_name, temp_remove_name)  
}

variable_info <-
  variable_info %>% 
  dplyr::filter(!variable_id %in% remove_name)

expression_data <- 
  expression_data[variable_info$variable_id,]

save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(expression_data, file = "expression_data")

save(sample_info, file = "metabolites/sample_info")
save(variable_info, file = "metabolites/variable_info")
save(expression_data, file = "metabolites/expression_data")
















