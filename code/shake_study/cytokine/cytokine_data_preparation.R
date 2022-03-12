##
no_function()


sxtTools::setwd_project()
setwd("data/")
library(tidyverse)
rm(list = ls())
load("shake_omics.RData")

sxtTools::setwd_project()
setwd("data/cytokine_data_analysis/data_preparation/")

##metabolomics data
data <- 
  omics %>% 
  dplyr::filter(Class == "Cytokine")

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

expression_data$mol_name == variable_info$mol_name

variable_info <- 
  variable_info %>% 
  dplyr::mutate(variable_id = paste("cytokine", 1:nrow(variable_info), sep = "_")) %>% 
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

variable_info


save(sample_info, file = "sample_info")

save(variable_info, file = "variable_info")

save(expression_data, file = "expression_data")
















