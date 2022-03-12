##
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/food_log/data_preparation/")

library(tidyverse)
rm(list = ls())

load("../../all_omes_wear_food")

dim(all_omes_wear_food)

data = 
  all_omes_wear_food %>% 
  dplyr::filter(MolClass == "Food")

data %>% 
  dplyr::filter(is.na(Intensity))

data$SampleID = as.character(data$DT)

###check duplicated samples for 
library(plyr)
data %>% 
  plyr::dlply(.variables = .(MolName, SampleID)) %>% 
  purrr::map(function(x){
    if(nrow(x) == 1){
      return(NULL)
    }else{
      x
    }
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

sample_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                accurate_time = DT,
                day = day,
                time = tod,
                hour = hour_of_day) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = "Mike1") %>%
  dplyr::select(subject_id, sample_id, everything())


variable_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                class = MolClass,
                subclass = MolSubclass) %>%
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
  dplyr::mutate(variable_id = paste("food", 1:nrow(variable_info), sep = "_")) %>% 
  dplyr::select(variable_id, everything())

rownames(expression_data) <- variable_info$variable_id

expression_data <- 
  expression_data %>% 
  dplyr::select(-mol_name) %>% 
  as.data.frame()

# expression_data %>% 
#   apply(2, unlist) 

dim(expression_data)
dim(variable_info)
dim(sample_info)

colnames(expression_data) == sample_info$sample_id

save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(expression_data, file = "expression_data")



