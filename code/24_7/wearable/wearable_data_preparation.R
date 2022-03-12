##
##steps
no_function()

library(tidyverse)

###sleep
sxtTools::setwd_project()
setwd("data/7_24_mike/sleep/data_preparation/")
load("sleep")

sleep %>% 
  dplyr::filter(Day == "2019-04-29") %>% 
  ggplot(aes(DT, 1)) +
  geom_tile(aes(fill = Intensity))

sleep$sample_id =
   as.character(sleep$DT)

sample_info <-
  sleep %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id,
                accurate_time = DT,
                day = Day,
                time = Time,
                time_window,
                seconds) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = "Mike1") %>%
  dplyr::select(subject_id, sample_id, everything())

sample_info$hour = 
  sample_info$time %>% 
  as.character() %>% 
  stringr::str_split(pattern = "\\:") %>% 
  purrr::map(function(x){
    x[1]
  }) %>% 
  unlist() %>% 
  as.numeric()

variable_info <-
  sleep %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                class = MolClass) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::mutate(variable_id = "sleep1") %>% 
  dplyr::select(variable_id, everything())

expression_data <-
  sleep %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id,
                mol_name = MolName,
                intensity = Intensity) %>%
  tidyr::pivot_wider(names_from = sample_id,
                     values_from = intensity) %>% 
  as.data.frame()

expression_data <- 
  expression_data %>% 
  dplyr::select(-mol_name) %>% 
  as.data.frame()

rownames(expression_data) = variable_info$variable_id

colnames(expression_data) == sample_info$sample_id

sample_info = 
sample_info %>% 
  dplyr::arrange(accurate_time)

expression_data = 
  expression_data[,sample_info$sample_id]

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")



################################################################################
######steps
sxtTools::setwd_project()
setwd("data/7_24_mike/steps/data_preparation/")

library(tidyverse)
rm(list = ls())

load("../../all_omes_wear_food")

dim(all_omes_wear_food)

data = 
  all_omes_wear_food %>% 
  dplyr::filter(MolClass == "Steps")

data$SampleID = as.character(data$DT)

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
  dplyr::mutate(variable_id = paste("step", 1:nrow(variable_info), sep = "_")) %>% 
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

save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(expression_data, file = "expression_data")



################################################################################
##HR heart rate
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/hr/data_preparation")

library(tidyverse)
rm(list = ls())

load("../../all_omes_wear_food")

dim(all_omes_wear_food)

data = 
  all_omes_wear_food %>% 
  dplyr::filter(MolClass == "HR")

data$SampleID = as.character(data$DT)

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
  dplyr::mutate(variable_id = paste("hr", 1:nrow(variable_info), sep = "_")) %>% 
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

save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(expression_data, file = "expression_data")





################################################################################
##CGM continuous glucose monitoring
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/cgm/data_preparation/")

library(tidyverse)
rm(list = ls())

load("../../all_omes_wear_food")

dim(all_omes_wear_food)

data = 
  all_omes_wear_food %>% 
  dplyr::filter(MolClass == "CGM")

data$SampleID = as.character(data$DT)

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
  dplyr::mutate(variable_id = paste("cgm", 1:nrow(variable_info), sep = "_")) %>% 
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

save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(expression_data, file = "expression_data")















