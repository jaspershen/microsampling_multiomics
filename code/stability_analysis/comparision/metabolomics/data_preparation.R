no_function()

masstools::setwd_project()
setwd("data/stability_analysis/metabolomics/")
rm(list = ls())
library(tidyverse)

data <-
  readr::read_csv("Ryan_mitra_plasma_metabolomics_all_Anno (1).csv")

data_hilic_pos <- readr::read_csv("210113_RyanMitra_HILIC_pos.csv")
data_hilic_neg <- readr::read_csv("210113_RyanMitra_HILIC_neg.csv")
data_rplc_pos <- readr::read_csv("2100113_RyanMitra_RPLC_pos.csv")
data_rplc_neg <- readr::read_csv("210113_RyanMitra_RPLC_neg.csv")

variable_info <-
  data %>%
  dplyr::select(`...1`:"Minimum.CV.", "Accepted.Compound.ID":"Pathway")

expression_data <-
  data %>%
  dplyr::select(-c(`...1`:"Minimum.CV.", "Accepted.Compound.ID":"Pathway"))


sample_info <-
  data.frame(sample_id = colnames(expression_data)) %>%
  dplyr::mutate(subject_id = stringr::str_replace(sample_id, "P|M", "")) %>%
  dplyr::mutate(class = stringr::str_extract(sample_id, "P|M"))

variable_info <-
  variable_info %>%
  dplyr::select(-c(`...1`, "X")) %>%
  dplyr::rename(variable_id = Compound,
                mz = m.z,
                rt = Retention.time..min.) %>%
  dplyr::mutate(rt = rt * 60) %>%
  dplyr::select(variable_id, mz, rt, everything())

variable_info$variable_id <-
  paste(variable_info$variable_id, variable_info$Mode, sep = "_")

rownames(expression_data) <-
  variable_info$variable_id

expression_data <-
  as.data.frame(expression_data)

variable_info <-
  as.data.frame(variable_info)

sample_info <-
  as.data.frame(sample_info)

library(massdataset)

metabolomics_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

masstools::setwd_project()
setwd("data/stability_analysis/metabolomics/data_preparation")

metabolomics_data <-
  metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  arrange(subject_id)

expression_data <- extract_expression_data(metabolomics_data)
sample_info <- extract_sample_info(metabolomics_data)
variable_info <- extract_variable_info(metabolomics_data)

# save(metabolomics_data, file = "metabolomics_data")
# save(expression_data, file = "expression_data")
# save(sample_info, file = "sample_info")
# save(variable_info, file = "variable_info")

metabolomics_data <- 
metabolomics_data %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  dplyr::select(1:3)

massdataset::export_mass_dataset(object = metabolomics_data,
                                 file_type = "xlsx", 
                                 path = ".")
