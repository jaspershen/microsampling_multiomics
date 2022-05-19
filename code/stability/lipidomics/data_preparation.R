no_function()

masstools::setwd_project()
setwd("data/stability/lipidomics/")
rm(list = ls())
library(tidyverse)

data1 <- readxl::read_xlsx("Results - Lipidyzer Batch - 20210805135007.xlsx", sheet = 1)
data2 <- readxl::read_xlsx("Results - Lipidyzer Batch - 20210809141328.xlsx", sheet = 1)

data1$Name
data2$Name

data1 <- 
  data1 %>% 
  tibble::column_to_rownames(var = "Name") %>% 
  t() %>% 
  as.data.frame()

data2 <- 
  data2 %>% 
  tibble::column_to_rownames(var = "Name") %>% 
  t() %>% 
  as.data.frame()

colnames(data1)
colnames(data2)

intersect_id <- 
  intersect(rownames(data1),
            rownames(data2))

data <- cbind(data1[intersect_id,],
              data2[intersect_id,])

variable_info <- 
  data.frame(
    variable_id = paste("lipid", 1:nrow(data), sep = "_"),
    lipid_name = rownames(data))

rownames(data) <- variable_info$variable_id

expression_data <- 
  data %>% 
  dplyr::select(-contains("Q")) %>% 
  dplyr::select(-contains("B"))

sample_info <-
  data.frame(sample_id = colnames(expression_data)) %>% 
  dplyr::mutate(subject_id = stringr::str_replace(sample_id, "P|M", "")) %>% 
  dplyr::mutate(class = stringr::str_extract(sample_id, "P|M")) 

expression_data <-
  as.data.frame(expression_data)

variable_info <-
  as.data.frame(variable_info)

sample_info <-
  as.data.frame(sample_info)

library(massdataset)

lipidomics_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

masstools::setwd_project()
dir.create("data/stability/lipidomics/data_preparation", recursive = TRUE)
setwd("data/stability/lipidomics/data_preparation")

lipidomics_data <- 
lipidomics_data %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  arrange(subject_id)

expression_data <- extract_expression_data(lipidomics_data)
sample_info <- extract_expression_data(lipidomics_data)
variable_info <- extract_expression_data(lipidomics_data)

save(lipidomics_data, file = "lipidomics_data")
save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
