##
no_function()

masstools::setwd_project()
rm(list = ls())
lipid_info <-
  read.table(
    "data/ensure_shake_study/lipidomics_data_analysis/DEG/Lipomat05.txt",
    header = TRUE,
    sep = "\t"
  )

setwd("data/24_7_study/lipidomics/")

library(tidyverse)

load("MOmics_2020-Apr-24.RData")
ls()
lipids = DF247_lipids2[, c(
  "LipidSpecies",
  "LipidClass",
  "log_sample_nmol_per_g_concentration",
  "SampleID",
  "CollectionTime_format",
  "collector"
)]
lipids$SampleID = as.numeric(gsub("24_7_prepsample", "", lipids$SampleID))
lipids = lipids[!is.na(lipids$SampleID), ]
lipids$MolClass = "Lipid"
lipids =
  dplyr::rename(
    lipids,
    MolName = LipidSpecies,
    MolSubclass = LipidClass,
    Intensity = log_sample_nmol_per_g_concentration,
    DT = CollectionTime_format
  )
lipids = dplyr::select(lipids, -collector)
lipids$SampleID = as.character(lipids$SampleID)

lipids %>%
  dplyr::filter(MolName == "CE.12.0.") %>%
  dplyr::pull(Intensity) %>%
  plot()

data <-
  lipids %>%
  dplyr::filter(MolClass == "Lipid")

sample_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                time = DT) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = "Mike1") %>%
  dplyr::select(subject_id, sample_id, everything())

variable_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
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
  dplyr::mutate(variable_id = paste("lipid", 1:nrow(variable_info), sep = "_")) %>%
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

expression_data[variable_info$mol_name == "CE.12.0.", ] %>%
  as.numeric() %>%
  plot()

setwd("data_preparation/")

variable_info
lipid_info$Lipid_ID =
  stringr::str_trim(lipid_info$Lipid_ID, side = "both")

variable_info =
  variable_info %>%
  dplyr::left_join(
    lipid_info %>% dplyr::distinct(Lipid_ID, .keep_all = TRUE),
    by = c("mol_name" = "Lipid_ID")
  )

dim(variable_info)
dim(expression_data)

rownames(expression_data) == variable_info$variable_id

sample_info =
  sample_info %>%
  dplyr::arrange(time)

expression_data = expression_data[, sample_info$sample_id]

sample_info =
  sample_info %>%
  dplyr::rename(accurate_time = time)

# save(sample_info, file = "sample_info")
# save(variable_info, file = "variable_info")
# save(expression_data, file = "expression_data")

write.csv(expression_data, "expression_data.csv", row.names = FALSE)
write.csv(sample_info, "sample_info.csv", row.names = FALSE)
write.csv(variable_info, "variable_info.csv", row.names = FALSE)
