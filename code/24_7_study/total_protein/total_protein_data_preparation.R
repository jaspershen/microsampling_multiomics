##
##steps
no_function()

masstools::setwd_project()
library(tidyverse)
setwd("data/24_7_study/total_protein/data_preparation/")

# metabolic protein panel
MetaData = read.csv("Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv")

MetaData = MetaData[MetaData$Study == "F0" &
                      MetaData$Original_SampleID %in% 1:97, ]
MetaData$Name = paste("DH-", MetaData$Original_SampleID, sep = "")

# total protein
tp = read.csv("27-4CytokineTotalProtein.csv")
tp = left_join(tp, MetaData[, c("Name", "CollectionTime", "PrepIndex")])
tp$DT = as.POSIXct(tp$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
tp = dplyr::rename(tp, SampleID = PrepIndex, Intensity = Total_protein)
tpL = tp # gather(tp,"Cytokine","MFI",5)
# tpL$Intensity = zscore(tpL$Intensity)
#tpL$Intensity = log10(as.numeric(tpL$MFI))
tpL$MolClass = "BCA"
tpL$MolSubclass = NA
tpL$MolName = "TotalProtein"
tpL = select(tpL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)

data = tpL

sample_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                accurate_time = DT) %>%
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
  dplyr::mutate(variable_id = paste("total_protein", 1:nrow(variable_info), sep = "_")) %>%
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

sample_info =
  sample_info %>%
  dplyr::arrange(accurate_time)

expression_data =
  expression_data[, sample_info$sample_id]

# save(sample_info, file = "sample_info")
# save(variable_info, file = "variable_info")
# save(expression_data, file = "expression_data")

write.csv(expression_data,
          "expression_data.csv",
          row.names = FALSE)
write.csv(sample_info, "sample_info.csv", row.names = FALSE)
write.csv(variable_info, "variable_info.csv", row.names = FALSE)
