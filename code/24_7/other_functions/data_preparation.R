##
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/")

library(tidyverse)
rm(list = ls())

load("raw_data_from_box/MOmics_01.RData")

time_anno = read.csv("raw_data_from_box/TimeAnno.csv")
time_anno = time_anno %>%
  dplyr::select(-time, -date) %>%
  dplyr::rename(DT = date_time)
events = gather(time_anno,"Type","Details",5:6)
events = events[!is.na(events$Details),]

library(lubridate)

start = as_datetime("2019-04-29 03:00:00", tz="America/Los_Angeles")
events$DT = as.POSIXct(events$DT, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
events = events %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))



####lipidomics
lipids = DF247_lipids2[,c(1,2,4,5,10,13)]
lipids$SampleID = as.numeric(gsub("24_7_prepsample","",lipids$SampleID))
lipids = lipids[!is.na(lipids$SampleID),]
lipids$MolClass = "Lipid"
lipids <- 
  dplyr::rename(lipids, MolName = LipidSpecies, MolSubclass = LipidClass, Intensity = log_sample_nmol_per_g_concentration, DT = CollectionTime_format)
lipids <- 
  select(lipids, -collector)
lipids$SampleID = as.character(lipids$SampleID)

sample_info <-
  lipids %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                # subject_id = PID,
                DT = DT) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame()

variable_info <-
  lipids %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                mol_class = MolClass,
                mol_subclass = MolSubclass) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mol_name = as.character(mol_name))

expression_data <-
  lipids %>%
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

save(sample_info, file = "lipidomics/data_preparation/sample_info")
save(variable_info, file = "lipidomics/data_preparation/variable_info")
save(expression_data, file = "lipidomics/data_preparation/expression_data")



###### proteins - plasma
library(data.table)
proteins <- data.table(fread("raw_data_from_box/allproteins_overtime_annotated_Mike247plasma.csv"), sep='\t')
uniprot = read.csv("raw_data_from_box/uniprot_match.csv")
uniprot$Entry_name = gsub("_HUMAN","",uniprot$Entry_name)
proteinsL = gather(proteins[,c(1:313,321,325)], Sample, Intensity, 1:313)
proteinsL = left_join(proteinsL,select(uniprot,Entry,Entry_name), by=c("Sample"="Entry"))
proteinsL = subset(proteinsL,!is.na(Entry_name))
proteinsL = unite(data = proteinsL, col = "MolName", Sample, Entry_name)
proteinsL <-
  proteinsL %>%
  dplyr::rename(DT = date_time, SampleID = SampleIndex)
proteinsL$MolClass = "Protein"
proteinsL$MolSubclass = ""
proteinsL$SampleID = as.character(proteinsL$SampleID)
proteinsL$DT = as.POSIXct(proteinsL$DT, tz = "America/Los_Angeles", format = "%Y-%m-%d %H:%M:%S")

sample_info <-
  proteinsL %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                # subject_id = PID,
                DT = DT) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame()

variable_info <-
  proteinsL %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                mol_class = MolClass,
                mol_subclass = MolSubclass) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mol_name = as.character(mol_name))

expression_data <-
  proteinsL %>%
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
  dplyr::mutate(variable_id = paste("protein", 1:nrow(variable_info), sep = "_")) %>% 
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

variable_info <-
  variable_info %>% 
  dplyr::mutate(protein_name2 = 
                  stringr::str_split(mol_name, "_", n = 2) %>% 
                  purrr::map(.f = function(x)x[2]) %>% 
                  unlist())

variable_info <- 
variable_info %>% 
  dplyr::left_join(uniprot, by = c("protein_name2" = "Entry_name"))

save(sample_info, file = "proteomics/data_preparation/sample_info")
save(variable_info, file = "proteomics/data_preparation/variable_info")
save(expression_data, file = "proteomics/data_preparation/expression_data")


############## metabolomics----------------------------------------------------
met_data <-
  read.csv("raw_data_from_box/Microsampling_247_combined_all.csv")

met_data$name <- paste(met_data$Mode, met_data$Compounds_ID, sep = '_')

library(plyr)

met_data <- 
met_data %>% 
  plyr::dlply(.variables = .(Mode)) %>% 
  purrr::map(.f = function(x){
    x <- 
      x %>% 
      dplyr::distinct(name, .keep_all = TRUE)
    x 
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

### metID
rplc_pos_anno = read.csv(
  "metabolomics/metID_Microsampling_247/RPLC/pos/identification.table.new_pRPLC.csv"
) %>% 
  dplyr::distinct(name, .keep_all = TRUE)

rplc_neg_anno = read.csv(
  "metabolomics/metID_Microsampling_247/RPLC/neg/identification.table.new_nRPLC.csv"
) %>% 
  dplyr::distinct(name, .keep_all = TRUE)

hilic_pos_anno = read.csv(
  "metabolomics/metID_Microsampling_247/HILIC/pos/identification.table.new_pHILIC.csv"
) %>% 
  dplyr::distinct(name, .keep_all = TRUE)

hilic_neg_anno = read.csv(
  "metabolomics/metID_Microsampling_247/HILIC/neg/identification.table.new_nHILIC.csv"
) %>% 
  dplyr::distinct(name, .keep_all = TRUE)

rplc_pos_anno$mode <- "RPLC_pos"
rplc_neg_anno$mode <- "RPLC_neg"
hilic_pos_anno$mode <- "HILIC_pos"
hilic_neg_anno$mode <- "HILIC_neg"

met_tags = rbind(rplc_pos_anno,rplc_neg_anno,hilic_pos_anno,hilic_neg_anno)

met_tags$name <- 
  paste(met_tags$mode, met_tags$name, sep = "_")

expression_data <-
  met_data %>%
  dplyr::select(-c(mz, m.z, RT)) %>% 
  dplyr::left_join(met_tags, by = c("name"))

variable_info <- 
  expression_data[,-grep("X[0-9]{1,3}|QC", colnames(expression_data))]

expression_data <- 
  expression_data[,grep("X[0-9]{1,3}|QC", colnames(expression_data))]

variable_info <- 
  variable_info %>% 
  dplyr::mutate(name = paste(Mode, Compounds_ID, sep = "_")) %>% 
  dplyr::select(-c(X.1, X, Neutral.mass..Da., Charge, Retention.time..min., Identifications, Accepted.Compound.ID, Accepted.Description,
                   Score, Fragmentation.Score, Mass.Error..ppm., Isotope.Similarity, Retention.Time.Error..mins.,
                   Compound.Link, Mode, Compounds_ID)) %>% 
  dplyr::select(name, mz, rt, MS2.spectrum.name:Database, everything())
  

variable_info$name

rownames(expression_data) <- variable_info$name

variable_info <- 
  variable_info %>% 
  dplyr::rename(variable_id = name)

sample_info <- 
  data.frame(injection_order = colnames(expression_data)) %>% 
  dplyr::mutate(class = case_when(
    stringr::str_detect(injection_order, "QC") ~ "QC",
    TRUE ~ "Subject"
  ))

##injection order
inj_order <- readxl::read_xlsx("raw_data_from_box/M-Sampling-Brittany.xlsx")
inj_order <-
  dplyr::rename(inj_order, sample_id = "PrepIndex (SampleName)")

inj_order$AcqOrder <- paste0("X", inj_order$AcqOrder)

sample_info <-
  sample_info %>%
  dplyr::left_join(inj_order, by = c("injection_order" = "AcqOrder")) %>% 
  dplyr::filter(!is.na(sample_id))

expression_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$injection_order))

colnames(expression_data) == sample_info$injection_order

colnames(expression_data) = as.character(sample_info$sample_id)

sample_info <- 
  sample_info %>% 
  dplyr::select(sample_id, injection_order, everything())

dim(expression_data)
dim(sample_info)
dim(variable_info)

sample_info$injection_order <- 
  sample_info$injection_order %>% 
  stringr::str_replace("X", "") %>% 
  as.numeric()

sample_info <- 
  sample_info %>% 
  dplyr::mutate(
    class = case_when(
      stringr::str_detect(Filenames, "QC") ~ "QC",
      stringr::str_detect(Filenames, "blank") ~ "Blank",
      TRUE ~ "Subject"
    )
  )

save(expression_data, file = "metabolomics/data_preparation/peaks/expression_data")
save(sample_info, file = "metabolomics/data_preparation/peaks/sample_info")
save(variable_info, file = "metabolomics/data_preparation/peaks/variable_info")

###metabolites
variable_info <- 
  variable_info %>%
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2)

expression_data <-
  expression_data[match(variable_info$variable_id, rownames(expression_data)),]

rownames(expression_data) == variable_info$variable_id

save(expression_data, file = "metabolomics/data_preparation/metabolites/expression_data")
save(sample_info, file = "metabolomics/data_preparation/metabolites/sample_info")
save(variable_info, file = "metabolomics/data_preparation/metabolites/variable_info")


# cortisol 
MetaData=read.csv("raw_data_from_box/Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv") 

MetaData = MetaData[MetaData$Study == "F0" & MetaData$Original_SampleID %in% 1:97,]
MetaData$Name = paste("DH-",MetaData$Original_SampleID,sep="")

cort = read.csv("raw_data_from_box/RYAN CORTISOL 8-28-19.csv")
cort = cort[10:102,]

cort = left_join(cort,MetaData[,c("Name","CollectionTime","PrepIndex")])
cort$DT = as.POSIXct(cort$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
cort = 
  dplyr::rename(cort, SampleID=PrepIndex)

cortL = gather(cort,"Cytokine","MFI",5)
cortL$Intensity = log10(as.numeric(cortL$MFI))
cortL$MolClass = "CortisolEnzymatic"
cortL$MolSubclass = NA
cortL$MolName = "Cortisol"
cortL = select(cortL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)

#41-plex
cyto = read.csv("raw_data_from_box/24-7 omics-HumanLuminexMAG42plex-Mitra H41 Plex-1.csv")
cyto = cyto[10:105,]

cyto = left_join(cyto,MetaData[,c("Name","CollectionTime","PrepIndex")])
cyto$DT = as.POSIXct(cyto$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")

cytoL = gather(cyto,"Cytokine","MFI",5:49)
cytoL$Intensity = log10(as.numeric(cytoL$MFI))
cytoL$MolClass = "Cytokine_41Plex"
cytoL$MolSubclass = NA
cytoL = dplyr::rename(cytoL, SampleID = PrepIndex, MolName = Cytokine)
cytoL = select(cytoL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)


sample_info <- 
  cytoL %>% 
  dplyr::select(sample_id = SampleID, DT) %>% 
  dplyr::distinct(sample_id, .keep_all = TRUE)

variable_info <- 
  cytoL %>% 
  dplyr::select(mol_name = MolName, mol_class = MolClass, mol_subclass = MolSubclass) %>% 
  dplyr::distinct(mol_name, .keep_all = TRUE)

expression_data <- 
  cytoL %>% 
  dplyr::select(sample_id = SampleID, mol_name = MolName, Intensity) %>% 
  tidyr::pivot_wider(names_from = sample_id, values_from = Intensity) %>% 
  tibble::column_to_rownames(var = "mol_name")

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$mol_name  

variable_info <- 
  variable_info %>% 
  dplyr::mutate(variable_id = mol_name) %>% 
  dplyr::select(variable_id, everything())

save(expression_data, file = "cortisol/data_preparation/expression_data")
save(variable_info, file = "cortisol/data_preparation/variable_info")
save(sample_info, file = "cortisol/data_preparation/sample_info")


