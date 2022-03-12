##
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/proteomics/data_preparation/")

library(tidyverse)
library(data.table)
rm(list = ls())

proteins <- data.table(fread("allproteins_overtime_annotated_Mike247plasma.csv"), sep='\t')
uniprot = read.csv("uniprot_match.csv")
uniprot$Entry_name = gsub("_HUMAN","",uniprot$Entry_name)

###no duplicated proteins
sum(duplicated(colnames(proteins)))

##check duplicated samples
temp =
  proteins %>%
  dplyr::select(SampleID, SampleIndex) %>%
  dplyr::arrange(SampleIndex)

duplicated_sample = 
temp %>% 
  dplyr::group_by(SampleIndex) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::pull(SampleIndex)

temp =
  temp %>% 
  dplyr::filter(SampleIndex %in% duplicated_sample)

temp  

##so we need to remove the duplicated samples
proteins[proteins$SampleID %in% c(temp$SampleID[c(1,2)]),c(1:313)] %>% 
  t() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point()

proteins[proteins$SampleID %in% c(temp$SampleID[c(3,4)]),c(1:313)] %>% 
  t() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point()

proteins[proteins$SampleID %in% c(temp$SampleID[c(5,6)]),c(1:313)] %>% 
  t() %>% 
  as.data.frame() %>% 
  ggplot(aes(x = V1, y = V2)) +
  geom_point()

##for duplicated samples, use the mean value
temp
library(plyr)
proteins = 
proteins %>% 
  plyr::dlply(.variables = .(SampleIndex)) %>% 
  purrr::map(function(x){
    x = 
      x[,c(1:313, 321, 325)]
    x = 
    x %>% 
      dplyr::mutate_if(is.numeric, mean)
  
    x %>% 
      dplyr::distinct(SampleIndex, .keep_all = TRUE)
        
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

proteins$SampleIndex

proteinsL = tidyr::gather(data = proteins,
                          key = Sample,
                          value = Intensity,
                          1:313)


proteinsL = left_join(proteinsL,dplyr::select(uniprot,Entry,Entry_name, Protein_names, Gene_names), by=c("Sample"="Entry"))
proteinsL = subset(proteinsL,!is.na(Entry_name))
# proteinsL = unite(data = proteinsL, col = "MolName", Sample,Entry_name)
proteinsL = dplyr::mutate(.data = proteinsL, MolName = paste(Sample,Entry_name, sep = "_"))
proteinsL = proteinsL %>% dplyr::rename(DT = date_time, SampleID = SampleIndex)
proteinsL$MolClass = "Protein"
proteinsL$MolSubclass = ""
proteinsL$SampleID = as.character(proteinsL$SampleID)
proteinsL$DT = as.POSIXct(proteinsL$DT, tz = "America/Los_Angeles", format = "%Y-%m-%d %H:%M:%S")

head(proteinsL)
length(unique(proteinsL$SampleID))
length(unique(proteinsL$MolName))

proteinsL %>% 
  dplyr::filter(MolName == "B9A064_IGLL5") %>% 
  dplyr::pull(Intensity) %>% 
  plot()

data = proteinsL

sample_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(sample_id = SampleID,
                accurate_time = DT) %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = "Mike1") %>%
  dplyr::select(subject_id, sample_id, everything())

# sample_info %>% 
#   dplyr::mutate(day = )
# 
# 
# day = day
# time = tod
# hour = hour_of_day

variable_info <-
  data %>%
  dplyr::ungroup() %>%
  dplyr::select(mol_name = MolName,
                Entry_name,
                gene_name = Gene_names,
                protein_name = Protein_names,
                subclass = MolSubclass) %>%
  dplyr::distinct(.keep_all = TRUE) %>% 
  as.data.frame()

library(plyr)
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

dim(expression_data)

###order samples
sample_info = 
sample_info %>% 
  dplyr::arrange(accurate_time)

expression_data = 
  expression_data[,sample_info$sample_id]

variable_info$UNIPROT =
  variable_info$mol_name %>%
  stringr::str_split(pattern = "_") %>%
  purrr::map(function(x)x[1]) %>%
  unlist()

KEGG =
  clusterProfiler::bitr_kegg(
    geneID = as.character(variable_info$UNIPROT),
    fromType = "uniprot",
    toType = "kegg",
    drop = FALSE,
    organism = "hsa"
  ) %>%
  dplyr::distinct(uniprot, .keep_all = TRUE)

library(org.Hs.eg.db)

ENTREZID =
  clusterProfiler::bitr(
    geneID = as.character(variable_info$UNIPROT),
    fromType = "UNIPROT",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db,
    drop = TRUE
  ) %>%
  dplyr::distinct(UNIPROT, .keep_all = TRUE)


variable_info = 
variable_info %>% 
  dplyr::left_join(KEGG, by = c("UNIPROT" = "uniprot")) %>% 
  dplyr::left_join(ENTREZID, by = c("UNIPROT"))




sample_info$accurate_time = 
  sample_info$accurate_time %>% 
  as.character() %>% 
  lubridate::as_datetime(tz = "America/Los_Angeles")

save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")
save(expression_data, file = "expression_data")
















