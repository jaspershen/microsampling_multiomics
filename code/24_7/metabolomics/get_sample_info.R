##
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/")

library(tidyverse)
library(data.table)
rm(list = ls())

load("raw_data_from_box/MOmics_01.RData")

time_anno = read.csv("raw_data_from_box/TimeAnno.csv")
time_anno = time_anno %>%
  dplyr::select(-time, -date) %>%
  dplyr::rename(DT = date_time)
events = gather(time_anno,"Type","Details",5:6)
events = events[!is.na(events$Details),]

library(lubridate)

start = lubridate::as_datetime(x = "2019-04-29 03:00:00", tz = "America/Los_Angeles")
events$DT = as.POSIXct(events$DT, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
events = events %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))


############## metabolomics----------------------------------------------------
met_data <-
  data.table(fread(file = "raw_data_from_box/Microsampling_247_combined_all.csv"))
  # read.csv("raw_data_from_box/Microsampling_247_combined_all.csv")

sum(duplicated(met_data$Compounds_ID[met_data$Mode == "RPLC_pos"]))
sum(duplicated(met_data$Compounds_ID[met_data$Mode == "RPLC_neg"]))
sum(duplicated(met_data$Compounds_ID[met_data$Mode == "HILIC_pos"]))
sum(duplicated(met_data$Compounds_ID[met_data$Mode == "HILIC_neg"]))

met_data$name <- paste(met_data$Mode, met_data$Compounds, sep = '_')

library(plyr)

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


##remove the reduntant peaks in one mode
met_data <-
  met_data %>%
  plyr::dlply(.variables = .(name)) %>%
  purrr::map(.f = function(x){
    x <-
      x %>%
      dplyr::filter(Maximum.Abundance == max(Maximum.Abundance))
    x = x[1,,drop = FALSE]  
    x
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

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
    dplyr::select(
      -c(
        V1,
        X,
        Neutral.mass..Da.,
        Charge,
        Retention.time..min.,
        Identifications,
        Accepted.Compound.ID,
        Accepted.Description,
        Score,
        Fragmentation.Score,
        Mass.Error..ppm.,
        Isotope.Similarity,
        Retention.Time.Error..mins.,
        Compound.Link,
        Mode,
        Compounds_ID
      )
    ) %>%
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


###
load("all_omes_wear_food")
  data = 
  all_omes_wear_food %>% 
  dplyr::filter(MolClass == "Metabolite")

  sample_info2 <-
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
  
sample_info =
  sample_info %>% 
  dplyr::mutate(sample_id = as.character(sample_id)) %>% 
  dplyr::left_join(sample_info2, by = "sample_id")
    

sxtTools::setwd_project()
setwd("data/7_24_mike/metabolomics/data_preparation")
dir.create("peaks")
dir.create("metabolites")

dim(expression_data)
length(unique(variable_info$variable_id))

x = 
variable_info %>% 
  dplyr::filter(Metabolite_val != "0") %>% 
  dplyr::select(Level, Compound.name, Metabolite_val, MSMS) %>% 
  dplyr::filter(Level == "1")


save(sample_info, file = "sample_info")


# save(expression_data, file = "peaks/expression_data")
# save(sample_info, file = "peaks/sample_info")
# save(variable_info, file = "peaks/variable_info")
# 
# ###metabolites
# variable_info <- 
#   variable_info %>%
#   dplyr::filter(!is.na(Level)) %>% 
#   dplyr::filter(Level == 1 | Level == 2)
# 
# expression_data <-
#   expression_data[match(variable_info$variable_id, rownames(expression_data)),]
# 
# rownames(expression_data) == variable_info$variable_id
# 
# 
# # ###remove redundant name
# # library(plyr)
# # 
# # variable_info = 
# # variable_info %>% 
# #   plyr::dlply(.variables = .(Compound.name)) %>% 
# #   purrr::map(.f = function(x){
# #     if(nrow(x) == 1){
# #       return(x)
# #     }else{
# #       x = 
# #       x %>% 
# #         dplyr::filter(Level == min(Level)) %>% 
# #         dplyr::filter(SS == max(SS)) %>% 
# #         dplyr::distinct(Compound.name, .keep_all = TRUE)
# #         
# #     }
# #   }) %>% 
# #   do.call(rbind, .) %>% 
# #   as.data.frame()
# # 
# # expression_data =
# #   expression_data[match(variable_info$variable_id, rownames(expression_data)),]
# 
# save(expression_data, file = "metabolites/expression_data")
# save(sample_info, file = "metabolites/sample_info")
# save(variable_info, file = "metabolites/variable_info")

