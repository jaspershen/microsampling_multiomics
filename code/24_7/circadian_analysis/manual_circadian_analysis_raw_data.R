no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

####load data all the omics data
{
  ###
  load(here::here("data/7_24_mike/metabolomics/data_preparation/peaks/expression_data"))
  load(here::here("data/7_24_mike/metabolomics/data_preparation/peaks/sample_info"))
  load(here::here("data/7_24_mike/metabolomics/data_preparation/peaks/variable_info"))
}

dir.create("data/7_24_mike/circadian_analysis/manual_raw_data", recursive = TRUE)
setwd("data/7_24_mike/circadian_analysis/manual_raw_data")

######manual check
load(here::here("data/7_24_mike/summary_info/day_night_df"))

day_night_df = 
  day_night_df %>% 
  dplyr::filter(!as.character(day) %in% c("2019-05-07"))

###read metadata
metdata =
  readr::read_csv(here::here("data/7_24_mike/raw_data_from_box/sample_registration.csv"))
# metdata =
#   readxl::read_xlsx(here::here("data/7_24_mike/raw_data_from_box/sample registration.xlsx"))

milk_time =
  metdata %>%
  dplyr::filter(!is.na(food)) %>%
  dplyr::filter(stringr::str_detect(food, "milk|Milk")) %>%
  pull(date_time) %>% 
  lubridate::as_datetime(tz = "America/Los_Angeles")

milk_time[8] = milk_time[8] - 24*60*60

load(here::here("data/7_24_mike/all_omes"))

####metabolite name

###Salicylic acid
##Salicylic acid
##Salicylic acid
metabolite_name = "Salicylic acid"
grep("Salicylic acid$",variable_info$Compound.name, value = TRUE)
grep("Salicylic acid$",variable_info$Compound.name, value = FALSE)

grep("Salicylic acid$",all_omes$MolName, value = TRUE)
grep("Salicylic acid$",all_omes$MolName, value = FALSE)

temp_data =
  all_omes[grep("Salicylic acid$",all_omes$MolName, value = FALSE), ]

temp_data <-
  temp_data %>%
  dplyr::select(SampleID, Intensity) %>%
  dplyr::left_join(sample_info, by = c("SampleID" = "sample_id"))

plot2 =
  time_plot(
    x = temp_data$Intensity,
    time = temp_data$accurate_time,
    day_night_df = day_night_df,
    x_name = "Salicylic acid",
    add_point = TRUE,
    x_color = class_color["metabolomics"],
    y_name = "Scaled intensity"
  ) +
  geom_vline(xintercept = milk_time, color = "red")

plot2

ggsave(plot2, filename = paste0(metabolite_name, ".pdf"), width = 21, height = 5)

temp_data %>% 
  ggplot(aes(time, Intensity)) +
  geom_line(aes(group = day, color = as.character(day)))








#####coffee

coffee_time =
  metdata %>%
  dplyr::filter(!is.na(food)) %>%
  dplyr::filter(stringr::str_detect(food, "coffee|Coffee")) %>%
  pull(date_time) %>% 
  lubridate::as_datetime(tz = "America/Los_Angeles")

load(here::here("data/7_24_mike/all_omes"))

####metabolite name

###Coffee
metabolite_name = "Caffeine"
grep("Caffeine$",variable_info$Compound.name, value = TRUE)
grep("Caffeine$",variable_info$Compound.name, value = FALSE)

grep("Caffeine$",all_omes$MolName, value = TRUE)
grep("Caffeine$",all_omes$MolName, value = FALSE)

temp_data =
  all_omes[grep("Caffeine$",all_omes$MolName, value = FALSE), ]

temp_data <-
  temp_data %>%
  dplyr::select(SampleID, Intensity) %>%
  dplyr::left_join(sample_info, by = c("SampleID" = "sample_id"))

plot2 =
  time_plot(
    x = temp_data$Intensity,
    time = temp_data$accurate_time,
    day_night_df = day_night_df,
    x_name = "Caffeine",
    add_point = TRUE,
    x_color = class_color["metabolomics"],
    y_name = "Scaled intensity"
  ) +
  geom_vline(xintercept = coffee_time, color = "red")

plot2

ggsave(plot2, filename = paste0(metabolite_name, ".pdf"), width = 21, height = 5)

temp_data %>% 
  ggplot(aes(time, Intensity)) +
  geom_line(aes(group = day, color = as.character(day)))



######annotation confirm
load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/POS/NCE25/result_hilic_pos25"))
load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/POS/NCE35/result_hilic_pos35"))

load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/NEG/NCE25/result_hilic_neg25"))
load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/NEG/NCE35/result_hilic_neg35"))


load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/RPLC/POS/NCE25/result_rplc_pos25"))
load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/RPLC/POS/NCE50/result_rplc_pos50"))

load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/RPLC/NEG/NCE25/result_rplc_neg25"))
load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/RPLC/NEG/NCE50/result_rplc_neg50"))


load(here::here("data/7_24_mike/metabolomics/orbitrapDatabase0.0.1"))
load(here::here("data/7_24_mike/metabolomics/massbankDatabase0.0.2"))
load(here::here("data/7_24_mike/metabolomics/hmdbDatabase0.0.2"))
load(here::here("data/7_24_mike/metabolomics/monaDatabase0.0.2"))
load(here::here("data/7_24_mike/metabolomics/metlinDatabase0.0.2"))

#####Salicylic acid
grep("Salicylic acid$",variable_info$Compound.name, value = TRUE)
grep("Salicylic acid$",variable_info$Compound.name, value = FALSE)
variable_info[grep("Salicylic acid$",variable_info$Compound.name, value = FALSE),]

library(metid)

ms2plot(object = result_rplc_pos50[["hmdbDatabase0.0.2"]], 
        database = hmdbDatabase0.0.2, 
        which.peak = "0.53_139.0391m/z")

ms2plot(object = result_rplc_pos25[["hmdbDatabase0.0.2"]], 
        database = hmdbDatabase0.0.2, 
        which.peak = "0.53_139.0391m/z")

ms2plot(object = result_hilic_pos25[["monaDatabase0.0.2"]], 
        database = monaDatabase0.0.2, 
        which.peak = "4.66_121.0285m/z")

ms2plot(object = result_hilic_pos35[["monaDatabase0.0.2"]], 
        database = monaDatabase0.0.2, 
        which.peak = "4.66_121.0285m/z")


#####1,2,3-benzenetriol sulfate
grep("1,2,3-benzenetriol sulfate$",variable_info$Compound.name, value = TRUE)
grep("1,2,3-benzenetriol sulfate$",variable_info$Compound.name, value = FALSE)
variable_info[grep("1,2,3-benzenetriol sulfate$",variable_info$Compound.name, value = FALSE),]

library(metid)

ms2plot(object = result_hilic_neg25[["msDatabase_hilic0.0.2"]], 
        database = msDatabase_hilic0.0.2, 
        which.peak = "0.95_204.9812m/z")

ms2plot(object = result_hilic_neg35[["msDatabase_hilic0.0.2"]], 
        database = msDatabase_hilic0.0.2, 
        which.peak = "0.95_204.9812m/z")






#####Hydroxyphenyllactic acid
grep("Hydroxyphenyllactic acid$",variable_info$Compound.name, value = TRUE)
grep("Hydroxyphenyllactic acid$",variable_info$Compound.name, value = FALSE)
variable_info[grep("Hydroxyphenyllactic acid$",variable_info$Compound.name, value = FALSE),]

library(metid)

ms2plot(object = result_hilic_neg25[["msDatabase_hilic0.0.2"]], 
        database = msDatabase_hilic0.0.2, 
        which.peak = "6.12_182.0579n")

ms2plot(object = result_hilic_neg35[["msDatabase_hilic0.0.2"]], 
        database = msDatabase_hilic0.0.2, 
        which.peak = "6.12_182.0579n")

