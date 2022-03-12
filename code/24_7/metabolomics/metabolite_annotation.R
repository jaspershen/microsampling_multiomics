no_source()


###HILIC
###positive mode
sxtTools::setwd_project()
setwd("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/POS/")
library(tidyverse)
library(data.table)
library(metID)

peak_table = readr::read_csv("PQI_Microsampling_HILIC_pos.csv")

# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)

# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_hilic0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# result_hilic_pos25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE25/", pattern = "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
#   path = "NCE25/"
# )
# 
# result_hilic_pos35 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir(path = "NCE35/", pattern = "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE35/"
# )
# 
# save(result_hilic_pos25, file = "NCE25/result_hilic_pos25")
# save(result_hilic_pos35, file = "NCE35/result_hilic_pos35")

load("NCE25/result_hilic_pos25")
load("NCE35/result_hilic_pos35")

annotation_table =
  metID::get_identification_table_all(result_hilic_pos25,
                                      result_hilic_pos35,
                                      candidate.num = 1)
head(annotation_table)

plot =
  summary_annotation_result(object = annotation_table)
plot
ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)

peak_table = readr::read_csv("PQI_Microsampling_HILIC_pos.csv")

annotation_table$name == peak_table$Compound

annotation_table =
  cbind(annotation_table, peak_table)

write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)


###HILIC
###negative mode
sxtTools::setwd_project()
setwd("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/NEG/")
library(tidyverse)
library(data.table)
library(metID)

# peak_table = readr::read_csv("PQI_Microsampling_HILIC_neg.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_hilic0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# 
# result_hilic_neg25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE25/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
#   path = "NCE25/"
# )
# 
# result_hilic_neg35 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE35/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE35/"
# )
# 
# 
# save(result_hilic_neg25, file = "NCE25/result_hilic_neg25")
# save(result_hilic_neg35, file = "NCE35/result_hilic_neg35")

load("NCE25/result_hilic_neg25")
load("NCE35/result_hilic_neg35")

annotation_table =
  metID::get_identification_table_all(result_hilic_neg25,
                                      result_hilic_neg35,
                                      candidate.num = 1)
head(annotation_table)

plot =
  summary_annotation_result(object = annotation_table)
plot
ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)

peak_table = readr::read_csv("PQI_Microsampling_HILIC_neg.csv")

annotation_table$name == peak_table$Compound

annotation_table =
  cbind(annotation_table, peak_table)

write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)










######################-----------------------------------------------------------
###RPLC
###positive mode
sxtTools::setwd_project()
setwd("data/7_24_mike/metabolomics/metabolite_annotation/RPLC/POS/")
library(tidyverse)
library(data.table)
library(metID)

# peak_table = readr::read_csv("PQI_Microsampling_RPLC_pos.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_rplc0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# result_rplc_pos25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE25/", pattern = "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
#   path = "NCE25/"
# )
# 
# 
# result_rplc_pos50 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir(path = "NCE50/", pattern = "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE50/"
# )
# 
# save(result_rplc_pos25, file = "NCE25/result_rplc_pos25")
# save(result_rplc_pos50, file = "NCE50/result_rplc_pos50")

load("NCE25/result_rplc_pos25")
load("NCE50/result_rplc_pos50")


annotation_table =
  metID::get_identification_table_all(result_rplc_pos25,
                                      result_rplc_pos50,
                                      candidate.num = 1)
head(annotation_table)

plot =
  summary_annotation_result(object = annotation_table)
plot

ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)

peak_table = readr::read_csv("PQI_Microsampling_RPLC_pos.csv")

annotation_table$name == peak_table$Compound

annotation_table =
  cbind(annotation_table, peak_table)

write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)


###RPLC
###negative mode
sxtTools::setwd_project()
setwd("data/7_24_mike/metabolomics/metabolite_annotation/RPLC/NEG/")
library(tidyverse)
library(data.table)
library(metID)

# peak_table = readr::read_csv("PQI_Microsampling_RPLC_neg.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_rplc0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# result_rplc_neg25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE25/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7, param8),
#   path = "NCE25/"
# )
# 
# result_rplc_neg50 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE50/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE50/"
# )
# 
# save(result_rplc_neg25, file = "NCE25/result_rplc_neg25")
# save(result_rplc_neg50, file = "NCE50/result_rplc_neg50")

load("NCE25/result_rplc_neg25")
load("NCE50/result_rplc_neg50")

annotation_table =
  metID::get_identification_table_all(result_rplc_neg25,
                                      result_rplc_neg50,
                                      candidate.num = 1)
head(annotation_table)

plot =
  summary_annotation_result(object = annotation_table)

plot

ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)

peak_table = readr::read_csv("PQI_Microsampling_RPLC_neg.csv")

annotation_table$name == peak_table$Compound

annotation_table =
  cbind(annotation_table, peak_table)

write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)












###combine them together
##rename
sxtTools::setwd_project()
setwd("data/7_24_mike/metabolomics/metabolite_annotation/")
library(tidyverse)
library(data.table)
library(metID)

hilic_pos_data = data.frame(fread("HILIC/POS/annotation_table.csv"))
hilic_neg_data = data.frame(fread("HILIC/NEG/annotation_table.csv"))

rplc_pos_data = data.frame(fread("RPLC/POS/annotation_table.csv"))
rplc_neg_data = data.frame(fread("RPLC/NEG/annotation_table.csv"))

hilic_pos_data$name = paste("HILIC_POS", hilic_pos_data$name, sep = "_")
hilic_neg_data$name = paste("HILIC_NEG", hilic_neg_data$name, sep = "_")

rplc_pos_data$name = paste("RPLC_POS", rplc_pos_data$name, sep = "_")
rplc_neg_data$name = paste("RPLC_NEG", rplc_neg_data$name, sep = "_")

sum(colnames(hilic_pos_data) == colnames(hilic_neg_data))

setdiff(colnames(hilic_pos_data), colnames(hilic_neg_data))
setdiff(colnames(hilic_neg_data), colnames(hilic_pos_data))
colnames(hilic_neg_data)[colnames(hilic_neg_data) == "QC5_bad"] = "QC5"

hilic_data = rbind(hilic_pos_data, hilic_neg_data)

sum(colnames(rplc_pos_data) == colnames(rplc_neg_data))

setdiff(colnames(rplc_pos_data), colnames(rplc_neg_data))
setdiff(colnames(rplc_neg_data), colnames(rplc_pos_data))
colnames(rplc_neg_data)[colnames(rplc_neg_data) == "X46_190703124001"] = "X46"

rplc_data = rbind(rplc_pos_data, rplc_neg_data)

all_name = unique(c(
  colnames(rplc_data),
  colnames(hilic_data)
))

setdiff(all_name, colnames(hilic_data))
setdiff(all_name, colnames(rplc_data))

hilic_data = 
  hilic_data %>% 
  dplyr::select(-c(blk3, blk4, QC13, QC14, QC15))

sum(colnames(hilic_data) == colnames(rplc_data))

dim(hilic_data)
dim(rplc_data)

expression_data = 
  rbind(rplc_data,
        hilic_data)

variable_info = 
  expression_data %>% 
  dplyr::select(c(name:`Minimum.CV.`, `Accepted.Compound.ID`:`Compound.Link`))

variable_info =   
  variable_info %>% 
  dplyr::rename(variable_id = name)

variable_info$variable_id
unique(variable_info$variable_id)

expression_data = 
  expression_data %>% 
  dplyr::select(-c(name:`Minimum.CV.`, `Accepted.Compound.ID`:`Compound.Link`))

sample_info = data.frame(sample_id = colnames(expression_data))

sample_info = 
  sample_info %>% 
  dplyr::mutate(class = case_when(
    stringr::str_detect(sample_id, "QC") ~ "QC",
    stringr::str_detect(sample_id, "dln") ~ "QC_DL",
    stringr::str_detect(sample_id, "blank") ~ "Blank",
    TRUE ~ "Subject"
  ))

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) <- variable_info$variable_id

variable_info

save(expression_data, file = "expression_data")
save(sample_info, file = "sample_info")
save(variable_info, file = "variable_info")


####only remain metabolite
variable_info = 
  variable_info %>% 
  dplyr::filter(!is.na(Compound.name)) %>% 
  dplyr::filter(Level != 3)

variable_info = 
  variable_info %>% 
  plyr::dlply(.variables = .(Compound.name)) %>% 
  purrr::map(function(x){
    x = 
      x %>% 
      dplyr::filter(Level == min(Level))
    
    if(nrow(x) == 1){
      return(x)
    }
    
    x = 
      x %>% 
      dplyr::filter(SS == max(SS))
    x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(expression_data)

expression_data = 
  expression_data[variable_info$variable_id,]

rownames(expression_data) = variable_info$variable_id

dir.create("data_preparation/metabolite")
save(expression_data, file = "data_preparation/metabolite/expression_data")
save(sample_info, file = "data_preparation/metabolite/sample_info")
save(variable_info, file = "data_preparation/metabolite/variable_info")


dim(variable_info)
table(variable_info$Level)

expression_data[grep("hilic_neg", rownames(expression_data)),"QC7"]



