##
no_function()

sxtTools::setwd_project()
setwd("data/7_24_mike/metabolomics/data_preparation/")

library(tidyverse)
library(data.table)
rm(list = ls())

load("sample_info")

sample_info1 = sample_info

load("../metabolite_annotation/expression_data")
load("../metabolite_annotation/variable_info")
load("../metabolite_annotation/sample_info")

annotation_from_kevin = readxl::read_xlsx("metabolites/annotation_result_from_kevin.xlsx")

sample_info2 = sample_info

dim(sample_info1)
dim(sample_info2)

table(sample_info1$class)
table(sample_info2$class)

sample_info2$sample_id = 
  sample_info2$sample_id %>% 
  stringr::str_replace_all("X", "")

setdiff(sample_info1$sample_id, sample_info2$sample_id)
setdiff(sample_info2$sample_id, sample_info1$sample_id)

sample_info = 
  sample_info1 

colnames(expression_data) = 
  colnames(expression_data) %>% 
  stringr::str_replace_all("X", "")

dim(sample_info)

dim(expression_data)

sum(colnames(expression_data) %in% sample_info$sample_id)

expression_data = 
  expression_data %>% 
  dplyr::select(sample_info$sample_id)

colnames(expression_data) == sample_info$sample_id

rownames(expression_data) == variable_info$variable_id

##median normalization
expression_data1 = 
expression_data %>% 
  purrr::map(function(x){
    x = x/mean(x)
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

sum(is.na(expression_data))

plot(density(log(expression_data$`1`)))
plot(density(log(expression_data1$`1`)))

plot(expression_data$`1`, expression_data1$`1`)

expression_data = expression_data1

rownames(expression_data) = variable_info$variable_id

sample_info = 
sample_info %>% 
  dplyr::arrange(accurate_time)

expression_data = 
  expression_data[,sample_info$sample_id]

###remove the features with zero percentage > 20
zero_percentage = 
apply(expression_data, 1, function(x){
  sum(x == 0)/ncol(expression_data)
})

remain_idx = which(zero_percentage < 0.2)

expression_data = 
  expression_data[remain_idx,]

variable_info = 
  variable_info[remain_idx,]

###add KEGG and HMDB ID to it
kegg_id1 = variable_info$KEGG.ID
hmdb_id1 = variable_info$HMDB.ID

idx = grep("C", hmdb_id1)
kegg_id1[idx]
hmdb_id1[idx]

variable_info$KEGG.ID[idx] = hmdb_id1[idx]
variable_info$HMDB.ID[idx] = kegg_id1[idx]

kegg_id1 = variable_info$KEGG.ID
hmdb_id1 = variable_info$HMDB.ID

idx = grep("HMDB", kegg_id1)
kegg_id1[idx]
hmdb_id1[idx]

variable_info$KEGG.ID[idx] = hmdb_id1[idx]
variable_info$HMDB.ID[idx] = kegg_id1[idx]

variable_info$KEGG.ID[variable_info$KEGG.ID == "NA"] = NA
variable_info$HMDB.ID[variable_info$HMDB.ID == "NA"] = NA

variable_info$HMDB.ID = 
  variable_info$HMDB.ID %>%
  purrr::map(function(x) {
    if (is.na(x)) {
      return(NA)
    }
    
    if (nchar(x) == 9) {
      x = stringr::str_replace(x, "HMDB", "HMDB00")
      return(x)
    }
    
    if (nchar(x) == 11) {
      return(x)
    }
    
    if (nchar(x) > 13) {
      x = stringr::str_split(x, "\\|")[[1]][1]
      return(x)
    }
    return(x)
  }) %>% 
  unlist()

variable_info$KEGG.ID = 
  variable_info$KEGG.ID %>%
  purrr::map(function(x) {
    if (is.na(x)) {
      return(NA)
    }

    if (nchar(x) > 10) {
      x = stringr::str_split(x, "\\|")[[1]][1]
      return(x)
    }
    return(x)
  }) %>% 
  unlist()

dim(variable_info)

variable_info$variable_id

annotation_from_kevin$variable_id =
  paste(annotation_from_kevin$Mode, annotation_from_kevin$Compound, sep = "_") %>% 
  stringr::str_replace("pos", "POS") %>% 
  stringr::str_replace("neg", "NEG") 


idx =  match(annotation_from_kevin$variable_id, variable_info$variable_id)

variable_info %>%
  dplyr::left_join(annotation_from_kevin[, c("variable_id", "Metabolite_val", "HMDB_val", "KEGG_val")],
                   by = "variable_id") %>% 
  dplyr::filter(!is.na(Metabolite_val)) %>% 
  dplyr::select(Compound.name, Metabolite_val, Database, Level) %>% 
  View()

variable_info =
  variable_info %>%
  dplyr::left_join(annotation_from_kevin[, c("variable_id", "Metabolite_val", "HMDB_val", "KEGG_val")],
                   by = "variable_id")

variable_info %>% 
  dplyr::filter(Compound.name == "Caffeine")

variable_info %>% 
  dplyr::filter(Metabolite_val == "Caffeine")

grep("Caffeine", variable_info$Metabolite_val, value = FALSE)
grep("Caffeine", variable_info$Compound.name, value = FALSE)

variable_info[c(9980, 23714, 23902),]

plot(sample_info$accurate_time, as.numeric(expression_data[23902,]), type = "b")

grep("benzenetriol", variable_info$Metabolite_val, value = FALSE)
grep("benzenetriol", variable_info$Compound.name, value = FALSE)

variable_info[c(19441,18834),]

plot(sample_info$accurate_time, as.numeric(expression_data[18834,]), type = "b")

grep("Hydroxyphenyllactic acid", variable_info$Metabolite_val, value = FALSE)
grep("Hydroxyphenyllactic acid", variable_info$Compound.name, value = FALSE)

variable_info[c(23874,24336,22502),]

plot(sample_info$accurate_time, as.numeric(expression_data[22502,]), type = "b")

save(expression_data, file = "peaks/expression_data")
save(sample_info, file = "peaks/sample_info")
save(variable_info, file = "peaks/variable_info")

####only remain metabolite
load("peaks/expression_data")
load("peaks/sample_info")
load("peaks/variable_info")

variable_info = 
  variable_info %>% 
  dplyr::filter(!is.na(Compound.name)) %>% 
  dplyr::filter(Level != 3)

library(plyr)

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

rownames(expression_data) == variable_info$variable_id

kegg_id1 = variable_info$KEGG.ID
hmdb_id1 = variable_info$HMDB.ID
sum(is.na(kegg_id1))
sum(is.na(hmdb_id1))

# kegg_id2 =
#   variable_info$Compound.name %>%
#   purrr::map(function(x){
#     metflow2::transID(
#       query = stringr::str_trim(x, side = "both"),
#       from = "Chemical Name",
#       to = "KEGG",
#       top = 1
#     ) %>%
#       dplyr::pull(KEGG)
#   }) %>%
#   unlist()
#
# hmdb_id2 =
#   variable_info$Compound.name %>%
#   purrr::map(function(x){
#     metflow2::transID(
#       query = stringr::str_trim(x, side = "both"),
#       from = "Chemical Name",
#       to = "Human Metabolome Database",
#       top = 1
#     ) %>%
#       dplyr::pull(`Human Metabolome Database`)
#   }) %>%
#   unlist()
# #
# save(hmdb_id2, file = "hmdb_id2")
# save(kegg_id2, file = "kegg_id2")

load("hmdb_id2")
load("kegg_id2")

hmdb_id =
  data.frame(hmdb_id1, hmdb_id2) %>%
  apply(1, function(x) {
    if (all(is.na(x))) {
      return(NA)
    }
    if (!is.na(x[1])) {
      return(x[1])
    }
    
    if (is.na(x[1])) {
      return(x[2])
    }
  })

kegg_id =
  data.frame(kegg_id1, kegg_id2) %>%
  apply(1, function(x) {
    if (all(is.na(x))) {
      return(NA)
    }
    if (!is.na(x[1])) {
      return(x[1])
    }
    
    if (is.na(x[1])) {
      return(x[2])
    }
  })

kegg_id[kegg_id == ""] = NA
hmdb_id[hmdb_id == ""] = NA

variable_info$KEGG.ID = kegg_id
variable_info$HMDB.ID = hmdb_id

idx1 = which(is.na(variable_info$KEGG.ID) & !is.na(variable_info$HMDB.ID))

variable_info$KEGG.ID[idx1]
variable_info$HMDB.ID[idx1]

miss_kegg =
  variable_info$HMDB.ID[idx1] %>%
  purrr::map(function(x){
    metflow2::transID(
      query = stringr::str_trim(x, side = "both"),
      from = "Human Metabolome Database",
      to = "KEGG",
      top = 1
    ) %>%
      dplyr::pull(KEGG)
  }) %>%
  unlist()

variable_info$KEGG.ID[idx1] = miss_kegg

idx2 = which(!is.na(variable_info$KEGG.ID) & is.na(variable_info$HMDB.ID))

variable_info$KEGG.ID[idx2]
variable_info$HMDB.ID[idx2]

miss_hmdb =
  variable_info$KEGG.ID[idx2] %>%
  purrr::map(function(x){
    metflow2::transID(
      query = stringr::str_trim(x, side = "both"),
      from = "KEGG",
      to = "Human Metabolome Database",
      top = 1
    ) %>%
      dplyr::pull(`Human Metabolome Database`)
  }) %>%
  unlist()

variable_info$HMDB.ID[idx1] = miss_hmdb

no_id = 
  variable_info %>% 
  dplyr::select(Compound.name, KEGG.ID, HMDB.ID, Database) %>% 
  dplyr::filter(is.na(KEGG.ID) | is.na(HMDB.ID))

write.csv(no_id, "metabolites/no_id.csv", row.names = FALSE)

no_id = read_csv("metabolites/no_id_manual.csv")

idx = 
match(no_id$Compound.name, variable_info$Compound.name)


variable_info$KEGG.ID[idx] = no_id$KEGG.ID
variable_info$HMDB.ID[idx] = no_id$HMDB.ID

save(expression_data, file = "metabolites/expression_data")
save(sample_info, file = "metabolites/sample_info")
save(variable_info, file = "metabolites/variable_info")

