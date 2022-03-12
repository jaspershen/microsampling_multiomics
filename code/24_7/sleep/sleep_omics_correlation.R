#' ---
#' title: "sleep omics correlation"
#' author: 
#'   - name: "Xiaotao Shen" 
#'     url: https://www.shenxt.info/
#'     affiliation: Stanford School of Medicine
#' date: "`r Sys.Date()`"
#' site: distill::distill_website
#' output: 
#'   distill::distill_article:
#'     code_folding: false
#' ---

no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source(here::here("R/tools.R"))
source(here::here("R/modified_dtw.R"))
source(here::here("R/lagged_correlation.R"))

load(here::here("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info"))

metabolomics_variable_info_raw = variable_info

{
  ####load data
  ###sleep
  load("data/7_24_mike/sleep/data_preparation/sample_info")
  load("data/7_24_mike/sleep/data_preparation/variable_info")
  load("data/7_24_mike/sleep/data_preparation/expression_data")
  
  sleep_expression_data = expression_data
  sleep_sample_info = sample_info
  sleep_variable_info = variable_info
  
  load("data/7_24_mike/summary_info/day_night_df")
  
  ####this is for the day night time
  day_night_df =
    day_night_df %>%
    dplyr::mutate(
      start_time = as.POSIXct(hms::as_hms(start)),
      end_time = as.POSIXct(hms::as_hms(end)),
      week = format(day, "%a")
    ) %>% 
    dplyr::mutate(week = paste(
      week,
      lubridate::month(day),
      lubridate::day(day),
      sep = "-"
    )) %>% 
    dplyr::mutate(week = factor(week, unique(week)))  
}


# load("data/7_24_mike/combine_omics/data_preparation/expression_data")
# load("data/7_24_mike/combine_omics/data_preparation/sample_info")
# load("data/7_24_mike/combine_omics/data_preparation/variable_info")

####load omics data
{
  ###metabolic panel data
  {
    ###metabolic_panel
    load("data/7_24_mike/metabolic_panel/data_preparation/sample_info")
    load("data/7_24_mike/metabolic_panel/data_preparation/variable_info")
    load("data/7_24_mike/metabolic_panel/data_preparation/expression_data")
    metabolic_panel_sample_info = sample_info
    metabolic_panel_variable_info = variable_info
    metabolic_panel_expression_data = expression_data
    
    remain_idx = 
      which(
        apply(expression_data, 1, function(x){
          sum(x == 0)/ncol(expression_data)
        }) < 0.5
      )
    
    metabolic_panel_variable_info = 
      metabolic_panel_variable_info[remain_idx,]
    
    metabolic_panel_expression_data = 
      metabolic_panel_expression_data[metabolic_panel_variable_info$variable_id,]
    
    ###remove CHEX
    metabolic_panel_variable_info = 
      metabolic_panel_variable_info %>% 
      dplyr::filter(!stringr::str_detect(mol_name, "CHEX"))
    
    metabolic_panel_expression_data = 
      metabolic_panel_expression_data[metabolic_panel_variable_info$variable_id,]
    
    metabolic_panel_variable_info$mol_name
    
    load("data/7_24_mike/summary_info/day_night_df")
    
    ####this is for the day night time
    day_night_df =
      day_night_df %>%
      dplyr::mutate(
        start_time = as.POSIXct(hms::as_hms(start)),
        end_time = as.POSIXct(hms::as_hms(end)),
        week = format(day, "%a")
      ) %>% 
      dplyr::mutate(week = paste(
        week,
        lubridate::month(day),
        lubridate::day(day),
        sep = "-"
      )) %>% 
      dplyr::mutate(week = factor(week, unique(week)))
    
  }
  
  
  ####cytokine data
  {
    ###cytokine
    load("data/7_24_mike/cytokine/data_preparation/sample_info")
    load("data/7_24_mike/cytokine/data_preparation/variable_info")
    load("data/7_24_mike/cytokine/data_preparation/expression_data")
    cytokine_sample_info = sample_info
    cytokine_variable_info = variable_info
    cytokine_expression_data = expression_data
    
    remain_idx = 
      which(
        apply(cytokine_expression_data, 1, function(x){
          sum(x == 0)/ncol(cytokine_expression_data)
        }) < 0.5
      )
    
    cytokine_variable_info = 
      cytokine_variable_info[remain_idx,]
    
    cytokine_expression_data = 
      cytokine_expression_data[cytokine_variable_info$variable_id,]
    
    
    ###remove the controls
    remove_idx = grep("CHEX", cytokine_variable_info$mol_name)
    if(length(remove_idx) > 0){
      cytokine_variable_info = 
        cytokine_variable_info[-remove_idx,]
      cytokine_expression_data = 
        cytokine_expression_data[-remove_idx,]
    }
  }
  
  ##lipids
  {
    ###lipidomics
    load("data/7_24_mike/lipidomics/data_preparation/sample_info")
    load("data/7_24_mike/lipidomics/data_preparation/variable_info")
    load("data/7_24_mike/lipidomics/data_preparation/expression_data")
    lipidomics_sample_info = sample_info
    lipidomics_variable_info = variable_info
    lipidomics_expression_data = expression_data
    
    
    remain_idx = 
      which(
        apply(lipidomics_expression_data, 1, function(x){
          sum(x == 0)/ncol(lipidomics_expression_data)
        }) < 0.5
      )
    
    lipidomics_variable_info = 
      lipidomics_variable_info[remain_idx,]
    
    lipidomics_expression_data = 
      lipidomics_expression_data[lipidomics_variable_info$variable_id,]
    
  }
  
  
  ###metabolomics
  {
    load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
    load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
    load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
    metabolomics_sample_info = sample_info
    metabolomics_variable_info = variable_info
    metabolomics_expression_data = expression_data
    
    ###remove wrong metabolites
    metabolomics_variable_info = 
      metabolomics_variable_info %>% 
      dplyr::filter(Database %in% c("hmdbDatabase0.0.2", 
                                    "metlinDatabase0.0.2",
                                    "msDatabase_hilic0.0.2", 
                                    "msDatabase_rplc0.0.2",
                                    "nistDatabase0.0.2",
                                    "orbitrapDatabase0.0.1"))
    
    metabolomics_expression_data = 
      metabolomics_expression_data[metabolomics_variable_info$variable_id,]
    
    ###remove QC samples
    metabolomics_sample_info = 
      metabolomics_sample_info %>% 
      dplyr::filter(!is.na(accurate_time))
    
    metabolomics_expression_data = 
      metabolomics_expression_data[,metabolomics_sample_info$sample_id]  
    
  }
  
  
  #####proteomics
  {
    ###proteomics
    load("data/7_24_mike/proteomics/data_preparation/sample_info")
    load("data/7_24_mike/proteomics/data_preparation/variable_info")
    load("data/7_24_mike/proteomics/data_preparation/expression_data")
    proteomics_sample_info = sample_info
    proteomics_variable_info = variable_info
    proteomics_expression_data = expression_data
    
    load("data/7_24_mike/summary_info/day_night_df")
    
    remain_idx = 
      which(
        apply(proteomics_expression_data, 1, function(x){
          sum(x == 0)/ncol(proteomics_expression_data)
        }) < 0.5
      )
    
    proteomics_variable_info = 
      proteomics_variable_info[remain_idx,]
    
    proteomics_expression_data = 
      proteomics_expression_data[proteomics_variable_info$variable_id,]
  }
  
  
  #####cortisol
  {
    ###cortisol
    load("data/7_24_mike/cortisol/data_preparation/sample_info")
    load("data/7_24_mike/cortisol/data_preparation/variable_info")
    load("data/7_24_mike/cortisol/data_preparation/expression_data")
    cortisol_sample_info = sample_info
    cortisol_variable_info = variable_info
    cortisol_expression_data = expression_data
    
    load("data/7_24_mike/summary_info/day_night_df")
    
    remain_idx = 
      which(
        apply(cortisol_expression_data, 1, function(x){
          sum(x == 0)/ncol(cortisol_expression_data)
        }) < 0.5
      )
    
    cortisol_variable_info = 
      cortisol_variable_info[remain_idx,]
    
    cortisol_expression_data = 
      cortisol_expression_data[cortisol_variable_info$variable_id,]
  }
  
  
  
  #####total_protein
  {
    ###total_protein
    load("data/7_24_mike/total_protein/data_preparation/sample_info")
    load("data/7_24_mike/total_protein/data_preparation/variable_info")
    load("data/7_24_mike/total_protein/data_preparation/expression_data")
    total_protein_sample_info = sample_info
    total_protein_variable_info = variable_info
    total_protein_expression_data = expression_data
    
    load("data/7_24_mike/summary_info/day_night_df")
    
    remain_idx = 
      which(
        apply(total_protein_expression_data, 1, function(x){
          sum(x == 0)/ncol(total_protein_expression_data)
        }) < 0.5
      )
    
    total_protein_variable_info = 
      total_protein_variable_info[remain_idx,]
    
    total_protein_expression_data = 
      total_protein_expression_data[total_protein_variable_info$variable_id,]
  }
  
  
  ####combine all omics together
  intersect_time = 
    Reduce(f = intersect, x = list(
      as.character(lipidomics_sample_info$accurate_time),
      as.character(metabolomics_sample_info$accurate_time),
      as.character(cytokine_sample_info$accurate_time),
      as.character(total_protein_sample_info$accurate_time),
      as.character(cortisol_sample_info$accurate_time),
      as.character(metabolic_panel_sample_info$accurate_time),
      as.character(proteomics_sample_info$accurate_time)
    ))
  
  lipidomics_sample_info =
    lipidomics_sample_info %>%
    dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
    dplyr::arrange(accurate_time)
  
  lipidomics_expression_data = lipidomics_expression_data[, lipidomics_sample_info$sample_id]
  
  metabolomics_sample_info =
    metabolomics_sample_info %>%
    dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
    dplyr::arrange(accurate_time)
  
  metabolomics_expression_data = metabolomics_expression_data[, metabolomics_sample_info$sample_id]
  
  cytokine_sample_info =
    cytokine_sample_info %>%
    dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
    dplyr::arrange(accurate_time)
  
  cytokine_expression_data = cytokine_expression_data[, cytokine_sample_info$sample_id]
  
  total_protein_sample_info =
    total_protein_sample_info %>%
    dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
    dplyr::arrange(accurate_time)
  
  total_protein_expression_data = total_protein_expression_data[, total_protein_sample_info$sample_id]
  
  cortisol_sample_info =
    cortisol_sample_info %>%
    dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
    dplyr::arrange(accurate_time)
  
  cortisol_expression_data = cortisol_expression_data[, cortisol_sample_info$sample_id]
  
  metabolic_panel_sample_info =
    metabolic_panel_sample_info %>%
    dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
    dplyr::arrange(accurate_time)
  
  metabolic_panel_expression_data = metabolic_panel_expression_data[, metabolic_panel_sample_info$sample_id]
  
  proteomics_sample_info =
    proteomics_sample_info %>%
    dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
    dplyr::arrange(accurate_time)
  
  proteomics_expression_data = proteomics_expression_data[, proteomics_sample_info$sample_id]
  
  lipidomics_sample_info$accurate_time == metabolic_panel_sample_info$accurate_time
  lipidomics_sample_info$accurate_time == metabolomics_sample_info$accurate_time
  lipidomics_sample_info$accurate_time == proteomics_sample_info$accurate_time
  lipidomics_sample_info$accurate_time == cytokine_sample_info$accurate_time
  lipidomics_sample_info$accurate_time == cortisol_sample_info$accurate_time
  lipidomics_sample_info$accurate_time == total_protein_sample_info$accurate_time
  
  colnames(lipidomics_expression_data) =
    colnames(metabolomics_expression_data) =
    colnames(cytokine_expression_data) =
    colnames(total_protein_expression_data) =
    colnames(cortisol_expression_data) =
    colnames(metabolic_panel_expression_data) =
    colnames(proteomics_expression_data) =
    as.character(lipidomics_sample_info$accurate_time)
  
  expression_data = 
    rbind(lipidomics_expression_data,
          metabolomics_expression_data,
          cytokine_expression_data,
          total_protein_expression_data,
          cortisol_expression_data,
          metabolic_panel_expression_data,
          proteomics_expression_data)
  
  lipidomics_sample_info$sample_id = 
    as.character(lipidomics_sample_info$accurate_time)
  
  colnames(lipidomics_sample_info)
  
  metabolic_panel_sample_info$sample_id = 
    as.character(metabolic_panel_sample_info$accurate_time)
  
  colnames(metabolic_panel_sample_info)
  
  metabolomics_sample_info$sample_id = 
    as.character(metabolomics_sample_info$accurate_time)
  
  colnames(metabolomics_sample_info)
  
  proteomics_sample_info$sample_id = 
    as.character(proteomics_sample_info$accurate_time)
  
  colnames(proteomics_sample_info)
  
  cytokine_sample_info$sample_id = 
    as.character(cytokine_sample_info$accurate_time)
  
  colnames(cytokine_sample_info)
  
  cortisol_sample_info$sample_id = 
    as.character(cortisol_sample_info$accurate_time)
  
  colnames(cortisol_sample_info)
  
  total_protein_sample_info$sample_id = 
    as.character(total_protein_sample_info$accurate_time)
  
  colnames(total_protein_sample_info)
  
  sample_info = 
    lipidomics_sample_info %>% 
    dplyr::full_join(metabolic_panel_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time")) %>% 
    dplyr::full_join(metabolomics_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time", "day", "time", "hour")) %>% 
    dplyr::full_join(proteomics_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time"))  %>% 
    dplyr::full_join(cytokine_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time", "day", "time", "hour")) %>% 
    dplyr::full_join(cortisol_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time", "day", "time", "hour")) %>% 
    dplyr::full_join(total_protein_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time")) 
  
  ###variable_info
  metabolic_panel_variable_info =
    metabolic_panel_variable_info %>% 
    dplyr::mutate(data_type = "metabolic_panel")
  
  metabolomics_variable_info= 
    metabolomics_variable_info %>% 
    dplyr::mutate(data_type = "metabolomics") %>% 
    dplyr::mutate(mol_name = Compound.name)
  
  lipidomics_variable_info = 
    lipidomics_variable_info %>% 
    dplyr::mutate(data_type = "lipidomics") %>% 
    dplyr::mutate(mol_name = case_when(
      is.na(Lipid_Name) ~ mol_name,
      !is.na(Lipid_Name) ~ Lipid_Name
    ))
  
  proteomics_variable_info = 
    proteomics_variable_info %>% 
    dplyr::mutate(data_type = "proteomics")
  
  cytokine_variable_info = 
    cytokine_variable_info %>% 
    dplyr::mutate(data_type = "cytokine") %>% 
    dplyr::rename(cytokine_classification = classification)
  
  cortisol_variable_info = 
    cortisol_variable_info %>% 
    dplyr::mutate(data_type = "cytokine") %>% 
    dplyr::rename(cortisol_class = class)
  
  total_protein_variable_info = 
    total_protein_variable_info %>% 
    dplyr::mutate(data_type = "total_protein") %>% 
    dplyr::rename(total_protein_class = class)
  
  variable_info =
    rbind(lipidomics_variable_info[,c("variable_id", "mol_name", "data_type")],
          metabolomics_variable_info[,c("variable_id", "mol_name", "data_type")],
          cytokine_variable_info[,c("variable_id", "mol_name", "data_type")],
          total_protein_variable_info[,c("variable_id", "mol_name", "data_type")],
          cortisol_variable_info[,c("variable_id", "mol_name", "data_type")],
          metabolic_panel_variable_info[,c("variable_id", "mol_name", "data_type")],
          proteomics_variable_info[,c("variable_id", "mol_name", "data_type")])
  
  dim(variable_info)
  
  colnames(expression_data) == sample_info$sample_id
  rownames(expression_data) == variable_info$variable_id  
}




######sleep vs omics
dir.create("data/7_24_mike/sleep/sleep_omics")
setwd("data/7_24_mike/sleep/sleep_omics")

grep("Caf", variable_info$mol_name)
variable_info[1116,]

library(metid)

load(here::here("data/7_24_mike/metabolomics/orbitrapDatabase0.0.1"))

which(metabolomics_variable_info_raw$variable_id == variable_info[1116,]$variable_id)

metabolomics_variable_info_raw[388,]

load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/POS/NCE25/result_hilic_pos25"))
load(here::here("data/7_24_mike/metabolomics/metabolite_annotation/HILIC/POS/NCE35/result_hilic_pos35"))

ms2plot(object = result_hilic_pos25[[7]], database = orbitrapDatabase0.0.1, 
        which.peak = "0.90_195.0916m/z")

ms2plot(object = result_hilic_pos35[[7]], database = orbitrapDatabase0.0.1, 
        which.peak = "0.90_195.0916m/z")

plot = 
time_plot(
  x = (as.numeric(expression_data[1116, ]) - mean(as.numeric(expression_data[1116, ])))/sd(as.numeric(expression_data[1116, ])),
  time = sample_info$accurate_time,
  day_night_df = day_night_df, x_name = "Caffeine",
    add_point = TRUE, 
  x_color = class_color["metabolomics"], y_name = "Caffeine"
)
plot
# ggsave(plot, filename = "caffeine_plot.pdf", width = 20, height = 7)

day_night_df

sleep_variable_info
sleep_expression_data

plot(sleep_sample_info$seconds)

library(plyr)

day = unique(sleep_sample_info$day)

sleep_section =
  data.frame(
    section = c(1, 2, 3, 4, 5, 6, 7, 8),
    from_time = c(
      "2019-04-29 06:00:00 PDT",
      "2019-04-30 06:00:00 PDT",
      "2019-05-01 06:00:00 PDT",
      "2019-05-02 06:00:00 PDT",
      "2019-05-03 06:00:00 PDT",
      "2019-05-04 06:00:00 PDT",
      "2019-05-05 06:00:00 PDT",
      "2019-05-06 06:00:00 PDT"
    ),
    end_time = c(
      "2019-04-30 06:00:00 PDT",
      "2019-05-01 06:00:00 PDT",
      "2019-05-02 06:00:00 PDT",
      "2019-05-03 06:00:00 PDT",
      "2019-05-04 06:00:00 PDT",
      "2019-05-05 06:00:00 PDT",
      "2019-05-06 06:00:00 PDT",
      "2019-05-07 06:00:00 PDT"
    )
  ) %>% 
  dplyr::mutate(
    from_time = as.POSIXct(from_time),
    end_time = as.POSIXct(end_time)
  )

sleep_data =
  sleep_expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample_id") %>% 
  dplyr::mutate(sleep_value = case_when(
    sleep1 == "restless" ~ 0,
    sleep1 == "wake" ~ 0,
    sleep1 == "awake" ~ 0,
    sleep1 == "asleep" ~ 1,
    sleep1 == "light" ~ 1,
    sleep1 == "deep" ~ 1,
    sleep1 == "rem" ~ 1
  ))

sleep_section_data =
sleep_section %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    temp_data = 
    sleep_sample_info %>% 
      dplyr::filter(accurate_time >= as.POSIXct(x[2]) & accurate_time <= as.POSIXct(x[3])) %>% 
      dplyr::left_join(sleep_data,
                       by = "sample_id") %>% 
      dplyr::select(seconds, sleep1, sleep_value) %>% 
      dplyr::arrange(sleep_value) %>% 
      dplyr::group_by(sleep1, sleep_value) %>% 
      dplyr::summarise(sum_seconds = sum(seconds)) %>% 
      dplyr::ungroup() %>% 
      dplyr::arrange(sleep_value) %>% 
      dplyr::mutate(section = x[1])
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(sleep_section_data) = NULL

dim(expression_data)
colnames(expression_data)

#####
omics_section_data =
  sleep_section %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    temp_sample_id = 
      sample_info %>% 
      dplyr::filter(accurate_time >= as.POSIXct(x[2]) & accurate_time <= as.POSIXct(x[3])) %>% 
      dplyr::pull(sample_id)
    
    expression_data[,temp_sample_id] %>% 
      apply(1, mean)
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

colnames(omics_section_data) = 1:ncol(omics_section_data)

dim(omics_section_data)
dim(sleep_section_data)

library(plyr)
sleep_section_data2 = 
sleep_section_data %>% 
  plyr::dlply(.variables = .(section)) %>% 
  purrr::map(function(x){
    sleep_score = sum(x$sleep_value * x$sum_seconds)
    x %>% 
      dplyr::mutate(sleep_score = sleep_score) %>% 
      dplyr::select(section, sleep_score) %>% 
      dplyr::distinct(section, .keep_all = TRUE)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()
  
dim(sleep_section_data2)
dim(omics_section_data)

####calculate the correlation between sleep quality and metabolite
result = 
omics_section_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    test = 
      cor.test(x, sleep_section_data2$sleep_score, method = "spearman")
    c(cor = unname(test$estimate), p = unname(test$p.value))
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(p.adjust = p.adjust(p, method = "BH"))

grep("Caf",variable_info$mol_name)
grep("Caf",variable_info$mol_name, value = TRUE)
  
variable_info[1116,]

which(rownames(result) == "HILIC_POS_0.90_195.0916m/z")
result[1116,]

temp_data = 
  data.frame(sleep = sleep_section_data2$sleep_score,
             caffeine = as.numeric(omics_section_data[1116,]),
             day = 1:8)

temp_data %>% 
  ggplot(aes(caffeine, sleep)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = day)) +
  geom_smooth(method = "lm", se = FALSE) +
  base_theme

result[1116,]

plot(result$cor)

result = 
result %>% 
  dplyr::arrange(desc(abs(cor))) %>% 
  dplyr::filter(abs(cor) > 0.7)

dir.create("scatter_plot")
# ######output the plot about the sleep and omics
# for (i in 1:nrow(result)) {
#   cat(i, " ")
#   temp_data = 
#     data.frame(sleep = (sleep_section_data2$sleep_score - mean(sleep_section_data2$sleep_score))/sd(sleep_section_data2$sleep_score),
#                omics = (as.numeric(omics_section_data[rownames(result)[i],]) - mean(as.numeric(omics_section_data[rownames(result)[i],])))/sd(as.numeric(omics_section_data[rownames(result)[i],])),
#                day = 1:8)
#   
#   omics_name = variable_info$mol_name[variable_info$variable_id == rownames(result)[i]]
#   omics_type = variable_info$data_type[variable_info$variable_id == rownames(result)[i]]
#   plot = 
#   temp_data %>% 
#     ggplot(aes(omics, sleep)) +
#     geom_point(size = 5, color = class_color[omics_type]) +
#     ggrepel::geom_text_repel(aes(label = day), size = 5) +
#     geom_smooth(method = "lm", se = FALSE, color = "black") +
#     base_theme +
#     labs(y = "Sleep score", x = omics_name, 
#          title = paste("Spearman correlation: ",round(result$cor[i], 2),
#                        "; p value:", round(result$p[i], 3), sep = ""))
#   plot
#   name = rownames(result)[i] %>% 
#     stringr::str_replace_all("\\/", "_")
#   ggsave(plot, 
#          filename = file.path("scatter_plot", paste(name, ".pdf", sep = "")),
#          width = 7, height = 7)
# }

grep("Caf", variable_info$mol_name, value = TRUE)
grep("Caf", variable_info$mol_name)
variable_info[1116,]

temp_data = 
  result %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  dplyr::left_join(variable_info, by = "variable_id")


####remove the metabolites which are not from non-human database
temp_data = 
temp_data %>%
  dplyr::left_join(metabolomics_variable_info[, c("variable_id", "Database")],
                   by = "variable_id") %>%
  dplyr::filter(Database %in% c("hmdbDatabase0.0.2", 
                                "metlinDatabase0.0.2",
                                "msDatabase_hilic0.0.2", 
                                "msDatabase_rplc0.0.2",
                                "nistDatabase0.0.2") | mol_name == "Caffeine" | is.na(Database)) %>% 
  dplyr::arrange(cor) %>%
  dplyr::mutate(mol_name = factor(mol_name, mol_name))

plot = 
temp_data %>%
  dplyr::mutate(hjust = case_when(
    cor > 0 ~ 1,
    cor < 0 ~ 0
  )) %>% 
  ggplot(aes(cor, mol_name)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(
    x = 0,
    xend = cor,
    y = mol_name,
    yend = mol_name,
    color = data_type
  )) +
  geom_point(aes(size = -log(p, 10),
                 color = data_type)) +
  base_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(x = "Spearman correlation", y = "") +
  scale_color_manual(values = class_color) +
  geom_text(
    aes(
      x = 0,
      y = mol_name,
      label = mol_name,
      color = data_type,
      hjust = hjust
    ),
    # hjust = hjust,
    size = 2
  )
plot
ggsave(plot, filename = "sleep_vs_molecules.pdf", width = 7, height = 10)

####output the result
openxlsx::write.xlsx(temp_data,
                     "sleep_vs_omics.xlsx",
                     asTable = TRUE,
                     overwrite = TRUE)



































