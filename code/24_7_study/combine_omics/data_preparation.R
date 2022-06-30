no_function()
library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")
source("code/modified_dtw.R")
source("code/lagged_correlation.R")

{
  ###metabolic panel data
  {
    ###metabolic_panel
    load("data/24_7_study/metabolic_panel/data_preparation/sample_info")
    load("data/24_7_study/metabolic_panel/data_preparation/variable_info")
    load("data/24_7_study/metabolic_panel/data_preparation/expression_data")
    metabolic_panel_sample_info = sample_info
    metabolic_panel_variable_info = variable_info
    metabolic_panel_expression_data = expression_data
    
    remain_idx =
      which(apply(expression_data, 1, function(x) {
        sum(x == 0) / ncol(expression_data)
      }) < 0.5)
    
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
    
    load("data/24_7_study/summary_info/day_night_df")
    
    ####this is for the day night time
    day_night_df =
      day_night_df %>%
      dplyr::mutate(
        start_time = as.POSIXct(hms::as_hms(start)),
        end_time = as.POSIXct(hms::as_hms(end)),
        week = format(day, "%a")
      ) %>%
      dplyr::mutate(week = paste(week,
                                 lubridate::month(day),
                                 lubridate::day(day),
                                 sep = "-")) %>%
      dplyr::mutate(week = factor(week, unique(week)))
    
  }
  
  
  ####cytokine data
  {
    ###cytokine
    load("data/24_7_study/cytokine/data_preparation/sample_info")
    load("data/24_7_study/cytokine/data_preparation/variable_info")
    load("data/24_7_study/cytokine/data_preparation/expression_data")
    cytokine_sample_info = sample_info
    cytokine_variable_info = variable_info
    cytokine_expression_data = expression_data
    
    remain_idx =
      which(apply(cytokine_expression_data, 1, function(x) {
        sum(x == 0) / ncol(cytokine_expression_data)
      }) < 0.5)
    
    cytokine_variable_info =
      cytokine_variable_info[remain_idx,]
    
    cytokine_expression_data =
      cytokine_expression_data[cytokine_variable_info$variable_id,]
    
    
    ###remove the controls
    remove_idx = grep("CHEX", cytokine_variable_info$mol_name)
    if (length(remove_idx) > 0) {
      cytokine_variable_info =
        cytokine_variable_info[-remove_idx,]
      cytokine_expression_data =
        cytokine_expression_data[-remove_idx,]
    }
  }
  
  ##lipids
  {
    ###lipidomics
    load("data/24_7_study/lipidomics/data_preparation/sample_info")
    load("data/24_7_study/lipidomics/data_preparation/variable_info")
    load("data/24_7_study/lipidomics/data_preparation/expression_data")
    lipidomics_sample_info = sample_info
    lipidomics_variable_info = variable_info
    lipidomics_expression_data = expression_data
    
    
    remain_idx =
      which(apply(lipidomics_expression_data, 1, function(x) {
        sum(x == 0) / ncol(lipidomics_expression_data)
      }) < 0.5)
    
    lipidomics_variable_info =
      lipidomics_variable_info[remain_idx,]
    
    lipidomics_expression_data =
      lipidomics_expression_data[lipidomics_variable_info$variable_id,]
    
  }
  
  
  ###metabolomics
  {
    load("data/24_7_study/metabolomics/data_preparation/metabolites/sample_info")
    load("data/24_7_study/metabolomics/data_preparation/metabolites/variable_info")
    load(
      "data/24_7_study/metabolomics/data_preparation/metabolites/expression_data"
    )
    metabolomics_sample_info = sample_info
    metabolomics_variable_info = variable_info
    metabolomics_expression_data = expression_data
    
    ###remove wrong metabolites
    metabolomics_variable_info =
      metabolomics_variable_info %>%
      dplyr::filter(
        Database %in% c(
          "hmdbDatabase0.0.2",
          "metlinDatabase0.0.2",
          "msDatabase_hilic0.0.2",
          "msDatabase_rplc0.0.2",
          "nistDatabase0.0.2"
        )
      )
    
    metabolomics_expression_data =
      metabolomics_expression_data[metabolomics_variable_info$variable_id,]
    
    remain_idx =
      which(apply(metabolomics_expression_data, 1, function(x) {
        sum(x == 0) / ncol(metabolomics_expression_data)
      }) < 0.5)
    
    metabolomics_variable_info =
      metabolomics_variable_info[remain_idx,]
    
    metabolomics_expression_data =
      metabolomics_expression_data[metabolomics_variable_info$variable_id,]
    
    ###remove QC samples
    metabolomics_sample_info =
      metabolomics_sample_info %>%
      dplyr::filter(!is.na(accurate_time))
    
    metabolomics_expression_data =
      metabolomics_expression_data[, metabolomics_sample_info$sample_id]
    
  }
  
  
  #####proteomics
  {
    ###proteomics
    load("data/24_7_study/proteomics/data_preparation/sample_info")
    load("data/24_7_study/proteomics/data_preparation/variable_info")
    load("data/24_7_study/proteomics/data_preparation/expression_data")
    proteomics_sample_info = sample_info
    proteomics_variable_info = variable_info
    proteomics_expression_data = expression_data
    
    load("data/24_7_study/summary_info/day_night_df")
    
    remain_idx =
      which(apply(proteomics_expression_data, 1, function(x) {
        sum(x == 0) / ncol(proteomics_expression_data)
      }) < 0.5)
    
    proteomics_variable_info =
      proteomics_variable_info[remain_idx,]
    
    proteomics_expression_data =
      proteomics_expression_data[proteomics_variable_info$variable_id,]
  }
  
  
  #####cortisol
  {
    ###cortisol
    load("data/24_7_study/cortisol/data_preparation/sample_info")
    load("data/24_7_study/cortisol/data_preparation/variable_info")
    load("data/24_7_study/cortisol/data_preparation/expression_data")
    cortisol_sample_info = sample_info
    cortisol_variable_info = variable_info
    cortisol_expression_data = expression_data
    
    load("data/24_7_study/summary_info/day_night_df")
    
    remain_idx =
      which(apply(cortisol_expression_data, 1, function(x) {
        sum(x == 0) / ncol(cortisol_expression_data)
      }) < 0.5)
    
    cortisol_variable_info =
      cortisol_variable_info[remain_idx,]
    
    cortisol_expression_data =
      cortisol_expression_data[cortisol_variable_info$variable_id,]
  }
  
  
  
  #####total_protein
  {
    ###total_protein
    load("data/24_7_study/total_protein/data_preparation/sample_info")
    load("data/24_7_study/total_protein/data_preparation/variable_info")
    load("data/24_7_study/total_protein/data_preparation/expression_data")
    total_protein_sample_info = sample_info
    total_protein_variable_info = variable_info
    total_protein_expression_data = expression_data
    
    load("data/24_7_study/summary_info/day_night_df")
    
    remain_idx =
      which(apply(total_protein_expression_data, 1, function(x) {
        sum(x == 0) / ncol(total_protein_expression_data)
      }) < 0.5)
    
    total_protein_variable_info =
      total_protein_variable_info[remain_idx,]
    
    total_protein_expression_data =
      total_protein_expression_data[total_protein_variable_info$variable_id,]
  }
  
  
  ####combine all omics together
  intersect_time =
    Reduce(
      f = intersect,
      x = list(
        as.character(lipidomics_sample_info$accurate_time),
        as.character(metabolomics_sample_info$accurate_time),
        as.character(cytokine_sample_info$accurate_time),
        as.character(total_protein_sample_info$accurate_time),
        as.character(cortisol_sample_info$accurate_time),
        as.character(metabolic_panel_sample_info$accurate_time),
        as.character(proteomics_sample_info$accurate_time)
      )
    )
  
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
  
  lipidomics_sample_info$accurate_time == sample_info$accurate_time
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
    rbind(
      lipidomics_expression_data,
      metabolomics_expression_data,
      cytokine_expression_data,
      total_protein_expression_data,
      cortisol_expression_data,
      metabolic_panel_expression_data,
      proteomics_expression_data
    )
  
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
    dplyr::full_join(
      metabolomics_sample_info,
      by = c(
        "subject_id",
        "sample_id",
        "accurate_time",
        "day",
        "time",
        "hour"
      )
    ) %>%
    dplyr::full_join(proteomics_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time"))  %>%
    dplyr::full_join(
      cytokine_sample_info,
      by = c(
        "subject_id",
        "sample_id",
        "accurate_time",
        "day",
        "time",
        "hour"
      )
    ) %>%
    dplyr::full_join(
      cortisol_sample_info,
      by = c(
        "subject_id",
        "sample_id",
        "accurate_time",
        "day",
        "time",
        "hour"
      )
    ) %>%
    dplyr::full_join(total_protein_sample_info,
                     by = c("subject_id", "sample_id", "accurate_time"))
  
  ###variable_info
  metabolic_panel_variable_info =
    metabolic_panel_variable_info %>%
    dplyr::mutate(data_type = "metabolic_panel")
  
  metabolomics_variable_info =
    metabolomics_variable_info %>%
    dplyr::mutate(data_type = "metabolomics") %>%
    dplyr::mutate(mol_name = Compound.name)
  
  lipidomics_variable_info =
    lipidomics_variable_info %>%
    dplyr::mutate(data_type = "lipidomics") %>%
    dplyr::mutate(mol_name = case_when(is.na(Lipid_Name) ~ mol_name,
                                       !is.na(Lipid_Name) ~ Lipid_Name))
  
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
    rbind(
      lipidomics_variable_info[, c("variable_id", "mol_name", "data_type")],
      metabolomics_variable_info[, c("variable_id", "mol_name", "data_type")],
      cytokine_variable_info[, c("variable_id", "mol_name", "data_type")],
      total_protein_variable_info[, c("variable_id", "mol_name", "data_type")],
      cortisol_variable_info[, c("variable_id", "mol_name", "data_type")],
      metabolic_panel_variable_info[, c("variable_id", "mol_name", "data_type")],
      proteomics_variable_info[, c("variable_id", "mol_name", "data_type")]
    )
  
  dim(variable_info)
  
  colnames(expression_data) == sample_info$sample_id
  rownames(expression_data) == variable_info$variable_id
}
