#' ---
#' title: "Milk components in body"
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

#+ echo=TRUE
# no_source()

library(tidyverse)
library(data.table)
library(tidymass)

tinyTools::setwd_project()
rm(list=ls())
source("R/tools.R")

###load data from shake and ensure shake study
###load shake study plasma metabolome
load("data/shake_study/metabolome_data_analysis/data_preparation/peaks/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/peaks/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/peaks/variable_info")

plasma_expression_data = expression_data
plasma_sample_info = sample_info
plasma_variable_info = variable_info

###load milk metabolomics data
load("data/ensure_shake/metabolomics/data_preparation/peak/expression_data")
load("data/ensure_shake/metabolomics/data_preparation/peak/sample_info")
load("data/ensure_shake/metabolomics/data_preparation/peak/variable_info")

milk_expression_data = expression_data
milk_sample_info = sample_info
milk_variable_info = variable_info

#+ r, message = 1:10
milk_variable_info$variable_id
plasma_variable_info$Mode

milk_variable_info = 
milk_variable_info %>% 
  dplyr::mutate(mode = case_when(
    stringr::str_detect(variable_id, "rplc_pos") ~ "pRPLC",
    stringr::str_detect(variable_id, "rplc_neg") ~ "nRPLC",
    stringr::str_detect(variable_id, "hilic_pos") ~ "pHILIC",
    stringr::str_detect(variable_id, "hilic_neg") ~ "nHILIC"
  ))
  
####set work directory
tinyTools::setwd_project()
setwd("data/shake_study/milk_component_in_body")

####match peaks in milk and plasma
# ###rplc pos
# data1 = 
#   milk_variable_info %>% dplyr::filter(mode == "pRPLC")
# 
# data2 =
#   plasma_variable_info %>% dplyr::filter(Mode == "pRPLC")
# 
# temp_data = 
# data1 %>% 
#   dplyr::filter(!is.na(Level)) %>% 
#   dplyr::filter(Level == 1 | Level == 2) %>% 
#   dplyr::select(mz, rt, Compound.name, Adduct) %>% 
#   dplyr::left_join(data2[,c("mz", "rt", "Metabolite", "Adduct")], 
#                    by = c("Compound.name" = "Metabolite",
#                           "Adduct" = "Adduct")) %>% 
#   dplyr::filter(!is.na(mz.y)) %>% 
#   dplyr::mutate(mz.error = abs(mz.x - mz.y)*10^6/ifelse(mz.y < 400, 400, mz.y)) %>% 
#   dplyr::mutate(rt.error = abs(rt.x - rt.y))
# 
# plot = 
# temp_data %>% 
#   ggplot(aes(1, mz.error)) +
#   geom_boxplot(aes(x = 1, y = mz.error)) +
#   geom_jitter(aes(x = 1, y = mz.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,1.5)) +
#   labs(x = "", y = "m/z error",
#        title = paste("m/z error, mean: ", round(mean(temp_data$mz.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$mz.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "rplc_pos_mz_error.pdf", width = 3, height = 4)
# 
# plot = 
#   temp_data %>% 
#   ggplot(aes(1, rt.error)) +
#   geom_boxplot(aes(x = 1, y = rt.error)) +
#   geom_jitter(aes(x = 1, y = rt.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,150)) +
#   labs(x = "", y = "RT error",
#        title = paste("RT error, mean: ", round(mean(temp_data$rt.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$rt.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "rplc_pos_rt_error.pdf", width = 3, height = 4)
# 
# rplc_pos_match_idx = 
# tinyTools::mz_rt_match(data1 = data1 %>% dplyr::select(mz, rt),
#                        data2 = data2 %>% dplyr::select(mz, rt),
#                        mz.tol = 1.06+3.05, 
#                        rt.tol = 60.7+97.62, 
#                        rt.error.type = "abs") %>% 
#   as.data.frame()
# 
# rplc_pos_match_idx$milk_variable_id = data1$variable_id[rplc_pos_match_idx$Index1] 
# rplc_pos_match_idx$plasma_variable_id = data2$variable_id[rplc_pos_match_idx$Index2] 
# 
# ###remove duplicated
# rplc_pos_match_idx = 
# rplc_pos_match_idx %>% 
#   dplyr::group_by(milk_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::group_by(plasma_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup()
#   
# save(rplc_pos_match_idx, file = "rplc_pos_match_idx")

getwd()
load("rplc_pos_match_idx")

# ###rplc neg
# data1 = 
#   milk_variable_info %>% dplyr::filter(mode == "nRPLC")
# 
# data2 =
#   plasma_variable_info %>% dplyr::filter(Mode == "nRPLC")
# 
# temp_data = 
#   data1 %>% 
#   dplyr::filter(!is.na(Level)) %>% 
#   dplyr::filter(Level == 1 | Level == 2) %>% 
#   dplyr::select(mz, rt, Compound.name, Adduct) %>% 
#   dplyr::left_join(data2[,c("mz", "rt", "Metabolite", "Adduct")], 
#                    by = c("Compound.name" = "Metabolite",
#                           "Adduct" = "Adduct")) %>% 
#   dplyr::filter(!is.na(mz.y)) %>% 
#   dplyr::mutate(mz.error = abs(mz.x - mz.y)*10^6/ifelse(mz.y < 400, 400, mz.y)) %>% 
#   dplyr::mutate(rt.error = abs(rt.x - rt.y))
# 
# dim(temp_data)
# 
# temp_data
# 
# temp_data = 
#   temp_data %>% 
#   dplyr::filter(rt.error < 150)
# 
# plot = 
#   temp_data %>% 
#   ggplot(aes(1, mz.error)) +
#   geom_boxplot(aes(x = 1, y = mz.error)) +
#   geom_jitter(aes(x = 1, y = mz.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,1.5)) +
#   labs(x = "", y = "m/z error",
#        title = paste("m/z error, mean: ", round(mean(temp_data$mz.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$mz.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "rplc_neg_mz_error.pdf", width = 3, height = 4)
# 
# plot = 
#   temp_data %>% 
#   ggplot(aes(1, rt.error)) +
#   geom_boxplot(aes(x = 1, y = rt.error)) +
#   geom_jitter(aes(x = 1, y = rt.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,150)) +
#   labs(x = "", y = "RT error",
#        title = paste("RT error, mean: ", round(mean(temp_data$rt.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$rt.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "rplc_neg_rt_error.pdf", width = 3, height = 4)
# 
# rplc_neg_match_idx = 
#   tinyTools::mz_rt_match(data1 = data1 %>% dplyr::select(mz, rt),
#                          data2 = data2 %>% dplyr::select(mz, rt),
#                          mz.tol = 0.66+0.6, 
#                          rt.tol = 47.23+58.81, 
#                          rt.error.type = "abs") %>% 
#   as.data.frame()
# 
# rplc_neg_match_idx$milk_variable_id = data1$variable_id[rplc_neg_match_idx$Index1] 
# rplc_neg_match_idx$plasma_variable_id = data2$variable_id[rplc_neg_match_idx$Index2] 
# 
# ###remove duplicated
# rplc_neg_match_idx = 
#   rplc_neg_match_idx %>% 
#   dplyr::group_by(milk_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::group_by(plasma_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup()
# 
# rplc_neg_match_idx
# 
# save(rplc_neg_match_idx, file = "rplc_neg_match_idx")
getwd()
load("rplc_neg_match_idx")

# ###hilic pos
# data1 = 
#   milk_variable_info %>% dplyr::filter(mode == "pHILIC")
# 
# data2 =
#   plasma_variable_info %>% dplyr::filter(Mode == "pHILIC")
# 
# temp_data = 
#   data1 %>% 
#   dplyr::filter(!is.na(Level)) %>% 
#   dplyr::filter(Level == 1 | Level == 2) %>% 
#   dplyr::select(mz, rt, Compound.name, Adduct) %>% 
#   dplyr::left_join(data2[,c("mz", "rt", "Metabolite", "Adduct")], 
#                    by = c("Compound.name" = "Metabolite",
#                           "Adduct" = "Adduct")) %>% 
#   dplyr::filter(!is.na(mz.y)) %>% 
#   dplyr::mutate(mz.error = abs(mz.x - mz.y)*10^6/ifelse(mz.y < 400, 400, mz.y)) %>% 
#   dplyr::mutate(rt.error = abs(rt.x - rt.y))
# 
# temp_data =
#   temp_data %>% 
#   dplyr::filter(rt.error < 150) %>% 
#   dplyr::filter(mz.error < 25)
# 
# plot = 
#   temp_data %>% 
#   ggplot(aes(1, mz.error)) +
#   geom_boxplot(aes(x = 1, y = mz.error)) +
#   geom_jitter(aes(x = 1, y = mz.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,1.5)) +
#   labs(x = "", y = "m/z error",
#        title = paste("m/z error, mean: ", round(mean(temp_data$mz.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$mz.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "hilic_pos_mz_error.pdf", width = 3, height = 4)
# 
# plot = 
#   temp_data %>% 
#   ggplot(aes(1, rt.error)) +
#   geom_boxplot(aes(x = 1, y = rt.error)) +
#   geom_jitter(aes(x = 1, y = rt.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,150)) +
#   labs(x = "", y = "RT error",
#        title = paste("RT error, mean: ", round(mean(temp_data$rt.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$rt.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "hilic_pos_rt_error.pdf", width = 3, height = 4)
# 
# 
# hilic_pos_match_idx = 
#   tinyTools::mz_rt_match(data1 = data1 %>% dplyr::select(mz, rt),
#                          data2 = data2 %>% dplyr::select(mz, rt),
#                          mz.tol = 4.32+4.51, 
#                          rt.tol = 38.41+34.83, 
#                          rt.error.type = "abs") %>% 
#   as.data.frame()
# 
# hilic_pos_match_idx$milk_variable_id = data1$variable_id[hilic_pos_match_idx$Index1] 
# hilic_pos_match_idx$plasma_variable_id = data2$variable_id[hilic_pos_match_idx$Index2] 
# 
# ###remove duplicated
# hilic_pos_match_idx = 
#   hilic_pos_match_idx %>% 
#   dplyr::group_by(milk_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::group_by(plasma_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup()
# 
# save(hilic_pos_match_idx, file = "hilic_pos_match_idx")

load("hilic_pos_match_idx")

# ###hilic neg
# data1 = 
#   milk_variable_info %>% dplyr::filter(mode == "nHILIC")
# 
# data2 =
#   plasma_variable_info %>% dplyr::filter(Mode == "nHILIC")
# 
# temp_data = 
#   data1 %>% 
#   dplyr::filter(!is.na(Level)) %>% 
#   dplyr::filter(Level == 1 | Level == 2) %>% 
#   dplyr::select(mz, rt, Compound.name, Adduct) %>% 
#   dplyr::left_join(data2[,c("mz", "rt", "Metabolite", "Adduct")], 
#                    by = c("Compound.name" = "Metabolite",
#                           "Adduct" = "Adduct")) %>% 
#   dplyr::filter(!is.na(mz.y)) %>% 
#   dplyr::mutate(mz.error = abs(mz.x - mz.y)*10^6/ifelse(mz.y < 400, 400, mz.y)) %>% 
#   dplyr::mutate(rt.error = abs(rt.x - rt.y))
# 
# dim(temp_data)
# 
# temp_data
# 
# temp_data = 
#   temp_data %>% 
#   dplyr::filter(rt.error < 150) %>% 
#   dplyr::filter(mz.error < 25)
# 
# plot = 
#   temp_data %>% 
#   ggplot(aes(1, mz.error)) +
#   geom_boxplot(aes(x = 1, y = mz.error)) +
#   geom_jitter(aes(x = 1, y = mz.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,1.5)) +
#   labs(x = "", y = "m/z error",
#        title = paste("m/z error, mean: ", round(mean(temp_data$mz.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$mz.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "hilic_neg_mz_error.pdf", width = 3, height = 4)
# 
# plot = 
#   temp_data %>% 
#   ggplot(aes(1, rt.error)) +
#   geom_boxplot(aes(x = 1, y = rt.error)) +
#   geom_jitter(aes(x = 1, y = rt.error), 
#               size = 5, alpha = 0.7) +
#   scale_y_continuous(limits = c(0,150)) +
#   labs(x = "", y = "RT error",
#        title = paste("RT error, mean: ", round(mean(temp_data$rt.error), 2), 
#                      "/",
#                      "sd:",round(sd(temp_data$rt.error), 2), sep = "")) +
#   base_theme +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.grid = element_blank())
# 
# plot
# ggsave(plot, filename = "hilic_neg_rt_error.pdf", width = 3, height = 4)
# 
# hilic_neg_match_idx = 
#   tinyTools::mz_rt_match(data1 = data1 %>% dplyr::select(mz, rt),
#                          data2 = data2 %>% dplyr::select(mz, rt),
#                          mz.tol = 1.19+2.58, 
#                          rt.tol = 40.39 + 29.31, 
#                          rt.error.type = "abs") %>% 
#   as.data.frame()
# 
# hilic_neg_match_idx$milk_variable_id = data1$variable_id[hilic_neg_match_idx$Index1] 
# hilic_neg_match_idx$plasma_variable_id = data2$variable_id[hilic_neg_match_idx$Index2] 
# 
# ###remove duplicated
# hilic_neg_match_idx = 
#   hilic_neg_match_idx %>% 
#   dplyr::group_by(milk_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::group_by(plasma_variable_id) %>% 
#   dplyr::filter(`rt error` == min(`rt error`)) %>% 
#   dplyr::ungroup()
# 
# hilic_neg_match_idx
# 
# save(hilic_neg_match_idx, file = "hilic_neg_match_idx")

load("hilic_neg_match_idx")

# plasma_milk_match_table =
# rbind(
#   rplc_pos_match_idx,
#   rplc_neg_match_idx,
#   hilic_pos_match_idx,
#   hilic_neg_match_idx
# ) %>%
#   dplyr::select(milk_variable_id, plasma_variable_id, mz.error = `mz error`, rt.error = `rt error`) %>%
#   dplyr::group_by(milk_variable_id) %>%
#   dplyr::filter(rt.error == min(rt.error)) %>%
#   dplyr::filter(mz.error == min(mz.error)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(plasma_variable_id) %>%
#   dplyr::filter(rt.error == min(rt.error)) %>%
#   dplyr::filter(mz.error == min(mz.error)) %>%
#   dplyr::ungroup()
# 
# sum(duplicated(plasma_milk_match_table$milk_variable_id))
# sum(duplicated(plasma_milk_match_table$plasma_variable_id))
# 
# ratio_to_blank =
# milk_expression_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x){
#     x = as.numeric(x)
#     mean(x[1])/x[3]
#   }) %>%
#   unlist()
# 
# ratio_to_blank[is.infinite(ratio_to_blank)] = 100
# ratio_to_blank[is.na(ratio_to_blank)] = 0
# 
# sum(is.na(ratio_to_blank))
# sum(is.infinite(ratio_to_blank))
# 
# plasma_milk_match_table =
# plasma_milk_match_table %>%
#   dplyr::filter(milk_variable_id %in% milk_variable_info$variable_id[which(ratio_to_blank > 2)]) %>%
#   dplyr::left_join(milk_variable_info[,c("variable_id", "Compound.name")], by = c("milk_variable_id" = "variable_id")) %>%
#   dplyr::left_join(plasma_variable_info[,c("variable_id", "Metabolite")], by = c("plasma_variable_id" = "variable_id")) %>%
#   dplyr::rename(milk_annotation = Compound.name, plasma_annotation = Metabolite) %>%
#   dplyr::distinct(milk_variable_id, .keep_all = TRUE) %>%
#   dplyr::distinct(plasma_variable_id, .keep_all = TRUE)
# 
# 
# ####remove the features which are detected in baseline samples
T0_int =
plasma_expression_data %>%
  dplyr::select(contains("T0")) %>%
  apply(1, mean)

T30_int =
  plasma_expression_data %>%
  dplyr::select(contains("T30")) %>%
  apply(1, mean)

T60_int =
  plasma_expression_data %>%
  dplyr::select(contains("T60")) %>%
  apply(1, mean)

T120_int =
  plasma_expression_data %>%
  dplyr::select(contains("T120")) %>%
  apply(1, mean)

T240_int =
  plasma_expression_data %>%
  dplyr::select(contains("T240")) %>%
  apply(1, mean)

quantile(c(T30_int, T60_int, T120_int, T240_int), probs = 0.05)

plot(density(log(c(T30_int, T60_int, T120_int, T240_int) + 1, 2)))

quantile(log(c(T30_int, T60_int, T120_int, T240_int) + 1, 2), probs = 0.05)

not_in_baseline = names(which(T0_int < 352))
# 
# plasma_milk_match_table = 
# plasma_milk_match_table %>% 
#   dplyr::mutate(in_baseline = case_when(
#     plasma_variable_id %in% not_in_baseline ~ "No",
#     !plasma_variable_id %in% not_in_baseline ~ "Yes"
#   ))
# 
# table(plasma_milk_match_table$in_baseline)
# 
# 
# which(duplicated(plasma_milk_match_table$milk_variable_id))
# 
# sum(duplicated(plasma_milk_match_table$milk_variable_id))
# sum(duplicated(plasma_milk_match_table$plasma_variable_id))
# 
# save(plasma_milk_match_table, file = "plasma_milk_match_table")
getwd()
load("plasma_milk_match_table")

library(openxlsx)

# wb = createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "Plasma milk match table", gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = plasma_milk_match_table,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "plasma_milk_match_table.xlsx", overwrite = TRUE) ## save to working directory

plasma_milk_match_table
dim(plasma_milk_match_table)

####there are 558 features are matched in milk and plasma metabolomics data
plasma_milk_match_table %>% 
  dplyr::left_join(plasma_variable_info, 
                   by = c("plasma_variable_id" = "variable_id")) %>% 
  pull(Mode) %>% 
  table()

plot = 
plasma_milk_match_table %>% 
  ggplot() +
  geom_vline(xintercept = median(plasma_milk_match_table$mz.error)) +
  geom_vline(xintercept = mean(plasma_milk_match_table$mz.error)) +
  geom_density(aes(x = mz.error)) +
  labs(x = "Absolute mass to charge ratio (m/z) error (ppm)",
       y = "Density") +
  base_theme

range(plasma_milk_match_table$mz.error)

plot

# ggsave(plot, filename = "feature_matching_mz_error_distributation.pdf", width = 7, height = 7)

plot = 
plasma_milk_match_table %>% 
  ggplot() +
  geom_vline(xintercept = median(plasma_milk_match_table$rt.error)) +
  geom_vline(xintercept = mean(plasma_milk_match_table$rt.error)) +
  geom_density(aes(x = rt.error)) +
  labs(x = "Absolute retention time (RT) error (second)",
       y = "Density") +
  base_theme 

plot

# ggsave(plot, filename = "feature_matching_rt_error_distributation.pdf", width = 7, height = 7)

mean(plasma_milk_match_table$rt.error)
median(plasma_milk_match_table$rt.error)
range(plasma_milk_match_table$rt.error)

sum(!is.na(plasma_milk_match_table$milk_annotation))
sum(!is.na(plasma_milk_match_table$plasma_annotation))

####overlap of peaks with annotations in plasma and milk
library(ggvenn)

milk = 
plasma_milk_match_table %>% 
  dplyr::filter(!is.na(milk_annotation)) %>% 
  pull(milk_variable_id)

plasma = 
  plasma_milk_match_table %>% 
  dplyr::filter(!is.na(plasma_annotation)) %>% 
  pull(milk_variable_id)

plot = 
ggvenn::ggvenn(data = list(plasma = plasma, milk = milk))

plot

# ggsave(plot, filename = "annotation_milk_plasma_venn.pdf", width = 7, height = 7)

annotation_table = 
plasma_milk_match_table %>% 
  dplyr::filter(!is.na(plasma_annotation) & !is.na(milk_annotation))

# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)
getwd()
annotation_table = readr::read_csv("annotation_table_manual.csv")

temp_data = 
annotation_table$Check %>% 
  table() %>% 
  data.frame() 

colnames(temp_data) = c("Class", "Freq")
  
plot = 
temp_data %>% 
ggplot(aes(x = 2, y = Freq)) +
  geom_bar(aes(x = 2, y = Freq,
               fill = Class), 
           stat = "identity",
           color = "white",
           show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("Same" = ggsci::pal_aaas()(n=10)[2],
                               "Different" = "grey")) +
  theme_void()


plot
# ggsave(plot, filename = "annotation_consistence.pdf", width = 7, height = 7)

# idx = match(plasma_milk_match_table$plasma_variable_id, plasma_variable_info$variable_id)
# save(idx, file = "idx")
# sum(is.na(match(plasma_milk_match_table$plasma_variable_id, plasma_variable_info$variable_id)))

plasma_milk_match_table$in_baseline

##here we define if the peaks are in baseline samples or not.
##we use the peaks in all the samples after drink milk, 5% quantile as the cutoff. 

table(plasma_milk_match_table$in_baseline)

plasma_milk_match_table %>% 
  dplyr::filter(in_baseline == "No")









