#' ---
#' title: "CGM cortisol correlation"
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
tinytools::setwd_project()
library(here)
rm(list = ls())

source(here::here("R/tools.R"))
source(here::here("R/modified_dtw.R"))
source(here::here("R/lagged_correlation.R"))

####load data
{
  ###CGM
  load(here::here("data/7_24_mike/cgm/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/cgm/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/cgm/data_preparation/expression_data"))
  
  cgm_expression_data = expression_data
  cgm_sample_info = sample_info
  cgm_variable_info = variable_info
  
  
  ###HR
  load("data/7_24_mike/hr/data_preparation/sample_info")
  load("data/7_24_mike/hr/data_preparation/variable_info")
  load("data/7_24_mike/hr/data_preparation/expression_data")
  
  hr_expression_data = expression_data
  hr_sample_info = sample_info
  hr_variable_info = variable_info
  
  ###Step
  load("data/7_24_mike/steps/data_preparation/sample_info")
  load("data/7_24_mike/steps/data_preparation/variable_info")
  load("data/7_24_mike/steps/data_preparation/expression_data")
  steps_expression_data = expression_data
  steps_sample_info = sample_info
  steps_variable_info = variable_info
  
}


#####load omics data
load("data/7_24_mike/combine_omics/data_preparation/expression_data")
load("data/7_24_mike/combine_omics/data_preparation/sample_info")
load("data/7_24_mike/combine_omics/data_preparation/variable_info")

omics_expression_data = expression_data
omics_sample_info = sample_info
omics_variable_info = variable_info

dim(omics_expression_data)
dim(omics_sample_info)
dim(omics_variable_info)

colnames(omics_expression_data) == omics_sample_info$sample_id
rownames(omics_expression_data) == omics_variable_info$variable_id

table(omics_variable_info$data_type)

omics_variable_info %>% 
  dplyr::filter(data_type == "cytokine")

load(here::here("data/7_24_mike/summary_info/day_night_df"))

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

tinytools::setwd_project()
setwd("data/7_24_mike/wearable_predict_omics/raw_data/")

#####heart rate, step and scg predict omics
###lasso 
dim(omics_expression_data)

###match omics with heart rate
###time window is 30, 40, 50, 60, 90, 120 mins
####

library(moments)
time_windows = c(30, 40, 50, 60, 90, 120)

# new_hr_expression_data =
# purrr::map(time_windows, function(window) {
#   cat(window, " ")
#   temp_data =
#   omics_sample_info$accurate_time %>%
#     purrr::map(function(x) {
#       diff_time = difftime(x, hr_sample_info$accurate_time, unit = "mins")
#       idx = which(diff_time <= window & diff_time > 0)
#       mean_value = mean(as.numeric(hr_expression_data[1,idx]))
#       median_value = median(as.numeric(hr_expression_data[1,idx]))
#       sd_value = sd(as.numeric(hr_expression_data[1,idx]))
#       max_value = max(as.numeric(hr_expression_data[1,idx]))
#       min_value = min(as.numeric(hr_expression_data[1,idx]))
#       skewness_value = moments::skewness(as.numeric(hr_expression_data[1,idx]))
#       kurtosis_value = moments::kurtosis(as.numeric(hr_expression_data[1,idx]))
#       range_value = abs(diff(range(as.numeric(hr_expression_data[1,idx]))))
#       time_hour = lubridate::hour(x)
#       c(mean_value = mean_value,
#                  median_value = median_value,
#                  sd_value = sd_value,
#                  max_value = max_value,
#                  min_value = min_value,
#                  skewness_value = skewness_value,
#                  kurtosis_value = kurtosis_value,
#                  range_value = range_value,
#                  time_hour = time_hour)
#     }) %>%
#     do.call(cbind, .) %>%
#     as.data.frame()
# 
#   colnames(temp_data) = omics_sample_info$sample_id
#   temp_data
# })
# 
# 
# new_hr_expression_data = 
#   new_hr_expression_data %>% 
#   purrr::map(function(x){
#     remove_idx1 = which(apply(x, 2, function(x)sum(is.na(x))) > 0)
#     remove_idx2 = which(apply(x, 2, function(x)sum(is.infinite(x))) > 0)
#     remove_idx = unique(c(remove_idx1, remove_idx2))
#     if(length(remove_idx)){
#       x[,-remove_idx,drop = FALSE]
#     }else{
#       x
#     }
#   })
# 
# lapply(new_hr_expression_data, ncol) %>% 
#   unlist
# 
# ####add rownames for each dataset
# names(new_hr_expression_data) = 
#   time_windows
# 
# new_hr_expression_data = 
#   purrr::map2(time_windows, new_hr_expression_data, function(x, y){
#     rownames(y) = paste("hr", x, sep = "_", rownames(y))
#     y
#   })
# 
# save(new_hr_expression_data, file = "new_hr_expression_data")

load("new_hr_expression_data")


# new_steps_expression_data =
#   purrr::map(time_windows, function(window) {
#     cat(window, " ")
#     temp_data =
#       omics_sample_info$accurate_time %>%
#       purrr::map(function(x) {
#         diff_time = difftime(x, steps_sample_info$accurate_time, unit = "mins")
#         idx = which(diff_time <= window & diff_time > 0)
#         mean_value = mean(as.numeric(steps_expression_data[1, idx]))
#         median_value = median(as.numeric(steps_expression_data[1, idx]))
#         sd_value = sd(as.numeric(steps_expression_data[1, idx]))
#         max_value = max(as.numeric(steps_expression_data[1, idx]))
#         min_value = min(as.numeric(steps_expression_data[1, idx]))
#         skewness_value = moments::skewness(as.numeric(steps_expression_data[1, idx]))
#         kurtosis_value = moments::kurtosis(as.numeric(steps_expression_data[1, idx]))
#         range_value = abs(diff(range(as.numeric(
#           steps_expression_data[1, idx]
#         ))))
#         time_hour = lubridate::hour(x)
#         c(
#           mean_value = mean_value,
#           median_value = median_value,
#           sd_value = sd_value,
#           max_value = max_value,
#           min_value = min_value,
#           skewness_value = skewness_value,
#           kurtosis_value = kurtosis_value,
#           range_value = range_value,
#           time_hour = time_hour
#         )
#       }) %>%
#       do.call(cbind, .) %>%
#       as.data.frame()
#     
#     colnames(temp_data) = omics_sample_info$sample_id
#     temp_data
#   })
# 
# new_steps_expression_data =
#   new_steps_expression_data %>%
#   purrr::map(function(x) {
#     remove_idx1 = which(apply(x, 2, function(x)
#       sum(is.na(x))) > 0)
#     remove_idx2 = which(apply(x, 2, function(x)
#       sum(is.infinite(x))) > 0)
#     remove_idx = unique(c(remove_idx1, remove_idx2))
#     if (length(remove_idx)) {
#       x[, -remove_idx, drop = FALSE]
#     } else{
#       x
#     }
#   })
# 
# lapply(new_steps_expression_data, ncol) %>%
#   unlist
# 
# 
# names(new_steps_expression_data) =
#   time_windows
# 
# new_steps_expression_data =
#   purrr::map2(time_windows, new_steps_expression_data, function(x, y) {
#     rownames(y) = paste("steps", x, sep = "_", rownames(y))
#     y
#   })
# 
# save(new_steps_expression_data, file = "new_steps_expression_data")

load("new_steps_expression_data")

###match omics with heart rate
###time window is 1, 5, 10, 20, 30, 40, 50, 60 minutes
####


# new_cgm_expression_data =
# purrr::map(time_windows, function(window) {
#   cat(window, " ")
#   temp_data =
#   omics_sample_info$accurate_time %>%
#     purrr::map(function(x) {
#       diff_time = difftime(x, cgm_sample_info$accurate_time, unit = "mins")
#       idx = which(diff_time <= window & diff_time > 0)
#       mean_value = mean(as.numeric(cgm_expression_data[1,idx]))
#       median_value = median(as.numeric(cgm_expression_data[1,idx]))
#       sd_value = sd(as.numeric(cgm_expression_data[1,idx]))
#       max_value = max(as.numeric(cgm_expression_data[1,idx]))
#       min_value = min(as.numeric(cgm_expression_data[1,idx]))
#       skewness_value = moments::skewness(as.numeric(cgm_expression_data[1,idx]))
#       kurtosis_value = moments::kurtosis(as.numeric(cgm_expression_data[1,idx]))
#       range_value = abs(diff(range(as.numeric(cgm_expression_data[1,idx]))))
#       time_hour = lubridate::hour(x)
#       c(mean_value = mean_value,
#                  median_value = median_value,
#                  sd_value = sd_value,
#                  max_value = max_value,
#                  min_value = min_value,
#                  skewness_value = skewness_value,
#                  kurtosis_value = kurtosis_value,
#                  range_value = range_value,
#                  time_hour = time_hour)
#     }) %>%
#     do.call(cbind, .) %>%
#     as.data.frame()
# 
#   colnames(temp_data) = omics_sample_info$sample_id
#   temp_data
# })
# 
# new_cgm_expression_data = 
#   new_cgm_expression_data %>% 
#   purrr::map(function(x){
#     remove_idx1 = which(apply(x, 2, function(x)sum(is.na(x))) > 0)
#     remove_idx2 = which(apply(x, 2, function(x)sum(is.infinite(x))) > 0)
#     remove_idx = unique(c(remove_idx1, remove_idx2))
#     if(length(remove_idx)){
#       x[,-remove_idx,drop = FALSE]
#     }else{
#       x
#     }
#   })
# 
# lapply(new_cgm_expression_data, ncol) %>% 
#   unlist
# 
# ####add rownames for each dataset
#   names(new_cgm_expression_data) = 
#   time_windows
# 
# new_cgm_expression_data = 
#   purrr::map2(time_windows, new_cgm_expression_data, function(x, y){
#     rownames(y) = paste("cgm", x, sep = "_", rownames(y))
#     y
#   })
# 
# save(new_cgm_expression_data, file = "new_cgm_expression_data")

load("new_cgm_expression_data")


#######RF model
#######to predict using time windows is 1 min
library(caret)

dir.create("RF")

for(i in 1:length(new_hr_expression_data)){
  cat(i, "\n")
  temp_hr_data =   new_hr_expression_data[[i]]
  temp_steps_data =   new_steps_expression_data[[i]]
  temp_cgm_data =   new_cgm_expression_data[[i]]

  intersect_sample_id =
    Reduce(intersect, list(
      colnames(temp_hr_data),
      colnames(temp_steps_data),
      colnames(temp_cgm_data)
    ))

  temp_x =
    rbind(temp_hr_data[,intersect_sample_id],
          temp_steps_data[,intersect_sample_id],
          temp_cgm_data[,intersect_sample_id])

  prediction_result =
  purrr::map(1:nrow(omics_expression_data), function(j){
    cat(j, " ")
    temp_y =
    as.numeric(omics_expression_data[j,intersect_sample_id])

    temp_x_y =
      data.frame(data.frame(t(temp_x)), y = temp_y)

    ###7-fold cross validation
    train.control <- caret::trainControl(method = "cv", 
                                         number = 5)

    model <- train(y ~., data = temp_x_y,
                   method = "rf",
                   trControl = train.control)

    result =
    model$results %>%
      as.data.frame() %>%
      dplyr::filter(RMSE == min(RMSE)) %>%
      head(1) %>%
      dplyr::mutate(variable_id = omics_variable_info$variable_id[j]) %>%
      dplyr::select(variable_id, everything())
    result

  }) %>%
    dplyr::bind_rows()

  save(prediction_result, 
       file = paste("RF/rf_prediction_result", time_windows[i], "min", sep = "_"))
}

load("prediction_result_5_min")
prediction_result_5_min = prediction_result
plot(prediction_result_5_min$Rsquared)

load("prediction_result_10_min")
prediction_result_10_min = prediction_result
plot(prediction_result_10_min$Rsquared)

load("prediction_result_20_min")
prediction_result_20_min = prediction_result
plot(prediction_result_20_min$Rsquared)

load("prediction_result_30_min")
prediction_result_30_min = prediction_result
plot(prediction_result_30_min$Rsquared)

load("prediction_result_40_min")
prediction_result_40_min = prediction_result
plot(prediction_result_40_min$Rsquared)

load("prediction_result_50_min")
prediction_result_50_min = prediction_result
plot(prediction_result_50_min$Rsquared)

load("prediction_result_60_min")
prediction_result_60_min = prediction_result
plot(prediction_result_60_min$Rsquared)

load("prediction_result_90_min")
prediction_result_90_min = prediction_result
plot(prediction_result_90_min$Rsquared)

load("prediction_result_120_min")
prediction_result_120_min = prediction_result
plot(prediction_result_120_min$Rsquared)

temp_data =
  data.frame(min5 = prediction_result_5_min$Rsquared,
             min10 = prediction_result_10_min$Rsquared,
             min20 = prediction_result_20_min$Rsquared,
             min30 = prediction_result_30_min$Rsquared,
             min40 = prediction_result_40_min$Rsquared,
             min50 = prediction_result_50_min$Rsquared,
             min60 = prediction_result_60_min$Rsquared,
             min90 = prediction_result_90_min$Rsquared,
             min120 = prediction_result_120_min$Rsquared)
rownames(temp_data) = prediction_result_10_min$variable_id

library(ComplexHeatmap)

plot = 
Heatmap(temp_data,
        cluster_columns = FALSE,
        show_row_names = FALSE, 
        name = "R2", 
        column_names_rot = 45)

library(ggplotify)
plot = ggplotify::as.ggplot(plot)
ggsave(plot, filename = "rf_plot.pdf", width = 5, height = 10)


idx = 
temp_data %>% 
  apply(1, function(x){
    sum(x > 0.3)
  }) %>% 
  `>`(0) %>% 
  which()

plot = 
  Heatmap(temp_data[idx,],
          cluster_columns = FALSE,
          show_row_names = FALSE, 
          name = "R2", 
          column_names_rot = 45)

library(ggplotify)
plot = ggplotify::as.ggplot(plot)
plot
ggsave(plot, filename = "rf_plot2.pdf", width = 7, height = 7)

omics_variable_info %>% 
  dplyr::filter(variable_id %in% names(idx)) %>% 
  pull(data_type) %>% 
  table()

omics_variable_info %>% 
  dplyr::filter(variable_id %in% names(idx)) %>% 
  pull(mol_name)





