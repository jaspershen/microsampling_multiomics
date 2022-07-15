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
masstools::setwd_project()
library(here)
rm(list = ls())

source(here::here("code/tools.R"))
source(here::here("code/modified_dtw.R"))
source(here::here("code/lagged_correlation.R"))

####load data
{
  ###CGM
  load(here::here("data/24_7_study/cgm/data_preparation/sample_info"))
  load(here::here("data/24_7_study/cgm/data_preparation/variable_info"))
  load(here::here("data/24_7_study/cgm/data_preparation/expression_data"))
  
  cgm_expression_data = expression_data
  cgm_sample_info = sample_info
  cgm_variable_info = variable_info
  
  
  ###HR
  load("data/24_7_study/hr/data_preparation/sample_info")
  load("data/24_7_study/hr/data_preparation/variable_info")
  load("data/24_7_study/hr/data_preparation/expression_data")
  
  hr_expression_data = expression_data
  hr_sample_info = sample_info
  hr_variable_info = variable_info
  
  ###Step
  load("data/24_7_study/steps/data_preparation/sample_info")
  load("data/24_7_study/steps/data_preparation/variable_info")
  load("data/24_7_study/steps/data_preparation/expression_data")
  steps_expression_data = expression_data
  steps_sample_info = sample_info
  steps_variable_info = variable_info
  
}

#####load omics data
load("data/24_7_study/combine_omics/data_preparation/expression_data")
load("data/24_7_study/combine_omics/data_preparation/sample_info")
load("data/24_7_study/combine_omics/data_preparation/variable_info")

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

load(here::here("data/24_7_study/summary_info/day_night_df"))

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

masstools::setwd_project()
setwd("data/24_7_study/wearable_predict_omics/raw_data/")

#####heart rate, step and scg predict omics
###lasso
dim(omics_expression_data)

###match omics with heart rate
###time window is c(5, 10, 20, 30, 40, 50, 60, 90, 120) mins
####

library(moments)
time_windows = c(5, 10, 20, 30, 40, 50, 60, 90, 120)

# new_hr_expression_data =
# purrr::map(time_windows, function(window) {
#   cat(window, " ")
#   temp_data =
#   omics_sample_info$accurate_time %>%
#     purrr::map(function(x) {
#       diff_time = abs(difftime(x, hr_sample_info$accurate_time, unit = "mins"))
#       idx = which(diff_time <= window)
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

#
# new_steps_expression_data =
#   purrr::map(time_windows, function(window) {
#     cat(window, " ")
#     temp_data =
#       omics_sample_info$accurate_time %>%
#       purrr::map(function(x) {
#         diff_time = abs(difftime(x, steps_sample_info$accurate_time, unit = "mins"))
#         idx = which(diff_time <= window)
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

# ###match omics with heart rate
# ###time window is 1, 5, 10, 20, 30, 40, 50, 60 minutes
# ####
# new_cgm_expression_data =
# purrr::map(time_windows, function(window) {
#   cat(window, " ")
#   temp_data =
#   omics_sample_info$accurate_time %>%
#     purrr::map(function(x) {
#       diff_time = abs(difftime(x, cgm_sample_info$accurate_time, unit = "mins"))
#       idx = which(diff_time <= window)
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

# for(i in 1:length(new_hr_expression_data)){
#   cat(i, "\n")
#   temp_hr_data =   new_hr_expression_data[[i]]
#   temp_steps_data =   new_steps_expression_data[[i]]
#   temp_cgm_data =   new_cgm_expression_data[[i]]
#
#   intersect_sample_id =
#     Reduce(intersect, list(
#       colnames(temp_hr_data),
#       colnames(temp_steps_data),
#       colnames(temp_cgm_data)
#     ))
#
#   temp_x =
#     rbind(temp_hr_data[,intersect_sample_id],
#           temp_steps_data[,intersect_sample_id],
#           temp_cgm_data[,intersect_sample_id])
#
#   prediction_result =
#     purrr::map(1:nrow(omics_expression_data), function(j){
#       cat(j, " ")
#       temp_y =
#         as.numeric(omics_expression_data[j,intersect_sample_id])
#
#       temp_x_y =
#         data.frame(data.frame(t(temp_x)), y = temp_y)
#
#       ###7-fold cross validation
#       train.control <- caret::trainControl(method = "cv",
#                                            number = 5)
#
#       model <- train(y ~., data = temp_x_y,
#                      method = "rf",
#                      trControl = train.control)
#
#       result =
#         model$results %>%
#         as.data.frame() %>%
#         dplyr::filter(RMSE == min(RMSE)) %>%
#         head(1) %>%
#         dplyr::mutate(variable_id = omics_variable_info$variable_id[j]) %>%
#         dplyr::select(variable_id, everything())
#       result
#
#     }) %>%
#     dplyr::bind_rows()
#
#   save(prediction_result,
#        file = paste("RF/rf_prediction_result", time_windows[i], "min", sep = "_"))
# }

load("RF/rf_prediction_result_5_min")
rf_prediction_result_5_min = prediction_result
plot(rf_prediction_result_5_min$Rsquared)

load("RF/rf_prediction_result_10_min")
rf_prediction_result_10_min = prediction_result
plot(rf_prediction_result_10_min$Rsquared)

load("RF/rf_prediction_result_20_min")
rf_prediction_result_20_min = prediction_result
plot(rf_prediction_result_20_min$Rsquared)

load("RF/rf_prediction_result_30_min")
rf_prediction_result_30_min = prediction_result
plot(rf_prediction_result_30_min$Rsquared)

load("RF/rf_prediction_result_40_min")
rf_prediction_result_40_min = prediction_result
plot(rf_prediction_result_40_min$Rsquared)

load("RF/rf_prediction_result_50_min")
rf_prediction_result_50_min = prediction_result
plot(rf_prediction_result_50_min$Rsquared)

load("RF/rf_prediction_result_60_min")
rf_prediction_result_60_min = prediction_result
plot(rf_prediction_result_60_min$Rsquared)

load("RF/rf_prediction_result_90_min")
rf_prediction_result_90_min = prediction_result
plot(rf_prediction_result_90_min$Rsquared)

load("RF/rf_prediction_result_120_min")
rf_prediction_result_120_min = prediction_result
plot(rf_prediction_result_120_min$Rsquared)


temp_data =
  data.frame(
    min5 = rf_prediction_result_5_min$Rsquared,
    min10 = rf_prediction_result_10_min$Rsquared,
    min20 = rf_prediction_result_20_min$Rsquared,
    min30 = rf_prediction_result_30_min$Rsquared,
    min40 = rf_prediction_result_40_min$Rsquared,
    min50 = rf_prediction_result_50_min$Rsquared,
    min60 = rf_prediction_result_60_min$Rsquared,
    min90 = rf_prediction_result_90_min$Rsquared,
    min120 = rf_prediction_result_120_min$Rsquared
  )

rownames(temp_data) = rf_prediction_result_30_min$variable_id

library(ComplexHeatmap)

plot =
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    name = "R2",
    column_names_rot = 45,
    border = TRUE
  )

library(ggplotify)
plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "rf_plot.pdf", width = 5, height = 10)

idx =
  temp_data %>%
  apply(1, function(x) {
    sum(x > 0.3)
  }) %>%
  `>`(0) %>%
  which()


#####add annotation
temp_data2 =
  temp_data[idx,]

class =
  variable_info$data_type[match(rownames(temp_data2),
                                variable_info$variable_id)]
non_lipid_idx = which(class != "lipidomics")
non_lipid_name = variable_info$mol_name[match(rownames(temp_data2)[non_lipid_idx],
                                              variable_info$variable_id)]

lipid_name = variable_info$mol_name[match(rownames(temp_data2)[-non_lipid_idx],
                                          variable_info$variable_id)]

row_ha =
  rowAnnotation(
    class = class,
    r2 = anno_boxplot(
      x = as.matrix(temp_data2),
      height = unit(4, "cm"),
      outline = FALSE,
      box_width = 0
    ),
    col = list(class =
                 class_color),
    mark = anno_mark(
      at = non_lipid_idx,
      labels = non_lipid_name,
      labels_gp = gpar(col = class_color[class[non_lipid_idx]],
                       cex = 0.8)
    )
  )

plot =
  Heatmap(
    temp_data2,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    name = "R2",
    column_names_rot = 45,
    border = TRUE,
    right_annotation = row_ha,
    row_split = factor(class, levels = c(
      "metabolic_panel",
      "lipidomics",
      "proteomics"
    ))
  )
plot

row_order =
  ComplexHeatmap::row_order(plot)

align_to = row_order[["lipidomics"]][1:250]

term =
  variable_info$mol_name[match(rownames(temp_data2)[which(class == "lipidomics")],
                               variable_info$variable_id)] %>%
  stringr::str_replace_all("\\.FA", replacement = "_") %>%
  stringr::str_split(pattern = "\\_") %>%
  unlist() %>%
  stringr::str_replace_all("\\.", replacement = "\\:") %>%
  stringr::str_split(pattern = "\\(") %>%
  unlist() %>%
  stringr::str_replace_all("\\)", replacement = "")

library(simplifyEnrichment)

align_to = list(align_to)
term = list(term)
names(align_to) =
  names(term) =
  "test"


plot =
  plot +
  rowAnnotation(
    foo = simplifyEnrichment::anno_word_cloud(
      align_to = align_to,
      term = term,
      max_words = 50,
      fontsize_range = c(5, 15),
      exclude_words = NULL,
      bg_gp = gpar(
        col = class_color["lipidomics"],
        # fill = grey(0.6, 0.3)),
        count_words_param = list(
          remove_numbers = FALSE,
          remove_punctuation = FALSE,
          transform_case = as.character
        ),
        add_new_line = FALSE,
        word_cloud_grob_param = list(max_width = unit(20, "mm"))
      )
    )
  )

plot

# pdf(file = "rf_plot2.pdf", width = 10, height = 7)
# plot
# dev.off()


omics_variable_info %>%
  dplyr::filter(variable_id %in% names(idx)) %>%
  pull(data_type) %>%
  table()

omics_variable_info %>%
  dplyr::filter(variable_id %in% names(idx)) %>%
  pull(mol_name)



#######For lipidomics what are important predictor

temp_hr_data =   new_hr_expression_data[[9]]
temp_steps_data =   new_steps_expression_data[[9]]
temp_cgm_data =   new_cgm_expression_data[[9]]

intersect_sample_id =
  Reduce(intersect, list(
    colnames(temp_hr_data),
    colnames(temp_steps_data),
    colnames(temp_cgm_data)
  ))

temp_x =
  rbind(temp_hr_data[, intersect_sample_id],
        temp_steps_data[, intersect_sample_id],
        temp_cgm_data[, intersect_sample_id])



temp_omics_expression_data =
  omics_expression_data[variable_info$variable_id[match(lipid_name, variable_info$mol_name)],]
#
# lipid_rf_feature_importance =
#   purrr::map(1:nrow(temp_omics_expression_data), function(j){
#     cat(j, " ")
#     temp_y =
#       as.numeric(temp_omics_expression_data[j,intersect_sample_id])
#
#     temp_x_y =
#       data.frame(data.frame(t(temp_x)), y = temp_y)
#
#     library(randomForest)
#     #Random Forest Modelling
#     model <- randomForest(y ~ ., data = temp_x_y, importance = TRUE)
#
#     importance =
#       importance(model) %>%
#       as.data.frame()
#     importance$mean =
#       apply(importance, 1, mean)
#     importance =
#       importance %>%
#       dplyr::arrange(desc(mean))
#     importance
#   })
#
# names(lipid_rf_feature_importance) =
#   rownames(temp_omics_expression_data)
#
# save(lipid_rf_feature_importance,
#      file = "RF/lipid_rf_feature_importance")


load("RF/lipid_rf_feature_importance")
temp_score =
  lipid_rf_feature_importance %>%
  purrr::map(function(x) {
    x =
      x %>%
      tibble::rownames_to_column(var = "feature") %>%
      dplyr::select(feature, mean) %>%
      dplyr::arrange(desc(mean)) %>%
      dplyr::filter(!stringr::str_detect(feature, "hour")) %>%
      head(5) %>%
      pull(feature)
  })  %>%
  unlist() %>%
  unname()

####word cloud

library(ggwordcloud)

temp_score =
  data.frame(temp_score) %>%
  dplyr::count(temp_score) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "hr", "HR")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "steps", "Step")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "cgm", "CGM")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "120_", "")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "_value", "")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "range", "Range")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "max", "Maximum")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "min", "Minumum")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "sd", "SD")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "skewness", "Skewness")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "kurtosis", "Kurtosis")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "mean", "Mean")) %>%
  dplyr::mutate(temp_score =
                  stringr::str_replace(temp_score, "median", "Median")) %>%
  dplyr::mutate(
    class = case_when(
      stringr::str_detect(temp_score, "Step") ~ "step",
      stringr::str_detect(temp_score, "HR") ~ "hr",
      stringr::str_detect(temp_score, "CGM") ~ "cgm"
    )
  )

temp_score$label = paste0(temp_score$temp_score, " (", temp_score$n, ")")

plot =
  temp_score %>%
  ggplot(aes(
    label = label,
    size = n,
    color = class
  )) +
  geom_text_wordcloud_area(shape = "circle", aes(color = class)) +
  # scale_size_area(trans = power_trans(1/.7)) +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  scale_color_manual(values = wearable_color)
plot

# ggsave(plot, filename = "lipid_importance_wordcloud.pdf", width = 10, height = 7)



#####lipid minion for analysis
dir.create("lipidminion")

lipidomics_variable_info =
  variable_info %>%
  dplyr::filter(data_type == "lipidomics") %>%
  dplyr::mutate(Lipid_Name = mol_name)

lipidomics_variable_info$Lipid_Name =
  lipidomics_variable_info$Lipid_Name %>%
  stringr::str_replace_all("\\_", "\\/")

lipidomics_variable_info$Lipid_Name[grep("TAG", lipidomics_variable_info$Lipid_Name)] =
  lipidomics_variable_info$Lipid_Name[grep("TAG", lipidomics_variable_info$Lipid_Name)] %>%
  purrr::map(function(x) {
    # main =
    stringr::str_extract(x, "TAG[0-9]{1,3}\\.[0-9]{1,2}") %>%
      stringr::str_replace("TAG", "") %>%
      stringr::str_split("\\.") %>%
      `[[`(1) %>%
      paste(collapse = ":") %>%
      paste("TAG(", ., ")", sep = "")
  }) %>%
  unlist()

###TAG to TG and DAG to TG
lipidomics_variable_info$Lipid_Name =
  lipidomics_variable_info$Lipid_Name %>%
  stringr::str_replace_all("TAG", "TG") %>%
  stringr::str_replace_all("DAG", "DG")


important_lipid =
  lipidomics_variable_info$Lipid_Name[match(rownames(temp_omics_expression_data),
                                            lipidomics_variable_info$variable_id)] %>%
  data.frame(lipid = .)


# write.table(
#   important_lipid,
#   file = "lipidminion/important_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )


universe_lipid =
  lipidomics_variable_info %>%
  dplyr::select(Lipid_Name) %>%
  dplyr::filter(!is.na(Lipid_Name)) %>%
  dplyr::rename(lipid = Lipid_Name) %>%
  dplyr::distinct(lipid)

# write.table(
#   universe_lipid,
#   file = "lipidminion/universe_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )


result =
  readr::read_delim("lipidminion/Fisher output table (3 pvals _ 0.05).txt",
                    delim = " ")

network_node =
  readr::read_delim("lipidminion/Network_nodes.txt", delim = "\t")
network_edge =
  readr::read_delim("lipidminion/Network_edges.txt", delim = "\t")

network_edge_attr =
  readr::read_delim("lipidminion/Network_edge_attributes.txt", delim = "\t")

library(igraph)
library(ggraph)
library(tidygraph)

network_node

node =
  network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class"))

edge =
  network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

netwrok =
  tidygraph::tbl_graph(nodes = node,
                       edges = edge)

plot =
  ggraph(netwrok,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "black"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 8,
                               "class" = 15)) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(plot, filename = "lipidminion/network.pdf", width = 7, height = 7)

temp_class =
  readr::read_delim("lipidminion/Fisher output table (9 pvals _ 0.05) (1).txt", delim = " ")

temp_class =
  temp_class %>%
  dplyr::filter(Test.performed == "Main class") %>%
  select(Classifier, `%.query`, `%.universe`, `FDR.q-value`) %>%
  dplyr::rename(query = `%.query`,
                all = `%.universe`,
                qvalue = `FDR.q-value`)

library(ggalluvial)

temp_class =
  rbind(temp_class,
        data.frame(
          Classifier = c("PE", "PI", "SM", "CE"),
          query = c(0, 0, 0, 0),
          all = c(23.68, 2.87, 2.87, 6.22),
          qvalue = c(1, 1, 1, 1)
        ))

plot <-
  temp_class %>%
  as.data.frame() %>%
  tidyr::pivot_longer(
    cols = -c(Classifier, qvalue),
    names_to = "sample",
    values_to = "value"
  ) %>%
  ggplot(aes(
    y = value,
    x = sample,
    stratum = Classifier,
    alluvium = Classifier,
    fill = Classifier
  )) +
  geom_alluvium(width = 0.3, alpha = 0.5) +
  geom_stratum(width = 0.5) +
  labs(x = "", y = "Proportion") +
  theme_bw() +
  ggsci::scale_fill_nejm() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

plot

# ggsave(plot, filename = "lipidminion/lipid_main_class.pdf", width = 11, height = 7)


####for protein and metabolic, optimize the random forest model
non_lipid_name

temp_hr_data =   new_hr_expression_data[[9]]
temp_steps_data =   new_steps_expression_data[[9]]
temp_cgm_data =   new_cgm_expression_data[[9]]

intersect_sample_id =
  Reduce(intersect, list(
    colnames(temp_hr_data),
    colnames(temp_steps_data),
    colnames(temp_cgm_data)
  ))

temp_x =
  rbind(temp_hr_data[, intersect_sample_id],
        temp_steps_data[, intersect_sample_id],
        temp_cgm_data[, intersect_sample_id])

temp_omics_expression_data =
  omics_expression_data[variable_info$variable_id[match(non_lipid_name, variable_info$mol_name)],]



non_lipid_prediction_result =
  purrr::map(1:nrow(temp_omics_expression_data), function(j) {
    cat(j, " ")
    temp_y =
      as.numeric(temp_omics_expression_data[j, intersect_sample_id])
    
    temp_x_y =
      data.frame(data.frame(t(temp_x)), y = temp_y)
    
    library(randomForest)
    
    #Random Forest Modelling
    library(Boruta)
    boruta_test <-
      Boruta(y ~ .,
             data = temp_x_y,
             doTrace = 3,
             holdHistory = TRUE)
    
    plot(boruta_test)
    
    marker_rf <-
      boruta_test$finalDecision[boruta_test$finalDecision == "Confirmed"] %>%
      names() %>%
      sort()
    
    no_infi <- function(x)
      all(!is.infinite(x))
    
    plot =
      boruta_test$ImpHistory[1:length(marker_rf), ] %>%
      t() %>%
      as_tibble() %>%
      select_if(., no_infi) %>%
      dplyr::transmute(
        .,
        mean = apply(., 1, mean),
        sd = apply(., 1, sd),
        ymax = mean + sd,
        ymin = mean - sd
      ) %>%
      mutate(name = colnames(boruta_test$ImpHistory)) %>%
      dplyr::filter(!stringr::str_detect(name, "hour")) %>%
      dplyr::filter(!stringr::str_detect(name, "shadowMax")) %>%
      arrange(., mean) %>%
      tail(5) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "hr", "HR")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "steps", "Step")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "cgm", "CGM")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "120_", "")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "_value", "")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "range", "Range")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "max", "Maximum")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "min", "Minumum")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "sd", "SD")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "skewness", "Skewness")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "kurtosis", "Kurtosis")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "mean", "Mean")) %>%
      dplyr::mutate(name =
                      stringr::str_replace(name, "median", "Median")) %>%
      ggplot(., aes(x = factor(name, name), y = mean)) +
      labs(x = "", y = "Importance") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax),
                    colour = "#155F83FF",
                    width = 0) +
      geom_point(size = 2, colour = "#FFA319FF") +
      theme_bw() +
      coord_flip() +
      theme(
        axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 10)
      )
    
    
    ggsave(
      plot,
      filename = file.path("RF/", paste0(
        rownames(temp_omics_expression_data)[j], "_marker_rf.pdf"
      )),
      width = 4,
      height = 7
    )
    
    
    marker_rf <-
      boruta_test$ImpHistory[1:length(marker_rf), ] %>%
      t() %>%
      as_tibble() %>%
      select_if(., no_infi) %>%
      dplyr::transmute(
        .,
        mean = apply(., 1, mean),
        sd = apply(., 1, sd),
        ymax = mean + sd,
        ymin = mean - sd
      ) %>%
      mutate(name = colnames(boruta_test$ImpHistory)) %>%
      filter(name %in% marker_rf) %>%
      dplyr::select(name, everything())
    
    colnames(marker_rf)[-1] <-
      colnames(marker_rf)[-1] %>%
      paste(., "importance", sep = "_")
    
    write.csv(marker_rf,
              file = file.path("RF/", paste0(
                rownames(temp_omics_expression_data)[j], "_marker_rf.csv"
              )),
              row.names = FALSE)
    
    ####parameter tunning
    temp_x_y <-
      temp_x_y %>%
      as.data.frame() %>%
      dplyr::select(y, one_of(marker_rf$name)) %>%
      as.matrix()
    
    ##7-fold cross validation
    train.control <- caret::trainControl(method = "cv",
                                         number = 5)
    
    model <- train(y ~ .,
                   data = temp_x_y,
                   method = "rf",
                   trControl = train.control)
    
    result =
      model$results %>%
      as.data.frame() %>%
      dplyr::filter(RMSE == min(RMSE)) %>%
      head(1) %>%
      dplyr::mutate(variable_id = rownames(temp_omics_expression_data)[j]) %>%
      dplyr::select(variable_id, everything())
    result
  })


variable_info$mol_name[match(rownames(temp_omics_expression_data),
                             variable_info$variable_id)]
