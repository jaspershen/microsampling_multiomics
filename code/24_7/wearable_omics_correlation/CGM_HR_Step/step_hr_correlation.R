no_function()
library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

{
  ####load data
  ###Step
  load("data/7_24_mike/steps/data_preparation/sample_info")
  load("data/7_24_mike/steps/data_preparation/variable_info")
  load("data/7_24_mike/steps/data_preparation/expression_data")
  steps_expression_data = expression_data
  steps_sample_info = sample_info
  steps_variable_info = variable_info
  
  ###HR
  load("data/7_24_mike/hr/data_preparation/sample_info")
  load("data/7_24_mike/hr/data_preparation/variable_info")
  load("data/7_24_mike/hr/data_preparation/expression_data")
  hr_expression_data = expression_data
  hr_sample_info = sample_info
  hr_variable_info = variable_info  
  
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

######steps vs heart rate
dir.create("data/7_24_mike/wearable_omics_correlation/steps_hr_correlation")
setwd("data/7_24_mike/wearable_omics_correlation/steps_hr_correlation")

# ####---------------------------------------------------------------------------
# ####---------------------------------------------------------------------------
# ###correlation between steps and hrs
# #global correlation
# 
# dir.create("lagged_correlation")
# 
# lagged_cor = rep(NA, nrow(hr_expression_data))
# global_cor = rep(NA, nrow(hr_expression_data))
# 
# lagged_result = vector(mode = "list", length = nrow(hr_expression_data))
# 
# for(i in 1:nrow(hr_expression_data)){
#   cat(i, " ")
#   x = as.numeric(steps_expression_data[i, ])
#   time1 = steps_sample_info$accurate_time
#   y = as.numeric(hr_expression_data[1, ])
#   time2 = hr_sample_info$accurate_time
#   
#   result = lagged_correlation(
#     x = x,
#     y = y,
#     time1 = time1,
#     time2 = time2,
#     time_tol = 60/60,
#     step = 1/60
#   )
#   lagged_result[[i]] = result
# }
# names(lagged_result) = rownames(hr_expression_data)
# save(lagged_result, file = "lagged_correlation/lagged_result")

load("lagged_correlation/lagged_result")

lagged_cor = 
  lagged_result %>% 
  purrr::map(function(x){
    x$max_cor
  }) %>% 
  unlist()

global_cor = 
  lagged_result %>% 
  purrr::map(function(x){
    x$global_cor
  }) %>% 
  unlist()

shift_time = 
  lagged_result %>% 
  purrr::map(function(x){
    x$shift_time[x$which_max_idx] %>% 
      stringr::str_replace("\\(", "") %>% 
      stringr::str_replace("\\]", "") %>%
      stringr::str_split(",") %>% 
      `[[`(1) %>% 
      as.numeric() %>% 
      mean()
      
  }) %>% 
  unlist()
  
names(lagged_cor) = names(global_cor) = 
  hr_variable_info$variable_id

cor_data =
  data.frame(wearable = "Step",
             hr_variable_info,
             global_cor = global_cor,
             lagged_cor = lagged_cor,
             shift_time = shift_time) 
  # dplyr::filter(abs(lagged_cor) > 0.2)

p_value = 
  cor_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    # cat(x[2], " ")
    x[!is.na(x)] = stringr::str_trim(x[!is.na(x)], side = "both")
    result = lagged_result[[x[2]]]
    
    ###lagged correlation p value
    x_value = result$x
    y_value = result$y
    
    y_value = 
      result$max_idx %>% 
      purrr::map(function(idx){
        mean(y_value[idx])
      }) %>% 
      unlist()
    
    x_value = x_value[!is.na(y_value)]
    y_value = y_value[!is.na(y_value)]
    lagged_cor_p = 
      cor.test(x = x_value, y = y_value, method = "pearson")$p.value
    
    ###global correlation p value
    x_value = result$x
    y_value = result$y
    
    y_value = 
      result$global_idx %>% 
      purrr::map(function(idx){
        mean(y_value[idx])
      }) %>% 
      unlist()
    
    x_value = x_value[!is.na(y_value)]
    y_value = y_value[!is.na(y_value)]
    global_cor_p = 
      cor.test(x = x_value, y = y_value, method = "pearson")$p.value
    
    c(global_cor_p = global_cor_p,
      lagged_cor_p = lagged_cor_p)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

cor_data = 
  data.frame(cor_data, p_value)

cor_data$global_cor_p_adjust = p.adjust(cor_data$global_cor_p, method = "BH")
cor_data$lagged_cor_p_adjust = p.adjust(cor_data$lagged_cor_p, method = "BH")

library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "Step hr global cor",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = cor_data,
#                colNames = TRUE, rowNames = TRUE)
# saveWorkbook(wb, "lagged_correlation/cor_data.xlsx", overwrite = TRUE)

##output the top 10 negative and top 100 positive
pos_top_10 =
cor_data %>% 
  dplyr::arrange(lagged_cor) %>% 
  dplyr::filter(lagged_cor > 0) %>% 
  tail(10)

neg_top_10 =
  cor_data %>% 
  dplyr::arrange(lagged_cor) %>% 
  dplyr::filter(lagged_cor < 0) %>% 
  head(10)

dir.create("cor_plot")

temp = 
  rbind(neg_top_10,
        pos_top_10)

# for (i in 1:nrow(temp)) {
#   cat(i, " ")
#   plot1 =
#   lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                         day_night_df = day_night_df,
#                         internal_omics_color = wearable_color["hr"],
#                         wearable_color = wearable_color["step"],
#                         internal_omics_name = temp$mol_name[i],
#                         warable_name = "Step",
#                         which = "max",
#                         x_limit = c(1,1000000),
#                         non_matched_point_size = 0.1,
#                         wearable_point_size = 0.1,
#                         internal_omics_point_size = 0.1,
#                         integrated = FALSE,
#                         add_connect_line = FALSE, add_point = FALSE)
# 
#   plot2 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = wearable_color["hr"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
#                           which = "max",
#                           x_limit = c(1,1500),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.5,
#                           internal_omics_point_size = 0.5,
#                           integrated = FALSE,
#                           add_connect_line = FALSE)
# 
#   plot3 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = wearable_color["hr"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
#                           which = "max",
#                           x_limit = c(1,1000000),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.1,
#                           internal_omics_point_size = 0.1,
#                           integrated = TRUE,
#                           add_connect_line = FALSE)
# 
# 
#   plot4 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = wearable_color["hr"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
#                           which = "max",
#                           x_limit = c(1,1500),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.5,
#                           internal_omics_point_size = 0.5,
#                           integrated = TRUE,
#                           add_connect_line = FALSE)
# 
#   name = paste("Step vs",temp$mol_name[i])
#   ggsave(plot1,
#          filename = file.path("cor_plot", paste(name, "plot1.pdf", sep = "")),
#          width = 20, height = 7)
# 
#   ggsave(plot2,
#          filename = file.path("cor_plot", paste(name, "plot2.pdf", sep = "")),
#          width = 20, height = 7)
# 
#   ggsave(plot3,
#          filename = file.path("cor_plot", paste(name, "plot3.pdf", sep = "")),
#          width = 20, height = 7)
# 
#   ggsave(plot4,
#          filename = file.path("cor_plot", paste(name, "plot4.pdf", sep = "")),
#          width = 20, height = 7)
# }


##Pathway analysis for the negative or positive correlation analysis
##because the correlation is too small, so here we just use the cutoff from
##0.75
##Pathway analysis for the negative or positive correlation analysis
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("pathway_enrichment")

important_step_hr =
  rbind(
    cor_data %>%
      dplyr::filter(lagged_cor > 0) %>%
      # dplyr::filter(lagged_cor > quantile(lagged_cor, 0.75)) %>%
      dplyr::mutate(class1 = "positive correlation"),
    cor_data %>%
      dplyr::filter(lagged_cor < 0) %>%
      # dplyr::filter(lagged_cor < quantile(lagged_cor, 0.25)) %>%
      dplyr::mutate(class1 = "negative correlation")
  )

save(important_step_hr, file = "important_step_hr")






####output the shift time vs lagged cor plot for the important proteins
dir.create("shift_time_vs_cor")

for(i in 1:nrow(temp)){
  cat(i, "")
  
  result = lagged_result[[temp$variable_id[i]]]
  
  temp_data =
    result[c("shift_time", "all_cor", "all_cor_p")] %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    dplyr::mutate(shift_time = stringr::str_replace(shift_time, "\\(", "")) %>%
    dplyr::mutate(shift_time = stringr::str_replace(shift_time, "\\]", "")) %>%
    dplyr::mutate(shift_time = stringr::str_split(shift_time, ",")) %>%
    dplyr::mutate(all_cor = round(as.numeric(all_cor), 4))
  
  temp_data$shift_time =
    temp_data$shift_time %>%
    purrr::map(function(x){
      mean(as.numeric(x))
    }) %>%
    unlist()
  
  temp_data$all_cor_p = as.numeric(temp_data$all_cor_p)
  temp_data$log_p = -log(temp_data$all_cor_p, 10)
  temp_data$log_p[is.infinite(temp_data$log_p)] = max(temp_data$log_p[!is.infinite(temp_data$log_p)])
  temp_data = 
    temp_data %>% 
    dplyr::mutate(yesornot = 
                    case_when(
                      all_cor_p < 0.05 ~ "yes",
                      all_cor_p >= 0.05 ~ "no"
                    ))
  
  plot =
    temp_data %>%
    ggplot(aes(x = shift_time, y = all_cor)) +
    geom_vline(xintercept = temp_data$shift_time[result$which_max_idx],
               color = "red") +
    geom_hline(yintercept = 0) +
    ggplot2::annotate(geom = "text",
                      x = temp_data$shift_time[result$which_max_idx],
                      y = result$max_cor,
                      label = round(result$max_cor, 4)) +
    ggplot2::annotate(geom = "text",
                      x = temp_data$shift_time[result$which_global_idx],
                      y = result$global_cor,
                      label = round(result$global_cor, 4)) +
    scale_color_manual(values = c("yes" = "red", "no" = "grey")) +
    geom_line(aes(group = 1)) +
    geom_point(aes(size = log_p, color = yesornot), show.legend = FALSE) +
    base_theme +
    labs(x = "Shift time (Omics - CGM, min)",
         y = "Pearsom correlation") +
    theme()
  
  name =
   temp$mol_name[i]
  
  ggsave(
    plot,
    file = file.path("shift_time_vs_cor", paste(name, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
}

