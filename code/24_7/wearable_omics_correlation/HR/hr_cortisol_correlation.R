no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

####load data
###HR
load("data/7_24_mike/hr/data_preparation/sample_info")
load("data/7_24_mike/hr/data_preparation/variable_info")
load("data/7_24_mike/hr/data_preparation/expression_data")

hr_expression_data = expression_data
hr_sample_info = sample_info
hr_variable_info = variable_info

###cortisol
load("data/7_24_mike/cortisol/data_preparation/sample_info")
load("data/7_24_mike/cortisol/data_preparation/variable_info")
load("data/7_24_mike/cortisol/data_preparation/expression_data")

cortisol_sample_info = sample_info
cortisol_variable_info = variable_info
cortisol_expression_data = expression_data

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

######hr vs cortisol
dir.create("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_cortisol")
setwd("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_cortisol")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between hr and cortisols
#global correlation

dir.create("lagged_correlation")

lagged_cor = rep(NA, nrow(cortisol_expression_data))
global_cor = rep(NA, nrow(cortisol_expression_data))

# lagged_result = vector(mode = "list", length = nrow(cortisol_expression_data))
# 
# for(i in 1:nrow(cortisol_expression_data)){
#   cat(i, " ")
#   x = as.numeric(cortisol_expression_data[i, ])
#   time1 = cortisol_sample_info$accurate_time
#   y = as.numeric(hr_expression_data[1, ])
#   time2 = hr_sample_info$accurate_time
# 
#   result = lagged_correlation(
#     x = x,
#     y = y,
#     time1 = time1,
#     time2 = time2,
#     time_tol = 60/60,
#     step = 5/60
#   )
#   lagged_result[[i]] = result
# }
# names(lagged_result) = rownames(cortisol_expression_data)
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
  cortisol_variable_info$variable_id

cor_data =
  data.frame(wearable = "HR",
             cortisol_variable_info,
             global_cor = global_cor,
             lagged_cor = lagged_cor,
             shift_time = shift_time) %>% 
  dplyr::filter(abs(lagged_cor) > 0.2)


p_value = 
  cor_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    x = stringr::str_trim(x, side = "both")
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
# addWorksheet(wb, sheetName = "HR cortisol global cor",
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
#                         internal_omics_color = class_color["cortisol"],
#                         wearable_color = wearable_color["hr"],
#                         internal_omics_name = temp$mol_name[i],
#                         warable_name = "HR",
#                         which = "max",
#                         x_limit = c(1,1000),
#                         non_matched_point_size = 0.1,
#                         wearable_point_size = 0.5,
#                         internal_omics_point_size = 2,
#                         integrated = FALSE)
# 
#   plot2 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["cortisol"],
#                           wearable_color = wearable_color["hr"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "HR",
#                           which = "max",
#                           x_limit = c(1,10),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.5,
#                           internal_omics_point_size = 2,
#                           integrated = FALSE)
# 
#   plot3 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["cortisol"],
#                           wearable_color = wearable_color["hr"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "HR",
#                           which = "max",
#                           x_limit = c(1,1000),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.5,
#                           internal_omics_point_size = 2,
#                           integrated = TRUE)
# 
# 
#   plot4 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["cortisol"],
#                           wearable_color = wearable_color["hr"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "HR",
#                           which = "max",
#                           x_limit = c(1,30),
#                           non_matched_point_size = 3,
#                           wearable_point_size = 3,
#                           internal_omics_point_size = 3,
#                           integrated = TRUE)
# 
#   name = paste("HR vs",temp$mol_name[i])
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

cor_data %>% 
  ggplot(aes(global_cor, lagged_cor)) +
  geom_point()

##cortisol and HR
## output the correlation changes according to shift time
dir.create("shift_time_vs_cor")
# for(i in 1:length(lagged_result)){
#   cat(i, "")
#   result = lagged_result[[i]]
#   x = result$x
#   y = result$max_idx %>%
#     lapply(function(x) {
#       mean(result$y[x])
#     }) %>%
#     unlist()
# 
#   temp_data =
#     result[c("shift_time", "all_cor")] %>%
#     do.call(cbind, .) %>%
#     as.data.frame() %>%
#     dplyr::mutate(shift_time = stringr::str_replace(shift_time, "\\(", "")) %>%
#     dplyr::mutate(shift_time = stringr::str_replace(shift_time, "\\]", "")) %>%
#     dplyr::mutate(shift_time = stringr::str_split(shift_time, ",")) %>%
#     dplyr::mutate(all_cor = round(as.numeric(all_cor), 4))
# 
#   temp_data$shift_time =
#     temp_data$shift_time %>%
#     purrr::map(function(x){
#       mean(as.numeric(x))
#     }) %>%
#     unlist()
# 
#   plot =
#   temp_data %>%
#     ggplot(aes(x = shift_time, y = all_cor)) +
#     geom_vline(xintercept = temp_data$shift_time[result$which_max_idx],
#                color = "red") +
#     geom_hline(yintercept = 0) +
#     annotate(geom = "text",
#              x = temp_data$shift_time[result$which_max_idx],
#              y = result$max_cor, label = result$max_cor) +
#     annotate(geom = "text",
#              x = temp_data$shift_time[result$which_global_idx],
#              y = result$global_cor, label = result$global_cor) +
#     geom_point() +
#     geom_line(aes(group = 1)) +
#     base_theme +
#     labs(x = "Shift time (Omics - HR, min)",
#          y = "Pearsom correlation") +
#     theme()
#   name =
#     cortisol_variable_info$mol_name[match(names(lagged_result)[i], cortisol_variable_info$variable_id)]
# 
#   ggsave(
#     plot,
#     file = file.path("shift_time_vs_cor", paste(name, ".pdf", sep = "")),
#     width = 8,
#     height = 7
#   )
# }

data.frame(time = result$time1,
           x, y) %>% 
  dplyr::mutate(date = as.character(lubridate::as_date(result$time1))) %>% 
  dplyr::filter(!is.na(y)) %>% 
  ggplot(aes(x, y)) +
  geom_point(aes(color = date)) +
  geom_smooth(aes(color = date), method = "lm", se = FALSE)

important_cortisol = 
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

save(important_cortisol, file = "important_cortisol")
