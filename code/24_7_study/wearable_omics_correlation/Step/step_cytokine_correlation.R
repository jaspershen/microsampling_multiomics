no_function()

library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")
source("code/modified_dtw.R")
source("code/lagged_correlation.R")

####load data
###Step
load("data/24_7_study/steps/data_preparation/sample_info")
load("data/24_7_study/steps/data_preparation/variable_info")
load("data/24_7_study/steps/data_preparation/expression_data")
steps_expression_data = expression_data
steps_sample_info = sample_info
steps_variable_info = variable_info

###cytokine
load("data/24_7_study/cytokine/data_preparation/sample_info")
load("data/24_7_study/cytokine/data_preparation/variable_info")
load("data/24_7_study/cytokine/data_preparation/expression_data")
cytokine_sample_info = sample_info
cytokine_variable_info = variable_info
cytokine_expression_data = expression_data

load("data/24_7_study/summary_info/day_night_df")

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

######steps vs cytokine
dir.create("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_cytokine")
setwd("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_cytokine")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between steps and cytokines
#global correlation

dir.create("lagged_correlation")

lagged_cor = rep(NA, nrow(cytokine_expression_data))
global_cor = rep(NA, nrow(cytokine_expression_data))

lagged_result = vector(mode = "list", length = nrow(cytokine_expression_data))

# for(i in 1:nrow(cytokine_expression_data)){
#   cat(i, " ")
#   x = as.numeric(cytokine_expression_data[i, ])
#   time1 = cytokine_sample_info$accurate_time
#   y = as.numeric(steps_expression_data[1, ])
#   time2 = steps_sample_info$accurate_time
# 
#   result = lagged_correlation(
#     x = x,
#     y = y,
#     time1 = time1,
#     time2 = time2,
#     time_tol = 60/60,
#     step = 5/60
#   )
# 
#   lagged_cor[i] = result$max_cor
#   global_cor[i] = result$global_cor
#   lagged_result[[i]] = result
# }
# names(lagged_result) = rownames(cytokine_expression_data)
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
  cytokine_variable_info$variable_id

cor_data =
  data.frame(wearable = "Step",
             cytokine_variable_info,
             global_cor = global_cor,
             lagged_cor = lagged_cor,
             shift_time = shift_time) 
  # dplyr::filter(abs(lagged_cor) > 0.2)


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
# addWorksheet(wb, sheetName = "Step cytokine global cor",
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
#                         internal_omics_color = class_color["cytokine"],
#                         wearable_color = wearable_color["step"],
#                         internal_omics_name = temp$mol_name[i],
#                         warable_name = "Step",
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
#                           internal_omics_color = class_color["cytokine"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
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
#                           internal_omics_color = class_color["cytokine"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
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
#                           internal_omics_color = class_color["cytokine"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
#                           which = "max",
#                           x_limit = c(1,30),
#                           non_matched_point_size = 3,
#                           wearable_point_size = 3,
#                           internal_omics_point_size = 3,
#                           integrated = TRUE)
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



##remove CHEX
cor_data = 
  cor_data %>% 
  dplyr::filter(!stringr::str_detect(mol_name, "CHEX"))

plot =
  cor_data %>%
  dplyr::mutate(direction = case_when(
    shift_time > 0 ~ "After",
    shift_time < 0 ~ "Before",
    shift_time == 0 ~ "Synchronization"
  )) %>%
  ggplot(aes(global_cor, lagged_cor)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(fill = direction,
                 size = -log(lagged_cor_p_adjust, 10)),
             shape = 21) +
  labs(x = "Global Pearson correlation",
       y = "Lagged correlation") +
  scale_fill_manual(values = c(
    "After" = ggsci::pal_aaas()(n = 10)[1],
    "Before" = ggsci::pal_aaas()(n = 10)[2],
    "Synchronization" = "grey"
  )) +
  guides(size = guide_legend(title = "-log10(p.adjust)")) +
  scale_size_continuous(range = c(1,10)) +
  shadowtext::geom_shadowtext(aes(label = ifelse(lagged_cor_p_adjust < 0.05,
                                                 mol_name, NA)),
                              check_overlap = TRUE,
                              bg.colour='white',
                              color = "black") +
  base_theme

plot

# ggsave(plot, filename = "global_lagged_correlation.pdf", width = 9, height = 7)

plot = 
  cor_data %>%
  dplyr::mutate(direction = case_when(
    shift_time > 0 ~ "After",
    shift_time < 0 ~ "Before",
    shift_time == 0 ~ "Synchronization"
  )) %>%
  ggplot(aes(shift_time, lagged_cor)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(fill = direction,
                 size = abs(lagged_cor)),
             shape = 21, 
             show.legend = TRUE) +
  labs(x = "Shift time",
       y = "Lagged correlation") +
  scale_fill_manual(values = c(
    "After" = ggsci::pal_aaas()(n = 10)[1],
    "Before" = ggsci::pal_aaas()(n = 10)[2],
    "Synchronization" = "grey"
  )) +
  ggrepel::geom_text_repel(aes(label = ifelse(lagged_cor > quantile(cor_data$lagged_cor[cor_data$lagged_cor > 0], 0.75) |
                                                lagged_cor < quantile(cor_data$lagged_cor[cor_data$lagged_cor < 0], 0.25),
                                              mol_name, NA))) +
  base_theme

plot

# ggsave(plot, filename = "shift_lagged_correlation.pdf", width = 9, height = 7)

cor_data %>% 
  dplyr::filter(lagged_cor == max(lagged_cor))

cor_data %>% 
  dplyr::filter(lagged_cor == min(lagged_cor))

cor_data %>% 
  dplyr::select(mol_name, classification, lagged_cor, shift_time) %>% 
  dplyr::arrange(lagged_cor)

cor_data %>% 
  dplyr::arrange(lagged_cor) %>% 
  dplyr::mutate(mol_name = factor(mol_name, levels = mol_name)) %>% 
  ggplot(aes(lagged_cor, mol_name)) +
  geom_point() +
  labs(x = "Lagged correlation", y = "") +
  base_theme


####output the shift time vs lagged cor plot for the important cytokines
important_cytokine =
  cor_data %>%
  dplyr::filter(lagged_cor > quantile(cor_data$lagged_cor[cor_data$lagged_cor > 0], 0.75) |
                  lagged_cor < quantile(cor_data$lagged_cor[cor_data$lagged_cor < 0], 0.25))

dir.create("shift_time_vs_cor")

# for(i in 1:nrow(important_cytokine)) {
#   cat(i, "")
#   result = lagged_result[[important_cytokine$variable_id[i]]]
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
#     purrr::map(function(x) {
#       mean(as.numeric(x))
#     }) %>%
#     unlist()
# 
#   plot =
#     temp_data %>%
#     ggplot(aes(x = shift_time, y = all_cor)) +
#     geom_vline(xintercept = temp_data$shift_time[result$which_max_idx],
#                color = "red") +
#     geom_hline(yintercept = 0) +
#     annotate(
#       geom = "text",
#       x = temp_data$shift_time[result$which_max_idx],
#       y = result$max_cor,
#       label = round(result$max_cor, 4)
#     ) +
#     annotate(
#       geom = "text",
#       x = temp_data$shift_time[result$which_global_idx],
#       y = result$global_cor,
#       label = round(result$global_cor, 4)
#     ) +
#     geom_point() +
#     geom_line(aes(group = 1)) +
#     base_theme +
#     labs(x = "Shift time (Omics - Step, min)",
#          y = "Pearsom correlation") +
#     theme()
# 
#   name =
#     cytokine_variable_info$mol_name[match(names(lagged_result)[i],
#                                               cytokine_variable_info$variable_id)]
# 
#   ggsave(
#     plot,
#     file = file.path("shift_time_vs_cor", paste(name, ".pdf", sep = "")),
#     width = 8,
#     height = 7
#   )
# }


##Pathway analysis for the negative or positive correlation analysis
###here we used the cytokineminion tool
###out the positive cytokines

library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("pathway_enrichment")

important_cytokine =
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

save(important_cytokine, file = "important_cytokine")
