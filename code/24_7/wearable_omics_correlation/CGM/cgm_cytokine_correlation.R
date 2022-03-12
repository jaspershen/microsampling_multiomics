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

#+ fig.fullwidth=TRUE

no_function()

library(tidyverse)
library(here)

rm(list = ls())

source(here::here("R/tools.R"))
source(here::here("R/modified_dtw.R"))
source(here::here("R/lagged_correlation.R"))

{
####load data
###CGM
load(here::here("data/7_24_mike/cgm/data_preparation/sample_info"))
load(here::here("data/7_24_mike/cgm/data_preparation/variable_info"))
load(here::here("data/7_24_mike/cgm/data_preparation/expression_data"))

cgm_expression_data = expression_data
cgm_sample_info = sample_info
cgm_variable_info = variable_info

###cytokine
load(here::here("data/7_24_mike/cytokine/data_preparation/sample_info"))
load(here::here("data/7_24_mike/cytokine/data_preparation/variable_info"))
load(here::here("data/7_24_mike/cytokine/data_preparation/expression_data"))

cytokine_sample_info = sample_info
cytokine_variable_info = variable_info
cytokine_expression_data = expression_data

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
}

######cgm vs cytokine
###work directory
###data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine
sxtTools::setwd_project()
setwd("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between cgm and cytokines
#global correlation
# 
# lagged_cor = rep(NA, nrow(cytokine_expression_data))
# global_cor = rep(NA, nrow(cytokine_expression_data))
# 
# lagged_result = vector(mode = "list", length = nrow(cytokine_expression_data))
# 
# for(i in 1:nrow(cytokine_expression_data)){
#   cat(i, " ")
#   x = as.numeric(cytokine_expression_data[i, ])
#   time1 = cytokine_sample_info$accurate_time
#   y = as.numeric(cgm_expression_data[1, ])
#   time2 = cgm_sample_info$accurate_time
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
# 
# names(lagged_result) = rownames(cytokine_expression_data)
# save(
#   lagged_result,
#   file = here::here(
#     "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
#     "lagged_correlation/lagged_result"
#   )
# )

load(
  here::here(
    "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
    "lagged_correlation/lagged_result"
  )
)

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
  data.frame(wearable = "CGM",
             cytokine_variable_info,
             global_cor = global_cor,
             lagged_cor = lagged_cor,
             shift_time = shift_time)
  # dplyr::filter(abs(lagged_cor) > 0.2)

cor_data$wearable
cor_data$variable_id

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

sum(cor_data$global_cor_p_adjust < 0.05)

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "CGM cytokine global cor",
#              gridLines = TRUE)
# freezePane(wb,
#            sheet = 1,
#            firstRow = TRUE,
#            firstCol = TRUE)
# writeDataTable(
#   wb,
#   sheet = 1,
#   x = cor_data,
#   colNames = TRUE,
#   rowNames = TRUE
# )
# saveWorkbook(
#   wb,
#   file = here::here(
#     "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
#     "lagged_correlation/cor_data.xlsx"
#   ),
#   overwrite = TRUE
# )

cor_data =
  readxl::read_xlsx(
    here::here(
      "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
      "lagged_correlation/cor_data.xlsx"
    )
  )

cor_data = 
  cor_data %>% 
  dplyr::filter(!stringr::str_detect(mol_name, "CHEX"))

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

temp = 
  rbind(neg_top_10,
        pos_top_10) %>% 
  dplyr::filter(lagged_cor_p_adjust < 0.05)

# for (i in 1:nrow(temp)) {
#   cat(i, " ")
#   plot1 =
#   lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                         day_night_df = day_night_df,
#                         internal_omics_color = class_color["cytokine"],
#                         wearable_color = wearable_color["cgm"],
#                         internal_omics_name = temp$mol_name[i],
#                         warable_name = "CGM",
#                         which = "max",
#                         x_limit = c(1,1000),
#                         non_matched_point_size = 0.5,
#                         wearable_point_size = 2,
#                         internal_omics_point_size = 2,
#                         integrated = FALSE,
#                         add_connect_line = FALSE)
# 
#   plot2 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["cytokine"],
#                           wearable_color = wearable_color["cgm"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "CGM",
#                           which = "max",
#                           x_limit = c(1,10),
#                           non_matched_point_size = 0.5,
#                           wearable_point_size = 2,
#                           internal_omics_point_size = 2,
#                           integrated = FALSE)
# 
#   plot3 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["cytokine"],
#                           wearable_color = wearable_color["cgm"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "CGM",
#                           which = "max",
#                           x_limit = c(1,1000),
#                           non_matched_point_size = 0.5,
#                           wearable_point_size = 2,
#                           internal_omics_point_size = 2,
#                           integrated = TRUE,
#                           add_connect_line = FALSE)
# 
#   plot4 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["cytokine"],
#                           wearable_color = wearable_color["cgm"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "CGM",
#                           which = "max",
#                           x_limit = c(1,30),
#                           non_matched_point_size = 3,
#                           wearable_point_size = 3,
#                           internal_omics_point_size = 3,
#                           integrated = TRUE)
# 
#   name = paste("CGM vs",temp$mol_name[i])
# 
#   ggsave(
#     plot1,
#     filename =
#       here::here(
#         "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
#         file.path("cor_plot", paste(name, "plot1.pdf", sep = ""))
#       ),
#     width = 20,
#     height = 7
#   )
# 
#   ggsave(
#     plot2,
#     filename =
#       here::here(
#         "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
#         file.path("cor_plot", paste(name, "plot2.pdf", sep = ""))
#       ),
#     width = 20,
#     height = 7
#   )
# 
#   ggsave(
#     plot3,
#     filename =
#       here::here(
#         "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
#         file.path("cor_plot", paste(name, "plot3.pdf", sep = ""))
#       ),
#     width = 20,
#     height = 7
#   )
# 
#   ggsave(
#     plot4,
#     filename =
#       here::here(
#         "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
#         file.path("cor_plot", paste(name, "plot4.pdf", sep = ""))
#       ),
#     width = 20,
#     height = 7
#   )
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
  shadowtext::geom_shadowtext(aes(label = ifelse(lagged_cor_p_adjust < 0.05 & abs(lagged_cor) > 0.3,
                                                 mol_name, NA)),
                              check_overlap = TRUE,
                              bg.colour='white',
                              color = "black") +
  base_theme

plot

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
#     "global_lagged_correlation.pdf"
#   ),
#   width = 9,
#   height = 7
# )

# plot = 
#   cor_data %>%
#   dplyr::mutate(direction = case_when(
#     shift_time > 0 ~ "After",
#     shift_time < 0 ~ "Before",
#     shift_time == 0 ~ "Synchronization"
#   )) %>%
#   ggplot(aes(shift_time, lagged_cor)) +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   geom_abline(slope = 1, intercept = 0) +
#   geom_point(aes(fill = direction,
#                  size = abs(lagged_cor)),
#              shape = 21, 
#              show.legend = TRUE) +
#   labs(x = "Shift time",
#        y = "Lagged correlation") +
#   scale_fill_manual(values = c(
#     "After" = ggsci::pal_aaas()(n = 10)[1],
#     "Before" = ggsci::pal_aaas()(n = 10)[2],
#     "Synchronization" = "grey"
#   )) +
#   ggrepel::geom_text_repel(aes(label = ifelse(lagged_cor > quantile(cor_data$lagged_cor[cor_data$lagged_cor > 0], 0.75) |
#                                                 lagged_cor < quantile(cor_data$lagged_cor[cor_data$lagged_cor < 0], 0.25),
#                                               mol_name, NA))) +
#   base_theme
# 
# plot
# 
# # ggsave(
# #   plot,
# #   filename = here::here(
# #     "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
# #     "shift_lagged_correlation.pdf"
# #   ),
# #   width = 9,
# #   height = 7
# # )

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
  dplyr::filter(lagged_cor_p_adjust < 0.05)
  # dplyr::filter(
  #   lagged_cor > quantile(cor_data$lagged_cor[cor_data$lagged_cor > 0], 0.75) |
  #     lagged_cor < quantile(cor_data$lagged_cor[cor_data$lagged_cor < 0], 0.25)
  # )

for(i in 1:nrow(important_cytokine)) {
  cat(i, "")
  result = lagged_result[[important_cytokine$variable_id[i]]]

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

  name = important_cytokine$mol_name[i]

  ggsave(
    plot,
    file = here::here(
      "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine",
      file.path("shift_time_vs_cor", paste(name, ".pdf", sep = ""))
    ),
    width = 8,
    height = 7
  )
}



#######peak quality of 
evaluate_peak_quality(object = lagged_result[[4]])

dir.create("peak")

for(i in 1:nrow(important_cytokine)){
  cat(i, " ")
  result =
    evaluate_peak_quality(object = lagged_result[[which(names(lagged_result) == important_cytokine$variable_id[i])]], plot = TRUE)
  name = 
    paste(round(result$score, 4),"_", important_cytokine$mol_name[i], ".pdf",sep = "") %>% 
    stringr::str_replace_all("/", "_")
  ggsave(result$plot, 
         filename = 
           here::here(file.path("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine/peak_quality", name)), 
                      width = 9, height = 7)
}


##Pathway analysis for the negative or positive correlation analysis
###here we used the cytokine minion tool
###out the positive cytokines

library(clusterProfiler)
library(org.Hs.eg.db)

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

# save(
#   important_cytokine,
#   file = here::here(
#     "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine/",
#     "important_cytokine"
#   )
# )

load(
  here::here(
    "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine/",
    "important_cytokine"
  )
)

temp_data = 
  important_cytokine %>% 
  dplyr::mutate(significant = case_when(
    lagged_cor_p_adjust < 0.05 & abs(lagged_cor) > 0.3 ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  dplyr::arrange(lagged_cor) %>% 
  dplyr::mutate(mol_name = factor(mol_name, levels = mol_name))
  
plot = 
  temp_data %>% 
  ggplot(aes(lagged_cor, mol_name)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(size = -log(lagged_cor_p_adjust, 10),
                 color = significant)) +
  geom_text(aes(lagged_cor, 
                mol_name, 
                label = mol_name), 
            hjust = -0.5) +
  labs(x = "Lagged correlation", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  base_theme +
  scale_color_manual(values = c("Yes" = ggsci::pal_aaas()(n=10)[2],
                                "No" = "black")) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
  plot
  
  # ggsave(
  #   plot,
  #   file = here::here(
  #     "data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine/cytokine_order.pdf"
  #   ),
  #   width = 8,
  #   height = 7
  # )

  important_cytokine %>%
    dplyr::filter(lagged_cor_p_adjust < 0.05 &
                    abs(lagged_cor) > 0.3)
    # dplyr::select(mol_name, classification)
  

  
  #####-----------------------------------------------------------------------------
# knitr::spin(
#   hair = here::here(
#     "R/24_7_mike/wearable_omics_correlation/CGM/cgm_cytokine_correlation.R"
#   ),
#   knit = FALSE
# )
# 
# file.copy(from = here::here("R/24_7_mike/wearable_omics_correlation/CGM/cgm_cytokine_correlation.Rmd"),
#           to = here::here("website_files/"), overwrite = TRUE)
