#' ---
#' title: "CGM protein correlation"
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

# no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())

source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

{
  ####load data
  ###CGM
  load("data/7_24_mike/cgm/data_preparation/sample_info")
  load("data/7_24_mike/cgm/data_preparation/variable_info")
  load("data/7_24_mike/cgm/data_preparation/expression_data")
  
  cgm_expression_data = expression_data
  cgm_sample_info = sample_info
  cgm_variable_info = variable_info
  
  ###proteomics
  load("data/7_24_mike/proteomics/data_preparation/sample_info")
  load("data/7_24_mike/proteomics/data_preparation/variable_info")
  load("data/7_24_mike/proteomics/data_preparation/expression_data")
  
  proteomics_sample_info = sample_info
  proteomics_variable_info = variable_info
  proteomics_expression_data = expression_data
  
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


######cgm vs proteomics
tinyTools::setwd_project()
dir.create("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_proteomics",
           showWarnings = FALSE)
setwd("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_proteomics")

#####---------------------------------------------------------------------------
###correlation between cgm and proteomics
#global correlation

dir.create("lagged_correlation")
lagged_cor = rep(NA, nrow(proteomics_expression_data))
global_cor = rep(NA, nrow(proteomics_expression_data))
lagged_result = vector(mode = "list", length = nrow(proteomics_expression_data))

for(i in 1:nrow(proteomics_expression_data)){
  cat(i, " ")
  x = as.numeric(proteomics_expression_data[i, ])
  time1 = proteomics_sample_info$accurate_time
  y = as.numeric(cgm_expression_data[1, ])
  time2 = cgm_sample_info$accurate_time

  result = lagged_correlation(
    x = x,
    y = y,
    time1 = time1,
    time2 = time2,
    time_tol = 60/60,
    step = 5/60
  )
  lagged_result[[i]] = result
}
names(lagged_result) = rownames(proteomics_expression_data)
save(lagged_result, file = "lagged_correlation/lagged_result")

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
  proteomics_variable_info$variable_id

cor_data =
  data.frame(wearable = "CGM",
             proteomics_variable_info,
             global_cor = global_cor,
             lagged_cor = lagged_cor,
             shift_time = shift_time) %>%
  dplyr::filter(abs(lagged_cor) > 0.2)

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
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
addWorksheet(wb, sheetName = "CGM proteomics global cor",
             gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
writeDataTable(wb, sheet = 1, x = cor_data,
               colNames = TRUE, rowNames = TRUE)
saveWorkbook(wb, "lagged_correlation/cor_data.xlsx", overwrite = TRUE)

cor_data = readxl::read_xlsx("lagged_correlation/cor_data.xlsx")

cor_data %>%
  ggplot(aes(global_cor, lagged_cor)) +
  geom_point()

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
  scale_size_continuous(range = c(0.1,5)) +
  shadowtext::geom_shadowtext(aes(label = ifelse(lagged_cor_p_adjust < 0.05,
                                                 mol_name, NA)),
                              check_overlap = TRUE,
                              bg.colour='white',
                              color = "black") +
  base_theme

plot

# ggsave(plot, filename = "global_lagged_correlation.pdf", width = 9, height = 7)

##output the plot with abs(cor) > 0.3 and p_adjust < 0.05
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
  cor_data %>%
  dplyr::arrange(lagged_cor) %>%
  dplyr::filter(abs(lagged_cor) > 0.3 & lagged_cor_p_adjust < 0.05)

for (i in 1:nrow(temp)) {
  cat(i, " ")
  plot1 =
  lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
                        day_night_df = day_night_df,
                        internal_omics_color = class_color["proteomics"],
                        wearable_color = wearable_color["cgm"],
                        internal_omics_name = temp$mol_name[i],
                        warable_name = "CGM",
                        which = "max",
                        x_limit = c(1,1000),
                        non_matched_point_size = 0.1,
                        wearable_point_size = 0.5,
                        internal_omics_point_size = 2,
                        integrated = FALSE,
                        add_connect_line = TRUE)

  plot2 =
    lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
                          day_night_df = day_night_df,
                          internal_omics_color = class_color["proteomics"],
                          wearable_color = wearable_color["cgm"],
                          internal_omics_name = temp$mol_name[i],
                          warable_name = "CGM",
                          which = "max",
                          x_limit = c(1,10),
                          non_matched_point_size = 0.1,
                          wearable_point_size = 0.5,
                          internal_omics_point_size = 2,
                          integrated = FALSE)

  plot3 =
    lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
                          day_night_df = day_night_df,
                          internal_omics_color = class_color["proteomics"],
                          wearable_color = wearable_color["cgm"],
                          internal_omics_name = temp$mol_name[i],
                          warable_name = "CGM",
                          which = "max",
                          x_limit = c(1,1000),
                          non_matched_point_size = 0.1,
                          wearable_point_size = 0.5,
                          internal_omics_point_size = 2,
                          integrated = TRUE,
                          add_connect_line = FALSE)


  plot4 =
    lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
                          day_night_df = day_night_df,
                          internal_omics_color = class_color["proteomics"],
                          wearable_color = wearable_color["cgm"],
                          internal_omics_name = temp$mol_name[i],
                          warable_name = "CGM",
                          which = "max",
                          x_limit = c(1,30),
                          non_matched_point_size = 3,
                          wearable_point_size = 3,
                          internal_omics_point_size = 3,
                          integrated = TRUE)

  name = paste("CGM vs",temp$mol_name[i])
  ggsave(plot1,
         filename = file.path("cor_plot", paste(name, "_plot1.pdf", sep = "")),
         width = 20, height = 7)

  ggsave(plot2,
         filename = file.path("cor_plot", paste(name, "_plot2.pdf", sep = "")),
         width = 20, height = 7)

  ggsave(plot3,
         filename = file.path("cor_plot", paste(name, "_plot3.pdf", sep = "")),
         width = 20, height = 7)

  ggsave(plot4,
         filename = file.path("cor_plot", paste(name, "_plot4.pdf", sep = "")),
         width = 20, height = 7)
}



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
    proteomics_variable_info$Entry_name[match(temp$variable_id[i], 
                                              proteomics_variable_info$variable_id)]

  ggsave(
    plot,
    file = file.path("shift_time_vs_cor", paste(name, ".pdf", sep = "")),
    width = 8,
    height = 7
  )
}


##Pathway analysis for the negative or positive correlation analysis
##because the correlation is too small, so here we just use the cutoff from
##0.75
##Pathway analysis for the negative or positive correlation analysis
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("pathway_enrichment")

important_protein =
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

# save(important_protein, file = "important_protein")

load("important_protein")


###here we only use the protein with p < 0.05
important_protein = 
  important_protein %>% 
  dplyr::filter(lagged_cor_p_adjust < 0.05)


important_protein
