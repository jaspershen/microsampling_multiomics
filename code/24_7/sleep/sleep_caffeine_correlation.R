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

#####load the raw data from Ryan
load("data/7_24_mike/all_omes")

temp_data =
  all_omes[which(all_omes$MolName == "Caffeine"), ]

source(here::here("R/tools.R"))

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
    dplyr::mutate(week = paste(week,
                               lubridate::month(day),
                               lubridate::day(day),
                               sep = "-")) %>%
    dplyr::mutate(week = factor(week, unique(week)))
}

####load omics data

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
    which(apply(expression_data, 1, function(x) {
      sum(x == 0) / ncol(expression_data)
    }) < 0.5)
  
  metabolic_panel_variable_info =
    metabolic_panel_variable_info[remain_idx, ]
  
  metabolic_panel_expression_data =
    metabolic_panel_expression_data[metabolic_panel_variable_info$variable_id, ]
  
  ###remove CHEX
  metabolic_panel_variable_info =
    metabolic_panel_variable_info %>%
    dplyr::filter(!stringr::str_detect(mol_name, "CHEX"))
  
  metabolic_panel_expression_data =
    metabolic_panel_expression_data[metabolic_panel_variable_info$variable_id, ]
  
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
    dplyr::mutate(week = paste(week,
                               lubridate::month(day),
                               lubridate::day(day),
                               sep = "-")) %>%
    dplyr::mutate(week = factor(week, unique(week)))
  
}


temp_data =
  temp_data %>%
  dplyr::select(SampleID, Intensity) %>%
  dplyr::left_join(sample_info, by = c("SampleID" = "sample_id"))

temp_data %>%
  ggplot(aes(accurate_time, Intensity)) +
  geom_point() +
  geom_line()


######sleep vs omics
dir.create("data/7_24_mike/sleep/sleep_omics", recursive = TRUE)
setwd("data/7_24_mike/sleep/sleep_omics")

library(metid)

load(here::here("data/7_24_mike/metabolomics/orbitrapDatabase0.0.1"))

plot =
  time_plot(
    x = temp_data$Intensity,
    time = temp_data$accurate_time,
    day_night_df = day_night_df,
    x_name = "Caffeine",
    add_point = TRUE,
    x_color = class_color["metabolomics"],
    y_name = "Caffeine"
  )

plot

####read metadata
metdata =
  readr::read_csv(here::here("data/7_24_mike/raw_data_from_box/sample_registration.csv"))

coffee_time =
  metdata %>%
  dplyr::filter(!is.na(food)) %>%
  dplyr::filter(stringr::str_detect(food, "coffee|Coffee")) %>%
  pull(date_time)

plot =
  plot +
  geom_vline(xintercept = lubridate::as_datetime(coffee_time, tz = "PDT"))

plot

# ggsave(plot, filename = "caffeine_plot.pdf", width = 20, height = 7)

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
  dplyr::mutate(from_time = as.POSIXct(from_time),
                end_time = as.POSIXct(end_time))

sleep_data =
  sleep_expression_data %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::mutate(
    sleep_value = case_when(
      sleep1 == "restless" ~ 0,
      sleep1 == "wake" ~ 1,
      sleep1 == "awake" ~ 2,
      sleep1 == "asleep" ~ 3,
      sleep1 == "light" ~ 4,
      sleep1 == "deep" ~ 5,
      sleep1 == "rem" ~ 6
    )
  )

sleep_section_data =
  sleep_section %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    temp_data =
      sleep_sample_info %>%
      dplyr::filter(accurate_time >= as.POSIXct(x[2]) &
                      accurate_time <= as.POSIXct(x[3])) %>%
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

library(plyr)
sleep_section_data2 =
  sleep_section_data %>%
  plyr::dlply(.variables = .(section)) %>%
  purrr::map(function(x) {
    sleep_score = sum(x$sleep_value * x$sum_seconds)
    x %>%
      dplyr::mutate(sleep_score = sleep_score) %>%
      dplyr::select(section, sleep_score) %>%
      dplyr::distinct(section, .keep_all = TRUE)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

#####
omics_section_data =
  sleep_section %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    temp_data %>%
      dplyr::filter(accurate_time >= as.POSIXct(x[2]) &
                      accurate_time <= as.POSIXct(x[3])) %>%
      pull(Intensity) %>%
      median()
  }) %>%
  unlist()

plot(omics_section_data, sleep_section_data2$sleep_score)

cor.test(omics_section_data, sleep_section_data2$sleep_score, method = "spearman")

####calculate the correlation between sleep quality and metabolite
omics_section_data

temp_data =
  data.frame(
    sleep = sleep_section_data2$sleep_score,
    caffeine = as.numeric(omics_section_data),
    day = 1:8
  )

temp_data %>%
  ggplot(aes(caffeine, sleep)) +
  geom_point(size = 5, color = class_color["metabolomics"]) +
  ggrepel::geom_text_repel(aes(label = day)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  base_theme +
  labs(y = "Sleep score",
       x = "Caffenine")
