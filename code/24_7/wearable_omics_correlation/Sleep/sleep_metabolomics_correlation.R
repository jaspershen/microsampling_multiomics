#' ---
#' title: "sleep metabolomics correlation"
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
source(here::here("R/tools.R"))
source(here::here("R/modified_dtw.R"))
source(here::here("R/lagged_correlation.R"))

{
  ####load data
  ###sleep
  load("data/7_24_mike/sleep/data_preparation/sample_info")
  load("data/7_24_mike/sleep/data_preparation/variable_info")
  load("data/7_24_mike/sleep/data_preparation/expression_data")
  
  sleep_expression_data = expression_data
  sleep_sample_info = sample_info
  sleep_variable_info = variable_info
  
  ###metabolomics
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
  
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  metabolomics_expression_data = expression_data
  
  metabolomics_sample_info = 
    metabolomics_sample_info %>% 
    dplyr::filter(!is.na(accurate_time))
  
  metabolomics_expression_data =
    metabolomics_expression_data[,metabolomics_sample_info$sample_id]
  
  plot(density(as.numeric(metabolomics_expression_data[1,])))
  
  ###remove some metabolites
  grep("Caf",metabolomics_variable_info$Compound.name)
  unique(metabolomics_variable_info$Database)
  
  metabolomics_variable_info = 
  metabolomics_variable_info %>% 
    dplyr::filter(Database %in% c("msDatabase_hilic0.0.2", "orbitrapDatabase0.0.1",
                                  "nistDatabase0.0.2", "hmdbDatabase0.0.2", "msDatabase_rplc0.0.2",
                                  "metlinDatabase0.0.2"))
  
  metabolomics_expression_data =
    metabolomics_expression_data[metabolomics_variable_info$variable_id,]
  
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

######sleep vs metabolomics
dir.create("data/7_24_mike/wearable_omics_correlation/sleep_omics_correlation/sleep_metabolomics")
setwd("data/7_24_mike/wearable_omics_correlation/sleep_omics_correlation/sleep_metabolomics")

grep("Caf", metabolomics_variable_info$Compound.name)

plot = 
time_plot(
  x = as.numeric(metabolomics_expression_data[305, ]),
  time = metabolomics_sample_info$accurate_time,
  day_night_df = day_night_df, x_name = "Caffeine",
    add_point = TRUE, x_color = class_color["metabolomics"]
)
plot
ggsave(plot, filename = "caffeine_plot.pdf", width = 20, height = 7)

day_night_df

sleep_variable_info
sleep_expression_data

plot(sleep_sample_info$seconds)

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
  dplyr::mutate(
    from_time = as.POSIXct(from_time),
    end_time = as.POSIXct(end_time)
  )

sleep_data =
  sleep_expression_data %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample_id") %>% 
  dplyr::mutate(sleep_value = case_when(
    sleep1 == "restless" ~ 0,
    sleep1 == "wake" ~ 0,
    sleep1 == "awake" ~ 0,
    sleep1 == "asleep" ~ 1,
    sleep1 == "light" ~ 1,
    sleep1 == "deep" ~ 1,
    sleep1 == "rem" ~ 1
  ))

sleep_section_data =
sleep_section %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    temp_data = 
    sleep_sample_info %>% 
      dplyr::filter(accurate_time >= as.POSIXct(x[2]) & accurate_time <= as.POSIXct(x[3])) %>% 
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

dim(metabolomics_expression_data)
colnames(metabolomics_expression_data)

#####
metabolomics_section_data =
  sleep_section %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    temp_sample_id = 
      metabolomics_sample_info %>% 
      dplyr::filter(accurate_time >= as.POSIXct(x[2]) & accurate_time <= as.POSIXct(x[3])) %>% 
      dplyr::pull(sample_id)
    
    metabolomics_expression_data[,temp_sample_id] %>% 
      apply(1, mean)
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

colnames(metabolomics_section_data) = 1:ncol(metabolomics_section_data)

dim(metabolomics_section_data)
dim(sleep_section_data)

library(plyr)
sleep_section_data2 = 
sleep_section_data %>% 
  plyr::dlply(.variables = .(section)) %>% 
  purrr::map(function(x){
    sleep_score = sum(x$sleep_value * x$sum_seconds)
    x %>% 
      dplyr::mutate(sleep_score = sleep_score) %>% 
      dplyr::select(section, sleep_score) %>% 
      dplyr::distinct(section, .keep_all = TRUE)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()
  
dim(sleep_section_data2)
dim(metabolomics_section_data)

####calculate the correlation between sleep quality and metabolite
result = 
metabolomics_section_data %>% 
  t() %>% 
  as.data.frame() %>% 
  purrr::map(function(x){
    test = 
      cor.test(x, sleep_section_data2$sleep_score, method = "pearson")
    c(cor = unname(test$estimate), p = unname(test$p.value)  )
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(p.adjust = p.adjust(p, method = "BH"))

grep("Caf",metabolomics_variable_info$Compound.name)
grep("Caf",metabolomics_variable_info$Compound.name, value = TRUE)
  
metabolomics_variable_info[305,]

which(rownames(result) == "HILIC_POS_0.90_195.0916m/z")
result[305,]

temp_data = 
  data.frame(sleep = sleep_section_data2$sleep_score,
             caffeine = as.numeric(metabolomics_section_data[305,]),
             day = 1:8)

temp_data %>% 
  ggplot(aes(caffeine, sleep)) +
  geom_point(size = 5) +
  ggrepel::geom_text_repel(aes(label = day)) +
  geom_smooth(method = "lm", se = FALSE) +
  base_theme

result[305,]
