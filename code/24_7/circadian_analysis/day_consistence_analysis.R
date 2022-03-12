no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

#######we first use the raw data to get the consistence score
####load data
{
  ###
  load(here::here(
    "data/7_24_mike/combine_omics/data_preparation/expression_data"
  ))
  load(here::here(
    "data/7_24_mike/combine_omics/data_preparation/sample_info"
  ))
  load(here::here(
    "data/7_24_mike/combine_omics/data_preparation/variable_info"
  ))
}

library(MetaCycle)

setwd("data/7_24_mike/circadian_analysis/all_omics/day_consistence")

###reference https://abego.cn/2019/05/31/the-rule-of-gene-expression-in-the-day-and-nigth/
###https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html


data.frame(x = as.numeric(expression_data["HILIC_POS_2.87_127.0727m/z", ]), sample_info) %>%
  dplyr::filter(day != "2019-05-07") %>%
  ggplot(aes(time, x)) +
  geom_point() +
  geom_line() +
  facet_grid(rows = vars(day))

sample_info %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(n = n())

seperate_sample_info =
  sample_info %>%
  plyr::dlply(.variables = .(day))


####match_info
temp =   
purrr::map(
    1:(length(seperate_sample_info) - 1),
    .f = function(idx1) {
      cat(idx1, " ")
      purrr::map((idx1 + 1):length(seperate_sample_info), function(idx2) {
        time1 = 
          seperate_sample_info[[idx1]][, c("sample_id", "time")]
        time2 = 
          seperate_sample_info[[idx2]][, c("sample_id", "time")]
        
        matched_idx2 = 
          purrr::map(time1$time, function(x) {
            diff_time = abs(as.numeric(difftime(x, time2$time, units = "min")))
            temp_idx = which(diff_time <= 30)
            temp_idx = temp_idx[which.min(diff_time[temp_idx])]
            if(length(temp_idx) == 0){
              temp_idx = NA
            }
            temp_idx
          }) %>% 
          unlist()
        
        match_info =
          data.frame(matched_idx1 = 1:length(matched_idx2),
                     matched_idx2 = matched_idx2) %>%
          dplyr::mutate(sample_id.x = time1$sample_id[matched_idx1],
                        sample_id.y = time2$sample_id[matched_idx2]) %>%
          dplyr::filter(!is.na(matched_idx2))
        
        match_info
    }
  )
    })



variable_cor_data =
  purrr::map(
    1:(length(seperate_sample_info) - 1),
    .f = function(idx1) {
      cat(idx1, " ")
      purrr::map((idx1 + 1):length(seperate_sample_info), function(idx2) {
        time1 = 
          seperate_sample_info[[idx1]][, c("sample_id", "time")]
        time2 = 
          seperate_sample_info[[idx2]][, c("sample_id", "time")]
        
        matched_idx2 = 
        purrr::map(time1$time, function(x) {
          diff_time = abs(as.numeric(difftime(x, time2$time, units = "min")))
          temp_idx = which(diff_time <= 30)
          temp_idx = temp_idx[which.min(diff_time[temp_idx])]
          if(length(temp_idx) == 0){
            temp_idx = NA
          }
          temp_idx
        }) %>% 
          unlist()
        
        match_info =
          data.frame(matched_idx1 = 1:length(matched_idx2),
                     matched_idx2 = matched_idx2) %>%
          dplyr::mutate(sample_id.x = time1$sample_id[matched_idx1],
                        sample_id.y = time2$sample_id[matched_idx2]) %>%
          dplyr::filter(!is.na(matched_idx2))
        
        # match_info =
        #   seperate_sample_info[[idx1]][, c("sample_id", "hour")] %>%
        #   dplyr::left_join(seperate_sample_info[[idx2]][, c("sample_id", "hour")],
        #                    by = "hour") %>%
        #   dplyr::filter(!is.na(sample_id.x) & !is.na(sample_id.y))
        
        if (nrow(match_info) < 5) {
          return(NULL)
        }
        
        temp_data1 =
          expression_data[, match_info$sample_id.x]
        
        temp_data2 =
          expression_data[, match_info$sample_id.y]
        
        cor_result =
          purrr::map(1:nrow(temp_data1), function(idx) {
            result =
              cor.test(as.numeric(temp_data1[idx, ]),
                       as.numeric(temp_data2[idx, ]),
                       method = "spearman")
            c(cor = unname(result$estimate),
              p = unname(result$p.value))
          })  %>%
          dplyr::bind_rows()
        
        cor_result =
          data.frame(
            variable_id = rownames(expression_data),
            day1 = seperate_sample_info[[idx1]]$day[1],
            day2 = seperate_sample_info[[idx2]]$day[1],
            cor_result
          )
      }) %>%
        dplyr::bind_rows()
    }
  ) %>%
  dplyr::bind_rows()


unique(variable_cor_data$variable_id)

library(plyr)

consistence_score =
  variable_cor_data %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(function(x) {
    median(x$cor)
  }) %>%
  unlist()

plot(consistence_score)

which(consistence_score > 0.6)

unique(c(variable_cor_data$day1, variable_cor_data$day2))

data.frame(x = as.numeric(expression_data["HILIC_POS_2.87_127.0727m/z", ]), sample_info) %>%
  dplyr::filter(day != "2019-05-07") %>%
  ggplot(aes(time, x)) +
  geom_point() +
  geom_line() +
  facet_grid(rows = vars(day))

consistence_score =
  data.frame(
    variable_id = names(consistence_score),
    consistence_score = unname(consistence_score)
  )


#####used the loess predicted value to get the consistence socre
####load data
{
  ###
  load(here::here(
    "data/7_24_mike/combine_omics/data_preparation/new_expression_data"
  ))
  load(here::here(
    "data/7_24_mike/combine_omics/data_preparation/new_sample_info"
  ))
  load(here::here(
    "data/7_24_mike/combine_omics/data_preparation/new_variable_info"
  ))
}

library(MetaCycle)

###reference https://abego.cn/2019/05/31/the-rule-of-gene-expression-in-the-day-and-nigth/
###https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html

new_sample_info %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(n = n())

new_seperate_sample_info =
  new_sample_info %>%
  plyr::dlply(.variables = .(day))

new_variable_cor_data =
  purrr::map(
    1:(length(new_seperate_sample_info) - 1),
    .f = function(idx1) {
      cat(idx1, " ")
      purrr::map((idx1 + 1):length(new_seperate_sample_info), function(idx2) {
        match_info =
          new_seperate_sample_info[[idx1]][, c("sample_id", "time")] %>%
          dplyr::left_join(new_seperate_sample_info[[idx2]][, c("sample_id", "time")],
                           by = "time") %>%
          dplyr::filter(!is.na(sample_id.x) & !is.na(sample_id.y))
        
        if (nrow(match_info) < 15) {
          return(NULL)
        }
        
        temp_data1 =
          new_expression_data[, match_info$sample_id.x]
        
        temp_data2 =
          new_expression_data[, match_info$sample_id.y]
        
        cor_result =
          purrr::map(1:nrow(temp_data1), function(idx) {
            result =
              cor.test(as.numeric(temp_data1[idx, ]),
                       as.numeric(temp_data2[idx, ]),
                       method = "spearman")
            c(cor = unname(result$estimate),
              p = unname(result$p.value))
          })  %>%
          dplyr::bind_rows()
        
        cor_result =
          data.frame(
            variable_id = rownames(new_expression_data),
            day1 = new_seperate_sample_info[[idx1]]$day[1],
            day2 = new_seperate_sample_info[[idx2]]$day[1],
            cor_result
          )
      }) %>%
        dplyr::bind_rows()
    }
  ) %>%
  dplyr::bind_rows()


unique(new_variable_cor_data$variable_id)

library(plyr)

new_consistence_score =
  new_variable_cor_data %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(function(x) {
    median(x$cor)
  }) %>%
  unlist()

plot(new_consistence_score)

which(new_consistence_score > 0.7)

unique(c(new_variable_cor_data$day1, new_variable_cor_data$day2))

data.frame(x = as.numeric(new_expression_data["lipid_135",]), new_sample_info) %>%
  dplyr::filter(day != "2019-05-07") %>%
  ggplot(aes(time, x)) +
  geom_point() +
  geom_line() +
  facet_grid(rows = vars(day))

new_consistence_score = 
  data.frame(
    variable_id = names(new_consistence_score),
    new_consistence_score = unname(new_consistence_score)
  )

dim(consistence_score)
dim(new_consistence_score)

consistence_score = 
consistence_score %>% 
  dplyr::left_join(new_consistence_score, by = "variable_id")

plot =
  consistence_score %>%
  ggplot(aes(consistence_score, new_consistence_score)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = "red") +
  geom_point() +
  base_theme +
  labs(x = "Consistence score (raw data)",
       y = "Consistence score (predicted data)")

ggsave(plot, filename = "consistence_plot.pdf", width = 7, height = 7)

consistence_score$consistence_score[is.na(consistence_score$consistence_score)] = 0

which(is.na(consistence_score$consistence_score))
which(is.na(consistence_score$new_consistence_scoreconsistence_score))

cor(consistence_score$consistence_score, 
    consistence_score$new_consistence_score,
    method = "pearson")

save(consistence_score, file = "consistence_score")







