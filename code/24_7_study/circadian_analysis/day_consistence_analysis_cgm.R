no_function()

library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")

#######we first use the raw data to get the consistence score
####load data
{
  ###
  load(here::here("data/24_7_study/cgm/data_preparation/expression_data"))
  load(here::here("data/24_7_study/cgm/data_preparation/sample_info"))
  load(here::here("data/24_7_study/cgm/data_preparation/variable_info"))
}

library(MetaCycle)

setwd("data/24_7_study/circadian_analysis/cgm/day_consistence")

###reference https://abego.cn/2019/05/31/the-rule-of-gene-expression-in-the-day-and-nigth/
###https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html

sample_info %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(n = n())

library(plyr)
seperate_sample_info =
  sample_info %>%
  plyr::dlply(.variables = .(day))

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
        
        match_info = 
        time1$time %>% 
          purrr::map(function(x){
          temp_idx = which(abs(x - time2$time)/60 < 1)
          temp_idx = temp_idx[which.min((abs(x - time2$time)/60)[temp_idx])]
          if(length(temp_idx) == 0){
            temp_idx = NA
          }
          temp_idx
          }) %>% 
          unlist()
                           
        match_info =
          data.frame(sample_id.x = time1$sample_id,
                     sample_id.y = time2$sample_id[match_info]) %>% 
          dplyr::filter(!is.na(sample_id.x) & !is.na(sample_id.y)) %>% 
          dplyr::distinct(sample_id.x, .keep_all = TRUE) %>% 
          dplyr::distinct(sample_id.y, .keep_all = TRUE)
        
        if (nrow(match_info) < 8) {
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

data.frame(x = as.numeric(expression_data["cgm_1", ]), sample_info) %>%
  ggplot(aes(time, x)) +
  geom_point() +
  geom_line() +
  facet_grid(rows = vars(day))

consistence_score = 
data.frame(
  variable_id = names(consistence_score),
  consistence_score = unname(consistence_score)
)

save(consistence_score, file = "consistence_score")







