no_function()
library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")
source("code/modified_dtw.R")
source("code/lagged_correlation.R")

###data loading
####load all omics data loess data
load(here::here("data/24_7_study/combine_omics/data_preparation/new_expression_data"))
load(here::here("data/24_7_study/combine_omics/data_preparation/new_sample_info"))
load(here::here("data/24_7_study/combine_omics/data_preparation/new_variable_info"))

expression_data = new_expression_data
sample_info = new_sample_info
variable_info = new_variable_info

######
dir.create("data/24_7_study/inter_omics_correlation")
dir.create("data/24_7_study/inter_omics_correlation/inter_all_omics_loess_data")
setwd("data/24_7_study/inter_omics_correlation/inter_all_omics_loess_data")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between metabolic_panels
#global correlation

###get the matched index
time1 = sample_info$accurate_time
time2 = sample_info$accurate_time
x = as.numeric(expression_data[1,])
y = as.numeric(expression_data[1,])
time_tol = 300 / 60
step = 30 / 60

time_window1 = 
  seq(from = step/2, to = time_tol, by = step)

time_window2 = 
  -rev(seq(from = step/2, to = time_tol, by = step))

time_window = sort(c(time_window2, time_window1))

shift_time = 
  paste("(",
        paste(round(time_window[-length(time_window)] * 60, 2),
              round(time_window[-1] * 60, 2), sep = ','),
        "]", sep = ""
  )


temp_fun = 
  function(temp_idx,
           time_window,
           x,
           y,
           time1,
           time2
  ){
    idx =
      time1 %>%
      purrr::map(function(x) {
        diff_time =
          difftime(x, time2, units = "hours")
        which(diff_time > time_window[temp_idx] &
                diff_time <= time_window[temp_idx + 1])
      })
  }

bpparam = BiocParallel::MulticoreParam(workers = 10, 
                                       progressbar = TRUE)

all_idx = 
    BiocParallel::bplapply(X = 1:(length(time_window) - 1), 
                           FUN = temp_fun, 
                           time_window = time_window,
                           x = x,
                           y = y,
                           time1 = time1,
                           time2 = time2, 
                           BPPARAM = bpparam)      

all_idx

shift_time = 
  lapply(1:(length(time_window) - 1), FUN = function(idx){
    mean(time_window[c(idx, idx+1)])
  }) %>% 
  unlist() %>% 
  `*`(60)

all_idx %>%
  lapply(function(x) {
    length(unique(unlist(x)))
  }) %>%
  unlist() %>%
  data.frame(number = .) %>%
  dplyr::mutate(index = shift_time) %>%
  ggplot(aes(index, number)) +
  geom_point(size = 3) +
  base_theme +
  geom_text(aes(x = index, y = number, label = number),
            nudge_x = 12,
            nudge_y = 1) +
  labs(x = "Shift time (min)", y = "Matched sample number")

temp_data = 
purrr::map2(.x = all_idx, 
            .y = shift_time,
              .f = function(idx,
                            t){
  purrr::map2(
    time1,
    idx,
    .f = function(x, y) {
      difftime(time1 = x, time2 = time1[y], units = "hour")
    }
  ) %>% 
    unlist() %>% 
    data.frame(time = ., shift_time = t) %>% 
                  dplyr::mutate(time = time * 60)
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

library(gghalves)

plot = 
temp_data %>% 
  dplyr::mutate(shift_time = as.character(shift_time)) %>% 
  dplyr::mutate(shift_time = factor(shift_time, 
                                    levels = as.character(unique(temp_data$shift_time)))) %>% 
  ggplot(aes(x=shift_time, time)) +
  geom_half_boxplot(aes(x=shift_time, time),
                    fill = "transparent",
                    show.legend = FALSE, side = "l") +
  geom_half_violin(aes(x=shift_time, time, color = shift_time), 
                   show.legend = FALSE, side = "r") +
  geom_half_point(aes(x=shift_time, time, color = shift_time), 
                  show.legend = FALSE, side = "l", alpha = 0.7) +
  ggsci::scale_color_lancet() +
  base_theme +
  labs(y = "Different time (min)", x = "Shift time (min)")

plot

# ggsave(plot, filename = "match_plot.pdf", width = 9, height = 7)
# 
# dir.create("lagged_correlation")
# total_number = nrow(expression_data)
# 
# lagged_result =
#   purrr::map(1:(total_number - 1), function(i) {
#     cat(i, " ")
#     temp_result =
#       purrr::map((i + 1):total_number,
#                  .f = function(j) {
#                    cat(paste(i, j, sep = "-"), " ")
#                    x = as.numeric(expression_data[i, ])
#                    time1 = sample_info$accurate_time
#                    y = as.numeric(expression_data[j, ])
#                    time2 = sample_info$accurate_time
# 
#                    result = lagged_correlation(
#                      x = x,
#                      y = y,
#                      time1 = time1,
#                      time2 = time2,
#                      time_tol = 300 / 60,
#                      step = 30 / 60,
#                      min_matched_sample = 10,
#                      all_idx = all_idx,
#                      progressbar = FALSE
#                    )
# 
#                    p = result$all_cor_p
#                    p = p.adjust(p,
#                                 method = "BH",
#                                 n = nrow(variable_info) * (nrow(variable_info) - 1) / 2)
#                    if (all(p > 0.05, na.rm = TRUE)) {
#                      result = NULL
#                    }
#                    result
#                  }
#       )
#     names(temp_result) = rownames(expression_data)[(i + 1):total_number]
#     temp_result
#   })
# 
# names(lagged_result) =  rownames(expression_data)[1:(total_number - 1)]
# save(lagged_result, file = "lagged_correlation/lagged_result")

load("lagged_correlation/lagged_result")

# for(i in 1:length(lagged_result)){
#   cat(i, " ")
#   for(j in 1:length(lagged_result[[i]])){
#     lagged_result[[i]][[j]]$shift_time = shift_time
#   }
# }


# #####evaluate peak quality
# 

evaluate_peak_quality(object = lagged_result[[2]][[15]])

###output the matched plot
dir.create("sample_matching_plot")

idx = 
lagged_result[[2]] %>% lapply(class) %>% unlist() %>% 
  `==`("list") %>% 
  which()


load(here::here("data/24_7_study/summary_info/day_night_df"))




# cor_data =
#   purrr::map(1:length(lagged_result), function(i) {
#     cat(i, " ")
#     x = lagged_result[[i]]
#     purrr::map(1:length(x), function(j) {
#       # cat(j, " ")
#       y = x[[j]]
#       if (is.null(y) | length(y) == 1) {
#         return(NULL)
#       }
#       temp_data =
#       data.frame(
#         from = names(lagged_result)[i],
#         to = names(x)[j],
#         shift_time = y$shift_time,
#         cor = y$all_cor,
#         p = y$all_cor_p
#       )
# 
#       result = evaluate_peak_quality(object = y, plot = FALSE)
#       temp_data$score = result$score
#       temp_data
#     }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# save(cor_data, file = "lagged_correlation/cor_data", compress = "xz")
load("lagged_correlation/cor_data")
  
##output the shift time vs correlation plots
temp_data =
  cor_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::group_by(from, to) %>%
  dplyr::filter(abs(cor) == max(abs(cor))) %>%
  dplyr::ungroup()

temp_data$p_adjust = p.adjust(temp_data$p, method = "bonferroni",
                               n = nrow(expression_data) * (nrow(expression_data) - 1)/2)

temp_data =
  temp_data %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::filter(shift_time != "(-15,15]")

dir.create("shift_time_vs_cor")

# for(i in 1:nrow(temp_data)){
#   cat(i, " ")
#   from = temp_data$from[i]
#   to = temp_data$to[i]
#   idx1 = which(names(lagged_result) == from)
#   idx2 = which(names(lagged_result[[idx1]]) == to)
#   result =
#     evaluate_peak_quality(object = lagged_result[[idx1]][[idx2]], plot = TRUE)
#   name =
#     paste(round(result$score, 4), "_",from, "_",to, ".pdf",sep = "") %>%
#     stringr::str_replace_all("/", "_")
#   ggsave(result$plot, filename = file.path("shift_time_vs_cor", name), width = 9, height = 7)
# }

###check the shift time vs correlation, and set the score cutoff as 0.5

lagged_cor =
  cor_data %>%
  dplyr::filter(!is.na(cor))

lagged_cor$cor[lagged_cor$shift_time != "(-15,15]" &
                 lagged_cor$score < 0.5] = 0

lagged_cor =
  lagged_cor %>%
  dplyr::group_by(from, to) %>%
  dplyr::filter(abs(cor) == max(abs(cor))) %>%
  dplyr::ungroup() 

lagged_cor$p_adjust = p.adjust(lagged_cor$p, method = "bonferroni", 
                               n = nrow(expression_data) * (nrow(expression_data) - 1)/2)

lagged_cor = 
  lagged_cor %>% 
  dplyr::filter(p_adjust < 0.05)

global_cor =
  cor_data %>%
  dplyr::filter(shift_time == "(-15,15]")

global_cor$p_adjust = p.adjust(global_cor$p, method = "bonferroni",
                               nrow(expression_data) * (nrow(expression_data) - 1)/2)

global_cor = 
  global_cor %>% 
  dplyr::filter(p_adjust < 0.05)

dim(lagged_cor)
dim(global_cor)

table(lagged_cor$shift_time)

lagged_cor =
  lagged_cor %>%
  dplyr::mutate(name = paste(from, to, sep = "_")) %>%
  dplyr::select(name, dplyr::everything()) %>%
  dplyr::arrange(name)

global_cor = 
  global_cor %>% 
  dplyr::mutate(name = paste(from, to, sep = "_")) %>% 
  dplyr::select(name, dplyr::everything()) %>% 
  dplyr::arrange(name)

global_cor %>% 
  dplyr::arrange(abs(cor)) %>% 
  tail()

plot(as.numeric(expression_data["lipid_465",]),
     as.numeric(expression_data["lipid_468",]))

dim(global_cor)
dim(lagged_cor)

library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "lagged correlation",
#              gridLines = TRUE)
# addWorksheet(wb, sheetName = "global correlation",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = lagged_cor,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 2, x = global_cor,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "lagged_correlation/cor_data.xlsx", overwrite = TRUE)


















