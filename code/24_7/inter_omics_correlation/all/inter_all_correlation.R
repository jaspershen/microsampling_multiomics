no_function()
library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

###data loading
{
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
    which(
      apply(metabolic_panel_expression_data, 1, function(x){
        sum(x == 0)/ncol(metabolic_panel_expression_data)
      }) < 0.5
    )
  
  metabolic_panel_variable_info = 
    metabolic_panel_variable_info[remain_idx,]
  
  metabolic_panel_expression_data = 
    metabolic_panel_expression_data[metabolic_panel_variable_info$variable_id,]
  
  ###remove CHEX
  metabolic_panel_variable_info = 
  metabolic_panel_variable_info %>% 
    dplyr::filter(!stringr::str_detect(mol_name, "CHEX"))
  
  metabolic_panel_expression_data = 
    metabolic_panel_expression_data[metabolic_panel_variable_info$variable_id,]
  
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
    dplyr::mutate(week = paste(
      week,
      lubridate::month(day),
      lubridate::day(day),
      sep = "-"
    )) %>% 
    dplyr::mutate(week = factor(week, unique(week)))
  
}


####cytokine data
{
  ###cytokine
  load("data/7_24_mike/cytokine/data_preparation/sample_info")
  load("data/7_24_mike/cytokine/data_preparation/variable_info")
  load("data/7_24_mike/cytokine/data_preparation/expression_data")
  cytokine_sample_info = sample_info
  cytokine_variable_info = variable_info
  cytokine_expression_data = expression_data

  remain_idx = 
    which(
      apply(cytokine_expression_data, 1, function(x){
        sum(x == 0)/ncol(cytokine_expression_data)
      }) < 0.5
    )
  
  cytokine_variable_info = 
    cytokine_variable_info[remain_idx,]
  
  cytokine_expression_data = 
    cytokine_expression_data[cytokine_variable_info$variable_id,]
  
  
  ###remove the controls
  remove_idx = grep("CHEX", cytokine_variable_info$mol_name)
  if(length(remove_idx) > 0){
    cytokine_variable_info = 
      cytokine_variable_info[-remove_idx,]
    cytokine_expression_data = 
      cytokine_expression_data[-remove_idx,]
  }
}

##lipids
{
  ###lipidomics
  load("data/7_24_mike/lipidomics/data_preparation/sample_info")
  load("data/7_24_mike/lipidomics/data_preparation/variable_info")
  load("data/7_24_mike/lipidomics/data_preparation/expression_data")
  lipidomics_sample_info = sample_info
  lipidomics_variable_info = variable_info
  lipidomics_expression_data = expression_data
  
  
  remain_idx = 
    which(
      apply(lipidomics_expression_data, 1, function(x){
        sum(x == 0)/ncol(lipidomics_expression_data)
      }) < 0.5
    )
  
  lipidomics_variable_info = 
    lipidomics_variable_info[remain_idx,]
  
  lipidomics_expression_data = 
    lipidomics_expression_data[lipidomics_variable_info$variable_id,]
  
}


###metabolomics
{
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  metabolomics_expression_data = expression_data

  ###remove wrong metabolites
  metabolomics_variable_info = 
    metabolomics_variable_info %>% 
    dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                                  "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                                  "nistDatabase0.0.2"))
  
  metabolomics_expression_data = 
    metabolomics_expression_data[metabolomics_variable_info$variable_id,]
  
  remain_idx = 
    which(
      apply(metabolomics_expression_data, 1, function(x){
        sum(x == 0)/ncol(metabolomics_expression_data)
      }) < 0.5
    )
  
  metabolomics_variable_info = 
    metabolomics_variable_info[remain_idx,]
  
  metabolomics_expression_data = 
    metabolomics_expression_data[metabolomics_variable_info$variable_id,]
  
  ###remove QC samples
  metabolomics_sample_info = 
    metabolomics_sample_info %>% 
    dplyr::filter(!is.na(accurate_time))
  
  metabolomics_expression_data = 
    metabolomics_expression_data[,metabolomics_sample_info$sample_id]  
  
}


#####proteomics
{
  ###proteomics
  load("data/7_24_mike/proteomics/data_preparation/sample_info")
  load("data/7_24_mike/proteomics/data_preparation/variable_info")
  load("data/7_24_mike/proteomics/data_preparation/expression_data")
  proteomics_sample_info = sample_info
  proteomics_variable_info = variable_info
  proteomics_expression_data = expression_data
  
  load("data/7_24_mike/summary_info/day_night_df")
  
  remain_idx = 
    which(
      apply(proteomics_expression_data, 1, function(x){
        sum(x == 0)/ncol(proteomics_expression_data)
      }) < 0.5
    )
  
  proteomics_variable_info = 
    proteomics_variable_info[remain_idx,]
  
  proteomics_expression_data = 
    proteomics_expression_data[proteomics_variable_info$variable_id,]
}


#####cortisol
{
  ###cortisol
  load("data/7_24_mike/cortisol/data_preparation/sample_info")
  load("data/7_24_mike/cortisol/data_preparation/variable_info")
  load("data/7_24_mike/cortisol/data_preparation/expression_data")
  cortisol_sample_info = sample_info
  cortisol_variable_info = variable_info
  cortisol_expression_data = expression_data
  
  load("data/7_24_mike/summary_info/day_night_df")
  
  remain_idx = 
    which(
      apply(cortisol_expression_data, 1, function(x){
        sum(x == 0)/ncol(cortisol_expression_data)
      }) < 0.5
    )
  
  cortisol_variable_info = 
    cortisol_variable_info[remain_idx,]
  
  cortisol_expression_data = 
    cortisol_expression_data[cortisol_variable_info$variable_id,]
}



#####total_protein
{
  ###total_protein
  load("data/7_24_mike/total_protein/data_preparation/sample_info")
  load("data/7_24_mike/total_protein/data_preparation/variable_info")
  load("data/7_24_mike/total_protein/data_preparation/expression_data")
  total_protein_sample_info = sample_info
  total_protein_variable_info = variable_info
  total_protein_expression_data = expression_data
  
  load("data/7_24_mike/summary_info/day_night_df")
  
  remain_idx = 
    which(
      apply(total_protein_expression_data, 1, function(x){
        sum(x == 0)/ncol(total_protein_expression_data)
      }) < 0.5
    )
  
  total_protein_variable_info = 
    total_protein_variable_info[remain_idx,]
  
  total_protein_expression_data = 
    total_protein_expression_data[total_protein_variable_info$variable_id,]
}


####combine all omics together
intersect_time = 
Reduce(f = intersect, x = list(
  as.character(lipidomics_sample_info$accurate_time),
  as.character(metabolomics_sample_info$accurate_time),
  as.character(cytokine_sample_info$accurate_time),
  as.character(total_protein_sample_info$accurate_time),
  as.character(cortisol_sample_info$accurate_time),
  as.character(metabolic_panel_sample_info$accurate_time),
  as.character(proteomics_sample_info$accurate_time)
))

lipidomics_sample_info =
  lipidomics_sample_info %>%
  dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
  dplyr::arrange(accurate_time)

lipidomics_expression_data = lipidomics_expression_data[, lipidomics_sample_info$sample_id]

metabolomics_sample_info =
  metabolomics_sample_info %>%
  dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
  dplyr::arrange(accurate_time)

metabolomics_expression_data = metabolomics_expression_data[, metabolomics_sample_info$sample_id]
  
cytokine_sample_info =
  cytokine_sample_info %>%
  dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
  dplyr::arrange(accurate_time)

cytokine_expression_data = cytokine_expression_data[, cytokine_sample_info$sample_id]

total_protein_sample_info =
  total_protein_sample_info %>%
  dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
  dplyr::arrange(accurate_time)

total_protein_expression_data = total_protein_expression_data[, total_protein_sample_info$sample_id]

cortisol_sample_info =
  cortisol_sample_info %>%
  dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
  dplyr::arrange(accurate_time)

cortisol_expression_data = cortisol_expression_data[, cortisol_sample_info$sample_id]

metabolic_panel_sample_info =
  metabolic_panel_sample_info %>%
  dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
  dplyr::arrange(accurate_time)

metabolic_panel_expression_data = metabolic_panel_expression_data[, metabolic_panel_sample_info$sample_id]

proteomics_sample_info =
  proteomics_sample_info %>%
  dplyr::filter(as.character(accurate_time) %in% intersect_time) %>% 
  dplyr::arrange(accurate_time)

proteomics_expression_data = proteomics_expression_data[, proteomics_sample_info$sample_id]

lipidomics_sample_info$accurate_time == sample_info$accurate_time
lipidomics_sample_info$accurate_time == metabolomics_sample_info$accurate_time
lipidomics_sample_info$accurate_time == proteomics_sample_info$accurate_time
lipidomics_sample_info$accurate_time == cytokine_sample_info$accurate_time
lipidomics_sample_info$accurate_time == cortisol_sample_info$accurate_time
lipidomics_sample_info$accurate_time == total_protein_sample_info$accurate_time

colnames(lipidomics_expression_data) =
  colnames(metabolomics_expression_data) =
  colnames(cytokine_expression_data) =
  colnames(total_protein_expression_data) =
  colnames(cortisol_expression_data) =
  colnames(metabolic_panel_expression_data) =
  colnames(proteomics_expression_data) =
  as.character(lipidomics_sample_info$accurate_time)

expression_data = 
  rbind(lipidomics_expression_data,
        metabolomics_expression_data,
        cytokine_expression_data,
        total_protein_expression_data,
        cortisol_expression_data,
        metabolic_panel_expression_data,
        proteomics_expression_data)

lipidomics_sample_info$sample_id = 
  as.character(lipidomics_sample_info$accurate_time)

colnames(lipidomics_sample_info)

metabolic_panel_sample_info$sample_id = 
  as.character(metabolic_panel_sample_info$accurate_time)

colnames(metabolic_panel_sample_info)

metabolomics_sample_info$sample_id = 
  as.character(metabolomics_sample_info$accurate_time)

colnames(metabolomics_sample_info)

proteomics_sample_info$sample_id = 
  as.character(proteomics_sample_info$accurate_time)

colnames(proteomics_sample_info)

cytokine_sample_info$sample_id = 
  as.character(cytokine_sample_info$accurate_time)

colnames(cytokine_sample_info)

cortisol_sample_info$sample_id = 
  as.character(cortisol_sample_info$accurate_time)

colnames(cortisol_sample_info)

total_protein_sample_info$sample_id = 
  as.character(total_protein_sample_info$accurate_time)

colnames(total_protein_sample_info)

sample_info = 
  lipidomics_sample_info %>% 
  dplyr::full_join(metabolic_panel_sample_info,
                   by = c("subject_id", "sample_id", "accurate_time")) %>% 
  dplyr::full_join(metabolomics_sample_info,
                   by = c("subject_id", "sample_id", "accurate_time", "day", "time", "hour")) %>% 
  dplyr::full_join(proteomics_sample_info,
                   by = c("subject_id", "sample_id", "accurate_time"))  %>% 
  dplyr::full_join(cytokine_sample_info,
                   by = c("subject_id", "sample_id", "accurate_time", "day", "time", "hour")) %>% 
  dplyr::full_join(cortisol_sample_info,
                   by = c("subject_id", "sample_id", "accurate_time", "day", "time", "hour")) %>% 
  dplyr::full_join(total_protein_sample_info,
                   by = c("subject_id", "sample_id", "accurate_time")) 

###variable_info
metabolic_panel_variable_info =
metabolic_panel_variable_info %>% 
  dplyr::mutate(data_type = "metabolic_panel")

metabolomics_variable_info= 
metabolomics_variable_info %>% 
  dplyr::mutate(data_type = "metabolomics") %>% 
  dplyr::mutate(mol_name = Compound.name)

lipidomics_variable_info = 
lipidomics_variable_info %>% 
  dplyr::mutate(data_type = "lipidomics") %>% 
  dplyr::mutate(mol_name = case_when(
    is.na(Lipid_Name) ~ mol_name,
    !is.na(Lipid_Name) ~ Lipid_Name
  ))

proteomics_variable_info = 
proteomics_variable_info %>% 
  dplyr::mutate(data_type = "proteomics")

cytokine_variable_info = 
cytokine_variable_info %>% 
  dplyr::mutate(data_type = "cytokine") %>% 
  dplyr::rename(cytokine_classification = classification)

cortisol_variable_info = 
cortisol_variable_info %>% 
  dplyr::mutate(data_type = "cytokine") %>% 
  dplyr::rename(cortisol_class = class)

total_protein_variable_info = 
total_protein_variable_info %>% 
  dplyr::mutate(data_type = "total_protein") %>% 
  dplyr::rename(total_protein_class = class)

variable_info =
  rbind(lipidomics_variable_info[,c("variable_id", "mol_name", "data_type")],
        metabolomics_variable_info[,c("variable_id", "mol_name", "data_type")],
        cytokine_variable_info[,c("variable_id", "mol_name", "data_type")],
        total_protein_variable_info[,c("variable_id", "mol_name", "data_type")],
        cortisol_variable_info[,c("variable_id", "mol_name", "data_type")],
        metabolic_panel_variable_info[,c("variable_id", "mol_name", "data_type")],
        proteomics_variable_info[,c("variable_id", "mol_name", "data_type")])

dim(variable_info)

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

}

######
dir.create("data/7_24_mike/inter_omics_correlation")
dir.create("data/7_24_mike/inter_omics_correlation/inter_all_omics")
setwd("data/7_24_mike/inter_omics_correlation/inter_all_omics")

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
step = 60 / 60

time_window1 = 
  seq(from = step/2, to = time_tol, by = step)

time_window2 = 
  -rev(seq(from = step/2, to = time_tol, by = step))

time_window = sort(c(time_window2, time_window1))

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

dir.create("lagged_correlation")
total_number = nrow(expression_data)

# lagged_result =
#   purrr::map(1:(total_number - 1), function(i) {
#     cat(i, " ")
#     temp_result =
#       purrr::map((i + 1):total_number,
#                  .f = function(j) {
#                    cat(paste(i, j, sep = "-"), " ")
#                    x = as.numeric(expression_data[i,])
#                    time1 = sample_info$accurate_time
#                    y = as.numeric(expression_data[j,])
#                    time2 = sample_info$accurate_time
# 
#                    result = lagged_correlation(
#                      x = x,
#                      y = y,
#                      time1 = time1,
#                      time2 = time2,
#                      time_tol = 300 / 60,
#                      step = 60 / 60,
#                      min_matched_sample = 10,
#                      all_idx = all_idx,
#                      progressbar = FALSE
#                    )
# 
#                    p = result$all_cor_p
#                    p = p.adjust(p, method = "BH", n = nrow(variable_info) * (nrow(variable_info)-1)/2)
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

# #####evaluate peak quality
# 

evaluate_peak_quality(object = lagged_result[[2]][[15]])

###output the matched plot
dir.create("sample_matching_plot")

idx = 
lagged_result[[2]] %>% lapply(class) %>% unlist() %>% 
  `==`("list") %>% 
  which()

# for(i in 1:length(lagged_result[[2]][[idx[1]]])){
#   cat(i, " ")
#   plot =
#     show_sample_matching(object = lagged_result[[2]][[idx[1]]],
#                          index = i,
#                          only_remain_matched = TRUE,
#                          day_night_df = day_night_df,
#                          add_text = TRUE)
#   ggsave(plot,
#          filename = file.path("sample_matching_plot", paste("matching_plot", i, ".pdf", sep = "")),
#          width = 20, height = 7)
# 
# }

# cor_data =
#   purrr::map(1:length(lagged_result), function(i) {
#     cat(i, " ")
#     x = lagged_result[[i]]
#     purrr::map(1:length(x), function(j) {
#       y = x[[j]]
#       if (is.null(y)) {
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
  
###output the shift time vs correlation plots
# temp_data =
#   cor_data %>%
#   dplyr::filter(!is.na(cor)) %>%
#   dplyr::group_by(from, to) %>%
#   dplyr::filter(abs(cor) == max(abs(cor))) %>%
#   dplyr::ungroup()
#   
# temp_data$p_adjust = p.adjust(temp_data$p, method = "bonferroni", 
#                                n = nrow(expression_data) * (nrow(expression_data) - 1)/2)
# 
# temp_data = 
#   temp_data %>% 
#   dplyr::filter(p_adjust < 0.05) %>% 
#   dplyr::filter(shift_time != "(-30,30]")
# 
# dir.create("shift_time_vs_cor")

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

# lagged_cor$cor[lagged_cor$shift_time != "(-30,30]" &
#                  lagged_cor$score < 0.5] = 0

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
  dplyr::filter(shift_time == "(-30,30]")

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

# library(openxlsx)
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



















