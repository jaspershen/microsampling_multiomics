no_function()
library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")
# source("code/modified_dtw.R")
source("code/lagged_correlation.R")

####we need the lipidomics module information
lipidomics_variable_info = 
  readxl::read_xlsx("data/24_7_study/lipidomics/k_means_clustering/final_info_manual.xlsx")

load("data/24_7_study/lipidomics/k_means_clustering/new_variable_info")
load("data/24_7_study/lipidomics/k_means_clustering/new_expression_data")

lipidomics_variable_info = 
  new_variable_info

####load all omics data loess data
load(here::here("data/24_7_study/combine_omics/data_preparation/new_expression_data"))
load(here::here("data/24_7_study/combine_omics/data_preparation/new_sample_info"))
load(here::here("data/24_7_study/combine_omics/data_preparation/new_variable_info"))

variable_info1 =
  new_variable_info %>% 
  dplyr::filter(data_type == "lipidomics")

expression_data1 =
  new_expression_data[variable_info1$variable_id,]

variable_info2 =
  new_variable_info %>% 
  dplyr::filter(data_type != "lipidomics")

expression_data2 =
  new_expression_data[variable_info2$variable_id,]

expression_data1 = 
lipidomics_variable_info$old_variable_id %>% 
  purrr::map(function(x){
    x = stringr::str_split(x, ";")[[1]]
    expression_data1[x,,drop = FALSE] %>% 
      colMeans()
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(expression_data1) = lipidomics_variable_info$variable_id

variable_info1 = 
  data.frame(lipidomics_variable_info[,c("variable_id", "mol_name")],
             data_type = "lipidomics")

variable_info = 
  rbind(variable_info1,
        variable_info2)

expression_data =
  rbind(expression_data1,
        expression_data2)

sample_info =
  new_sample_info

rownames(expression_data) == variable_info$variable_id
colnames(expression_data) == sample_info$sample_id

####this is for the day night time
load("data/24_7_study/summary_info/day_night_df")

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

######
setwd("data/24_7_study/inter_omics_correlation/inter_all_omics_loess_data")
load("lagged_correlation/lagged_result")

# lagged_cor =
#   readxl::read_xlsx("lagged_correlation/cor_data.xlsx",1)
# 
# table(lagged_cor$shift_time)
# 
# lagged_cor %>%
#   dplyr::filter(shift_time != "(-15,15]") %>%
#   pull(cor) %>%
#   plot()
# 
# unique(c(lagged_cor$from,
#          lagged_cor$to))
# 
# dim(lagged_cor)
# 
# lagged_cor =
#   lagged_cor %>%
#   dplyr::left_join(variable_info, by = c("from" = "variable_id")) %>%
#   dplyr::rename(from_mol_name = mol_name,
#                 from_data_type = data_type) %>%
#   dplyr::left_join(variable_info, by = c("to" = "variable_id")) %>%
#   dplyr::rename(to_mol_name = mol_name,
#                 to_data_type = data_type)
# 
# 
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "lagged correlation",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = lagged_cor,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "lagged_correlation/lagged_cor.xlsx", overwrite = TRUE)

lagged_cor =
  readxl::read_xlsx("lagged_correlation/lagged_cor.xlsx")

table(paste(lagged_cor$from_data_type, lagged_cor$to_data_type, sep = "_"))

temp = lagged_result[[2]][[2]]

temp$max_idx
temp$max_cor
temp$global_cor

x = temp$x
y = purrr::map(.x = temp$max_idx, function(idx){
  mean(temp$y[idx])
}) %>% 
  unlist()

cor(x,y, method = "spearman", use = "complete.obs")
cor(x,y, method = "pearson", use = "complete.obs")
temp$max_cor

####
lagged_cor1 = 
  lagged_cor %>% 
  dplyr::filter(shift_time != "(-15,15]")

lagged_cor2 = 
  lagged_cor %>% 
  dplyr::filter(shift_time == "(-15,15]")

dim(lagged_cor1)
dim(lagged_cor2)

table(lagged_cor1$shift_time)

paste(lagged_cor1$from_data_type, 
      lagged_cor1$to_data_type, sep = "_") %>% 
  table()

table((c(lagged_cor1$from_data_type,
  lagged_cor1$to_data_type)))

#####lagged correlation network
####output the shift_time_vs_cor
dim(lagged_cor1)
length(unique(c(lagged_cor1$from,lagged_cor1$to)))

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
addWorksheet(wb, sheetName = "lagged correlation",
             gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
writeDataTable(wb, sheet = 1, x = lagged_cor1,
               colNames = TRUE, rowNames = FALSE)
saveWorkbook(wb, "lagged_cor_network/lagged_cor1.xlsx", overwrite = TRUE)


####shift time vs correlation
##only output 100 plots
idx = sample(1:nrow(lagged_cor1), 100)
# for(i in idx){
#   cat(i, " ")
# 
#   from = lagged_cor1$from[i]
#   to = lagged_cor1$to[i]
#   idx1 = which(names(lagged_result) == from)
#   idx2 = which(names(lagged_result[[idx1]]) == to)
#   result =
#     evaluate_peak_quality(object = lagged_result[[idx1]][[idx2]], plot = TRUE)
#   name = paste(round(result$score, 4), "_",from, "_",to, ".pdf",sep = "") %>%
#     stringr::str_replace_all("/", '_') %>% 
#     stringr::str_replace_all("\\(", '_') %>% 
#     stringr::str_replace_all("\\)", '_') %>% 
#     stringr::str_replace_all("\\:", '_')
#   ggsave(result$plot, filename = file.path("shift_time_vs_cor", name), width = 9, height = 7)
# }

##output cor plot

dir.create("cor_plot")

name = 
lagged_cor1 %>% 
  dplyr::filter(from_data_type != "lipidomics" | to_data_type != "lipidomics") %>% 
  dplyr::arrange(desc(abs(cor))) %>% 
  head(100) %>% 
  dplyr::pull(name)

idx = match(name, lagged_cor1$name)

# for (i in idx) {
#   cat(i, " ")
# 
#   plot1 =
#     lagged_alignment_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       day_night_df = day_night_df,
#       internal_omics_color = class_color[lagged_cor1$to_data_type[i]],
#       wearable_color = class_color[lagged_cor1$from_data_type[i]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "max",
#       x_limit = c(1, 1000),
#       non_matched_point_size = 2,
#       wearable_point_size = 2,
#       internal_omics_point_size = 2,
#       integrated = FALSE
#     )
# 
#   name = lagged_cor1$name[i] %>%
#     stringr::str_replace_all("/", "_") %>% 
#     stringr::str_replace_all("\\(", '_') %>% 
#     stringr::str_replace_all("\\)", '_') %>% 
#     stringr::str_replace_all("\\:", '_')
#   ggsave(plot1,
#          filename = file.path("cor_plot", paste(name, "plot1.pdf", sep = "")),
#          width = 20, height = 7)
# 
# 
#   plot2 =
#     lagged_alignment_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       day_night_df = day_night_df,
#       internal_omics_color = class_color[lagged_cor1$to_data_type[i]],
#       wearable_color = class_color[lagged_cor1$from_data_type[i]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "max",
#       x_limit = c(1, 1000),
#       non_matched_point_size = 2,
#       wearable_point_size = 2,
#       internal_omics_point_size = 2,
#       add_connect_line = FALSE,
#       integrated = FALSE
#     )
# 
#   ggsave(plot2,
#          filename = file.path("cor_plot", paste(name, "plot2.pdf", sep = "")),
#          width = 20, height = 7)
# }


##output scatter plot
dir.create("scatter_plot")
# for (i in idx) {
#   cat(i, " ")
# 
#   plot1 =
#     lagged_sactter_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "global"
#     )
# 
#   plot2 =
#     lagged_sactter_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "max"
#     )
# 
#   library(patchwork)
# 
#   plot =
#   plot1 + plot2 + patchwork::plot_layout(ncol = 2)
# 
#   name = lagged_cor1$name[i] %>%
#     stringr::str_replace_all("/", "_")
# 
#   ggsave(plot,
#          filename = file.path("scatter_plot", paste(name, "plot.pdf", sep = "")),
#          width = 14, height = 7)
# 
# }

