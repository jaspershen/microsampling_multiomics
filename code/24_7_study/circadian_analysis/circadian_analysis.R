no_function()

library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")

####load data all the omics data
{
  ###
  load(here::here("data/24_7_study/combine_omics/data_preparation/expression_data"))
  load(here::here("data/24_7_study/combine_omics/data_preparation/sample_info"))
  load(here::here("data/24_7_study/combine_omics/data_preparation/variable_info"))
}

####load consistence score for each molecules
load(here::here("data/24_7_study/circadian_analysis/all_omics/day_consistence/consistence_score"))

consistence_score$consistence_score1 = 
  consistence_score$consistence_score * 0.5 + consistence_score$new_consistence_score * 0.5

variable_info =
  variable_info %>%
  dplyr::left_join(consistence_score, by = "variable_id")

library(MetaCycle)

setwd("data/24_7_study/circadian_analysis/all_omics")

###scale for omics
dim(expression_data)
expression_data =
expression_data %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

###reference https://abego.cn/2019/05/31/the-rule-of-gene-expression-in-the-day-and-nigth/
###https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html

sample_info$time

data.frame(
  time = sample_info$time,
  value = as.numeric(expression_data["lipid_135", ]),
  day = as.character(sample_info$day)
) %>%
  ggplot(aes(time, value)) +
  geom_line(aes(group = day, color = day)) +
  facet_wrap(facets = vars(day), ncol = 1) +
  geom_point(aes(color = day))

##temp_data is the expression data, one row is one variable, and one column is
##one sample (time point)
temp_data =
  expression_data %>%
  tibble::rownames_to_column(var = "variable_id")

time_point = sample_info$time
sample_info$day

# ####remove day 2019-05-03 and 2019-05-07
# sample_info = 
# sample_info %>% 
#   dplyr::filter(!as.character(day) %in% c("2019-05-03", "2019-05-07"))
# 
# expression_data = 
#   expression_data[,sample_info$sample_id]

###as,numeric transfer time to seconds of this day.
time_point = 
  as.numeric(sample_info$time) / 60 / 60 
  # as.numeric((sample_info$day - sample_info$day[1]) * 24)

plot(time_point)

temp_data =
  expression_data %>% 
  tibble::rownames_to_column(var = "variable_id")

colnames(temp_data)[-1] = time_point

# write.csv(temp_data, "temp_data.csv", row.names = FALSE)

# meta2d(
#     infile = "temp_data.csv",
#     filestyle = "csv",
#     outdir = ".",
#     timepoints = time_point,
#     # cycMethod = "LS",
#     outputFile = TRUE,
#     outRawData = FALSE
#   )

####read result
result1 = readr::read_csv("LSresult_temp_data.csv")
result2 = readr::read_csv("meta2d_temp_data.csv")

plot = 
result1 %>%
  dplyr::left_join(variable_info, by = c("CycID" = "variable_id")) %>%
  ggplot(aes(-log(BH.Q, 10), consistence_score1)) +
  geom_point() +
  geom_vline(xintercept = 0.69897, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  base_theme +
  labs(y = "Consistence score") +
  ggrepel::geom_label_repel(aes(label = ifelse(consistence_score1 > 0.6, mol_name, NA)))

plot

ggsave(plot, filename = "circadian_p_vs_consistence_score.pdf", width = 9, height = 7)

result1 %>% 
  dplyr::filter(BH.Q < 0.001)

result1 %>% 
  dplyr::left_join(variable_info, by = c("CycID" = "variable_id")) %>% 
  dplyr::filter(consistence_score1 < 0 & BH.Q < 0.2)

result1$BH.Q

####so the molecules with consistence score < 0 will not be considered for circiadian
####analysis.
new_result =
  result1 %>% 
  dplyr::left_join(variable_info[,c("variable_id", "consistence_score1")],
                   by = c("CycID" = "variable_id")) %>% 
  dplyr::filter(consistence_score1 > 0)

new_result$BH.Q = p.adjust(new_result$p, method = "BH")

idx = which(new_result$BH.Q < 0.2)
idx

load(here::here("data/24_7_study/summary_info/day_night_df"))

day_night_df = 
data.frame(start = hms::hms(seconds = 0,minutes = 0,hours = 6),
           end = hms::hms(seconds = 0,minutes = 0,hours = 18))

# #####output plot
# dir.create("plot")
# for (i in idx) {
#   p = new_result$p[i]
#   bh = new_result$BH.Q[i]
#   mol_id = new_result$CycID[i]
#   mol_name = variable_info$mol_name[match(mol_id, variable_info$variable_id)]
#   
#   value = as.numeric(expression_data[mol_id,])
# 
#   temp_data =
#   data.frame(
#     time = sample_info$time,
#     value = value,
#     day = as.character(sample_info$day)
#   )
# 
#   plot =
#   ggplot() +
#     geom_rect(
#       mapping = aes(
#         xmin = start,
#         xmax = end,
#         ymin = -Inf,
#         ymax = Inf
#       ),
#       fill = "lightyellow",
#       data = day_night_df,
#       show.legend = FALSE
#     ) +
#     geom_point(aes(x = time,
#                    y = value,
#                    color = day),
#                alpha = 1,
#                size = 2,
#                data = temp_data) +
#     geom_smooth(aes(x = time,
#                     y = value),
#                 color = "black",
#                 size = 2,
#                 se = FALSE,
#                 data = temp_data) +
#     ggsci::scale_color_lancet() +
#     labs(y = "Continuous glucose monitoring", x = "",
#          title = paste(mol_name, ",p:",round(p, 5), ",BH", round(bh, 5))) +
#     base_theme
# 
#   if(bh < 0.2){
#   plot =
#     plot +
#     stat_smooth(geom='line',
#                 aes(x = time,
#                     y = value,
#                     color = day),
#                 alpha=0.3,
#                 se=FALSE,
#                 data = temp_data)
#   }else{
#     plot =
#       plot +
#       stat_smooth(geom='line',
#                   aes(x = time,
#                       y = value,
#                       color = day),
#                   alpha=0.3,
#                   se=FALSE,
#                   data = temp_data)
#   }
# 
#   name =
#     paste(round(bh,5),mol_id, "pdf", sep =".")
#   ggsave(plot, filename = file.path("plot", name),
#          width = 10, height = 7)
# }


####heatmap
temp_data = 
expression_data[new_result$CycID[idx],]

rownames(temp_data) = variable_info$mol_name[match(rownames(temp_data),variable_info$variable_id)]

sample_info = 
  sample_info %>% 
  dplyr::arrange(time)

temp_data =
  temp_data[,sample_info$sample_id]

temp_data =
  sort(unique(sample_info$hour)) %>%
  purrr::map(function(x) {
    temp_idx = which(sample_info$hour == x)
    temp_data[, temp_idx, drop = FALSE] %>%
      apply(1, median)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()
  
colnames(temp_data) = sort(unique(sample_info$hour))

temp_data =
  temp_data %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>% 
  t()

library(ComplexHeatmap)
library(circlize)

range(temp_data)

col_fun = colorRamp2(
  breaks = c(-3, 0, 3),
  colors =
    c("#366A9FFF", "white", "red"),
  transparency = 0
)

p = -log(new_result$BH.Q[idx], 10)
bar_col = ifelse(p > 1.30103, "red", "grey")

text_cor = 
  class_color[variable_info$data_type[idx]]

sample_number = as.numeric(table(sample_info$hour))

plot = 
Heatmap(
  as.matrix(temp_data),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  name = "z-score",
  border = TRUE,
  col = col_fun,
  column_names_rot = 45,
  clustering_method_rows = "ward.D",
  clustering_distance_rows = "euclidean",
  rect_gp = gpar(col = "white"),
  row_names_gp = gpar(cex = 0.8, col = text_cor),
  column_names_gp = gpar(cex = 0.8),
  row_km = 4,
  show_row_names = TRUE,
  row_names_side = "left", 
  top_annotation =  
    columnAnnotation("Sample number" = anno_barplot(
    x = sample_number,
    gp = gpar(col = "black", 
              fill = ggsci::pal_lancet()(n=10)[4])
  ))
) +
  rowAnnotation("-log(10, BH)" = anno_barplot(
    x = -log(new_result$BH.Q[idx], 10),
    gp = gpar(col = "black", 
              fill = bar_col)
  ))
plot
library(ggplotify)

pdf(file = "heatmap.pdf", width = 7, height = 10)
plot
dev.off()

module = 
ComplexHeatmap::row_order(plot) %>% 
  purrr::map(function(x){
    rownames(temp_data)[x]
  }) 

module =
  purrr::map2(module, names(module), function(x, y) {
    data.frame(mol_name = x, module = y)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

plot =
  temp_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "mol_name") %>%
  tidyr::pivot_longer(cols = -mol_name,
                      names_to = "time",
                      values_to = "value") %>%
  dplyr::mutate(time = as.numeric(time)) %>%
  dplyr::left_join(module, by = "mol_name") %>%
  ggplot(aes(time, value, group = mol_name)) +
  geom_hline(yintercept = 0) +
  # geom_line(aes(group = mol_name),
  #           show.legend = FALSE) +
  geom_smooth(aes(group = mol_name,
                  col = module),
              se = FALSE,
              show.legend = FALSE) +
  labs(x = "Hour of one day", y = "z-score") +
  ggsci::scale_colour_npg() +
  base_theme +
  facet_grid(rows = vars(module))
plot
ggsave(plot, filename = "module.pdf", width = 7, height = 14)









