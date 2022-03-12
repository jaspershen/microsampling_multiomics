no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

####load data
{
  ###
  load(here::here("data/7_24_mike/metabolomics/data_preparation/peaks/expression_data"))
  load(here::here("data/7_24_mike/metabolomics/data_preparation/peaks/sample_info"))
  load(here::here("data/7_24_mike/metabolomics/data_preparation/peaks/variable_info"))

  ###remove QC and blanks
  sample_info = 
  sample_info %>% 
    dplyr::filter(!is.na(accurate_time))
  
  expression_data = 
    expression_data[,sample_info$sample_id]
  }

library(MetaCycle)

setwd("data/7_24_mike/circadian_analysis/metabolic_feature")

###reference https://abego.cn/2019/05/31/the-rule-of-gene-expression-in-the-day-and-nigth/
###https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html

sample_info$time

data.frame(
  time = sample_info$time,
  value = as.numeric(expression_data[1, ]),
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

result1$BH.Q

idx = which(result1$p < 0.05)
idx
plot(result1$BH.Q[idx])

load(here::here("data/7_24_mike/summary_info/day_night_df"))

day_night_df = 
data.frame(start = hms::hms(seconds = 0,minutes = 0,hours = 6),
           end = hms::hms(seconds = 0,minutes = 0,hours = 18))

#####output plot
dir.create("plot")
for (i in idx) {
  p = result1$p[i]
  bh = result1$BH.Q[i]
  mol_name = variable_info$mol_name[i]
  value = as.numeric(expression_data[i,])

  temp_data =
  data.frame(
    time = sample_info$time,
    value = value,
    day = as.character(sample_info$day)
  )

  plot =
  ggplot() +
    geom_rect(
      mapping = aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "lightyellow",
      data = day_night_df,
      show.legend = FALSE
    ) +
    geom_point(aes(x = time,
                   y = value,
                   color = day),
               alpha = 1,
               size = 2,
               data = temp_data) +
    geom_smooth(aes(x = time,
                    y = value),
                color = "black",
                size = 2,
                se = FALSE,
                data = temp_data) +
    ggsci::scale_color_lancet() +
    labs(y = "Continuous glucose monitoring", x = "",
         title = paste(mol_name, ",p:",round(p, 5), ",BH", round(bh, 5))) +
    base_theme

  if(bh < 0.2){
  plot =
    plot +
    stat_smooth(geom='line',
                aes(x = time,
                    y = value,
                    color = day),
                alpha=0.3,
                se=FALSE,
                data = temp_data)
  }else{
    plot =
      plot +
      stat_smooth(geom='line',
                  aes(x = time,
                      y = value,
                      color = day),
                  alpha=0.3,
                  se=FALSE,
                  data = temp_data)
  }

  name =
    paste(round(bh,5),variable_info$variable_id[i], "pdf", sep =".")
  ggsave(plot, filename = file.path("plot", name),
         width = 10, height = 7)
}


####heatmap
temp_data = 
expression_data[idx,]

rownames(temp_data) = variable_info$mol_name[idx]

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

p = -log(result1$BH.Q[idx], 10)
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
    x = -log(result1$BH.Q[idx], 10),
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




