no_function()
library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

#####load food log data
{
  load(here::here("data/7_24_mike/food_log/data_preparation/expression_data"))
  load(here::here("data/7_24_mike/food_log/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/food_log/data_preparation/variable_info"))
  food_expression_data = expression_data
  food_sample_info = sample_info
  food_variable_info = variable_info
}

food_expression_data[is.na(food_expression_data)] = 0

food_expression_data = 
  food_expression_data %>% 
  apply(1, function(x){
    x/max(x)
  }) %>% 
  t() %>% 
  as.data.frame()

library(ComplexHeatmap)
library(circlize)

col_fun = 
  circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

Heatmap(food_expression_data, col = col_fun, cluster_columns = FALSE)

{
  load(here::here("data/7_24_mike/summary_info/day_night_df"))
  
  ####load data (combined omics data)
  load("data/7_24_mike/cgm/data_preparation/expression_data")
  load("data/7_24_mike/cgm/data_preparation/sample_info")
  load("data/7_24_mike/cgm/data_preparation/variable_info")
}

setwd("data/7_24_mike/cgm/cgm_vs_food")

table(variable_info$data_type)

dim(expression_data)

sample_info %>% 
  ggplot(aes(accurate_time, time)) +
  geom_point() +
  base_theme

  temp =
    expression_data %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id),
      names_to = "sample_id",
      values_to = "value"
    ) %>%
    dplyr::left_join(sample_info[, c("sample_id", "accurate_time")], by = "sample_id") 

  
  plot <-
    ggplot() +
    geom_rect(
      mapping = aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "lightyellow",
      data = day_night_df %>% dplyr::filter(as.character(day)
                                            %in% as.character(sample_info$day)),
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0) +
    geom_line(aes(accurate_time,
                  value),
              color = wearable_color["cgm"],
              data = temp,
              show.legend = FALSE) +
    scale_x_datetime(
      breaks = scales::date_breaks("12 hour"),
      date_labels = "%a %H:%M",
      timezone = "America/Los_Angeles"
    ) +
    theme_bw() +
    theme(
      panel.spacing = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 12
      ),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(
      x = "",
      y = "Z-score"
    ) +
    theme(panel.background = element_rect(fill = alpha("grey", 0.2)))

  plot

  ###food plot
  temp_data =
    food_expression_data %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(cols = -variable_id,
                        names_to = "sample_id",
                        values_to = "value") %>% 
    dplyr::left_join(food_sample_info, by = "sample_id") %>% 
    dplyr::left_join(food_variable_info, by = "variable_id")  %>% 
    dplyr::filter(value != 0)
  
  plot2 = 
  ggplot() +
    geom_rect(
      mapping = aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "lightyellow",
      data = day_night_df %>% dplyr::filter(as.character(day)
                                            %in% as.character(sample_info$day)),
      show.legend = FALSE
    ) +
    geom_point(aes(accurate_time, mol_name,
                   size = value),
               color = "black",
               alpha = 0.8,
               data = temp_data) +
    scale_size_continuous(range = c(0.01, 3)) +
    # scale_color_gradient(low = "white", high = "red") +
    scale_x_datetime(limits = c(min(temp$accurate_time), max(temp$accurate_time)),
      breaks = scales::date_breaks("12 hour"),
      date_labels = "%a %H:%M",
      timezone = "America/Los_Angeles"
    ) +
    base_theme +
    theme(
      panel.spacing = margin(0, 0, 0, 0),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.justification = c(0, 1),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = alpha("grey", 0.2)),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(x = "", y = "")
  
  library(patchwork)

  final_plot =   
  plot2 + plot + patchwork::plot_layout(ncol = 1)
  
  ggsave(
    final_plot,
    filename = "cgm_vs_food.pdf",
    width = 20,
    height = 7
  )

####Heatmap to show the clusters
library(ComplexHeatmap)
dim(expression_data)
rownames(expression_data) == cluster_info$variable_id

temp_data = expression_data

temp_cluster_info =
  cluster_info[match(rownames(expression_data), cluster_info$variable_id),]

rownames(expression_data) == temp_cluster_info$variable_id

cluster = temp_cluster_info$cluster

library(lubridate)

temp_data = 
  temp_data %>% 
  apply(1, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>% 
  t() %>% 
  as.data.frame()

range(temp_data, na.rm = TRUE)
temp_data[temp_data > 3] = 3
temp_data[temp_data < -3] = -3

dim(temp_data)

dim(sample_info)

x_labels = 
  colnames(temp_data) %>% 
  stringr::str_replace("\\:00", "") %>% 
  stringr::str_replace("2019-04-29", "Mon") %>% 
  stringr::str_replace("2019-04-30", "Tue") %>% 
  stringr::str_replace("2019-05-01", "Wed") %>% 
  stringr::str_replace("2019-05-02", "Thu") %>% 
  stringr::str_replace("2019-05-03", "Fri") %>% 
  stringr::str_replace("2019-05-04", "Sat") %>% 
  stringr::str_replace("2019-05-05", "Sun") %>% 
  stringr::str_replace("2019-05-06", "Mon") %>% 
  stringr::str_replace("2019-05-07", "Tue") 

x_labels[-seq(1, length(x_labels), by = 5)] = ""

table(cluster)

library(circlize)

col_fun = colorRamp2(
  breaks = c(-3, 0, 3),
  colors =
    c("#366A9FFF", "white", "red"),
  transparency = 0
)

table(cluster)

plot = 
  Heatmap(
    temp_data,
    col = col_fun,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    name = "Z-score",
    border = TRUE,
    column_names_gp = gpar(fontsize = 10, angle = 45),
    column_names_rot = 45,
    column_labels = x_labels,
    row_split = cluster,
    row_title = rep("", length(unique(cluster))),
    row_title_rot = 0,
    na_col = "grey",
    right_annotation = rowAnnotation(foo = anno_block(
      # gp = gpar(fill = omics_color),
      # labels = c("Lipidomics", "Targeted assay", "Proteomics", "Metabolomics"),
      labels_gp = gpar(col = "white", fontsize = 10)
    ))
  )
rownames(temp_data)[unlist(row_order(plot))]
rownames(temp_data)[unlist(row_order(plot))] == temp_cluster_info$variable_id

plot

temp_idx = 
row_order(plot) %>% 
  purrr::map(function(x){
    rownames(temp_data)[x]
  })

temp_data =
  temp_data[unlist(temp_idx),]

#####use ggplot2 to get the heatmap
colnames(temp_data)

# plot1 = ggplotify::as.ggplot(plot)
# 
# plot1
# 
# ggsave(plot1, filename = "heatmap.dent.pdf", width = 7, height = 7)

# time = lubridate::as_datetime(colnames(temp_data), tz = "America/Los_Angeles")
time = sample_info$accurate_time

min_time = lubridate::as_datetime("2019-04-29 00:00:01", tz = "America/Los_Angeles")
max_time = lubridate::as_datetime("2019-05-07 23:59:00", tz = "America/Los_Angeles")

# time[length(time)] = time[length(time)] + 30*60

colnames(temp_data)

time1 = c(time[1] - 30*60, time[-length(time)])
time2 = c(time)

names(time1) = 
  names(time2) = 
  colnames(temp_data)

time = 
  data.frame(sample_name = colnames(temp_data),
             time1, 
             time2)

####for the time when there are no sampling, just set it as NA
time = 
  time %>% 
  dplyr::mutate(width = as.numeric(difftime(time2, time1, units = "min")))

library(plyr)

time = 
time %>%
  plyr::dlply(.variables = .(sample_name)) %>%
  purrr::map(function(x) {
    if(x$width == 30){
      return(x)
    }else{
      x1 = x
      x2 = x
      x1$time2 = x1$time1 + 30*60
      x1$width = as.numeric(difftime(x1$time2, x1$time1, units = "min"))
      x2$time1 = x1$time2
      x2$width = as.numeric(difftime(x2$time2, x2$time1, units = "min"))
      x2$sample_name = paste(x2$sample_name, 1, sep = "_")
      rbind(x1, x2)
    }
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

temp_data1 = 
  temp_data

add_data = matrix(NA,
                  nrow = nrow(temp_data1),
                  ncol = length(time$sample_name[which(time$width > 30)])) %>%
  as.data.frame()

colnames(add_data) = time$sample_name[which(time$width > 30)]

temp_data1 = 
  cbind(temp_data1, add_data)

temp_data1 =
  temp_data1 %>%
  tibble::rownames_to_column(var = "variable_name") %>%
  dplyr::left_join(cluster_info[,c("variable_id", "cluster")], 
                   by = c("variable_name" = "variable_id")) %>% 
  tibble::rowid_to_column(var = "variable_id") %>%
  dplyr::select(-variable_name) %>% 
  dplyr::select(variable_id, cluster, everything()) %>% 
  dplyr::mutate(variable_id = as.numeric(variable_id)) %>%
  tidyr::pivot_longer(cols = -c(variable_id, cluster),
                      names_to = "sample_name",
                      values_to = "fill") %>%
  dplyr::left_join(time, by = "sample_name")

sum(is.na(temp_data1$fill))

range(temp_data1$fill, na.rm = TRUE)

plot = 
  ggplot(temp_data1, aes(
    xmin = time1,
    xmax = time2,
    ymin = variable_id,
    ymax = variable_id + 1
  )) +
  geom_rect(aes(fill = fill), colour = NA) +
  scale_fill_gradient2(low = "#366A9FFF", 
                       mid = "white", 
                       high = "red", 
                       na.value = alpha(RColorBrewer::brewer.pal(n = 10, name = "RdBu")[6]), 0.1) +
  scale_x_datetime(
    breaks = scales::date_breaks("12 hour"),
    date_labels = "%a %H:%M",
    limits = c(min_time,
               max_time),
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
        panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        panel.spacing = unit(0.1, "lines")) +
  facet_grid(rows = vars(cluster), 
             scales = "free_y",
             space="free")

plot

final_plot = 
plot2 + plot + patchwork::plot_layout(ncol = 1, heights = c(3, 8))

ggsave(final_plot, filename = "molecular_heatmap.pdf", width = 15, height = 10)


#####
###output the node information
dim(variable_info)

variable_info =
  variable_info %>%
  dplyr::left_join(node_info[, c("node", "cluster","module")], 
                   by = c("variable_id" = "node"))

library(openxlsx)
openxlsx::write.xlsx(
  variable_info,
  "variable_info_cluster_info.csv",
  asTable = TRUE,
  overwrite = TRUE
)


variable_info %>%
  plyr::dlply(.variables = .(cluster)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::group_by(data_type)  %>%
      dplyr::summarise(n = n())
  })









