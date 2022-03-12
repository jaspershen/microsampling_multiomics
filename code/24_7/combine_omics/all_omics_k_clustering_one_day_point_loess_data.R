no_function()
library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

#####load food log data
{
  load("data/7_24_mike/food_log/data_preparation/expression_data")
  load("data/7_24_mike/food_log/data_preparation/sample_info")
  load("data/7_24_mike/food_log/data_preparation/variable_info")
  food_expression_data = expression_data
  food_variable_info = sample_info
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
  load("data/7_24_mike/combine_omics/data_preparation/new_sample_info")
  load("data/7_24_mike/combine_omics/data_preparation/new_variable_info")
  load("data/7_24_mike/combine_omics/data_preparation/new_expression_data")
  
  expression_data = new_expression_data
  sample_info = new_sample_info
  variable_info = new_variable_info
  
  load("data/7_24_mike/combine_omics/data_preparation/cortisol_variable_info")
  load("data/7_24_mike/combine_omics/data_preparation/cytokine_variable_info")
  load("data/7_24_mike/combine_omics/data_preparation/lipidomics_variable_info")
  load("data/7_24_mike/combine_omics/data_preparation/metabolomics_variable_info")
  load("data/7_24_mike/combine_omics/data_preparation/proteomics_variable_info")
  load("data/7_24_mike/combine_omics/data_preparation/total_protein_variable_info")
  load("data/7_24_mike/combine_omics/data_preparation/metabolic_panel_variable_info")  
}

#####load consistence score
load("data/7_24_mike/circadian_analysis/all_omics/day_consistence/consistence_score")

##load circadian result
load("data/7_24_mike/circadian_analysis/all_omics_loess_data/new_result")

setwd("data/7_24_mike/combine_omics/k_means_clustering_one_day_loess_data")

###only remain the variables with consistence score > 0 and circadian BH < 0.05

cutoff = as.numeric(quantile(consistence_score$new_consistence_score, probs = 0.75))

variable_info = 
variable_info %>% 
  dplyr::left_join(consistence_score, by = "variable_id") %>% 
  dplyr::left_join(new_result %>% dplyr::select(CycID:BH.Q), by = c("variable_id" = "CycID")) %>% 
  dplyr::filter(new_consistence_score > cutoff & BH.Q < 0.05) 

expression_data = 
  expression_data[variable_info$variable_id,]

table(variable_info$data_type)

dim(expression_data)

###scale data for each day
expression_data = 
unique(sample_info$day) %>%
  purrr::map(function(day) {
    temp_expression_data =
      expression_data[, which(sample_info$day == day)]
    
    temp_expression_data =
      temp_expression_data %>% 
      apply(1, function(x){
        (x - mean(x)) / sd(x)
      }) %>% 
      t() %>% 
      as.data.frame()
    
    temp_expression_data
    
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

expression_data = 
  expression_data[,sample_info$sample_id]

####combine samples 
expression_data = 
sample_info %>%
  dplyr::mutate(hour = as.numeric(lubridate::hour(time))) %>%
  plyr::dlply(.variables = .(time)) %>%
  purrr::map(function(x) {
  rowMeans(expression_data[,x$sample_id,drop = FALSE])
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

sample_info =
  data.frame(sample_id = colnames(expression_data)
             # hour = as.numeric(colnames(expression_data))
             ) %>% 
  dplyr::mutate(time = hms::as_hms(colnames(expression_data)))

library(Mfuzz)

temp_data <- 
  expression_data

rownames(temp_data)

time <- c(1:ncol(temp_data))

temp_data <- rbind(time, temp_data)

temp_data2 =
  temp_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

row.names(temp_data)[1] <- "time"

# write.table(
#   temp_data,
#   file = "temp_data.txt",
#   sep = '\t',
#   quote = FALSE,
#   col.names = NA
# )

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

# plot <-
# Dmin(
#   data.s,
#   m = m1,
#   crange = seq(2, 40, 1),
#   repeats = 3,
#   visu = TRUE
# )
# 
# plot <-
# plot %>%
#   data.frame(distance = plot,
#              k = seq(2,40,1)) %>%
#   ggplot(aes(k, distance)) +
#   geom_point(shape = 21, size = 4, fill = "black") +
#   geom_smooth() +
#   geom_segment(aes(x = k, y = 0, xend = k, yend = distance)) +
#   theme_bw() +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) +
#   labs(
#     x = "Cluster number",
#     y = "Min. centroid distance"
#   ) +
#   scale_y_continuous(expand = expansion(mult = c(0,0.1)))
# 
# plot
# 
# ggsave(plot, filename = "distance_k_number.pdf", width = 7, height = 7)

clust = 3

c <- mfuzz(data.s, c = clust, m = m1)
# mfuzz.plot(eset = data.s,
#            min.mem = 0.8,
#            cl = c,
#            mfrow=c(3,4),
#            time.labels = time,
#            new.window = FALSE)
# 
# names(c$cluster) <- rownames(temp_data2)[-1]
# rownames(c$membership) <- rownames(temp_data2)[-1]
# save(c, file = "c")
load("c")

####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
center <- c$centers
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

acore <- acore(data.s,c,min.acore=0)
acore

centers <- c$centers
names(c$cluster) == rownames(c$membership)

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

openxlsx::write.xlsx(x = cluster_info,
                     file = "cluster_info.xlsx",
                     asTable = TRUE, overwrite = TRUE)

####output the expression data of different clusters
##plot for each cluster

day_night_df = 
  data.frame(start = hms::hms(seconds = 0,minutes = 0,hours = 6),
             end = hms::hms(seconds = 0,minutes = 0,hours = 18))

# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
#   dir.create(paste("cluster", cluster_idx, sep = "_"))
#   cluster_data <-
#     cluster_info %>%
#     dplyr::filter(cluster_idx == cluster_idx) %>%
#     dplyr::select(1, 1 + cluster_idx)
#   
#   colnames(cluster_data) <- c("variable_id", "membership")
#   
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.9)
#   
#   openxlsx::write.xlsx(
#     x = cluster_data,
#     file = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".xlsx", sep = "")
#     ),
#     asTable = TRUE,
#     overwrite = TRUE
#   )
#   
#   ###cluster plot
#   
#   temp =
#     temp_data2[cluster_data$variable_id,] %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership),
#       names_to = "sample_id",
#       values_to = "value"
#     ) %>%
#     dplyr::left_join(sample_info[, c("sample_id", "time")], by = "sample_id")
#   
#   plot <-
#     ggplot() +
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
#     geom_hline(yintercept = 0) +
#     geom_line(aes(x = time, y = value,
#                   group = variable_id,
#                   color = membership),
#               data = temp) +
#     # scale_x_datetime(
#     #   breaks = scales::date_breaks("12 hour"),
#     #   date_labels = "%a %H:%M",
#     #   timezone = "America/Los_Angeles"
#     # ) +
#     scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1,
#         vjust = 1,
#         size = 12
#       ),
#       panel.background = element_rect(fill = alpha("grey", 0.2)),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste(
#         "Cluster ",
#         cluster_idx,
#         " (",
#         nrow(cluster_data),
#         " molecules)",
#         sep = ""
#       )
#     )
#   
#   plot
#   
#   ggsave(
#     plot,
#     filename = file.path(
#       paste("cluster", cluster_idx, sep = "_"),
#       paste("cluster", cluster_idx, ".pdf", sep = "")
#     ),
#     width = 9,
#     height = 7
#   )
# }


## (4) feature number
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx") 
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx") 
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx") 

###annotation for each cluster
cluster1 = data.frame(cluster1, cluster = "1")
cluster2 = data.frame(cluster2, cluster = "2")
cluster3 = data.frame(cluster3, cluster = "3")

cluster = 
  rbind(cluster1,
        cluster2,
        cluster3)

variable_info = 
  variable_info %>% 
  dplyr::left_join(cluster, by = "variable_id")

variable_info$cluster[is.na(variable_info$cluster)] = "Other"

variable_info$mol_name[!is.na(variable_info$Lipid_Name)] =
  variable_info$Lipid_Name[!is.na(variable_info$Lipid_Name)]

#load global correlation
###network for each cluster
# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
# 
#   temp_variable_id =
#   variable_info %>%
#     dplyr::filter(cluster == cluster_idx) %>%
#     dplyr::pull(variable_id)
# 
#   ###calculate correlation
#   library(corrr)
#   cor_data =
#   corrr::correlate(x = t(expression_data[temp_variable_id, ]),
#                    method = "spearman",
#                    quiet = TRUE)
# 
#   cor_data =
#   cor_data %>%
#     shave() %>%
#     tidyr::pivot_longer(cols = -term, names_to = "to", values_to = "cor") %>%
#     dplyr::filter(!is.na(cor)) %>%
#     dplyr::rename(from = term)
# 
# 
#   edge_data =
#     cor_data %>%
#     dplyr::filter(abs(cor) > 0.7)
# 
#   node_data =
#     variable_info %>%
#     dplyr::rename(node = variable_id) %>%
#     dplyr::filter(node %in% c(edge_data$from, edge_data$to)) %>%
#     dplyr::filter(cluster == cluster_idx)
# 
#   edge_data =
#     edge_data %>%
#     dplyr::filter(from %in% node_data$node & to %in% node_data$node)
# 
#   node_data =
#     node_data %>%
#     dplyr::filter(node %in% c(edge_data$from, edge_data$to))
# 
#   if (nrow(node_data) > 0) {
#     library(ggraph)
#     library(igraph)
#     library(tidygraph)
# 
#     graph <-
#       tidygraph::tbl_graph(nodes = node_data,
#                            edges = edge_data,
#                            directed = FALSE) %>%
#       dplyr::mutate(Degree = centrality_degree(mode = 'all'))
# 
#     subnetworks <-
#       igraph::cluster_fast_greedy(graph = graph,
#                                        weights = abs(edge_attr(graph,
#                                                                "cor")))
#     plot =
#       modularity_plot(subnetworks = subnetworks)
# 
#     ggsave(
#       plot,
#       filename = file.path(
#         paste("cluster", cluster_idx, sep = "_"),
#         "all_modularity.pdf"
#       ),
#       width = 9,
#       height = 7
#     )
# 
#     if (max(subnetworks$modularity) < 0.4) {
#       module = paste(cluster_idx, rep(1, length(
#         membership(communities = subnetworks)
#       )), sep = "_")
#     } else{
#       module = paste(cluster_idx, membership(communities = subnetworks), sep = "_")
#     }
# 
#     graph =
#       graph %>%
#       tidygraph::activate(what = "nodes") %>%
#       dplyr::mutate(module = module)
# 
#     node =
#       vertex_attr(graph) %>%
#       dplyr::bind_cols() %>%
#       as.data.frame()
# 
#     save(node, file = file.path(paste("cluster", cluster_idx, sep = "_"), "node"))
# 
#     plot <-
#       ggraph(graph, layout = "fr",
#              circular = FALSE) +
#       geom_edge_link(aes(color = cor),
#                      alpha = 1,
#                      show.legend = TRUE) +
#       geom_node_point(
#         aes(size = Degree,
#             fill = data_type),
#         shape = 21,
#         alpha = 0.7,
#         # fill = class_color["lipidomics"],
#         show.legend = FALSE
#       ) +
#       scale_fill_manual(values = class_color) +
#       ggnewscale::new_scale_fill() +
#       ggforce::geom_mark_hull(
#         aes(
#           x = x,
#           y = y,
#           group = module,
#           fill = module
#         ),
#         show.legend = FALSE,
#         alpha = 0.15
#       ) +
#       scale_size_continuous(range = c(2, 8)) +
#       guides(
#         linetype = "none",
#         color = guide_colorbar(
#           title = "Lagged correlation",
#           override.aes = list(linetype = "none")
#         ),
#         size = guide_legend(
#           title = "Degree",
#           override.aes = list(
#             linetype = NA,
#             fill = "transparent",
#             shape = 21,
#             color = "black"
#           )
#         ),
#         fill = guide_legend(
#           title = "Class",
#           override.aes = list(
#             shape = 21,
#             size = 3,
#             alpha = 1
#           )
#         )
#       ) +
#       scale_edge_width_continuous(range = c(0.3, 2)) +
#       ggraph::scale_edge_color_gradient2(
#         low = alpha("#3B4992FF", 0.7),
#         mid = "white",
#         high = alpha("#EE0000FF", 0.7)
#       ) +
#       ggraph::theme_graph() +
#       theme(
#         plot.background = element_rect(fill = "transparent", color = NA),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         legend.position = "right",
#         legend.background = element_rect(fill = "transparent", color = NA)
#       )
# 
#     if (nrow(node_data) < 150) {
#       plot =
#         plot +
#         geom_node_text(aes(x = x,
#                            y = y,
#                            label = mol_name),
#                        size = 2,
#                        check_overlap = TRUE)
#     }
# 
#     # extrafont::loadfonts()
# 
#     plot
#     ggsave(
#       plot,
#       filename = file.path(paste("cluster", cluster_idx, sep = "_"), "network.pdf"),
#       width = 9,
#       height = 7
#     )
#   } else{
#     node = node_data
#     save(node, file = file.path(paste("cluster", cluster_idx, sep = "_"), "node"))
#   }
# }

#####load the new node information
load("cluster_1/node")
node1 = node
load("cluster_2/node")
node2 = node
load("cluster_3/node")
node3 = node

node_info = 
  rbind(node1, 
        node2,
        node3)

###plot for each module
# for (module_idx in 1:length(unique(node_info$module))) {
#   cat(module_idx, " ")
#   module = unique(node_info$module)[module_idx]
#   cluster_idx = stringr::str_split(module, pattern = "_")[[1]][1]
#   ###cluster plot
# 
#   temp =
#     temp_data2[node_info$node[node_info$module == module], ] %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id),
#       names_to = "sample_id",
#       values_to = "value"
#     ) %>%
#     dplyr::left_join(sample_info[, c("sample_id", "time")], by = "sample_id") %>%
#     dplyr::left_join(variable_info[, c("variable_id", "data_type")], by = "variable_id")
# 
#   temp_data_type =
#   temp %>%
#     dplyr::distinct(variable_id, .keep_all = TRUE) %>%
#     pull(data_type)
# 
#   plot <-
#     ggplot() +
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
#     geom_hline(yintercept = 0) +
#     geom_line(aes(time,
#                   value,
#                   group = variable_id,
#                   color = data_type),
#               data = temp,
#               show.legend = FALSE) +
#     # scale_x_datetime(
#     #   breaks = scales::date_breaks("12 hour"),
#     #   date_labels = "%a %H:%M",
#     #   timezone = "America/Los_Angeles"
#     # ) +
#     scale_color_manual(values = class_color) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1,
#         vjust = 1,
#         size = 12
#       ),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste("Module ",
#                     module,
#                     " (",
#                     length(unique(temp$variable_id)),
#                     " molecules:",
#                     paste(paste(names(table(temp_data_type)),
#                                 table(temp_data_type)),
#                           collapse = ','),
#                     ")",
#                     sep = "")
#     ) +
#     theme(panel.background = element_rect(fill = alpha("grey", 0.2)))
# 
#   plot
# 
#   ggsave(
#     plot,
#     filename = file.path(paste("cluster", cluster_idx, sep = "_"),
#                          paste("module", module, ".pdf", sep = "")),
#     width = 9,
#     height = 7
#   )
# }

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
time = sample_info$time

min_time = hms::hms(seconds = 0, minutes = 0, hours = 0)
max_time = hms::hms(seconds = 0, minutes = 0, hours = 24)

# time[length(time)] = time[length(time)] + 30*60

colnames(temp_data)

time1 = c(hms::as_hms(time[1] - 30*60), time[-length(time)])
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
  # scale_x_datetime(
  #   breaks = scales::date_breaks("12 hour"),
  #   date_labels = "%a %H:%M",
  #   limits = c(min_time,
  #              max_time),
  #   timezone = "America/Los_Angeles"
  # ) +
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

# ggsave(plot, filename = "molecular_heatmap.pdf", width = 15, height = 5)


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


node_info %>% 
  dplyr::filter(cluster == 3) %>% 
  dplyr::pull(node) %>% 
  purrr::map(function(x){
    name = 
    grep(pattern = paste(x,".pdf",sep=""), 
         dir("../../circadian_analysis/all_omics_loess_data/plot/"), 
         value = TRUE)
    file.copy(file.path("../../circadian_analysis/all_omics_loess_data/plot/",name),
                to = "cluster_3/plot")
  })












