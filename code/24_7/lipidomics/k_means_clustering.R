no_function()
library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

{
  ####this is for the day night time
  load("data/7_24_mike/summary_info/day_night_df")
  day_night_df =
    day_night_df %>%
    dplyr::mutate(
      start_time = as.POSIXct(hms::as_hms(start)),
      end_time = as.POSIXct(hms::as_hms(end)),
      week = format(day, "%a")
    ) %>%
    dplyr::mutate(week = paste(week,
                               lubridate::month(day),
                               lubridate::day(day),
                               sep = "-")) %>%
    dplyr::mutate(week = factor(week, unique(week)))
  
  ###lipidomics
  load("data/7_24_mike/lipidomics/data_preparation/sample_info")
  load("data/7_24_mike/lipidomics/data_preparation/variable_info")
  load("data/7_24_mike/lipidomics/data_preparation/expression_data")
  
  remain_idx =
    which(apply(expression_data, 1, function(x) {
      sum(x == 0) / ncol(expression_data)
    }) < 0.5)
  
  variable_info =
    variable_info[remain_idx, ]
  
  expression_data =
    expression_data[variable_info$variable_id, ]
}

######
setwd("data/7_24_mike/lipidomics/k_means_clustering")

dim(expression_data)

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

clust = 9

# c <- mfuzz(data.s, c = clust, m = m1)
# 
# mfuzz.plot(eset = data.s,
#            min.mem = 0.8,
#            cl = c,
#            mfrow=c(3,3),
#            time.labels = time,
#            new.window = FALSE)
# 
# names(c$cluster) <- rownames(temp_data2)
# rownames(c$membership) <- rownames(temp_data2)
# save(c, file = "c")
load("c")

####any two clusters with correlation > 0.8 should be considered as one
# library(corrplot)
# layout(1)
# center <- c$centers
# rownames(center) <- paste("Cluster", rownames(center), sep = ' ')
# 
# cor_data = get_cor_matrix(
#   data = temp_data,
#   c = c,
#   scale = TRUE,
#   method = "spearman",
#   which = "median"
# )
# 
# head(cor_data)
# 
# cor_data_median = 
# cor_data %>% 
#   dplyr::group_by(from, to) %>% 
#   dplyr::summarise(median = median(cor)) %>% 
#   dplyr::ungroup()
# 
# library(gghalves)
# plot = 
# cor_data %>%
#   dplyr::group_by(from, to) %>%
#   dplyr::mutate(median = median(cor)) %>%
#   ggplot(aes(x = "1", cor)) +
#   geom_half_boxplot(aes(x = 1, y = cor, fill = median),
#                     outlier.shape = NA,
#                     side = "l") +
#   geom_half_violin(aes(x = 1, y = cor, fill = median),
#                    outlier.shape = NA,
#                    side = "r") +
#   scale_fill_gradient2(low = ggsci::pal_aaas()(n = 10)[2],
#                        mid = "white",
#                        high = "red") +
#   labs(x = "", y = "Spearman correlation") +
#   base_theme +
#   facet_grid(rows = vars(from), cols = vars(to)) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.grid = element_blank(),
#     axis.text.y = element_text(size = 10)
#   ) 
# 
# plot
# 
# cor_data_median2 = 
#   cor_data_median
# 
# cor_data_median2$from = 
#   cor_data_median$to
# 
# cor_data_median2$to = 
#   cor_data_median$from
# 
# plot = 
# plot +
#   # geom_point(
#   #   data = cor_data_median2,
#   #   mapping = aes(
#   #     x = 1,
#   #     y = median(range(cor_data$cor)),
#   #     size = abs(median),
#   #     fill = median
#   #   ),
#   #   shape = 21,
#   #   show.legend = FALSE
#   # ) +
#   # ggnewscale::new_scale_fill() +
#   # scale_size_continuous(range = c(1, 20)) +
#   # scale_fill_gradient2(low = ggsci::pal_aaas()(n = 10)[2],
#   #                      mid = "white",
#   #                      high = "red") +
#   geom_text(data = cor_data_median2,
#             mapping = aes(
#               x = 1,
#               y = median(range(cor_data$cor)),
#               label = round(median, 2)
#             ))
# 
# plot


# ggsave(plot, filename = "cor_plot.pdf", width = 8, height = 7)

# corrplot::corrplot(
#   corr = cor(t(center)),
#   type = "full",
#   diag = TRUE,
#   order = "hclust",
#   hclust.method = "ward.D",
#   # addrect = 5,
#   col = colorRampPalette(colors = rev(
#     RColorBrewer::brewer.pal(n = 11, name = "Spectral")
#   ))(n = 100),
#   number.cex = .7,
#   addCoef.col = "black"
# )

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

# openxlsx::write.xlsx(x = cluster_info,
#                      file = "cluster_info.xlsx", asTable = TRUE)


#####output the expression data of different clusters

##plot for each cluster
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
#   openxlsx::write.xlsx(x = cluster_data,
#                        file = file.path(
#                          paste("cluster", cluster_idx, sep = "_"),
#                          paste("cluster", cluster_idx, ".xlsx", sep = "")
#                        ),
#                        asTable = TRUE, overwrite = TRUE)
# 
# ###cluster plot
#   
#   temp =
#     temp_data2[cluster_data$variable_id, ] %>%
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
#     dplyr::left_join(sample_info[, c("sample_id", "accurate_time")], by = "sample_id")
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
#     geom_line(aes(accurate_time, value, group = variable_id, color = membership),
#               data = temp) +
#     scale_x_datetime(
#       breaks = scales::date_breaks("12 hour"),
#       date_labels = "%a %H:%M",
#       timezone = "America/Los_Angeles"
#     ) +
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
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste("Cluster ",
#                     cluster_idx,
#                     " (",
#                     nrow(cluster_data),
#                     " lipids)",
#                     sep = "")
#     ) 
#   
#   plot
# 
#   ggsave(
#     plot,
#     filename = file.path(paste("cluster", cluster_idx, sep = "_"),
#                          paste("cluster", cluster_idx, ".pdf", sep = "")),
#     width = 20,
#     height = 7
#   )
# }

## (4) feature number
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx") 
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx") 
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx") 
cluster4 <- readxl::read_xlsx("cluster_4/cluster4.xlsx") 
cluster5 <- readxl::read_xlsx("cluster_5/cluster5.xlsx") 
cluster6 <- readxl::read_xlsx("cluster_6/cluster6.xlsx") 
cluster7 <- readxl::read_xlsx("cluster_7/cluster7.xlsx") 
cluster8 <- readxl::read_xlsx("cluster_8/cluster8.xlsx") 
cluster9 <- readxl::read_xlsx("cluster_9/cluster9.xlsx") 

###annotation for each cluster
cluster1 = data.frame(cluster1, cluster = "1")
cluster2 = data.frame(cluster2, cluster = "2")
cluster3 = data.frame(cluster3, cluster = "3")
cluster4 = data.frame(cluster4, cluster = "4")
cluster5 = data.frame(cluster5, cluster = "5")
cluster6 = data.frame(cluster6, cluster = "6")
cluster7 = data.frame(cluster7, cluster = "7")
cluster8 = data.frame(cluster8, cluster = "8")
cluster9 = data.frame(cluster9, cluster = "9")

cluster = 
  rbind(cluster1,
        cluster2,
        cluster3,
        cluster4,
        cluster5,
        cluster6,
        cluster7,
        cluster8,
        cluster9)

variable_info = 
  variable_info %>% 
  dplyr::left_join(cluster, by = "variable_id")

variable_info$cluster[is.na(variable_info$cluster)] = "Other"

variable_info$mol_name[!is.na(variable_info$Lipid_Name)] =
  variable_info$Lipid_Name[!is.na(variable_info$Lipid_Name)]

####load global correlation
load(here::here("data/7_24_mike/inter_omics_correlation/inter_all_omics/global correlation network/global_cor"))
global_cor

##network for each cluster
# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
#   edge_data =
#     global_cor %>%
#     dplyr::filter(from_data_type == "lipidomics" &
#                     to_data_type == "lipidomics") %>%
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
#       igraph::cluster_edge_betweenness(graph = graph,
#                                        weights = abs(edge_attr(graph,
#                                                                "cor")))
#     
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
#       geom_edge_link(aes(color = cor,
#                          width = -log(p_adjust, 10)),
#                      alpha = 1,
#                      show.legend = TRUE) +
#       geom_node_point(
#         aes(size = Degree,
#             fill = module),
#         shape = 21,
#         alpha = 0.7,
#         # fill = class_color["lipidomics"],
#         show.legend = FALSE
#       ) +
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
load("cluster_4/node")
node4 = node
load("cluster_5/node")
node5 = node
load("cluster_6/node")
node6 = node
load("cluster_7/node")
node7 = node
load("cluster_8/node")
node8 = node
load("cluster_9/node")
node9 = node

node_info = 
  rbind(node1, 
        node2,
        node3,
        node4,
        node5,
        node6,
        node7,
        node8,
        node9)

# #####plot for each module
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
#     dplyr::left_join(sample_info[, c("sample_id", "accurate_time")], by = "sample_id")
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
#     geom_line(aes(accurate_time, value, group = variable_id),
#               data = temp) +
#     scale_x_datetime(
#       breaks = scales::date_breaks("12 hour"),
#       date_labels = "%a %H:%M",
#       timezone = "America/Los_Angeles"
#     ) +
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
#                     " lipids)",
#                     sep = "")
#     )
# 
#   plot
# 
#   ggsave(
#     plot,
#     filename = file.path(paste("cluster", cluster_idx, sep = "_"),
#                          paste("module", module, ".pdf", sep = "")),
#     width = 20,
#     height = 7
#   )
# }




#####output the combined lipidomics data
dim(node_info)
dim(variable_info)

final_info =
  variable_info %>% 
  dplyr::left_join(node_info[,c("node", "module")], by = c("variable_id" = "node"))

final_info$module[is.na(final_info$module)] = "Other"

final_info %>% 
  dplyr::filter(module == "1_1") %>% 
  dplyr::pull(subclass)

final_info %>% 
  dplyr::filter(module == "1_2") %>% 
  dplyr::pull(subclass)

library(plyr)
final_info %>% 
  plyr::dlply(.(module)) %>% 
  purrr::map(function(x){
    table(x$subclass)
  })

###output
library(openxlsx)
final_info = 
  final_info %>% 
  dplyr::arrange(module, subclass)

openxlsx::write.xlsx(final_info, file = "final_info.xlsx", asTable = TRUE, overwrite = TRUE)  

final_info =
  readxl::read_xlsx("final_info_manual.xlsx")

table(final_info$module)

#####annotation for each module
library(plyr)

final_info %>% 
  plyr::dlply(.variables = .(module)) %>% 
  purrr::map(function(x){
    x[,c("mol_name", "subclass")]
  })

dim(expression_data)
dim(final_info)
dim(sample_info)

rownames(expression_data) == final_info$variable_id

new_variable_info = 
final_info %>%
  plyr::dlply(.variables = .(module)) %>%
  purrr::map(function(x) {
    cat(x$variable_id[1], " ")
    if(x$module[1] == "Other"){
      x$annotation = x$mol_name
      return(x)
    }
    if (nrow(x) == 1) {
      x$annotation = x$mol_name
    } else{
      x$variable_id = paste(x$variable_id, collapse = ";")
      x$annotation = paste(c("Lipidomics module", unique(x$module), unique(x$subclass)), collapse = "_")
      x$mol_name = paste(x$mol_name, collapse = ";")
      x$subclass = paste(x$subclass, collapse = ";")
      x$Lipid_Name = paste(x$Lipid_Name, collapse = ";")
      x$Lipid_Name_Tot = paste(x$Lipid_Name_Tot, collapse = ";")
      x$Lipid_Name_LION = paste(x$Lipid_Name_LION, collapse = ";")
      x$Lipid_Class_Detailed = paste(x$Lipid_Class_Detailed, collapse = ";")
      x$Lipid_Class = paste(x$Lipid_Class, collapse = ";")
      x$FA_1 = paste(x$FA_1, collapse = ";")
      x$FA_2 = paste(x$FA_2, collapse = ";")
      x$FA_All = paste(x$FA_All, collapse = ";")
      x$Sat_status_FA1 = paste(x$Sat_status_FA1, collapse = ";")
      x$Sat_status_FA2 = paste(x$Sat_status_FA2, collapse = ";")
      x$Saturation_All = paste(x$Saturation_All, collapse = ";")
      x$OddEven_FA1 = paste(x$OddEven_FA1, collapse = ";")
      x$OddEven_FA2 = paste(x$OddEven_FA2, collapse = ";")
      x$OddEven_All = paste(x$OddEven_All, collapse = ";")
      x$Omega_FA1 = paste(x$Omega_FA1, collapse = ";")
      x$Omega_FA2 = paste(x$Omega_FA2, collapse = ";")
      x$Omega_All = paste(x$Omega_All, collapse = ";")
      x$HMDB_ID = paste(x$HMDB_ID, collapse = ";")
      x$KEGG_ID = paste(x$KEGG_ID, collapse = ";")
      x$LION_ID = paste(x$LION_ID, collapse = ";")
      x$cellular_component = paste(x$cellular_component, collapse = ";")
      x$Function = paste(x$Function, collapse = ";")
      x$Total_Carb = paste(x$Total_Carb, collapse = ";")
      x$Total_Unsat = paste(x$Total_Unsat, collapse = ";")
      x$Carbon_numb_FA1 = paste(x$Carbon_numb_FA1, collapse = ";")
      x$Carbon_numb_FA2 = paste(x$Carbon_numb_FA2, collapse = ";")
      x$CarbonFA_All = paste(x$CarbonFA_All, collapse = ";")
      x$Unsat_numb_FA1 = paste(x$Unsat_numb_FA1, collapse = ";")
      x$Unsat_numb_FA2 = paste(x$Unsat_numb_FA2, collapse = ";")
      x$UnsatFA_All = paste(x$UnsatFA_All, collapse = ";")
      x$PhysicoChemical = paste(x$PhysicoChemical, collapse = ";")
      x$PathwaysKegg = paste(x$PathwaysKegg, collapse = ";")
      x$membership = paste(x$membership, collapse = ";")
      x$cluster = paste(x$cluster, collapse = ";")
      x$module = paste(x$module, collapse = ";")
      x =
        x %>%
        dplyr::distinct(variable_id, .keep_all = TRUE)
    }
    x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

new_expression_data = 
new_variable_info$variable_id %>% 
  purrr::map(function(x){
    x = stringr::str_split(x, pattern = ";")[[1]]
    expression_data[x,] %>% 
      colSums()
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(new_expression_data) = new_variable_info$annotation

new_variable_info$old_variable_id = new_variable_info$variable_id
new_variable_info$variable_id = new_variable_info$annotation

new_sample_info = sample_info

colnames(new_expression_data) == new_sample_info$sample_id
rownames(new_expression_data) == new_variable_info$variable_id

save(new_expression_data, file = "new_expression_data")
save(new_sample_info, file = "new_sample_info")
save(new_variable_info, file = "new_variable_info")









