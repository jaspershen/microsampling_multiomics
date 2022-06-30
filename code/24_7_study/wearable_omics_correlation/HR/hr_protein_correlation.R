no_function()

library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")
source("code/modified_dtw.R")
source("code/lagged_correlation.R")

{
  ####load data
  ###HR
  load("data/24_7_study/hr/data_preparation/sample_info")
  load("data/24_7_study/hr/data_preparation/variable_info")
  load("data/24_7_study/hr/data_preparation/expression_data")
  hr_expression_data = expression_data
  hr_sample_info = sample_info
  hr_variable_info = variable_info
  
  ###proteomics
  load("data/24_7_study/proteomics/data_preparation/sample_info")
  load("data/24_7_study/proteomics/data_preparation/variable_info")
  load("data/24_7_study/proteomics/data_preparation/expression_data")
  
  proteomics_sample_info = sample_info
  proteomics_variable_info = variable_info
  proteomics_expression_data = expression_data
  
  load("data/24_7_study/summary_info/day_night_df")
  
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

######hr vs proteomics
dir.create("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_proteomics")
setwd("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_proteomics")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between hr and proteomicss
#global correlation

dir.create("lagged_correlation")

lagged_cor = rep(NA, nrow(proteomics_expression_data))
global_cor = rep(NA, nrow(proteomics_expression_data))

# lagged_result = vector(mode = "list", length = nrow(proteomics_expression_data))
# 
# for(i in 1:nrow(proteomics_expression_data)){
#   cat(i, " ")
#   x = as.numeric(proteomics_expression_data[i, ])
#   time1 = proteomics_sample_info$accurate_time
#   y = as.numeric(hr_expression_data[1, ])
#   time2 = hr_sample_info$accurate_time
# 
#   result = lagged_correlation(
#     x = x,
#     y = y,
#     time1 = time1,
#     time2 = time2,
#     time_tol = 60/60,
#     step = 5/60
#   )
#   lagged_result[[i]] = result
# }
# names(lagged_result) = rownames(proteomics_expression_data)
# save(lagged_result, file = "lagged_correlation/lagged_result")

load("lagged_correlation/lagged_result")

# lagged_cor =
#   lagged_result %>%
#   purrr::map(function(x){
#     x$max_cor
#   }) %>%
#   unlist()
# 
# global_cor =
#   lagged_result %>%
#   purrr::map(function(x){
#     x$global_cor
#   }) %>%
#   unlist()
# 
# shift_time =
#   lagged_result %>%
#   purrr::map(function(x){
#     x$shift_time[x$which_max_idx] %>%
#       stringr::str_replace("\\(", "") %>%
#       stringr::str_replace("\\]", "") %>%
#       stringr::str_split(",") %>%
#       `[[`(1) %>%
#       as.numeric() %>%
#       mean()
# 
#   }) %>%
#   unlist()
# 
# names(lagged_cor) = names(global_cor) =
#   proteomics_variable_info$variable_id
# 
# cor_data =
#   data.frame(wearable = "HR",
#              proteomics_variable_info,
#              global_cor = global_cor,
#              lagged_cor = lagged_cor,
#              shift_time = shift_time) %>% 
#   dplyr::filter(abs(lagged_cor) > 0.2)
# 
# p_value =
#   cor_data %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x){
#     # cat(x[2], " ")
#     x[!is.na(x)] = stringr::str_trim(x[!is.na(x)], side = "both")
#     result = lagged_result[[x[2]]]
# 
#     ###lagged correlation p value
#     x_value = result$x
#     y_value = result$y
# 
#     y_value =
#       result$max_idx %>%
#       purrr::map(function(idx){
#         mean(y_value[idx])
#       }) %>%
#       unlist()
# 
#     x_value = x_value[!is.na(y_value)]
#     y_value = y_value[!is.na(y_value)]
#     lagged_cor_p =
#       cor.test(x = x_value, y = y_value, method = "pearson")$p.value
# 
#     ###global correlation p value
#     x_value = result$x
#     y_value = result$y
# 
#     y_value =
#       result$global_idx %>%
#       purrr::map(function(idx){
#         mean(y_value[idx])
#       }) %>%
#       unlist()
# 
#     x_value = x_value[!is.na(y_value)]
#     y_value = y_value[!is.na(y_value)]
#     global_cor_p =
#       cor.test(x = x_value, y = y_value, method = "pearson")$p.value
# 
#     c(global_cor_p = global_cor_p,
#       lagged_cor_p = lagged_cor_p)
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# cor_data =
#   data.frame(cor_data, p_value)
# 
# cor_data$global_cor_p_adjust = p.adjust(cor_data$global_cor_p, method = "BH")
# cor_data$lagged_cor_p_adjust = p.adjust(cor_data$lagged_cor_p, method = "BH")
# 
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "HR proteomics global cor",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = cor_data,
#                colNames = TRUE, rowNames = TRUE)
# saveWorkbook(wb, "lagged_correlation/cor_data.xlsx", overwrite = TRUE)

cor_data = readxl::read_xlsx("lagged_correlation/cor_data.xlsx")

##output the top 10 negative and top 100 positive
pos_top_10 =
cor_data %>% 
  dplyr::arrange(lagged_cor) %>% 
  dplyr::filter(lagged_cor > 0.3) %>% 
  tail(10)

neg_top_10 =
  cor_data %>% 
  dplyr::arrange(lagged_cor) %>% 
  dplyr::filter(lagged_cor < -0.3) %>% 
  head(10)

dir.create("cor_plot")

temp = 
  rbind(neg_top_10,
        pos_top_10)

# for (i in 1:nrow(temp)) {
#   cat(i, " ")
#   plot1 =
#   lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                         day_night_df = day_night_df,
#                         internal_omics_color = class_color["proteomics"],
#                         wearable_color = wearable_color["hr"],
#                         internal_omics_name = temp$mol_name[i],
#                         warable_name = "HR",
#                         which = "max",
#                         x_limit = c(1,1000),
#                         non_matched_point_size = 0.1,
#                         wearable_point_size = 0.5,
#                         internal_omics_point_size = 2,
#                         integrated = FALSE, add_connect_line = FALSE)
# 
#   plot2 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["proteomics"],
#                           wearable_color = wearable_color["hr"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "HR",
#                           which = "max",
#                           x_limit = c(1,10),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.5,
#                           internal_omics_point_size = 2,
#                           integrated = FALSE)
# 
#   plot3 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["proteomics"],
#                           wearable_color = wearable_color["hr"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "HR",
#                           which = "max",
#                           x_limit = c(1,1000),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.5,
#                           internal_omics_point_size = 2,
#                           integrated = TRUE, add_connect_line = FALSE)
# 
# 
#   plot4 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["proteomics"],
#                           wearable_color = wearable_color["hr"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "HR",
#                           which = "max",
#                           x_limit = c(1,30),
#                           non_matched_point_size = 3,
#                           wearable_point_size = 3,
#                           internal_omics_point_size = 3,
#                           integrated = TRUE)
# 
#   name = paste("HR vs",temp$mol_name[i])
#   ggsave(plot1,
#          filename = file.path("cor_plot", paste(name, "plot1.pdf", sep = "")),
#          width = 20, height = 7)
# 
#   ggsave(plot2,
#          filename = file.path("cor_plot", paste(name, "plot2.pdf", sep = "")),
#          width = 20, height = 7)
# 
#   ggsave(plot3,
#          filename = file.path("cor_plot", paste(name, "plot3.pdf", sep = "")),
#          width = 20, height = 7)
# 
#   ggsave(plot4,
#          filename = file.path("cor_plot", paste(name, "plot4.pdf", sep = "")),
#          width = 20, height = 7)
# }

cor_data %>% 
  ggplot(aes(global_cor, lagged_cor)) +
  geom_point()

##Pathway analysis for the negative or positive correlation analysis
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("pathway_enrichment")

important_protein =
  rbind(
    cor_data %>%
      dplyr::filter(lagged_cor > 0) %>%
      # dplyr::filter(lagged_cor > quantile(lagged_cor, 0.75)) %>%
      dplyr::mutate(class1 = "positive correlation"),
    cor_data %>%
      dplyr::filter(lagged_cor < 0) %>%
      # dplyr::filter(lagged_cor < quantile(lagged_cor, 0.25)) %>%
      dplyr::mutate(class1 = "negative correlation")
  )

# save(important_protein, file = "important_protein")
load("important_protein")

###here we only use the protein with p < 0.05
important_protein = 
  important_protein %>% 
  dplyr::filter(lagged_cor_p_adjust < 0.05)



# ####BP
# negative_correlation_go_bp =
#   clusterProfiler::enrichGO(
#     gene = important_protein$ENTREZID[important_protein$class1 == "negative correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "negative correlation"])],
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     # universe = proteomics_variable_info$ENTREZID,
#     ont = "BP",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE
#   )
# 
# negative_correlation_go_bp =
# clusterProfiler::simplify(
#   negative_correlation_go_bp,
#   cutoff = 0.7,
#   by = "p.adjust",
#   select_fun = min,
#   measure = "Wang"
# )
# 
# ####MF
# negative_correlation_go_mf =
#   clusterProfiler::enrichGO(
#     gene = important_protein$ENTREZID[important_protein$class1 == "negative correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "negative correlation"])],
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "MF",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE
#   )
# 
# negative_correlation_go_mf =
#   clusterProfiler::simplify(
#     negative_correlation_go_mf,
#     cutoff = 0.7,
#     by = "p.adjust",
#     select_fun = min,
#     measure = "Wang"
#   )
# 
# 
# # ####CC
# # negative_correlation_go_cc =
# #   clusterProfiler::enrichGO(
# #     gene = important_protein$ENTREZID[important_protein$class1 == "negative correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "negative correlation"])],
# #     OrgDb = org.Hs.eg.db,
# #     keyType = "ENTREZID",
# #     ont = "CC",
# #     pvalueCutoff = 0.05,
# #     pAdjustMethod = "BH",
# #     qvalueCutoff = 0.05,
# #     readable = TRUE
# #   )
# # 
# # negative_correlation_go_cc =
# #   clusterProfiler::simplify(
# #     negative_correlation_go_cc,
# #     cutoff = 0.7,
# #     by = "p.adjust",
# #     select_fun = min,
# #     measure = "Wang"
# #   )
# 
# negative_correlation_go =
#   clusterProfiler::enrichGO(
#     gene = important_protein$ENTREZID[important_protein$class1 == "negative correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "negative correlation"])],
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE
#   )
# 
# result1 =
#   negative_correlation_go_bp@result %>%
#   data.frame(ONTOLOGY = 'BP') %>%
#   dplyr::select(ONTOLOGY, everything())
# 
# result2 =
#   negative_correlation_go_mf@result %>%
#   data.frame(ONTOLOGY = 'MF') %>%
#   dplyr::select(ONTOLOGY, everything())
# 
# negative_correlation_go@result =
#   rbind(result1,
#         result2) %>%
#   as.data.frame() %>%
#   dplyr::arrange(p.adjust)
# 
# save(negative_correlation_go,
#      file = "pathway_enrichment/negative_correlation_go")

load("pathway_enrichment/negative_correlation_go")


###only remain the terms with at lest 5 genes
negative_correlation_go@result =
  negative_correlation_go@result %>% 
  dplyr::filter(Count > 5)

# ###calculate the similarity of terms
# mf_sim_matrix <-
#   simplifyEnrichment::GO_similarity(go_id = negative_correlation_go@result$ID[negative_correlation_go@result$ONTOLOGY == "MF"],
#                                     ont = "MF",
#                                     measure = "Wang") %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "name1") %>%
#   tidyr::pivot_longer(cols = -name1,
#                       names_to = "name2",
#                       values_to = "sim") %>%
#   dplyr::filter(name1 != name2) %>%
#   dplyr::filter(sim > 0.3)
# 
# name <- apply(mf_sim_matrix, 1, function(x) {
#   paste(sort(x[1:2]), collapse = "_")
# })
# 
# mf_sim_matrix <-
#   mf_sim_matrix %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::arrange(name) %>%
#   dplyr::distinct(name, .keep_all = TRUE) %>%
#   dplyr::select(-name)
# 
# bp_sim_matrix <-
#   simplifyEnrichment::GO_similarity(go_id = negative_correlation_go@result$ID[negative_correlation_go@result$ONTOLOGY == "BP"],
#                                     ont = "BP",
#                                     measure = "Wang") %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "name1") %>%
#   tidyr::pivot_longer(cols = -name1,
#                       names_to = "name2",
#                       values_to = "sim") %>%
#   dplyr::filter(name1 != name2) %>%
#   dplyr::filter(sim > 0.3)
# 
# name <- apply(bp_sim_matrix, 1, function(x) {
#   paste(sort(x[1:2]), collapse = "_")
# })
# 
# bp_sim_matrix <-
#   bp_sim_matrix %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::arrange(name) %>%
#   dplyr::distinct(name, .keep_all = TRUE) %>%
#   dplyr::select(-name)
# 
# # cc_sim_matrix <-
# #   simplifyEnrichment::GO_similarity(go_id = nausea_down_pregnancy_up_go@result$ID[nausea_down_pregnancy_up_go@result$ONTOLOGY == "CC"],
# #                                     ont = "CC",
# #                                     measure = "Wang") %>%
# #   as.data.frame() %>%
# #   tibble::rownames_to_column(var = "name1") %>%
# #   tidyr::pivot_longer(cols = -name1,
# #                       names_to = "name2",
# #                       values_to = "sim") %>%
# #   dplyr::filter(name1 != name2) %>%
# #   dplyr::filter(sim > 0.3)
# #
# # name <- apply(cc_sim_matrix, 1, function(x) {
# #   paste(sort(x[1:2]), collapse = "_")
# # })
# #
# # cc_sim_matrix <-
# #   cc_sim_matrix %>%
# #   dplyr::mutate(name = name) %>%
# #   dplyr::arrange(name) %>%
# #   dplyr::distinct(name, .keep_all = TRUE) %>%
# #   dplyr::select(-name)
# 
# go_sim_matrix <-
#   rbind(bp_sim_matrix, mf_sim_matrix) %>%
#   as.data.frame()
# 
# negative_go_sim_matrix =
#   go_sim_matrix
# 
# save(negative_go_sim_matrix, file = "negative_go_sim_matrix")

load("negative_go_sim_matrix")

edge_data <- 
  rbind(negative_go_sim_matrix) %>% 
  dplyr::rename(from = name1, to = name2)

node_data <-
  rbind(negative_correlation_go@result) %>%
  as.data.frame() %>% 
  dplyr::select(ID, everything()) %>% 
  dplyr::rename(node = ID) 

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)
library(igraph)
subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

cluster1 <-
  purrr::map(cluster, function(x) {
    paste("Cluster", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

cluster1[cluster1 == "Cluster NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(cluster = factor(cluster1, levels = stringr::str_sort(unique(cluster1), numeric = TRUE)))

# ###remove other clusters from the graph project
temp_graph =
  temp_graph %>%
  tidygraph::filter(cluster1 != "Other")

###clustered different GO terms
result_go <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(ONTOLOGY, cluster, p.adjust)

result_go =
  result_go %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::ungroup()

node_for_label = 
  result_go %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::filter(p.adjust == min(p.adjust)) %>% 
  dplyr::slice_head(n = 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::pull(node)

negative_node_for_label = node_for_label

negative_core_term = 
result_go %>% 
  dplyr::filter(node %in% node_for_label)

save(negative_core_term, file = "negative_core_term")

library(ggforce)

plot <-
  temp_graph %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = cluster,
        size = -log(p.adjust, 10)
    ),
    alpha = 1,
    shape = 21,
    show.legend = TRUE
  ) +
  ggforce::geom_mark_hull(aes(x = x, y = y, 
                              group = cluster, 
                              fill = cluster), 
                          show.legend = FALSE,
                          alpha = 0.15) +
  shadowtext::geom_shadowtext(aes(x = x, y = y,
                                  label = ifelse(node %in% node_for_label, Description, NA),
                                  color = cluster), 
                              bg.colour = "white", 
                              check_overlap = FALSE, 
                              show.legend = FALSE) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 7)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) 

plot
extrafont::loadfonts()

# ggsave(
#   plot,
#   filename = file.path(
#     "pathway_enrichment/",
#     paste("negative_sim_plot.pdf", sep = "_")
#   ),
#   width = 9,
#   height = 7
# )




# ####positaive
# ####BP
# positive_correlation_go_bp =
#   clusterProfiler::enrichGO(
#     gene = important_protein$ENTREZID[important_protein$class1 == "positive correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "positive correlation"])],
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "BP",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE
#   )
# 
# positive_correlation_go_bp =
# clusterProfiler::simplify(
#   positive_correlation_go_bp,
#   cutoff = 0.7,
#   by = "p.adjust",
#   select_fun = min,
#   measure = "Wang"
# )
# 
# ####MF
# positive_correlation_go_mf =
#   clusterProfiler::enrichGO(
#     gene = important_protein$ENTREZID[important_protein$class1 == "positive correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "positive correlation"])],
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "MF",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE
#   )
# 
# positive_correlation_go_mf =
#   clusterProfiler::simplify(
#     positive_correlation_go_mf,
#     cutoff = 0.7,
#     by = "p.adjust",
#     select_fun = min,
#     measure = "Wang"
#   )
# 
# 
# # ####CC
# # positive_correlation_go_cc =
# #   clusterProfiler::enrichGO(
# #     gene = important_protein$ENTREZID[important_protein$class1 == "positive correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "positive correlation"])],
# #     OrgDb = org.Hs.eg.db,
# #     keyType = "ENTREZID",
# #     ont = "CC",
# #     pvalueCutoff = 0.05,
# #     pAdjustMethod = "BH",
# #     qvalueCutoff = 0.05,
# #     readable = TRUE
# #   )
# #
# # positive_correlation_go_cc =
# #   clusterProfiler::simplify(
# #     positive_correlation_go_cc,
# #     cutoff = 0.7,
# #     by = "p.adjust",
# #     select_fun = min,
# #     measure = "Wang"
# #   )
# 
# 
# positive_correlation_go =
#   clusterProfiler::enrichGO(
#     gene = important_protein$ENTREZID[important_protein$class1 == "positive correlation"][!is.na(important_protein$ENTREZID[important_protein$class1 == "positive correlation"])],
#     OrgDb = org.Hs.eg.db,
#     keyType = "ENTREZID",
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.05,
#     readable = TRUE
#   )
# 
# result1 =
#   positive_correlation_go_bp@result %>%
#   data.frame(ONTOLOGY = 'BP') %>%
#   dplyr::select(ONTOLOGY, everything())
# 
# result2 =
#   positive_correlation_go_mf@result %>%
#   data.frame(ONTOLOGY = 'MF') %>%
#   dplyr::select(ONTOLOGY, everything())
# 
# positive_correlation_go@result =
#   rbind(result1,
#         result2) %>%
#   as.data.frame() %>%
#   dplyr::arrange(p.adjust)
# 
# save(positive_correlation_go,
#      file = "pathway_enrichment/positive_correlation_go")

load("pathway_enrichment/positive_correlation_go")


###only remain the terms with at lest 5 genes
positive_correlation_go@result =
  positive_correlation_go@result %>% 
  dplyr::filter(Count > 5)

# ###calculate the similarity of terms
# mf_sim_matrix <-
#   simplifyEnrichment::GO_similarity(go_id = positive_correlation_go@result$ID[positive_correlation_go@result$ONTOLOGY == "MF"],
#                                     ont = "MF",
#                                     measure = "Wang") %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "name1") %>%
#   tidyr::pivot_longer(cols = -name1,
#                       names_to = "name2",
#                       values_to = "sim") %>%
#   dplyr::filter(name1 != name2) %>%
#   dplyr::filter(sim > 0.3)
# 
# name <- apply(mf_sim_matrix, 1, function(x) {
#   paste(sort(x[1:2]), collapse = "_")
# })
# 
# mf_sim_matrix <-
#   mf_sim_matrix %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::arrange(name) %>%
#   dplyr::distinct(name, .keep_all = TRUE) %>%
#   dplyr::select(-name)
# 
# bp_sim_matrix <-
#   simplifyEnrichment::GO_similarity(go_id = positive_correlation_go@result$ID[positive_correlation_go@result$ONTOLOGY == "BP"],
#                                     ont = "BP",
#                                     measure = "Wang") %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "name1") %>%
#   tidyr::pivot_longer(cols = -name1,
#                       names_to = "name2",
#                       values_to = "sim") %>%
#   dplyr::filter(name1 != name2) %>%
#   dplyr::filter(sim > 0.3)
# 
# name <- apply(bp_sim_matrix, 1, function(x) {
#   paste(sort(x[1:2]), collapse = "_")
# })
# 
# bp_sim_matrix <-
#   bp_sim_matrix %>%
#   dplyr::mutate(name = name) %>%
#   dplyr::arrange(name) %>%
#   dplyr::distinct(name, .keep_all = TRUE) %>%
#   dplyr::select(-name)
# 
# # cc_sim_matrix <-
# #   simplifyEnrichment::GO_similarity(go_id = nausea_down_pregnancy_up_go@result$ID[nausea_down_pregnancy_up_go@result$ONTOLOGY == "CC"],
# #                                     ont = "CC",
# #                                     measure = "Wang") %>%
# #   as.data.frame() %>%
# #   tibble::rownames_to_column(var = "name1") %>%
# #   tidyr::pivot_longer(cols = -name1,
# #                       names_to = "name2",
# #                       values_to = "sim") %>%
# #   dplyr::filter(name1 != name2) %>%
# #   dplyr::filter(sim > 0.3)
# #
# # name <- apply(cc_sim_matrix, 1, function(x) {
# #   paste(sort(x[1:2]), collapse = "_")
# # })
# #
# # cc_sim_matrix <-
# #   cc_sim_matrix %>%
# #   dplyr::mutate(name = name) %>%
# #   dplyr::arrange(name) %>%
# #   dplyr::distinct(name, .keep_all = TRUE) %>%
# #   dplyr::select(-name)
# 
# go_sim_matrix <-
#   rbind(bp_sim_matrix, mf_sim_matrix) %>%
#   as.data.frame()
# 
# 
# positive_go_sim_matrix =
#   go_sim_matrix
# 
# save(positive_go_sim_matrix, file = "positive_go_sim_matrix")

load("positive_go_sim_matrix")

edge_data <- 
  rbind(positive_go_sim_matrix) %>% 
  dplyr::rename(from = name1, to = name2)

node_data <-
  rbind(positive_correlation_go@result) %>%
  as.data.frame() %>% 
  dplyr::select(ID, everything()) %>% 
  dplyr::rename(node = ID) 

temp_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

library(ggraph)
library(igraph)
subnetwork <-
  igraph::cluster_edge_betweenness(graph = temp_graph,
                                   weights = abs(edge_attr(temp_graph,
                                                           "sim")))
cluster <-
  as.character(membership(subnetwork)) %>%
  purrr::map(function(x) {
    if (sum(x == as.character(membership(subnetwork))) == 1) {
      return("Other")
    } else{
      return(x)
    }
  }) %>%
  unlist()

cluster1 <-
  purrr::map(cluster, function(x) {
    paste("Cluster", match(x, unique(cluster)[unique(cluster) != "Other"]))
  }) %>%
  unlist()

cluster1[cluster1 == "Cluster NA"] <- "Other"

temp_graph <-
  temp_graph %>%
  tidygraph::mutate(cluster = factor(cluster1, levels = stringr::str_sort(unique(cluster1), numeric = TRUE)))

# ###remove other clusters from the graph project
temp_graph =
  temp_graph %>%
  tidygraph::filter(cluster1 != "Other")

###clusterd different GO terms
result_go <-
  igraph::vertex_attr(temp_graph) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
  dplyr::arrange(ONTOLOGY, cluster, p.adjust)

result_go =
  result_go %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::ungroup()

node_for_label = 
  result_go %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::filter(p.adjust == min(p.adjust)) %>% 
  dplyr::slice_head(n = 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::pull(node)

positive_core_term = 
  result_go %>% 
  dplyr::filter(node %in% node_for_label)

save(positive_core_term, file = "positive_core_term")

library(ggforce)

plot <-
  temp_graph %>%
  ggraph(layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    aes(width = sim),
    strength = 1,
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = cluster,
        size = -log(p.adjust, 10)
    ),
    alpha = 1,
    shape = 21,
    show.legend = TRUE
  ) +
  ggforce::geom_mark_hull(aes(x = x, y = y, 
                              group = cluster, 
                              fill = cluster), 
                          show.legend = FALSE,
                          alpha = 0.15) +
  shadowtext::geom_shadowtext(aes(x = x, y = y,
                                  label = ifelse(node %in% node_for_label, Description, NA),
                                  color = cluster), 
                              bg.colour = "white", 
                              check_overlap = FALSE, 
                              show.legend = FALSE) +
  guides(fill = guide_legend(ncol = 1,
                             override.aes = list(size = 3, 
                                                 shape = 21))) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(1, 7)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) 

plot

# ggsave(
#   plot,
#   filename = file.path("pathway_enrichment/",
#                        paste("positive_sim_plot.pdf", sep = "_")
#   ),
#   width = 9,
#   height = 7
# )

