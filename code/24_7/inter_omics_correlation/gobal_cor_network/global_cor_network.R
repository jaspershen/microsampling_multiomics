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
  load("data/7_24_mike/lipidomics/k_means_clustering/new_sample_info")
  load("data/7_24_mike/lipidomics/k_means_clustering/new_variable_info")
  load("data/7_24_mike/lipidomics/k_means_clustering/new_expression_data")
  lipidomics_sample_info = new_sample_info
  lipidomics_variable_info = new_variable_info
  lipidomics_expression_data = new_expression_data
  
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
  load("data/7_24_mike/metabolomics/k_means_clustering/new_sample_info")
  load("data/7_24_mike/metabolomics/k_means_clustering/new_variable_info")
  load("data/7_24_mike/metabolomics/k_means_clustering/new_expression_data")
  metabolomics_sample_info = new_sample_info
  metabolomics_variable_info = new_variable_info
  metabolomics_expression_data = new_expression_data

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
  load("data/7_24_mike/proteomics/k_means_clustering/new_sample_info")
  load("data/7_24_mike/proteomics/k_means_clustering/new_variable_info")
  load("data/7_24_mike/proteomics/k_means_clustering/new_expression_data")
  proteomics_sample_info = new_sample_info
  proteomics_variable_info = new_variable_info
  proteomics_expression_data = new_expression_data
  
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

lipidomics_sample_info$accurate_time == metabolic_panel_sample_info$accurate_time
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

metabolomics_variable_info = 
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
  dplyr::mutate(data_type = "cortisol") %>% 
  dplyr::rename(cortisol_class = class)

total_protein_variable_info = 
total_protein_variable_info %>% 
  dplyr::mutate(data_type = "total_protein") %>% 
  dplyr::rename(total_protein_class = class)

variable_info =
  rbind(
    lipidomics_variable_info[, c("variable_id", "mol_name", "data_type")],
    metabolomics_variable_info[, c("variable_id", "mol_name", "data_type")],
    cytokine_variable_info[, c("variable_id", "mol_name", "data_type")],
    total_protein_variable_info[, c("variable_id", "mol_name", "data_type")],
    cortisol_variable_info[, c("variable_id", "mol_name", "data_type")],
    metabolic_panel_variable_info[, c("variable_id", "mol_name", "data_type")],
    proteomics_variable_info[, c("variable_id", "mol_name", "data_type")]
  )

dim(variable_info)

table(variable_info$data_type)

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) == variable_info$variable_id

}

######
dir.create("data/7_24_mike/inter_omics_correlation/global_cor_network")
setwd("data/7_24_mike/inter_omics_correlation/global_cor_network")

table(variable_info$data_type)

#####get the correlation between each two variables
# cor_data =
#   purrr::map(
#     1:(nrow(expression_data) - 1),
#     .f = function(idx1) {
#       cat(idx1, " ")
#       purrr::map((idx1 + 1):nrow(expression_data),
#                  .f = function(idx2) {
#                    x1 = as.numeric(expression_data[idx1, ])
#                    x2 = as.numeric(expression_data[idx2, ])
#                    test =
#                      cor.test(x1, x2, method = "spearman")
#                    data.frame(
#                      from = variable_info$variable_id[idx1],
#                      to = variable_info$variable_id[idx2],
#                      cor = unname(test$estimate),
#                      p = unname(test$p.value)
#                    )
#                  }
#       ) %>%
#         do.call(rbind, .) %>%
#         as.data.frame()
#     }
#   ) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# dim(cor_data)
# 
# save(cor_data, file = "cor_data")

load("cor_data")

cor_data =
  cor_data %>%
  dplyr::left_join(variable_info, by = c("from" = "variable_id")) %>%
  dplyr::rename(from_mol_name = mol_name,
                from_data_type = data_type) %>%
  dplyr::left_join(variable_info, by = c("to" = "variable_id")) %>%
  dplyr::rename(to_mol_name = mol_name,
                to_data_type = data_type)

table(paste(cor_data$from_data_type, cor_data$to_data_type, sep = "_"))

####add p value
cor_data$p_adjust = p.adjust(cor_data$p, method = "BH")

sum(cor_data$p_adjust < 0.05)

cor_data =
  cor_data %>%
  dplyr::filter(p_adjust < 0.05)

dim(cor_data)
length(unique(c(cor_data$from, cor_data$to)))

table(paste(cor_data$from_data_type, cor_data$to_data_type, sep = "_"))


######if we should set the cutoff for cor to remove some correlations?
# cor_data =
#   cor_data %>%
#   dplyr::filter(abs(cor) > 0.5)

dim(cor_data)
table(paste(cor_data$from_data_type, cor_data$to_data_type, sep = "_"))

edge_data =
  cor_data

edge_data$p_adjust[edge_data$p_adjust == 0] = min(edge_data$p_adjust[edge_data$p_adjust != 0])

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(variable_info, by = c("node" = "variable_id")) %>% 
  dplyr::mutate(node_type = case_when(
    stringr::str_detect(node, "module") ~ "Module",
    TRUE ~ "Molecule"
  ))

dim(node_data)
dim(edge_data)

library(igraph)
library(tidygraph)

total_global_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

library(ggraph)

plot <-
  ggraph(total_global_graph, layout = "kk") +
  geom_edge_link(
    aes(
      color = cor,
      width = -log(p_adjust, 10)
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = data_type,
        size = Degree,
        shape = node_type),
    alpha = 0.7,
    show.legend = TRUE
  ) +
  scale_shape_manual(values = c("Module" = 24,
                                "Molecule" = 21)) +
  scale_size_continuous(range = c(0.5, 8)) +
  scale_fill_manual(values = c(class_color)) +
  scale_color_manual(values = c(class_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Global correlation",
                           override.aes = list(linetype = "none")),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = NA,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    fill = guide_legend(
      title = "Class",
      override.aes = list(shape = 21, size = 3, alpha = 1)
    )
  ) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7),
                                                 "white",
                                                 alpha("#EE0000FF", 0.7))) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

ggsave(
  plot,
  filename = ("total_global_network.pdf"),
  width = 8.6,
  height = 7
)








######example network
load("cor_data")

cor_data =
  cor_data %>%
  dplyr::left_join(variable_info, by = c("from" = "variable_id")) %>%
  dplyr::rename(from_mol_name = mol_name,
                from_data_type = data_type) %>%
  dplyr::left_join(variable_info, by = c("to" = "variable_id")) %>%
  dplyr::rename(to_mol_name = mol_name,
                to_data_type = data_type)

table(paste(cor_data$from_data_type, cor_data$to_data_type, sep = "_"))

####add p value
cor_data$p_adjust = p.adjust(cor_data$p, method = "BH")

sum(cor_data$p_adjust < 0.05)

cor_data =
  cor_data %>%
  dplyr::filter(p_adjust < 0.05)

dim(cor_data)
length(unique(c(cor_data$from, cor_data$to)))

table(paste(cor_data$from_data_type, cor_data$to_data_type, sep = "_"))


#####if we should set the cutoff for cor to remove some correlations?
cor_data =
  cor_data %>%
  dplyr::filter(abs(cor) > 0.5)

dim(cor_data)
table(paste(cor_data$from_data_type, cor_data$to_data_type, sep = "_"))

edge_data =
  cor_data

edge_data$p_adjust[edge_data$p_adjust == 0] = min(edge_data$p_adjust[edge_data$p_adjust != 0])

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(variable_info, by = c("node" = "variable_id")) %>% 
  dplyr::mutate(node_type = case_when(
    stringr::str_detect(node, "module") ~ "Module",
    TRUE ~ "Molecule"
  ))

dim(node_data)
dim(edge_data)

library(igraph)
library(tidygraph)

example_global_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

library(ggraph)

####subnetwork
subnetworks =
  igraph::cluster_edge_betweenness(
    graph = example_global_graph,
    weights = abs(edge_attr(example_global_graph, "cor")),
    membership = TRUE, merges = TRUE)

plot(as.hclust(subnetworks), label=F)

plot = modularity_plot(subnetworks = subnetworks)
plot
ggsave(
  plot,
  filename = "all_modularity.pdf",
  width = 9,
  height = 7
)

###get good cluster number
membership = 
igraph::cut_at(communities = subnetworks, no = 12)

module = paste("Module", membership,
               sep = "_")

example_global_graph =
  example_global_graph %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::mutate(cluster = as.character(module))

plot <-
  ggraph(example_global_graph, layout = "fr") +
  geom_edge_link(
    aes(
      color = cor,
      width = -log(p_adjust, 10)
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = data_type,
        size = Degree,
        shape = node_type),
    alpha = 0.9,
    show.legend = TRUE
  ) +
  scale_size_continuous(range = c(0.5, 5)) +
  scale_fill_manual(values = c(class_color)) +
  scale_color_manual(values = c(class_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Global correlation",
                           override.aes = list(linetype = "none")),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = NA,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    fill = guide_legend(
      title = "Class",
      override.aes = list(shape = 21, size = 3, alpha = 1)
    )
  ) +
  scale_edge_width_continuous(range = c(0.1, 0.5)) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7),
                                                 "white",
                                                 alpha("#EE0000FF", 0.7))) +
  ggnewscale::new_scale_fill() +
  ggforce::geom_mark_hull(aes(
    x = x,
    y = y,
    fill = cluster
  ),
  concavity = 20,
  show.legend = FALSE,
  alpha = 0.15) +
  scale_shape_manual(values = c("Module" = 24,
                                "Molecule" = 21)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# extrafont::loadfonts()

# ggsave(
#   plot,
#   filename = ("example_global_network.pdf"),
#   width = 8.6,
#   height = 7
# )



####output results
node_info =
  igraph::vertex_attr(example_global_graph) %>% 
  dplyr::bind_rows() %>% 
  as.data.frame()

# openxlsx::write.xlsx(node_info,
#                      "node_info.xlsx",
#                      asTable = TRUE,
#                      overwrite = TRUE)

###annotation for each cluster














