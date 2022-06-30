#to avoind source
no_exist_function()

masstools::setwd_project()
rm(list = ls())
library(tidyverse)
source("R/tools.R")

###load data
##metabolome
load("data/ensure_shake_study/metabolome_data_analysis/data_preparation/expression_data")
load("data/ensure_shake_study/metabolome_data_analysis/data_preparation/variable_info")
load("data/ensure_shake_study/metabolome_data_analysis/data_preparation/sample_info")
load("data/ensure_shake_study/metabolome_data_analysis/metabolites/DEG/anova_marker_name")

variable_info <- 
  variable_info %>% 
  dplyr::filter(!is.na(Metabolite)) %>% 
  dplyr::filter(variable_id %in% anova_marker_name)

expression_data <- 
  expression_data[match(variable_info$variable_id, rownames(expression_data)),] %>% 
  `+`(1) %>% 
  log(2)

metabolome_expression_data <- expression_data
metabolome_variable_info <- variable_info
metabolome_sample_info <- sample_info

##lipidomics
load("data/ensure_shake_study/lipidomics_data_analysis/data_preparation/expression_data")
load("data/ensure_shake_study/lipidomics_data_analysis/data_preparation/variable_info")
load("data/ensure_shake_study/lipidomics_data_analysis/data_preparation/sample_info")
load("data/ensure_shake_study/lipidomics_data_analysis/DEG/anova_marker_name")

variable_info <- 
  variable_info %>% 
  dplyr::filter(variable_id %in% anova_marker_name)

expression_data <- 
  expression_data[match(variable_info$variable_id, rownames(expression_data)),] +
  `+`(1) %>% 
  log(2)

lipidomics_expression_data <- expression_data
lipidomics_variable_info <- variable_info
lipidomics_sample_info <- sample_info

##cytokine
load("data/ensure_shake_study/cytokine_data_analysis/data_preparation/expression_data")
load("data/ensure_shake_study/cytokine_data_analysis/data_preparation/variable_info")
load("data/ensure_shake_study/cytokine_data_analysis/data_preparation/sample_info")
load("data/ensure_shake_study/cytokine_data_analysis/DEG/anova_marker_name")

variable_info <- 
  variable_info %>% 
  dplyr::filter(variable_id %in% anova_marker_name)

expression_data <- 
  expression_data[match(variable_info$variable_id, rownames(expression_data)),] +
  `+`(1) %>% 
  log(2)

cytokine_expression_data <- expression_data
cytokine_variable_info <- variable_info
cytokine_sample_info <- sample_info

masstools::setwd_project()

setwd(
  "data/ensure_shake_study/correlation_network/inter_omics"
)


# ####calculate correlation
# ###metabolome vs other omics data
# ##met vs other omics
# library(corrr)
# 
# ##met vs lipidomics
# data1 <- metabolome_expression_data
# data2 <- lipidomics_expression_data
# 
# intersect_name <- intersect(colnames(data1), colnames(data2))
# 
# data1 <- as.data.frame(data1)[,intersect_name]
# data2 <- as.data.frame(data2)[,intersect_name]
# 
# metabolome_lipidomics_cor <-
#   cor(x = t(data1),
#       y = t(data2),
#       method = "spearman", use = "pairwise.complete.obs")
# 
# metabolome_lipidomics_cor <-
#   metabolome_lipidomics_cor %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "x") %>%
#   tidyr::pivot_longer(cols = -x,
#                       names_to = "y",
#                       values_to = "r")
# 
# p <-
#   as.data.frame(t(metabolome_lipidomics_cor)) %>%
#   purrr::map(.f = function(x){
#     cor.test(as.numeric(data1[x[1],]),
#              as.numeric(data2[x[2],]),
#              method = "spearman"
#     )$p.value
#   }) %>%
#   unlist()
# 
# fdr <- p.adjust(p, method = "fdr")
# 
# fdr[fdr == 0] <- min(fdr[fdr != 0])
# 
# metabolome_lipidomics_cor <-
#   data.frame(metabolome_lipidomics_cor, p, fdr, stringsAsFactors = FALSE) %>%
#   dplyr::filter(fdr < 0.05)
# 
# dim(metabolome_lipidomics_cor)
# 
# save(metabolome_lipidomics_cor, file = "metabolome_lipidomics_cor")
# load("metabolome_lipidomics_cor")
# 
# 
# ##met vs cytokine
# data1 <- metabolome_expression_data
# data2 <- cytokine_expression_data
# 
# intersect_name <- intersect(colnames(data1), colnames(data2))
# 
# data1 <- as.data.frame(data1)[,intersect_name]
# data2 <- as.data.frame(data2)[,intersect_name]
# 
# metabolome_cytokine_cor <-
#   cor(x = t(data1),
#       y = t(data2),
#       method = "spearman", use = "pairwise.complete.obs")
# 
# metabolome_cytokine_cor <-
#   metabolome_cytokine_cor %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "x") %>%
#   tidyr::pivot_longer(cols = -x,
#                       names_to = "y",
#                       values_to = "r")
# 
# p <-
#   as.data.frame(t(metabolome_cytokine_cor)) %>%
#   purrr::map(.f = function(x){
#     cor.test(as.numeric(data1[x[1],]),
#              as.numeric(data2[x[2],]),
#              method = "spearman"
#     )$p.value
#   }) %>%
#   unlist()
# 
# fdr <- p.adjust(p, method = "fdr")
# 
# fdr[fdr == 0] <- min(fdr[fdr != 0])
# 
# metabolome_cytokine_cor <-
#   data.frame(metabolome_cytokine_cor, p, fdr, stringsAsFactors = FALSE) %>%
#   dplyr::filter(fdr < 0.05)
# 
# dim(metabolome_cytokine_cor)
# 
# save(metabolome_cytokine_cor, file = "metabolome_cytokine_cor")
# load("metabolome_cytokine_cor")
# 
# 
# 
# 
# 
# 
# ##met vs cytokine
# data1 <- lipidomics_expression_data
# data2 <- cytokine_expression_data
# 
# intersect_name <- intersect(colnames(data1), colnames(data2))
# 
# data1 <- as.data.frame(data1)[,intersect_name]
# data2 <- as.data.frame(data2)[,intersect_name]
# 
# lipidomics_cytokine_cor <-
#   cor(x = t(data1),
#       y = t(data2),
#       method = "spearman", use = "pairwise.complete.obs")
# 
# lipidomics_cytokine_cor <-
#   lipidomics_cytokine_cor %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "x") %>%
#   tidyr::pivot_longer(cols = -x,
#                       names_to = "y",
#                       values_to = "r")
# 
# p <-
#   as.data.frame(t(lipidomics_cytokine_cor)) %>%
#   purrr::map(.f = function(x){
#     cor.test(as.numeric(data1[x[1],]),
#              as.numeric(data2[x[2],]),
#              method = "spearman"
#     )$p.value
#   }) %>%
#   unlist()
# 
# fdr <- p.adjust(p, method = "fdr")
# 
# fdr[fdr == 0] <- min(fdr[fdr != 0])
# 
# lipidomics_cytokine_cor <-
#   data.frame(lipidomics_cytokine_cor, p, fdr, stringsAsFactors = FALSE) %>%
#   dplyr::filter(fdr < 0.05)
# 
# dim(lipidomics_cytokine_cor)
# 
# save(lipidomics_cytokine_cor, file = "lipidomics_cytokine_cor")
# load("lipidomics_cytokine_cor")

load("metabolome_cytokine_cor")
load("metabolome_lipidomics_cor")
load("lipidomics_cytokine_cor")

####example of network
##trans
metabolome_lipidomics_cor <-
  metabolome_lipidomics_cor %>%
  mutate(var1_class = rep("metabolome", nrow(metabolome_lipidomics_cor)),
         var2_class = rep("lipidomics", nrow(metabolome_lipidomics_cor)))

metabolome_cytokine_cor <-
  metabolome_cytokine_cor %>%
  mutate(var1_class = rep("metabolome", nrow(metabolome_cytokine_cor)),
         var2_class = rep("cytokine", nrow(metabolome_cytokine_cor)))

lipidomics_cytokine_cor <-
  lipidomics_cytokine_cor %>%
  mutate(var1_class = rep("lipidomics", nrow(lipidomics_cytokine_cor)),
         var2_class = rep("cytokine", nrow(lipidomics_cytokine_cor)))

top50 <-
  rbind(metabolome_lipidomics_cor,
        metabolome_cytokine_cor,
        lipidomics_cytokine_cor) %>%
  distinct() %>% 
  dplyr::arrange(desc(abs(r))) %>% 
  head(200)

m <- apply(top50[, c("x", "y")], 1, sort) %>%
  t()

remove_idx <- which(duplicated(m))
remove_idx
if(length(remove_idx) > 0){
  top50 <- 
    top50[-remove_idx,] 
}

# save(top50, file = "top50")
load("top50")

library(igraph)
library(ggraph)
library(tidygraph)

value <- 
  c("Lipid" = ggsci::pal_aaas()(10)[1],
    "Metabolite" = ggsci::pal_aaas()(10)[3],
    "Cytokine" = ggsci::pal_aaas()(10)[4]
  )

correlation_data <- top50

edge_data <-
  correlation_data %>%
  dplyr::mutate(
    from = paste(var1_class, x, sep = "_"),
    to = paste(var2_class, y, sep = "_"),
    fdr = -log(fdr, 10),
    cor = r
  ) %>%
  dplyr::select(from, to, cor, fdr)
  
node_data <-
  data.frame(node = unique(c(edge_data$from, edge_data$to)),
             stringsAsFactors = FALSE) %>%
  dplyr::mutate(class = stringr::str_extract(node,
                                             "metabolome_|lipidomics_|cytokine_")) %>%
  dplyr::mutate(
    class = case_when(
      class == "lipidomics_" ~ "Lipid",
      class == "metabolome_" ~ "Metabolite",
      class == "cytokine_" ~ "Cytokine"
    )
  )
  
variable_info <-
  rbind(
    metabolome_variable_info %>%
      dplyr::select(variable_id, true_name = Metabolite) %>%
      dplyr::mutate(variable_id = paste("metabolome", variable_id, sep = "_")),
    lipidomics_variable_info %>%
      dplyr::select(variable_id, true_name = mol_name) %>%
      dplyr::mutate(variable_id = paste("lipidomics", variable_id, sep = "_")),
    cytokine_variable_info %>%
      dplyr::select(variable_id, true_name = mol_name) %>%
      dplyr::mutate(variable_id = paste("cytokine", variable_id, sep = "_"))
  )
  
node_data <-
  node_data %>%
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))
    
node_data$class <-
  factor(node_data$class,
         levels = c("Metabolite", "Lipid", "Cytokine"))
  
temp_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))
 
plot <-
  ggraph(temp_data,
         layout = 'linear',
         circular = TRUE) +
  geom_edge_arc0(
    strength = 1,
    aes(edge_width = fdr,
        color = cor),
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(color = class, size = Degree),
    shape = 16,
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = true_name,
      hjust = 'outward',
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      size = 3,
      colour = class
    ),
    size = 2,
    alpha = 1,
    show.legend = FALSE
  ) +
  guides(
    linetype = "none",
    color = guide_legend(
      title = "Class",
      override.aes = list(size = 3, linetype = "none")
    ),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = NA,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    edge_width = guide_legend(title = "FDR adjusted P value", override.aes = list(shape = NA)),
    edge_color = guide_edge_colorbar(title = "Spearman correlation")
  ) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7),
                                                 "white",
                                                 alpha("#EE0000FF", 0.7))) +
  ggraph::scale_edge_width(range = c(0.1, 1.5)) +
  scale_size_continuous(range = c(1, 5)) +
  scale_color_manual(values = value) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_cartesian(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4))

plot

extrafont::loadfonts()

# ggsave(
#   plot,
#   filename = "top50_correlation.pdf",
#   width = 8.5,
#   height = 7,
#   bg = "transparent"
# )

######total correlation
total_cor <-
  rbind(
    metabolome_cytokine_cor,
    metabolome_lipidomics_cor,
    lipidomics_cytokine_cor
  )  

dim(total_cor)
# save(total_cor, file = "total_cor")

load("total_cor")
unique(c(total_cor$x,total_cor$y)) %>% length()

dim(total_cor)
length(unique(c(total_cor$x, total_cor$y)))
quantile(abs(total_cor$r))

library(igraph)
library(ggraph)
library(tidygraph)

# all_omics_graph <-
#   construct_graph(correlation_data = total_cor,
#                   metabolite_variable_info = variable_info)
# # 

correlation_data <- total_cor

edge_data <-
  correlation_data %>%
  dplyr::mutate(
    from = paste(var1_class, x, sep = "_"),
    to = paste(var2_class, y, sep = "_"),
    fdr = -log(fdr, 10),
    cor = r
  ) %>%
  dplyr::select(from, to, cor, fdr)

node_data <-
  data.frame(node = unique(c(edge_data$from, edge_data$to)),
             stringsAsFactors = FALSE) %>%
  dplyr::mutate(class = stringr::str_extract(node,
                                             "metabolome_|lipidomics_|cytokine_")) %>%
  dplyr::mutate(
    class = case_when(
      class == "lipidomics_" ~ "Lipid",
      class == "metabolome_" ~ "Metabolite",
      class == "cytokine_" ~ "Cytokine"
    )
  )

variable_info <-
  rbind(
    metabolome_variable_info %>% 
      dplyr::select(variable_id, true_name = Metabolite) %>% 
      dplyr::mutate(variable_id = paste("metabolome", variable_id, sep = "_")),
    lipidomics_variable_info %>% 
      dplyr::select(variable_id, true_name = mol_name) %>% 
      dplyr::mutate(variable_id = paste("lipidomics", variable_id, sep = "_")),
    cytokine_variable_info %>% 
      dplyr::select(variable_id, true_name = mol_name) %>% 
      dplyr::mutate(variable_id = paste("cytokine", variable_id, sep = "_"))
  )

node_data <- 
  node_data %>% 
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))


node_data <-
  node_data %>%
  dplyr::arrange(class)

node_data$class <-
  factor(node_data$class, 
         levels = c("Metabolite", "Lipid", "Cytokine"))

all_omics_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

# save(all_omics_graph, file = "all_omics_graph")
load("all_omics_graph")

####statistics of network
##edge information
node_data <- 
  igraph::vertex_attr(graph = all_omics_graph) %>% 
  purrr::map(.f = function(x){
    if(class(x) == "factor"){
      x <- as.character(x)
    }
    x
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

edge_data <- total_cor %>%
  dplyr::rename(from = x, to = y)

dim(node_data)
dim(edge_data)
table(node_data$class)

plot <-
  ggraph(all_omics_graph,
         layout = 'fr') +
  geom_edge_link(
    strength = 1,
    aes(edge_width = fdr,
        color = cor),
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  guides(
    color = guide_legend(
      title = "Class",
      override.aes = list(size = 3, linetype = NA)
    ),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = 0,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    edge_width = guide_legend(title = "FDR adjusted P value")
  ) +
  ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7),
                                                 "white",
                                                 alpha("#EE0000FF", 0.7))) +
  ggraph::scale_edge_width(range = c(0.05, 0.5)) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_manual(values = value) +
  scale_fill_manual(values = value) +
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
#   filename = "all_correlation.pdf",
#   width = 15,
#   height = 10,
#   bg = "transparent"
# )

####subnetwork
# all_subnetworks <-
#   igraph::cluster_edge_betweenness(graph = all_omics_graph,
#                                    weights = abs(edge_attr(all_omics_graph,
#                                                            "cor")))
# save(all_subnetworks, file = "all_subnetworks")

load("all_subnetworks")

table(membership(all_subnetworks))

plot <- 
  ggplot(
    data.frame(index = 1:length(all_subnetworks$modularity),
               modu = all_subnetworks$modularity, stringsAsFactors = FALSE),
    aes(index, modu) 
  ) +
  geom_vline(xintercept = which.max(all_subnetworks$modularity), 
             linetype = 2, colour = "#800000B2") + 
  labs(x = "Community analysis iteration", y = "Modularity") +
  geom_line(colour = "black") +
  # geom_point() +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank())

plot <-
  plot + 
  ggplot2::annotate(geom = "point", 
                    x = which.max(all_subnetworks$modularity),
                    y = max(all_subnetworks$modularity), 
                    size = 3, 
                    colour = "#FFA319FF") +
  annotate(geom = "text", 
           x = which.max(all_subnetworks$modularity),
           y = max(all_subnetworks$modularity), 
           label = paste("(",  which.max(all_subnetworks$modularity),
                         ",", 
                         max(all_subnetworks$modularity) %>% round(3),
                         ")"),
           size = 5,
           colour = "#FFA319FF"
  )

plot

# ggsave(plot, filename = "all_modularity.pdf", width = 7, height = 7)

table(membership(communities = all_subnetworks))

which(table(membership(communities = all_subnetworks)) >= 5)

which(table(membership(communities = all_subnetworks)) >= 5) %>% 
  length()

idx <- 1
table(membership(communities = all_subnetworks))[idx]

# for(idx in which(table(membership(communities = all_subnetworks)) >= 5)){
#   cat(idx, " ")
# 
#   subnetwork <-
#     igraph::induced_subgraph(graph = all_omics_graph,
#                              v = which(membership(all_subnetworks) == idx)) %>%
#     tidygraph::as_tbl_graph() %>%
#     dplyr::mutate(Degree = centrality_degree(mode = 'all'))
# 
#   dir.create("subnetworks", showWarnings = FALSE)
#   path <- paste("subnetworks/subnetwork", idx, sep = "_")
#   dir.create(path)
# 
#   save(subnetwork, file = file.path(path, "subnetwork"))
# 
#   degree1 <-
#     igraph::degree(subnetwork, mode = "all", normalized = FALSE)
# 
#   degree2 <-
#     igraph::degree(subnetwork, mode = "all", normalized = TRUE)
# 
#   betweenness1 <- igraph::betweenness(
#     graph = subnetwork,
#     weights = abs(igraph::edge_attr(graph = subnetwork, name = "cor")),
#     normalized = FALSE
#   )
# 
#   betweenness2 <- igraph::betweenness(
#     graph = subnetwork,
#     weights = abs(igraph::edge_attr(graph = subnetwork, name = "cor")),
#     normalized = TRUE
#   )
# 
#   closeness1 <-
#     igraph::closeness(
#       graph = subnetwork,
#       mode = "all",
#       weights = abs(igraph::edge_attr(graph = subnetwork, name = "cor")),
#       normalized = FALSE
#     )
# 
#   closeness2 <-
#     igraph::closeness(
#       graph = subnetwork,
#       mode = "all",
#       weights = abs(igraph::edge_attr(graph = subnetwork, name = "cor")),
#       normalized = TRUE
#     )
# 
#   importance <-
#     data.frame((degree1) / sd(degree1),
#                (closeness1) / sd(closeness1),
#                (betweenness1) / sd(betweenness1),
#                stringsAsFactors = FALSE
#     ) %>%
#     apply(1, mean)
# 
#   node_name <- igraph::vertex.attributes(subnetwork)$node
#   class <- igraph::vertex.attributes(subnetwork)$class
#   true_name <- igraph::vertex.attributes(subnetwork)$true_name
# 
#   node_info <-
#     data.frame(
#       node_name,
#       true_name,
#       class,
#       importance,
#       degree1,
#       degree2,
#       betweenness1,
#       betweenness2,
#       closeness1,
#       closeness2,
#       stringsAsFactors = FALSE
#     ) %>%
#     dplyr::arrange(desc(importance))
# 
#   write.csv(node_info,
#             file = file.path(path, paste("subnetwork", idx, "_node_info.csv", sep = "")),
#             row.names = FALSE)
# 
#   if(nrow(node_info) <= 5){
#     hub_gene1 <- node_info %>%
#       dplyr::pull(node_name)
# 
#     hub_gene2 <- node_info %>%
#       dplyr::pull(node_name)
# 
#     hub_gene3 <- node_info %>%
#       dplyr::pull(node_name)
# 
#     hub_gene4 <- node_info %>%
#       dplyr::pull(node_name)
#   }else{
#     hub_gene1 <- node_info %>%
#       dplyr::filter(importance > quantile(importance, 0.75)) %>%
#       dplyr::pull(node_name)
# 
#     hub_gene2 <- node_info %>%
#       dplyr::filter(degree1 > quantile(degree1, 0.75)) %>%
#       dplyr::pull(node_name)
# 
#     hub_gene3 <- node_info %>%
#       dplyr::filter(betweenness1 > quantile(betweenness1, 0.75)) %>%
#       dplyr::pull(node_name)
# 
#     hub_gene4 <- node_info %>%
#       dplyr::filter(closeness1 > quantile(closeness1, 0.75)) %>%
#       dplyr::pull(node_name)
#   }
# 
# 
#   hub_gene <-
#     c(hub_gene1,
#       hub_gene2,
#       hub_gene3,
#       hub_gene4) %>% unique()
# 
#   temp_data <-
#     node_info %>%
#     dplyr::filter(node_name %in% hub_gene) %>%
#     dplyr::select(true_name, importance, contains("1")) %>%
#     dplyr::arrange(importance)
# 
#   plot <-
#     temp_data %>%
#     dplyr::select(true_name, importance, contains("1")) %>%
#     dplyr::arrange(desc(importance)) %>%
#     head(10) %>%
#     dplyr::mutate(degree1 = degree1 / sd(degree1)) %>%
#     dplyr::mutate(betweenness1 = betweenness1 / sd(betweenness1)) %>%
#     dplyr::mutate(closeness1 = closeness1 / sd(closeness1)) %>%
#     tidyr::pivot_longer(-true_name, names_to = "class", values_to = "value") %>%
#     dplyr::mutate(true_name = factor(true_name, levels = rev(node_info$true_name))) %>%
#     dplyr::mutate(
#       class =
#         case_when(
#           class == "importance" ~ "Importance",
#           class == "degree1" ~ "Degree",
#           class == "betweenness1" ~ "Betweenness",
#           class == "closeness1" ~ "Closeness"
#         )
#     ) %>%
#     dplyr::mutate(class = factor(
#       class,
#       levels = c("Importance", "Degree", "Betweenness", "Closeness")
#     )) %>%
#     ggplot(aes(value, true_name)) +
#     geom_point(
#       aes(fill = class),
#       shape = 21,
#       show.legend = FALSE,
#       size = 3
#     ) +
#     geom_segment(aes(
#       x = 0,
#       y = true_name,
#       xend = value,
#       yend = true_name,
#       color = class
#     ),
#     show.legend = FALSE) +
#     scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
#     ggsci::scale_fill_futurama() +
#     labs(x = "", y = "") +
#     facet_wrap(facets = "class",
#                nrow = 1,
#                scales = "free_x") +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank())
# 
#   ggsave(plot, file = file.path(path, paste("subnetwork",idx,"_hub_genes.pdf", sep = "")),
#          width = 8, height = 7)
# 
#   plot <-
#     ggraph(subnetwork,
#            layout = 'linear',
#            circular = TRUE) +
#     geom_edge_arc(
#       strength = 1,
#       aes(edge_width = fdr,
#           color = cor),
#       alpha = 1,
#       show.legend = TRUE
#     ) +
#     geom_node_point(
#       aes(fill = class, size = Degree),
#       alpha = 1,
#       show.legend = TRUE,
#       shape = 21
#     ) +
#     geom_node_text(
#       aes(
#         x = x * 1.05,
#         y = y * 1.05,
#         label = true_name,
#         hjust = 'outward',
#         angle = -((-node_angle(x, y) + 90) %% 180) + 90,
#         size = 3,
#         colour = class
#       ),
#       size = 2,
#       alpha = 1,
#       show.legend = FALSE
#     ) +
#     guides(
#       linetype = "none",
#       color = guide_legend(
#         title = "Class",
#         override.aes = list(size = 3, linetype = "none")
#       ),
#       size = guide_legend(
#         title = "Degree",
#         override.aes = list(
#           linetype = NA,
#           fill = "transparent",
#           shape = 21,
#           color = "black"
#         )
#       ),
#       edge_width = guide_legend(title = "FDR adjusted P value", override.aes = list(shape = NA)),
#       edge_color = guide_edge_colorbar(title = "Spearman correlation")
#     ) +
#     ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 1),
#                                                    "white",
#                                                    alpha("#EE0000FF", 1))) +
#     ggraph::scale_edge_width(range = c(0.1, 0.7)) +
#     scale_size_continuous(range = c(1, 5)) +
#     scale_fill_manual(values = value) +
#     scale_color_manual(values = value) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "right",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     coord_cartesian(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4))
# 
#   ggsave(
#     plot,
#     filename = file.path(path, paste("subnetwork", idx, "_circle.pdf", sep = "")),
#     width = 8.5,
#     height = 7,
#     bg = "transparent"
#   )
# 
#   plot <-
#     ggraph(subnetwork,
#            layout = 'kk') +
#     geom_edge_arc(aes(edge_width = fdr,
#                        color = cor),
#                    alpha = 1,
#                    show.legend = TRUE) +
#     geom_node_point(
#       aes(fill = class, size = Degree),
#       alpha = 1,
#       show.legend = TRUE,
#       shape = 21
#     ) +
#     geom_node_text(aes(label = true_name), repel = TRUE) +
#     guides(
#       linetype = "none",
#       fill = guide_legend(
#         title = "Class",
#         override.aes = list(size = 3, linetype = "none")
#       ),
#       size = guide_legend(
#         title = "Degree",
#         override.aes = list(
#           linetype = NA,
#           fill = "transparent",
#           shape = 21,
#           color = "black"
#         )
#       ),
#       edge_width = guide_legend(title = "FDR adjusted P value", override.aes = list(shape = NA)),
#       edge_color = guide_edge_colorbar(title = "Spearman correlation")
#     ) +
#     ggraph::scale_edge_color_gradientn(colours = c(alpha("#3B4992FF", 0.7),
#                                                    "white",
#                                                    alpha("#EE0000FF", 0.7))) +
#     ggraph::scale_edge_width(range = c(0.1, 1)) +
#     scale_size_continuous(range = c(3, 10)) +
#     scale_fill_manual(values = value) +
#     scale_color_manual(values = value) +
#     ggraph::theme_graph() +
#     theme(
#       plot.background = element_rect(fill = "transparent", color = NA),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       legend.position = "right",
#       legend.background = element_rect(fill = "transparent", color = NA)
#     )
# 
#   ggsave(
#     plot,
#     filename = file.path(path, paste("subnetwork",idx,"_kk.pdf", sep = "")),
#     width = 8.5,
#     height = 7,
#     bg = "transparent"
#   )
# }


# ###output all subnetwork information
# idx <- 2
# idx <- which(table(membership(communities = all_subnetworks)) >= 5)
# subnetwork_info <- vector(mode = "list", length = max(idx))
# 
# for (i in idx) {
#   cat(i, " ")
# 
#   subnetwork <-
#     igraph::induced_subgraph(graph = all_omics_graph,
#                              v = which(membership(all_subnetworks) == i))
# 
#   temp_node <-
#     igraph::vertex_attr(subnetwork)
# 
#   temp_node <- lapply(temp_node, as.character) %>%
#     do.call(cbind, .) %>%
#     as.data.frame()
# 
#   temp_edge <- igraph::as_data_frame(subnetwork)
# 
#   temp_edge$from <- temp_node$node[temp_edge$from]
#   temp_edge$to <- temp_node$node[temp_edge$to]
# 
#   path <-
#     file.path("subnetworks", paste("subnetwork", i, sep = "_"))
# 
#   dir.create(path)
# 
#   write.csv(temp_node,
#             file = file.path(path,
#                              paste("subnetwork", i, "node.csv", sep = "_")),
#             row.names = FALSE)
# 
#   write.csv(temp_edge,
#             file = file.path(path,
#                              paste("subnetwork", i, "edge.csv", sep = "_")),
#             row.names = FALSE)
# 
#   node_num <- dim(temp_node)[1]
#   edge_num <- dim(temp_edge)[1]
#   node_num2 <- table(temp_node$class)
#   num <- c(node_num, edge_num, node_num2)
#   names(num)[1:2] <- c("node", "edge")
#   num2 <- rep(0, 5)
#   names(num2) <- c(
#     "node",
#     "edge",
#     "Metabolite",
#     'Lipid',
#     "Cytokine"
#   )
# 
#   num2[match(names(num), names(num2))] <- num
# 
#   subnetwork_info[[i]] <- c(subnetwork = i, num2)
# }
# 
# subnetwork_info <-
#   subnetwork_info %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# write.csv(subnetwork_info,
#           "subnetworks/subnetwork_info.csv",
#           row.names = FALSE)

subnetwork_info <- readr::read_csv("subnetworks/subnetwork_info.csv")

####plot
library(scatterpie)
subnetwork_info1 <- 
  subnetwork_info %>%
  # dplyr::mutate(node = log(node, 10),
  #               edge = log(edge, 10)) %>%
  dplyr::mutate(group = factor(1:nrow(subnetwork_info))) %>% 
  dplyr::select(subnetwork:edge,group, everything())

subnetwork_info2 <- subnetwork_info1

subnetwork_info2$node[3] <- subnetwork_info2$node[3] + 1
subnetwork_info2$edge[3] <- subnetwork_info2$edge[3] + 1

subnetwork_info2$node[5] <- subnetwork_info2$node[5] - 1
subnetwork_info2$edge[5] <- subnetwork_info2$edge[5] - 1

plot <-
  subnetwork_info2 %>% 
  ggplot() +
  geom_scatterpie(aes(x = node, 
                      y = edge, 
                      group = group,
                      r = ifelse(node > 30, 7, 2),
                      ),
                  data = subnetwork_info2 %>% dplyr::select(-subnetwork),
                  cols = names(value)[-7],
                  color = "black",
                  alpha = 0.7) + 
  scale_fill_manual(values = value) +
  # coord_fixed(ratio = 1, xlim = c(0.6,1.5), ylim = c(0.5,1.5)) +
  geom_text(aes(node, edge,
                label = subnetwork),
            hjust = 2.5) +
  theme_bw() +
  labs(x = "log10(Node number)", y = "log10(Edge number)") +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = c(0,1),
    legend.justification = c(0,1),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
ggforce::facet_zoom(xlim = c(0,25), 
                    ylim = c(0,75), horizontal = TRUE, zoom.size = 0.5)

plot

# ggsave(
#   plot,
#   filename = "subnetworks/subnetwork_plot2.pdf",
#   width = 10,
#   height = 7,
#   bg = "transparent"
# )



#######information of 
####subnetwork 2
load("subnetworks/subnetwork_2/subnetwork")
subnetwork2 <- subnetwork
node_data =
  igraph::vertex_attr(subnetwork2) %>% 
  purrr::map(function(x){
    as.character(x)
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

edge_data =
  igraph::as_data_frame(subnetwork2)

edge_data$from_name =
  node_data$true_name[edge_data$from]

edge_data$to_name =
  node_data$true_name[edge_data$to]

edge_data$from = 
  node_data$node[edge_data$from]

edge_data$to = 
  node_data$node[edge_data$to]

node_data %>%
  dplyr::filter(true_name == "LEPTIN")

temp_data = 
node_data %>%
  dplyr::filter(class != "Lipid") %>%
  dplyr::pull(true_name) %>% 
purrr::map(function(x) {
  x <- 
  edge_data %>% 
    dplyr::filter(from_name %in% x | to_name %in% x)   %>% 
    dplyr::filter(stringr::str_detect(from_name, 'AG') | stringr::str_detect(to_name, 'AG'))
  c(sum(x$cor > 0), sum(x$cor < 0), nrow(x))
}) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

colnames(temp_data) = c("pos", "neg", "total")
rownames(temp_data) = 
  node_data %>%
  dplyr::filter(class != "Lipid") %>%
  dplyr::pull(true_name)


temp_data = 
  temp_data %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  dplyr::filter(!variable_id %in% c("1,2-Dihexadecanoyl-3-(9Z-octadecenoyl)-sn-glycerol",
                                    "BETAINE", "L-Lactic acid", "LysoPE(18:2)",
                                    "Triacylglycerol 18:1-18:1-18:2")) %>% 
  dplyr::left_join(node_data[,c("true_name", "class")], by = c("variable_id" = "true_name"))

angle = ifelse(temp_data$pos == 0, 0, 90)

plot =
  temp_data %>%
  ggplot(aes(x = pos, y = neg)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(size = total, fill = class),
             shape = 21,
             show.legend = FALSE) +
  scale_size_continuous(range = c(5, 10)) +
  geom_text(aes(pos, neg, label = variable_id),
            angle = angle) +
  labs(x = "Number of positive correlation with TGA/DAG",
       y = "Number of negative correlation with TGA/DAG") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

plot

# ggsave(plot, filename = "subnetworks/subnetwork_2/number_plot.pdf", width = 7, height = 7)

plot =
  temp_data %>%
  dplyr::arrange(pos, desc(neg)) %>%
  dplyr::mutate(neg = -neg) %>%
  dplyr::mutate(variable_id = factor(variable_id, levels = variable_id)) %>%
  dplyr::mutate(pos = case_when(pos == 0 ~ "NA",
                                TRUE ~ as.character(pos))) %>%
  dplyr::mutate(pos = as.numeric(pos)) %>% 
  dplyr::mutate(neg = case_when(neg == 0 ~ "NA",
                                TRUE ~ as.character(neg))) %>%
  dplyr::mutate(neg = as.numeric(neg)) %>% 
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_segment(aes(y = variable_id, x = 0, yend = variable_id, xend = pos),
               color = ggsci::pal_aaas()(n=10)[2]) +
  geom_point(aes(y = variable_id, x = pos,
                 size = total), shape = 21, 
             fill = ggsci::pal_aaas()(n=10)[2]) +
  geom_segment(aes(y = variable_id, x = 0, yend = variable_id, xend = neg),
               color = ggsci::pal_aaas()(n=10)[1]) +
  geom_point(aes(y = variable_id, 
                 x = neg,
                 size = total), shape = 21, 
             fill = ggsci::pal_aaas()(n=10)[1]) +
  scale_size_continuous(range = c(5, 10)) +
  labs(x = "Number of correlation with TGA/DAG",
       y = "") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

plot

# ggsave(plot, filename = "subnetworks/subnetwork_2/number_plot2.pdf", width = 7, height = 10)


###subnetworks
