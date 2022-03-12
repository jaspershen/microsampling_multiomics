
#' ---
#' title: "Wearable omics correlation network"
#' author: 
#'   - name: "Xiaotao Shen" 
#'     url: https://www.shenxt.info/
#'     affiliation: Stanford School of Medicine
#' date: "`r Sys.Date()`"
#' site: distill::distill_website
#' output: 
#'   distill::distill_article:
#'     code_folding: false
#' ---

no_function()

library(plyr)
library(here)
library(tidyverse)
rm(list = ls())

{
  source(here::here("R/tools.R"))
  source(here::here("R/modified_dtw.R"))
  source(here::here("R/lagged_correlation.R"))
}

###load data
####wearbale and omics data
###CGM
{
  load(here::here("data/7_24_mike/cgm/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/cgm/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/cgm/data_preparation/expression_data"))  
  
  cgm_expression_data = expression_data
  cgm_sample_info = sample_info
  cgm_variable_info = variable_info
  
  ###HR
  load(here::here("data/7_24_mike/hr/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/hr/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/hr/data_preparation/expression_data"))
  
  hr_expression_data = expression_data
  hr_sample_info = sample_info
  hr_variable_info = variable_info
  
  ###Step
  load(here::here("data/7_24_mike/steps/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/steps/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/steps/data_preparation/expression_data"))
  
  steps_expression_data = expression_data
  steps_sample_info = sample_info
  steps_variable_info = variable_info
  
  ###cortisol
  load(here::here("data/7_24_mike/cortisol/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/cortisol/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/cortisol/data_preparation/expression_data"))
  
  cortisol_sample_info = sample_info
  cortisol_variable_info = variable_info
  cortisol_expression_data = expression_data
  
  ###cytokine
  load(here::here("data/7_24_mike/cytokine/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/cytokine/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/cytokine/data_preparation/expression_data"))
  
  cytokine_sample_info = sample_info
  cytokine_variable_info = variable_info
  cytokine_expression_data = expression_data
  
  ###lipidomics
  load(here::here("data/7_24_mike/lipidomics/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/lipidomics/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/lipidomics/data_preparation/expression_data"))
  
  lipidomics_sample_info = sample_info
  lipidomics_variable_info = variable_info
  lipidomics_expression_data = expression_data
  
  ###metabolic_panel
  load(here::here("data/7_24_mike/metabolic_panel/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/metabolic_panel/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/metabolic_panel/data_preparation/expression_data"))
  
  metabolic_panel_sample_info = sample_info
  metabolic_panel_variable_info = variable_info
  metabolic_panel_expression_data = expression_data
  
  ###metabolomics
  load(here::here("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info"))
  load(here::here("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info"))
  load(here::here("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data"))
  
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  metabolomics_expression_data = expression_data
  
  ###remove metablites which are from other databases
  metabolomics_variable_info = 
    metabolomics_variable_info %>% 
    dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                                  "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                                  "nistDatabase0.0.2"))
  
  metabolomics_expression_data = 
    metabolomics_expression_data[metabolomics_variable_info$variable_id,]
  
  ###proteomics
  load(here::here("data/7_24_mike/proteomics/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/proteomics/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/proteomics/data_preparation/expression_data"))
  
  proteomics_sample_info = sample_info
  proteomics_variable_info = variable_info
  proteomics_expression_data = expression_data
  
  ###total_protein
  load(here::here("data/7_24_mike/total_protein/data_preparation/sample_info"))
  load(here::here("data/7_24_mike/total_protein/data_preparation/variable_info"))
  load(here::here("data/7_24_mike/total_protein/data_preparation/expression_data"))
  
  total_protein_sample_info = sample_info
  total_protein_variable_info = variable_info
  total_protein_expression_data = expression_data
  
  ####load lagged correlation 
  #####CGM
  load(here::here("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cortisol/important_cortisol"))
  cgm_cortisol_cor = important_cortisol
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine/important_cytokine"))
  cgm_cytokine_cor = important_cytokine
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_lipidomics/important_lipid"))
  cgm_lipid_cor = important_lipid
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_metabolic_panel/important_metabolic_panel"))
  cgm_metabolic_panel_cor = important_metabolic_panel
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_metabolomics/important_metabolite"))
  cgm_metabolite_cor = important_metabolite
  cgm_metabolite_cor = 
    cgm_metabolite_cor %>% 
    dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                                  "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                                  "nistDatabase0.0.2"))
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_proteomics/important_protein"))
  cgm_protein_cor = important_protein
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/cgm_omics_correlation/cgm_total_protein/important_total_protein"))
  cgm_total_protein_cor = important_total_protein
  
  #####HR
  load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_cortisol/important_cortisol"))
  hr_cortisol_cor = important_cortisol
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_cytokine/important_cytokine"))
  hr_cytokine_cor = important_cytokine
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_lipidomics/important_lipid"))
  hr_lipid_cor = important_lipid
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_metabolic_panel/important_metabolic_panel"))
  hr_metabolic_panel_cor = important_metabolic_panel
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_metabolomics/important_metabolite"))
  hr_metabolite_cor = important_metabolite
  hr_metabolite_cor = 
    hr_metabolite_cor %>% 
    dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                                  "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                                  "nistDatabase0.0.2"))
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_proteomics/important_protein"))
  hr_protein_cor = important_protein
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_total_protein/important_total_protein"))
  hr_total_protein_cor = important_total_protein
  
  #####Step
  load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_cortisol/important_cortisol"))
  steps_cortisol_cor = important_cortisol
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_cytokine/important_cytokine"))
  steps_cytokine_cor = important_cytokine
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_lipidomics/important_lipid"))
  steps_lipid_cor = important_lipid
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_metabolic_panel/important_metabolic_panel"))
  steps_metabolic_panel_cor = important_metabolic_panel
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_metabolomics/important_metabolite"))
  steps_metabolite_cor = important_metabolite
  
  steps_metabolite_cor = 
    steps_metabolite_cor %>% 
    dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                                  "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                                  "nistDatabase0.0.2"))
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_proteomics/important_protein"))
  steps_protein_cor = important_protein
  
  load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_total_protein/important_total_protein"))
  steps_total_protein_cor = important_total_protein
}

####set work directory
###output directory is 
##("data/7_24_mike/wearable_omics_correlation/total_network")


# ###cortisol from metabolomics data and cortisol
# grep('Cortisol', metabolomics_variable_info$Compound.name)
# metabolomics_variable_info[387,]
# dim(cortisol_expression_data)
# dim(metabolomics_expression_data)
# 
# intersect_name = intersect(
#   colnames(metabolomics_expression_data),
#   colnames(cortisol_expression_data)
# )
# 
# x = as.numeric(cortisol_expression_data[1,intersect_name])
# y = as.numeric(metabolomics_expression_data[387, intersect_name])
# plot(density(x))
# plot(density(y))
# y = log(y+1, 2)
# x = scale(x) %>% as.numeric()
# y = scale(y) %>% as.numeric()
# 
# cortisol_sample_info[,c("sample_id", "accurate_time")] %>% 
#   dplyr::left_join(metabolomics_sample_info[,c("sample_id", "accurate_time")],
#                    by = "sample_id")
# range(x)
# range(y)
# data.frame(x, y) %>% 
#   # dplyr::filter(y < 0.51) %>% 
#   ggplot(aes(x, y)) +
#     geom_point()
# 
# cor.test(x, y)
# 
# plot(x, y) 

###the correlation between cortisol and cortisol from metabolomics data are not good

# cgm_cortisol_cor
# cgm_total_protein_cor
# 
# hr_cortisol_cor
# hr_total_protein_cor
# 
# steps_cortisol_cor
# steps_total_protein_cor
# 
# colnames(cgm_cytokine_cor)
# colnames(cgm_lipid_cor)
# colnames(cgm_metabolic_panel_cor)
# colnames(cgm_metabolite_cor)
# colnames(cgm_protein_cor)
# colnames(cgm_total_protein_cor)

intersect_name =
  Reduce(
    f = intersect,
    x = list(
      colnames(cgm_cytokine_cor),
      colnames(cgm_lipid_cor),
      colnames(cgm_metabolic_panel_cor),
      colnames(cgm_metabolite_cor),
      colnames(cgm_protein_cor),
      colnames(cgm_total_protein_cor)
    )
  )

class_color =
  c(
    "lipidomics" = ggsci::pal_aaas()(10)[1],
    "metabolomics" = ggsci::pal_aaas()(10)[3],
    "cytokine" = ggsci::pal_aaas()(10)[4],
    "total_protein" = ggsci::pal_aaas()(10)[5],
    "cortisol" = ggsci::pal_aaas()(10)[6],
    "metabolic_panel" = ggsci::pal_aaas()(10)[7],
    "proteomics" = ggsci::pal_aaas()(10)[8]
  )

temp_data = 
rbind(
  data.frame(cgm_cytokine_cor[, intersect_name], class = "cytokine"),
  data.frame(cgm_lipid_cor[, intersect_name], class = "lipidomics"),
  data.frame(cgm_metabolic_panel_cor[, intersect_name], class = "metabolic_panel"),
  data.frame(cgm_metabolite_cor[, intersect_name], class = "metabolomics"),
  data.frame(cgm_protein_cor[, intersect_name], class = "proteomics"),
  data.frame(cgm_total_protein_cor[, intersect_name], class = "total_protein"),
  data.frame(cgm_cortisol_cor[, intersect_name], class = "cortisol"),
  
  data.frame(hr_cytokine_cor[, intersect_name], class = "cytokine"),
  data.frame(hr_lipid_cor[, intersect_name], class = "lipidomics"),
  data.frame(hr_metabolic_panel_cor[, intersect_name], class = "metabolic_panel"),
  data.frame(hr_metabolite_cor[, intersect_name], class = "metabolomics"),
  data.frame(hr_protein_cor[, intersect_name], class = "proteomics"),
  data.frame(hr_total_protein_cor[, intersect_name], class = "total_protein"),
  data.frame(hr_cortisol_cor[, intersect_name], class = "cortisol"),
  
  data.frame(steps_cytokine_cor[, intersect_name], class = "cytokine"),
  data.frame(steps_lipid_cor[, intersect_name], class = "lipidomics"),
  data.frame(steps_metabolic_panel_cor[, intersect_name], class = "metabolic_panel"),
  data.frame(steps_metabolite_cor[, intersect_name], class = "metabolomics"),
  data.frame(steps_protein_cor[, intersect_name], class = "proteomics"),
  data.frame(steps_total_protein_cor[, intersect_name], class = "total_protein"),
  data.frame(steps_cortisol_cor[, intersect_name], class = "cortisol")
)

dim(temp_data)


##############total network

temp_data = 
  temp_data %>% 
  dplyr::filter(lagged_cor_p_adjust < 0.05)

dim(temp_data)

edge_data =
  temp_data %>%
  dplyr::rename(from = wearable,
                to = variable_id)
node_data = 
  edge_data %>% 
  dplyr::select(from, to, class) %>% 
  tidyr::pivot_longer(cols = -class, names_to = "node", values_to = "node2") %>% 
  dplyr::distinct(node2, .keep_all = TRUE) %>% 
  dplyr::mutate(class = 
                  case_when(node == "from" ~ node2,
                            node == "to" ~ class)) %>% 
  dplyr::select(-node) %>% 
  dplyr::rename(node = node2) %>% 
  dplyr::select(node, class) %>% 
  dplyr::mutate(class = 
                  case_when(
                    class == "CGM" ~ "cgm",
                    class == "HR" ~ "hr",
                    class == "Step" ~ "step",
                    TRUE ~ class
                  )) 

dim(node_data)
dim(edge_data)

edge_data %>% 
  dplyr::filter(from == "CGM") %>% 
  pull(class) %>% 
  table()

metabolomics_variable_info = 
  metabolomics_variable_info %>% 
  dplyr::mutate(mol_name = Compound.name)

variable_info =
  rbind(cortisol_variable_info[,c("variable_id", "mol_name")],
        cytokine_variable_info[,c("variable_id", "mol_name")],
        lipidomics_variable_info[,c("variable_id", "mol_name")],
        metabolic_panel_variable_info[,c("variable_id", "mol_name")],
        metabolomics_variable_info[,c("variable_id", "mol_name")],
        proteomics_variable_info[,c("variable_id", "mol_name")],
        total_protein_variable_info[,c("variable_id", "mol_name")]
  )

node_data = 
  node_data %>% 
  dplyr::left_join(variable_info, by = c("node" = "variable_id")) %>% 
  dplyr::mutate(mol_name = 
                  case_when(is.na(mol_name) ~ node,
                            !is.na(mol_name) ~ mol_name)) %>% 
  dplyr::mutate(class2 = 
                  case_when(
                    class %in% c("cgm", "hr", "step") ~ "wearable",
                    !class %in% c("cgm", "hr", "step") ~ "internal-omics"
                  ))

table(node_data$class)

library(igraph)
library(tidygraph)

node_data = 
  node_data %>% 
  dplyr::arrange(class) %>% 
  dplyr::mutate(node = factor(node, levels = node))

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


#####up-down
g <- total_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

coords$y[coords$class2 == "wearable"] = 0.3
coords$y[coords$class2 == "internal-omics"] = 1

table(coords$class)

coords$y[coords$class2 == "wearable"] = -0.5
coords$y[coords$class == "metabolomics"] = 1
coords$y[coords$class == "lipidomics"] = 0.5
coords$y[coords$class == "proteomics"] = 0
coords$y[coords$class == "cytokine"] = 0
coords$y[coords$class == "metabolic_panel"] = 0
coords$y[coords$class == "total_protein"] = 0

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )


plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(
      color = lagged_cor,
      width = -log(lagged_cor_p_adjust, 10)
    ),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),
      nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    color = "black",
    size = 2,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color, wearable_color)) +
  scale_color_manual(values = c(class_color, wearable_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
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
  scale_edge_width_continuous(range = c(0.3, 2)) +
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

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/total_network",
#     "total_network.pdf"
#   ),
#   width = 8.6,
#   height = 7
# )






##############example network
##only remained the correlation with adjusted p value < 0.05 and top 10 correlations
# library(plyr)
# 
# temp_data = 
#   temp_data %>% 
#   dplyr::filter(lagged_cor_p_adjust < 0.05) %>% 
#   dplyr::filter(abs(lagged_cor) > 0.3)

temp_data = 
temp_data %>% 
  dplyr::filter(lagged_cor_p_adjust < 0.05) %>% 
  plyr::dlply(.variables = .(wearable, class)) %>% 
  purrr::map(function(x){
    x %>% 
      dplyr::arrange(abs(lagged_cor)) %>% 
      head(50)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

dim(temp_data)

edge_data =
  temp_data %>%
  dplyr::rename(from = wearable,
                to = variable_id)

node_data = 
  edge_data %>% 
  dplyr::select(from, to, class) %>% 
  tidyr::pivot_longer(cols = -class, names_to = "node", values_to = "node2") %>% 
  dplyr::distinct(node2, .keep_all = TRUE) %>% 
  dplyr::mutate(class = 
                  case_when(node == "from" ~ node2,
                            node == "to" ~ class)) %>% 
  dplyr::select(-node) %>% 
  dplyr::rename(node = node2) %>% 
  dplyr::select(node, class) %>% 
  dplyr::mutate(class = 
                  case_when(
                    class == "CGM" ~ "cgm",
                    class == "HR" ~ "hr",
                    class == "Step" ~ "step",
                    TRUE ~ class
                  )) 

dim(node_data)
dim(edge_data)

edge_data %>% 
  dplyr::filter(from == "CGM") %>% 
  pull(class) %>% 
  table()

metabolomics_variable_info = 
metabolomics_variable_info %>% 
  dplyr::mutate(mol_name = Compound.name)

variable_info =
  rbind(cortisol_variable_info[,c("variable_id", "mol_name")],
        cytokine_variable_info[,c("variable_id", "mol_name")],
        lipidomics_variable_info[,c("variable_id", "mol_name")],
        metabolic_panel_variable_info[,c("variable_id", "mol_name")],
        metabolomics_variable_info[,c("variable_id", "mol_name")],
        proteomics_variable_info[,c("variable_id", "mol_name")],
        total_protein_variable_info[,c("variable_id", "mol_name")]
        )

node_data = 
node_data %>% 
  dplyr::left_join(variable_info, by = c("node" = "variable_id")) %>% 
  dplyr::mutate(mol_name = 
                  case_when(is.na(mol_name) ~ node,
                            !is.na(mol_name) ~ mol_name)) %>% 
  dplyr::mutate(class2 = 
                  case_when(
                    class %in% c("cgm", "hr", "step") ~ "wearable",
                    !class %in% c("cgm", "hr", "step") ~ "internal-omics"
                  ))

table(node_data$class)

library(igraph)
library(tidygraph)

node_data = 
node_data %>% 
  dplyr::arrange(class) %>% 
  dplyr::mutate(node = factor(node, levels = node))

total_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))


#####up-down
g <- total_graph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

coords$y[coords$class2 == "wearable"] = 0.3
coords$y[coords$class2 == "internal-omics"] = 1

table(coords$class)

coords$y[coords$class2 == "wearable"] = -0.5
coords$y[coords$class == "metabolomics"] = 1
coords$y[coords$class == "lipidomics"] = 0.5
coords$y[coords$class == "proteomics"] = 0
coords$y[coords$class == "cytokine"] = 0
coords$y[coords$class == "metabolic_panel"] = 0
coords$y[coords$class == "total_protein"] = 0

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )


plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(
      color = lagged_cor,
      width = -log(lagged_cor_p_adjust, 10)
    ),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),
      color = class,
      nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    bg.color = "white",
    size = 2,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color, wearable_color)) +
  scale_color_manual(values = c(class_color, wearable_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
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
  scale_edge_width_continuous(range = c(0.3, 2)) +
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

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/example_network",
#     "example_network.pdf"
#   ),
#   width = 9.1,
#   height = 7
# )


######subnetwork
cgm_edge_data =
  edge_data %>% 
  dplyr::filter(from == 'CGM')

hr_edge_data =
  edge_data %>% 
  dplyr::filter(from == 'HR')

step_edge_data =
  edge_data %>% 
  dplyr::filter(from == 'Step')

###CGM
cgm_node = unique(c(cgm_edge_data$from, cgm_edge_data$to))
cgm_sub_mygraph = 
total_graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::filter(node %in% cgm_node) %>% 
  tidygraph::arrange(class) %>% 
  dplyr::mutate(node = factor(node, levels = node))


#####up-down
g <- cgm_sub_mygraph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

coords = 
coords %>% 
  dplyr::left_join(edge_data %>% dplyr::filter(from == "CGM") %>% dplyr::select(to, lagged_cor),
                   by = c("node" = "to"))

coords$lagged_cor[is.na(coords$lagged_cor)] = 0

coords <-
  coords %>%
  dplyr::mutate(y1 = x) %>% 
  dplyr::mutate(x1 = case_when(
    y == 1 ~ 0,
    y == 0 ~ 1
  )) %>% 
  dplyr::select(-c(x, y)) %>% 
  dplyr::select(x = x1, y = y1, everything())
  
coords$x[coords$lagged_cor < 0] = -1

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot = 
ggraph(cgm_sub_mygraph, 
       layout = "manual", 
       x = coords$x,
       y = coords$y) +
  geom_edge_diagonal(aes(
    color = lagged_cor,
    width = -log(lagged_cor_p_adjust, 10)
  ),
  alpha = 0.5,
  show.legend = TRUE) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),
      color = class
      # nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      # angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    bg.color = "white",
    size = 3,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_size_continuous(range = c(5, 10)) +
  scale_fill_manual(values = c(class_color, wearable_color)) +
  scale_color_manual(values = c(class_color, wearable_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
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
      override.aes = list(
        shape = 21,
        size = 3,
        alpha = 1
      )
    )
  ) +
  scale_edge_width_continuous(range = c(1, 4)) +
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

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/example_network",
#     "cgm_subnetwork.pdf"
#   ),
#   width = 7,
#   height = 7
# )


###HR Step
hr_node = unique(c(hr_edge_data$from, hr_edge_data$to))
step_node = unique(c(step_edge_data$from, step_edge_data$to))
hr_step_sub_mygraph = 
  total_graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::filter(node %in% c(hr_node, step_node)) %>% 
  tidygraph::arrange(class) %>% 
  dplyr::mutate(node = factor(node, levels = node))

#####up-down
g <- hr_step_sub_mygraph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

coords$y[coords$class2 == "wearable"] = 0.3
coords$y[coords$class2 == "internal-omics"] = 1

table(coords$class)

coords$y[coords$class2 == "wearable"] = -0.5
coords$y[coords$class == "metabolomics"] = 1
coords$y[coords$class == "lipidomics"] = 0.5
coords$y[coords$class == "proteomics"] = 0
coords$y[coords$class == "cytokine"] = 0
coords$y[coords$class == "metabolic_panel"] = 0
coords$y[coords$class == "total_protein"] = 0

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(
      color = lagged_cor,
      width = -log(lagged_cor_p_adjust, 10)
    ),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),
      color = class,
      nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    bg.color = "white",
    size = 2,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color, wearable_color)) +
  scale_color_manual(values = c(class_color, wearable_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
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
  scale_edge_width_continuous(range = c(0.3, 2)) +
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

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/example_network",
#     "hr_step_correlation_network.pdf"
#   ),
#   width = 9.8,
#   height = 7
# )



###HR
hr_node = unique(c(hr_edge_data$from, hr_edge_data$to))
hr_sub_mygraph = 
  total_graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::filter(node %in% c(hr_node)) %>% 
  tidygraph::arrange(class) %>% 
  dplyr::mutate(node = factor(node, levels = node))

#####up-down
g <- hr_sub_mygraph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

coords = 
  coords %>% 
  dplyr::left_join(edge_data %>% dplyr::filter(from == "HR") %>% dplyr::select(to, lagged_cor),
                   by = c("node" = "to"))

table(coords$class)

coords$y[coords$class2 == "wearable"] = -0.8
coords$y[coords$class == "metabolomics"] = 1
coords$y[coords$class == "lipidomics"] = 0.5
coords$y[coords$class == "proteomics"] = 0.5
coords$y[coords$class == "cytokine"] = 0.75
coords$y[coords$class == "metabolic_panel"] = 0.75
coords$y[coords$class == "total_protein"] = 0.75
coords$y[coords$class == "cortisol"] = 0.75

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(
      color = lagged_cor,
      width = -log(lagged_cor_p_adjust, 10)
    ),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),
      color = class,
      nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    bg.color = "white",
    size = 2,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color, wearable_color)) +
  scale_color_manual(values = c(class_color, wearable_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
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
  scale_edge_width_continuous(range = c(0.3, 2)) +
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

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/total_network",
#     "hr_correlation_network.pdf"
#   ),
#   width = 9.8,
#   height = 7
# )



###step
step_node = unique(c(step_edge_data$from, step_edge_data$to))
step_sub_mygraph = 
  total_graph %>% 
  tidygraph::activate(what = "nodes") %>% 
  tidygraph::filter(node %in% c(step_node)) %>% 
  tidygraph::arrange(class) %>% 
  dplyr::mutate(node = factor(node, levels = node))

#####up-down
g <- step_sub_mygraph
library(igraph)
library(ggraph)

V(g)$type <- bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, class2, x, y)

coords = 
  coords %>% 
  dplyr::left_join(edge_data %>% dplyr::filter(from == "Step") %>% dplyr::select(to, lagged_cor),
                   by = c("node" = "to"))

coords$lagged_cor[is.na(coords$lagged_cor)] = 0

table(coords$class)

coords$y[coords$class2 == "wearable"] = -0.8
coords$y[coords$class == "metabolomics"] = 1
coords$y[coords$class == "lipidomics"] = 0.5
coords$y[coords$class == "proteomics"] = 0.5
coords$y[coords$class == "cytokine"] = 0.75
coords$y[coords$class == "metabolic_panel"] = 0.75
coords$y[coords$class == "total_protein"] = 0.75
coords$y[coords$class == "cortisol"] = 0.75

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
    # node.position = coords
  )

plot <-
  ggraph(my_graph,
         layout = 'bipartite') +
  geom_edge_diagonal(
    aes(
      color = lagged_cor,
      width = -log(lagged_cor_p_adjust, 10)
    ),
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = class,
        size = Degree),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", mol_name, NA),
      color = class
    ),
    bg.color = "white",
    size = 5,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class2 == "wearable", NA, mol_name),
      color = class,
      nudge_y = ifelse(class2 == "wearable", "inward", 'outward'),
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    bg.color = "white",
    size = 2,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color, wearable_color)) +
  scale_color_manual(values = c(class_color, wearable_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
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
  scale_edge_width_continuous(range = c(0.3, 2)) +
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

# extrafont::loadfonts()

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/example_network",
#     "step_correlation_network.pdf"
#   ),
#   width = 9.5,
#   height = 7
# )



##########overlap between HR and Step vs protein
load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_proteomics/positive_core_term"))
hr_positive_core_term = positive_core_term

load(here::here("data/7_24_mike/wearable_omics_correlation/hr_omics_correlation/hr_proteomics/negative_core_term"))
hr_negative_core_term = negative_core_term

load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_proteomics/positive_core_term"))
step_positive_core_term = positive_core_term

load(here::here("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_proteomics/negative_core_term"))
step_negative_core_term = negative_core_term


intersect(hr_positive_core_term$Description,hr_negative_core_term$Description)
intersect(step_positive_core_term$Description,step_negative_core_term$Description)

intersect(hr_positive_core_term$Description,step_positive_core_term$Description)
intersect(hr_negative_core_term$Description,step_negative_core_term$Description)

intersect(hr_positive_core_term$Description,step_negative_core_term$Description)
intersect(hr_negative_core_term$Description,step_positive_core_term$Description)

library(ComplexUpset)

temp_data = 
rbind(
  data.frame(node = hr_positive_core_term$node, 
             class = "HR positive",
             value = TRUE),
  data.frame(node = hr_negative_core_term$node, 
             class = "HR negative",
             value = TRUE),
  data.frame(node = step_positive_core_term$node, 
             class = "Step positive",
             value = TRUE),
  data.frame(node = step_negative_core_term$node, 
             class = "Step negative",
             value = TRUE)  
) %>% 
  tidyr::pivot_wider(names_from = class, values_from = "value")

plot = 
  upset(
    data = temp_data,
    intersect = c("HR positive", "HR negative", "Step positive", "Step negative")
  )
plot

# ggsave(
#   plot,
#   filename = here::here(
#     "data/7_24_mike/wearable_omics_correlation/total_network",
#     "step_hr_protein_go_overlap.pdf"
#   ),
#   width = 7,
#   height = 7
# )

