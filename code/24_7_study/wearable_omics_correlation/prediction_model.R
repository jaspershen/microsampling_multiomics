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
library(ggsankey)
library(plyr)
library(here)
library(tidyverse)
rm(list = ls())

{
  source(here::here("code/tools.R"))
  source(here::here("code/modified_dtw.R"))
  source(here::here("code/lagged_correlation.R"))
}

###load data
####wearbale and omics data
###CGM
{
  load(here::here("data/24_7_study/cgm/data_preparation/sample_info"))
  load(here::here("data/24_7_study/cgm/data_preparation/variable_info"))
  load(here::here("data/24_7_study/cgm/data_preparation/expression_data"))  

cgm_expression_data = expression_data
cgm_sample_info = sample_info
cgm_variable_info = variable_info

###HR
load(here::here("data/24_7_study/hr/data_preparation/sample_info"))
load(here::here("data/24_7_study/hr/data_preparation/variable_info"))
load(here::here("data/24_7_study/hr/data_preparation/expression_data"))

hr_expression_data = expression_data
hr_sample_info = sample_info
hr_variable_info = variable_info

###Step
load(here::here("data/24_7_study/steps/data_preparation/sample_info"))
load(here::here("data/24_7_study/steps/data_preparation/variable_info"))
load(here::here("data/24_7_study/steps/data_preparation/expression_data"))

steps_expression_data = expression_data
steps_sample_info = sample_info
steps_variable_info = variable_info

###cortisol
load(here::here("data/24_7_study/cortisol/data_preparation/sample_info"))
load(here::here("data/24_7_study/cortisol/data_preparation/variable_info"))
load(here::here("data/24_7_study/cortisol/data_preparation/expression_data"))

cortisol_sample_info = sample_info
cortisol_variable_info = variable_info
cortisol_expression_data = expression_data

###cytokine
load(here::here("data/24_7_study/cytokine/data_preparation/sample_info"))
load(here::here("data/24_7_study/cytokine/data_preparation/variable_info"))
load(here::here("data/24_7_study/cytokine/data_preparation/expression_data"))

cytokine_sample_info = sample_info
cytokine_variable_info = variable_info
cytokine_expression_data = expression_data

###lipidomics
load(here::here("data/24_7_study/lipidomics/data_preparation/sample_info"))
load(here::here("data/24_7_study/lipidomics/data_preparation/variable_info"))
load(here::here("data/24_7_study/lipidomics/data_preparation/expression_data"))

lipidomics_sample_info = sample_info
lipidomics_variable_info = variable_info
lipidomics_expression_data = expression_data

###metabolic_panel
load(here::here("data/24_7_study/metabolic_panel/data_preparation/sample_info"))
load(here::here("data/24_7_study/metabolic_panel/data_preparation/variable_info"))
load(here::here("data/24_7_study/metabolic_panel/data_preparation/expression_data"))

metabolic_panel_sample_info = sample_info
metabolic_panel_variable_info = variable_info
metabolic_panel_expression_data = expression_data

###metabolomics
load(here::here("data/24_7_study/metabolomics/data_preparation/metabolites/sample_info"))
load(here::here("data/24_7_study/metabolomics/data_preparation/metabolites/variable_info"))
load(here::here("data/24_7_study/metabolomics/data_preparation/metabolites/expression_data"))

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
load(here::here("data/24_7_study/proteomics/data_preparation/sample_info"))
load(here::here("data/24_7_study/proteomics/data_preparation/variable_info"))
load(here::here("data/24_7_study/proteomics/data_preparation/expression_data"))

proteomics_sample_info = sample_info
proteomics_variable_info = variable_info
proteomics_expression_data = expression_data

###total_protein
load(here::here("data/24_7_study/total_protein/data_preparation/sample_info"))
load(here::here("data/24_7_study/total_protein/data_preparation/variable_info"))
load(here::here("data/24_7_study/total_protein/data_preparation/expression_data"))

total_protein_sample_info = sample_info
total_protein_variable_info = variable_info
total_protein_expression_data = expression_data

####load lagged correlation 
#####CGM
load(here::here("data/24_7_study/wearable_omics_correlation/cgm_omics_correlation/cgm_cortisol/important_cortisol"))
cgm_cortisol_cor = important_cortisol

load(here::here("data/24_7_study/wearable_omics_correlation/cgm_omics_correlation/cgm_cytokine/important_cytokine"))
cgm_cytokine_cor = important_cytokine

load(here::here("data/24_7_study/wearable_omics_correlation/cgm_omics_correlation/cgm_lipidomics/important_lipid"))
cgm_lipid_cor = important_lipid

load(here::here("data/24_7_study/wearable_omics_correlation/cgm_omics_correlation/cgm_metabolic_panel/important_metabolic_panel"))
cgm_metabolic_panel_cor = important_metabolic_panel

load(here::here("data/24_7_study/wearable_omics_correlation/cgm_omics_correlation/cgm_metabolomics/important_metabolite"))
cgm_metabolite_cor = important_metabolite
cgm_metabolite_cor = 
  cgm_metabolite_cor %>% 
  dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                                "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                                "nistDatabase0.0.2"))

load(here::here("data/24_7_study/wearable_omics_correlation/cgm_omics_correlation/cgm_proteomics/important_protein"))
cgm_protein_cor = important_protein

load(here::here("data/24_7_study/wearable_omics_correlation/cgm_omics_correlation/cgm_total_protein/important_total_protein"))
cgm_total_protein_cor = important_total_protein

#####HR
load(here::here("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_cortisol/important_cortisol"))
hr_cortisol_cor = important_cortisol

load(here::here("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_cytokine/important_cytokine"))
hr_cytokine_cor = important_cytokine

load(here::here("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_lipidomics/important_lipid"))
hr_lipid_cor = important_lipid

load(here::here("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_metabolic_panel/important_metabolic_panel"))
hr_metabolic_panel_cor = important_metabolic_panel

load(here::here("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_metabolomics/important_metabolite"))
hr_metabolite_cor = important_metabolite
hr_metabolite_cor = 
  hr_metabolite_cor %>% 
dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                              "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                              "nistDatabase0.0.2"))

load(here::here("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_proteomics/important_protein"))
hr_protein_cor = important_protein

load(here::here("data/24_7_study/wearable_omics_correlation/hr_omics_correlation/hr_total_protein/important_total_protein"))
hr_total_protein_cor = important_total_protein

#####Step
load(here::here("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_cortisol/important_cortisol"))
steps_cortisol_cor = important_cortisol

load(here::here("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_cytokine/important_cytokine"))
steps_cytokine_cor = important_cytokine

load(here::here("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_lipidomics/important_lipid"))
steps_lipid_cor = important_lipid

load(here::here("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_metabolic_panel/important_metabolic_panel"))
steps_metabolic_panel_cor = important_metabolic_panel

load(here::here("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_metabolomics/important_metabolite"))
steps_metabolite_cor = important_metabolite

steps_metabolite_cor = 
  steps_metabolite_cor %>% 
  dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2",
                                "msDatabase_hilic0.0.2", "msDatabase_rplc0.0.2",
                                "nistDatabase0.0.2"))

load(here::here("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_proteomics/important_protein"))
steps_protein_cor = important_protein

load(here::here("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_total_protein/important_total_protein"))
steps_total_protein_cor = important_total_protein
}

####set work directory
###output directory is 
##("data/24_7_study/wearable_omics_correlation/total_network")

# ###cortisol from metabolomics data and cortisol
# grep('Cortisol', metabolomics_variable_info$Compound.name)
# metabolomics_variable_info[428,]
# dim(cortisol_expression_data)
# dim(metabolomics_expression_data)
#
# intersect_name = intersect(
#   colnames(metabolomics_expression_data),
#   colnames(cortisol_expression_data)
# )
#
# x = as.numeric(cortisol_expression_data[1,intersect_name])
# y = as.numeric(metabolomics_expression_data[428, intersect_name])
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
#
# data.frame(x, y) %>%
#   # dplyr::filter(y < 0.51) %>%
#   ggplot(aes(x, y)) +
#   geom_point() +
#   base_theme +
#   labs(x = "From cortisol data", y = "From metabolomics data") +
#   geom_smooth(method = "lm")
#
#
#
#
#
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


##############whole network
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

#####node and edge distribution
# Libraries
library(ggsankey)
temp_data = 
  edge_data %>% 
  dplyr::select(Wearable = from, Internal_omics = class) %>% 
  ggsankey::make_long(Wearable, Internal_omics)

color = c(class_color, wearable_color)
names(color)[9:11]= c("CGM", "HR", "Step")

plot = 
  ggplot(
    temp_data,
    aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = factor(node),
      label = node
    )
  ) +
  geom_sankey(flow.alpha = 0.6,
              node.color = "black", width = 0.3) +
  geom_sankey_label(size = 4, 
                    color = "white", 
                    fill = "gray40") +
  scale_fill_manual(values = color) +
  theme_sankey(base_size = 10) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
        plot.background = element_blank(),
        panel.background = element_blank())

plot

# ggsave(
#   plot,
#   filename = here::here(
#     "data/24_7_study/wearable_omics_correlation/total_network",
#     "edge_distributation.pdf"
#   ),
#   width = 12,
#   height = 7
# )


####pie to show for each wearable
library(ggpie)

temp_data = 
  edge_data %>% 
  dplyr::filter(from == 'CGM') %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup()

labs = 
  paste0(temp_data$class,
         " (",
         temp_data$n,
         ",",
         round(temp_data$n * 100 / sum(temp_data$n), 2),
         "%)")

plot = 
  ggpubr::ggpie(
    temp_data,
    x = "n",
    label = labs,
    fill = "class",
    color = "white",
    palette = class_color,
    legend = "none"
  )
plot

# ggsave(
#   plot,
#   filename = here::here(
#     "data/24_7_study/wearable_omics_correlation/total_network",
#     "cgm_edge_distributation.pdf"
#   ),
#   width = 9,
#   height = 9
# )

temp_data = 
  edge_data %>% 
  dplyr::filter(from == 'HR') %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup()

labs = 
  paste0(temp_data$class,
         " (",
         temp_data$n,
         ",",
         round(temp_data$n * 100 / sum(temp_data$n), 2),
         "%)")

plot = 
  ggpubr::ggpie(
    temp_data,
    x = "n",
    label = labs,
    fill = "class",
    color = "white",
    palette = class_color,
    legend = "none"
  )
plot

# ggsave(
#   plot,
#   filename = here::here(
#     "data/24_7_study/wearable_omics_correlation/total_network",
#     "hr_edge_distributation.pdf"
#   ),
#   width = 9,
#   height = 9
# )

temp_data = 
  edge_data %>% 
  dplyr::filter(from == 'Step') %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup()

labs = 
  paste0(temp_data$class,
         " (",
         temp_data$n,
         ",",
         round(temp_data$n * 100 / sum(temp_data$n), 2),
         "%)")

plot = 
  ggpubr::ggpie(
    temp_data,
    x = "n",
    label = labs,
    fill = "class",
    color = "white",
    palette = class_color,
    legend = "none"
  )
plot
# ggsave(
#   plot,
#   filename = here::here(
#     "data/24_7_study/wearable_omics_correlation/total_network",
#     "step_edge_distributation.pdf"
#   ),
#   width = 9,
#   height = 9
# )


###node distribution
table(node_data$class)

#draw a parliament diagram
temp_data =
  node_data %>% 
  dplyr::filter(!class %in% c("cgm", "hr", "step")) %>% 
  dplyr::select(class) %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(class = factor(class, levels = class))

library(ggpol)

color = class_color
color = color[match(temp_data$class, names(color))]
names(color) = NULL

plot = 
temp_data %>%
  ggplot() +
  # geom_parliament(aes(seats = n, fill = class),
  #                 color = "white") +
  geom_arcbar(aes(shares = n, 
                  r0 = 5, r1 = 10, fill = class),
              color = "black" ,
              sep = 0.05, show.legend = TRUE) +
  scale_fill_manual(values = color,
                    labels = temp_data$class) +
  # scale_color_manual(values = color,
  #                   labels = temp_data$class) +
  coord_fixed() +
  theme_void() +
  labs(title  = "282 internal moleculars")+
  theme(title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.caption = element_text(vjust = -3,hjust = 0.9),
        legend.position = 'bottom',
        legend.direction = "horizontal",
        legend.spacing.y = unit(0.1,"cm"),
        legend.spacing.x = unit(0.1,"cm"),
        legend.key.size = unit(0.8, 'lines'),
        legend.text = element_text(margin = margin(r = 1, unit = 'cm')),
        legend.text.align = 0)+
  guides(fill=guide_legend(nrow = 1,
                           byrow=TRUE,
                           reverse = TRUE,
                           title=NULL, 
                           label.theme = element_text(size = 12),
                           override.aes = list(size = 10, 
                                               shape = 21)))
plot
# ggsave(
#   plot,
#   filename = here::here(
#     "data/24_7_study/wearable_omics_correlation/total_network",
#     "node_distributation2.pdf"
#   ),
#   width = 9,
#   height = 9
# )

plot = 
parlDiag(
  Parties = temp_data$class,
  shares = temp_data$n,
  cols = color,
  repr = "absolute"
)

plot

# ggsave(
#   plot,
#   filename = here::here(
#     "data/24_7_study/wearable_omics_correlation/total_network",
#     "node_distributation.pdf"
#   ),
#   width = 9,
#   height = 9
# )


#####upset plot to show the overlap between three wearables with internal omics
library(ComplexUpset)
edge_data
cgm_edge_data =
  edge_data %>% 
  dplyr::filter(from == 'CGM')

hr_edge_data =
  edge_data %>% 
  dplyr::filter(from == 'HR')

step_edge_data =
  edge_data %>% 
  dplyr::filter(from == 'Step')

temp_data =
  rbind(
    data.frame(from = "CGM", to = unique(edge_data$to)),
    data.frame(from = "HR", to = unique(edge_data$to)),
    data.frame(from = "Step", to = unique(edge_data$to))
  ) %>% 
  dplyr::mutate(edge_id = paste(from, to, sep = "_")) %>% 
  dplyr::select(edge_id)

temp_data =
edge_data %>% 
  dplyr::select(from, to, class) %>% 
  dplyr::mutate(value = TRUE) %>% 
  tidyr::pivot_wider(names_from = from, values_from = value)

temp_data$CGM[is.na(temp_data$CGM)] = FALSE
temp_data$HR[is.na(temp_data$HR)] = FALSE
temp_data$Step[is.na(temp_data$Step)] = FALSE

plot = 
upset(
  data = temp_data,
  intersect = c("CGM", "HR", "Step"),
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts = FALSE,
      width = 0.8,
      mapping = aes(fill = class),
      show.legend = FALSE,
      text=list(check_overlap=TRUE)
    ) +
      scale_fill_manual(values = class_color)
  ),
  set_sizes=(
    upset_set_size(
      geom = geom_bar(
        aes(fill = class),
        width = 0.8
      ),
      position = 'right'
    ) +
      scale_fill_manual(values = class_color) +
      geom_text(aes(label=..count..), 
                hjust = 0,
                stat='count')
  ),
  annotations = list(
    'Intersection percent'=(
      ggplot(mapping=aes(fill = class))
      + geom_bar(stat='count', position='fill', 
                 show.legend = FALSE,
                 width = 0.8)
      + scale_y_continuous(labels=scales::percent_format())
      + scale_fill_manual(values=class_color)
      + ylab('Intersection percentg')
    )
  ),
  width_ratio = 0.2,
  height_ratio = 0.5,
  guides='over', 
  stripes = upset_stripes(geom=geom_segment(size=11))
)
plot

# ggsave(
#   plot,
#   filename = here::here(
#     "data/24_7_study/wearable_omics_correlation/total_network",
#     "wearable_edge_overlap.pdf"
#   ),
#   width = 7,
#   height = 7
# )

####let see the overlap have same correlation or not.
step_edge_data[,c("to", "class1")] %>% 
  dplyr::left_join(hr_edge_data[,c("to", "class1", "class")], by = "to") %>% 
  dplyr::filter(!is.na(class1.y)) %>% 
  dplyr::mutate(direction = case_when(
    class1.x == class1.y ~ "Same",
    class1.x != class1.y ~ "Different"
  )) %>% 
    dplyr::group_by(class) %>% 
  dplyr::summarise(n = c(sum(direction == "Same"), sum(direction == "Different")),
                   direction = c("same", "different"))

step_edge_data[,c("to", "class1")] %>% 
  dplyr::left_join(cgm_edge_data[,c("to", "class1", "class")], by = "to") %>% 
  dplyr::filter(!is.na(class1.y)) %>% 
  dplyr::mutate(direction = case_when(
    class1.x == class1.y ~ "Same",
    class1.x != class1.y ~ "Different"
  )) %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = c(sum(direction == "Same"), sum(direction == "Different")),
                   direction = c("same", "different"))

hr_edge_data[,c("to", "class1")] %>% 
  dplyr::left_join(cgm_edge_data[,c("to", "class1", "class")], by = "to") %>% 
  dplyr::filter(!is.na(class1.y)) %>% 
  dplyr::mutate(direction = case_when(
    class1.x == class1.y ~ "Same",
    class1.x != class1.y ~ "Different"
  )) %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = c(sum(direction == "Same"), sum(direction == "Different")),
                   direction = c("same", "different"))


####looks like that step and heart rate have the positive correlation
###so we will get the correlation plot between them
{
  load(
    here::here(
      "data/24_7_study/wearable_omics_correlation/hr_cgm_correlation/important_hr_cgm"
    )
  )
  load(
    here::here(
      "data/24_7_study/wearable_omics_correlation/steps_cgm_correlation/important_step_cgm"
    )
  )
  load(
    here::here(
      "data/24_7_study/wearable_omics_correlation/steps_hr_correlation/important_step_hr"
    )
  )
  
  }

important_hr_cgm
important_step_cgm
important_step_hr

load(here::here("data/24_7_study/wearable_omics_correlation/steps_hr_correlation/lagged_correlation/lagged_result"))

x = lagged_result$hr_1$x
x = (x - mean(x))/sd(x)
y = lagged_result$hr_1$y
idx = lagged_result$hr_1$max_idx

temp_data = 
purrr::map(idx, function(i){
  mean(y[i])
}) %>% 
  unlist() %>% 
  data.frame(step = x, hr = .) %>% 
  dplyr::filter(!is.na(step) & !is.na(hr))

library(viridis)
my_palette <- rev(magma(8))[c(-1,-8)]

plot = 
temp_data %>%
  # dplyr::filter(step != min(step)) %>% 
  ggplot(aes(step, hr)) +
  # geom_point() +
  geom_hex(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  geom_smooth(method = "lm") +
  labs(x = "Scaled step", y = "Scaled heart rate") +
  annotate(
    geom = "text",
    x = -Inf,
    y = Inf,
    hjust = 0,
    vjust = 1,
    label = paste("Lagged correlation: ", round(lagged_result$hr_1$max_cor, 2),
                  "\n", "p.adj: ", round(important_step_hr$lagged_cor_p_adjust, 4)),
    size = 5
  ) +
  base_theme 

plot

ggsave(
  plot,
  filename = here::here(
    "data/24_7_study/wearable_omics_correlation/total_network",
    "step_vs_hr.pdf"
  ),
  width = 9,
  height = 7
)

