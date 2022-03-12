no_function()
library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

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
        apply(expression_data, 1, function(x){
          sum(x == 0)/ncol(expression_data)
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
    load("data/7_24_mike/lipidomics/data_preparation/sample_info")
    load("data/7_24_mike/lipidomics/data_preparation/variable_info")
    load("data/7_24_mike/lipidomics/data_preparation/expression_data")
    lipidomics_sample_info = sample_info
    lipidomics_variable_info = variable_info
    lipidomics_expression_data = expression_data
    
    
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
    load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
    load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
    load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
    metabolomics_sample_info = sample_info
    metabolomics_variable_info = variable_info
    metabolomics_expression_data = expression_data
    
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
    load("data/7_24_mike/proteomics/data_preparation/sample_info")
    load("data/7_24_mike/proteomics/data_preparation/variable_info")
    load("data/7_24_mike/proteomics/data_preparation/expression_data")
    proteomics_sample_info = sample_info
    proteomics_variable_info = variable_info
    proteomics_expression_data = expression_data
    
    load("data/7_24_mike/summary_info/day_night_df")
    
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
  
  lipidomics_sample_info$accurate_time == sample_info$accurate_time
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
  
  metabolomics_variable_info= 
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
    dplyr::mutate(data_type = "cytokine") %>% 
    dplyr::rename(cortisol_class = class)
  
  total_protein_variable_info = 
    total_protein_variable_info %>% 
    dplyr::mutate(data_type = "total_protein") %>% 
    dplyr::rename(total_protein_class = class)
  
  variable_info =
    rbind(lipidomics_variable_info[,c("variable_id", "mol_name", "data_type")],
          metabolomics_variable_info[,c("variable_id", "mol_name", "data_type")],
          cytokine_variable_info[,c("variable_id", "mol_name", "data_type")],
          total_protein_variable_info[,c("variable_id", "mol_name", "data_type")],
          cortisol_variable_info[,c("variable_id", "mol_name", "data_type")],
          metabolic_panel_variable_info[,c("variable_id", "mol_name", "data_type")],
          proteomics_variable_info[,c("variable_id", "mol_name", "data_type")])
  
  dim(variable_info)
  
  colnames(expression_data) == sample_info$sample_id
  rownames(expression_data) == variable_info$variable_id  
}

######
setwd("data/7_24_mike/inter_omics_correlation/inter_all_omics/")
load("lagged_correlation/lagged_result")
load("lagged_correlation/cor_data")

# global_cor = 
#   cor_data %>% 
#   dplyr::filter(shift_time == "(-30,30]")
# 
# global_cor$p_adjust = p.adjust(global_cor$p, method = "bonferroni",
#                                nrow(expression_data) * (nrow(expression_data) - 1)/2)
# 
# global_cor = 
#   global_cor %>% 
#   dplyr::filter(p_adjust < 0.05)
# 
# plot(global_cor$cor)
# 
# global_cor = 
# global_cor %>% 
#   dplyr::left_join(variable_info, by = c("from" = "variable_id")) %>% 
#   dplyr::rename(from_mol_name = mol_name,
#                 from_data_type = data_type) %>% 
#   dplyr::left_join(variable_info, by = c("to" = "variable_id")) %>% 
#   dplyr::rename(to_mol_name = mol_name,
#                 to_data_type = data_type)
# 
# dim(global_cor)
# 
# table(paste(global_cor$from_data_type, global_cor$to_data_type, sep = "_"))

###set work directory
# dir.create("global correlation network")
setwd("global correlation network")
# save(global_cor, file = "global_cor")
load("global_cor")

dim(global_cor)
length(unique(c(global_cor$from, global_cor$to)))

library(igraph)
library(ggraph)
library(tidygraph)

edge_data = 
  global_cor
  # dplyr::filter(abs(cor) > 0.8 & p_adjust < 0.01) 

edge_data$edge_name =
  paste(edge_data$from, edge_data$to, sep = "_")

edge_data$edge_type =
  paste(edge_data$from_data_type, edge_data$to_data_type, sep = "_")

edge_data = 
edge_data %>% 
  dplyr::mutate(direction = case_when(
    cor > 0 ~"positive",
    cor < 0 ~"negative"
  ))

table(edge_data$edge_type)

##mosaic plot to show the distributation of edge type
library(ggmosaic)

unique(edge_data$edge_type)
color = 
colorRampPalette(colors = ggsci::pal_lancet()(n=9))(n=length(unique(edge_data$edge_type)))

edge_type_color = color
names(edge_type_color) = unique(edge_data$edge_type)

plot = 
edge_data %>%
  dplyr::mutate(direction = factor(direction, levels = c("positive", "negative"))) %>%
  ggplot() +
  geom_mosaic(aes(x = product(edge_type, direction), fill = edge_type),
              offset = 0.05) +
  theme_mosaic()+
  scale_fill_manual(values = edge_type_color)

plot
# ggsave(plot, filename = "edge_class.pdf", width = 9, height = 7)



###show the 
library(gghalves)
alpha_pos = 
edge_data %>% 
  dplyr::filter(cor > 0) %>% 
  dplyr::group_by(edge_type) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(alpha = case_when(
    n > 5000 ~ 0.008,
    n > 2000 & n <= 5000 ~ 0.05,
    n > 1000 & n <= 2000 ~ 0.1,
    n > 100 & n <= 1000 ~ 0.5,
    n > 50 & n <= 100 ~ 0.7,
    n > 0 & n <= 50 ~ 0.9
  ))
alpha = 
alpha_pos$alpha
names(alpha) = alpha_pos$edge_type

plot1 = 
edge_data %>% 
  dplyr::mutate(edge_type = factor(edge_type, 
                                   levels = unique(edge_type))) %>% 
  dplyr::filter(cor > 0) %>%
  ggplot() +
  geom_jitter(aes(x = edge_type, y = cor, 
                  color = edge_type,
                  alpha = edge_type),
              show.legend = FALSE) +
  scale_alpha_manual(values = alpha) +
  geom_half_boxplot(aes(x = edge_type, y = cor),
                    outlier.shape = NA,
                    show.legend = FALSE,
                    fill = "transparent") +
  geom_half_violin(aes(x = edge_type, 
                       y = cor),
                    outlier.shape = NA,
                    show.legend = FALSE, 
                   fill = "transparent",
                   side = "r") +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "", y = "Spearman correlation") +
  scale_x_discrete(breaks=factor(edge_data$edge_type), drop=FALSE) +
  scale_color_manual(values = edge_type_color)

plot1  

alpha_neg = 
  edge_data %>% 
  dplyr::filter(cor < 0) %>% 
  dplyr::group_by(edge_type) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(alpha = case_when(
    n > 5000 ~ 0.008,
    n > 2000 & n <= 5000 ~ 0.05,
    n > 1000 & n <= 2000 ~ 0.1,
    n > 100 & n <= 1000 ~ 0.5,
    n > 50 & n <= 100 ~ 0.7,
    n > 0 & n <= 50 ~ 0.9
  ))
alpha = 
  alpha_neg$alpha
names(alpha) = alpha_neg$edge_type

plot2 = 
  edge_data %>% 
  dplyr::mutate(edge_type = factor(edge_type, levels = unique(edge_type))) %>% 
  dplyr::filter(cor < 0) %>%
  ggplot() +
  geom_jitter(aes(x = edge_type, y = cor, 
                  color = edge_type,
                  alpha = edge_type),
              show.legend = FALSE) +
  scale_alpha_manual(values = alpha) +
  geom_half_boxplot(aes(x = edge_type, y = cor),
                    outlier.shape = NA,
                    show.legend = FALSE,
                    fill = "transparent") +
  geom_half_violin(aes(x = edge_type, y = cor),
                   outlier.shape = NA,
                   show.legend = FALSE, 
                   fill = "transparent",
                   side = "r") +
  base_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 10)) +
  labs(x = "", y = "Spearman correlation") +
  scale_x_discrete(breaks=factor(edge_data$edge_type), drop=FALSE) +
  scale_color_manual(values = edge_type_color)

plot2

library(patchwork)
plot =
  plot1 + plot2 + patchwork::plot_layout(ncol = 1)

ggsave(plot, filename = "cor_distributation.pdf", width = 10, height = 7)

node_data =
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))


temp_data <- 
  tidygraph::tbl_graph(nodes = node_data, 
                       edges = edge_data,
                       directed = FALSE) %>% 
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

angle <- 360 * (c(1:nrow(node_data)) - 0.5)/nrow(node_data)
hjust <- ifelse(angle > 180, 1, 0)
angle <- ifelse(angle > 180, 90 - angle + 180, 90 - angle)



#####so we first cluster for each omics data and then calculated the network for them
###lipidomics data
dir.create("lipidomics")
dir.create("metabolomics")
dir.create("proteomics")
dir.create("cytokine")
dim(lipidomics_expression_data)

library(Mfuzz)
temp_data <- 
  lipidomics_expression_data







