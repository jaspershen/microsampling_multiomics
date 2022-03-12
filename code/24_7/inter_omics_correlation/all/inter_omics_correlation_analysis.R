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
setwd("data/7_24_mike/inter_omics_correlation/inter_all_omics")
load("lagged_correlation/lagged_result")
load("lagged_correlation/cor_data")

# lagged_cor =
#   readxl::read_xlsx("lagged_correlation/cor_data.xlsx",1)
# 
# table(lagged_cor$shift_time)
# 
# lagged_cor %>%
#   dplyr::filter(shift_time != "(-30,30]") %>%
#   pull(cor) %>%
#   plot()
# 
# unique(c(lagged_cor$from,
#          lagged_cor$to))
# 
# dim(lagged_cor)
# 
# lagged_cor =
#   lagged_cor %>%
#   dplyr::left_join(variable_info, by = c("from" = "variable_id")) %>%
#   dplyr::rename(from_mol_name = mol_name,
#                 from_data_type = data_type) %>%
#   dplyr::left_join(variable_info, by = c("to" = "variable_id")) %>%
#   dplyr::rename(to_mol_name = mol_name,
#                 to_data_type = data_type)
# 
# 
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "lagged correlation",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = lagged_cor,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "lagged_correlation/lagged_cor.xlsx", overwrite = TRUE)

lagged_cor =
  readxl::read_xlsx("lagged_correlation/lagged_cor.xlsx")

table(paste(lagged_cor$from_data_type, lagged_cor$to_data_type, sep = "_"))


temp = lagged_result[[2]][[2]]

temp$max_idx
temp$max_cor
temp$global_cor

x = temp$x
y = purrr::map(.x = temp$max_idx, function(idx){
  mean(temp$y[idx])
}) %>% 
  unlist()

cor(x,y, method = "spearman")
cor(x,y, method = "pearson")
temp$max_cor

####
lagged_cor1 = 
  lagged_cor %>% 
  dplyr::filter(shift_time != "(-30,30]")

lagged_cor2 = 
  lagged_cor %>% 
  dplyr::filter(shift_time == "(-30,30]")

dim(lagged_cor1)

table(lagged_cor1$shift_time)

paste(lagged_cor1$from_data_type, 
      lagged_cor1$to_data_type, sep = "_") %>% 
  table()

table((c(lagged_cor1$from_data_type,
  lagged_cor1$to_data_type)))

#####lagged correlation network
####output the shift_time_vs_cor
dim(lagged_cor1)
length(unique(c(lagged_cor1$from,lagged_cor1$to)))

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "lagged correlation",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = lagged_cor1,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "lagged_cor_network//lagged_cor1.xlsx", overwrite = TRUE)


# #####shift time vs correlation
# for(i in 1:nrow(lagged_cor1)){
#   cat(i, " ")
#   from = lagged_cor1$from[i]
#   to = lagged_cor1$to[i]
#   idx1 = which(names(lagged_result) == from)
#   idx2 = which(names(lagged_result[[idx1]]) == to)
#   result =
#     evaluate_peak_quality(object = lagged_result[[idx1]][[idx2]], plot = TRUE)
#   name = paste(round(result$score, 4), "_",from, "_",to, ".pdf",sep = "") %>%
#     stringr::str_replace_all("/", '_')
#   ggsave(result$plot, filename = file.path("shift_time_vs_cor", name), width = 9, height = 7)
# }

# ##output cor plot
# for (i in 1:nrow(lagged_cor1)) {
#   cat(i, " ")
# 
#   plot1 =
#     lagged_alignment_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       day_night_df = day_night_df,
#       internal_omics_color = class_color[lagged_cor1$to_data_type[i]],
#       wearable_color = class_color[lagged_cor1$from_data_type[i]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "max",
#       x_limit = c(1, 1000),
#       non_matched_point_size = 2,
#       wearable_point_size = 2,
#       internal_omics_point_size = 2,
#       integrated = FALSE
#     )
# 
#   name = lagged_cor1$name[i] %>%
#     stringr::str_replace_all("/", "_")
#   ggsave(plot1,
#          filename = file.path("cor_plot", paste(name, "plot1.pdf", sep = "")),
#          width = 20, height = 7)
# 
# 
#   plot2 =
#     lagged_alignment_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       day_night_df = day_night_df,
#       internal_omics_color = class_color[lagged_cor1$to_data_type[i]],
#       wearable_color = class_color[lagged_cor1$from_data_type[i]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "max",
#       x_limit = c(1, 1000),
#       non_matched_point_size = 2,
#       wearable_point_size = 2,
#       internal_omics_point_size = 2,
#       add_connect_line = FALSE,
#       integrated = FALSE
#     )
# 
#   ggsave(plot2,
#          filename = file.path("cor_plot", paste(name, "plot2.pdf", sep = "")),
#          width = 20, height = 7)
# }


##output scatter plot
dir.create("scatter_plot")
# for (i in 1:nrow(lagged_cor1)) {
#   cat(i, " ")
# 
#   plot1 =
#     lagged_sactter_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "global"
#     )
# 
#   plot2 =
#     lagged_sactter_plot(
#       object = lagged_result[[lagged_cor1$from[i]]][[lagged_cor1$to[i]]],
#       internal_omics_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$to[i]],
#       warable_name = variable_info$mol_name[variable_info$variable_id == lagged_cor1$from[i]],
#       which = "max"
#     )
# 
#   library(patchwork)
# 
#   plot =
#   plot1 + plot2 + patchwork::plot_layout(ncol = 2)
# 
#   name = lagged_cor1$name[i] %>%
#     stringr::str_replace_all("/", "_")
# 
#   ggsave(plot,
#          filename = file.path("scatter_plot", paste(name, "plot.pdf", sep = "")),
#          width = 14, height = 7)
# 
# }



####construct the network for lagged correlation
##############total network
temp_data = 
  lagged_cor1 %>% 
  dplyr::filter(p_adjust < 0.05)

dim(temp_data)

edge_data =
  temp_data %>% 
  dplyr::select(from, to, shift_time, cor, p_adjust, dplyr::everything()) %>% 
  dplyr::mutate(time = 
                  stringr::str_replace_all(shift_time, "\\(|\\]", "") %>% 
                  stringr::str_split(",") %>% 
                  purrr::map(function(x){mean(as.numeric(x))}) %>% 
                  unlist()) %>% 
  dplyr::mutate(
    shift = case_when(
      time > 0 ~ "1",
      time < 0 ~ "-1"
    )
  )

edge_data$p_adjust[edge_data$p_adjust == 0] = 
  min(edge_data$p_adjust[edge_data$p_adjust != 0])

node_data = 
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>% 
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))
  
dim(node_data)
dim(edge_data)

library(igraph)
library(tidygraph)

# node_data = 
#   node_data %>% 
#   dplyr::arrange(class) %>% 
#   dplyr::mutate(node = factor(node, levels = node))

total_lagged_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

#####up-down
# g <- total_lagged_graph
# library(igraph)
# library(ggraph)
# 
# V(g)$type <- bipartite_mapping(g)$type
# 
# coords <-
#   create_layout(g, layout = "bipartite") %>%
#   dplyr::select(node, class, mol_name, class2, x, y)
# 
# coords$y[coords$class2 == "wearable"] = 0.3
# coords$y[coords$class2 == "internal-omics"] = 1
# 
# table(coords$class)
# 
# coords$y[coords$class2 == "wearable"] = -0.5
# coords$y[coords$class == "metabolomics"] = 1
# coords$y[coords$class == "lipidomics"] = 0.5
# coords$y[coords$class == "proteomics"] = 0
# coords$y[coords$class == "cytokine"] = 0
# coords$y[coords$class == "metabolic_panel"] = 0
# coords$y[coords$class == "total_protein"] = 0
# 
# coords <-
#   coords %>%
#   dplyr::select(x, y) %>%
#   dplyr::mutate(
#     theta = x / (max(x) + 1) * 2 * pi,
#     r = y + 1,
#     x = r * cos(theta),
#     y = r * sin(theta)
#   )
# 
# my_graph <-
#   create_layout(
#     graph = g,
#     layout = "manual",
#     x = coords$x,
#     y = coords$y
#     # node.position = coords
#   )

library(ggraph)
plot <-
  ggraph(total_lagged_graph, layout = "fr") +
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
        size = Degree),
    shape = 21,
    alpha = 0.7,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(aes(
    x = x,
    y = y,
    color = data_type,
    label = ifelse(data_type != "lipidomics", mol_name, NA)
  ), check_overlap = TRUE, bg.colour = "white") +
  shadowtext::geom_shadowtext(aes(
    x = x,
    y = y,
    color = data_type,
    label = ifelse(data_type == "lipidomics" & Degree > 10, mol_name, NA)
  ), check_overlap = TRUE, bg.colour = "white") +
  scale_size_continuous(range = c(2, 8)) +
  scale_fill_manual(values = c(class_color)) +
  scale_color_manual(values = c(class_color)) +
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
#   filename = (
#     "lagged_cor_network/total_lagged_network.pdf"
#   ),
#   width = 8.6,
#   height = 7
# )

####output the network data
node_data = 
igraph::vertex.attributes(graph = total_lagged_graph) %>% 
  dplyr::bind_cols() %>% 
  as.data.frame() %>% 
  dplyr::arrange(desc(Degree))

# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "Node data",
#              gridLines = TRUE)
# addWorksheet(wb, sheetName = "Edge data",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = node_data,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 2, x = edge_data,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "lagged_cor_network/network_data.xlsx", overwrite = TRUE)

node_data

