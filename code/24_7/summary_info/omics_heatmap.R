no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

{
  ###wearable
  
  ###internal omics data
  ###cortisol
  load("data/7_24_mike/cortisol/data_preparation/sample_info")
  load("data/7_24_mike/cortisol/data_preparation/variable_info")
  load("data/7_24_mike/cortisol/data_preparation/expression_data")
  cortisol_expression_data = expression_data
  cortisol_sample_info = sample_info
  cortisol_variable_info = variable_info
  
  ###cytokine
  load("data/7_24_mike/cytokine/data_preparation/sample_info")
  load("data/7_24_mike/cytokine/data_preparation/variable_info")
  load("data/7_24_mike/cytokine/data_preparation/expression_data")
  cytokine_sample_info = sample_info
  cytokine_variable_info = variable_info
  cytokine_expression_data = expression_data
  
  ###lipidomics
  load("data/7_24_mike/lipidomics/data_preparation/sample_info")
  load("data/7_24_mike/lipidomics/data_preparation/variable_info")
  load("data/7_24_mike/lipidomics/data_preparation/expression_data")
  lipidomics_sample_info = sample_info
  lipidomics_variable_info = variable_info
  lipidomics_expression_data = expression_data
  
  ###metabolic_panel
  load("data/7_24_mike/metabolic_panel/data_preparation/sample_info")
  load("data/7_24_mike/metabolic_panel/data_preparation/variable_info")
  load("data/7_24_mike/metabolic_panel/data_preparation/expression_data")
  metabolic_panel_sample_info = sample_info
  metabolic_panel_variable_info = variable_info
  metabolic_panel_expression_data = expression_data
  
  ###remove CHEX
  metabolic_panel_variable_info = 
    metabolic_panel_variable_info %>% 
    dplyr::filter(!stringr::str_detect(mol_name, "CHEX"))
  
  metabolic_panel_expression_data = 
    metabolic_panel_expression_data[metabolic_panel_variable_info$variable_id,]
  
  # ###metabolomics
  # load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
  # load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
  # load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
  # metabolomics_sample_info = sample_info
  # metabolomics_variable_info = variable_info
  # metabolomics_expression_data = expression_data
  # 
  # load("data/7_24_mike/metabolomics/data_preparation/peaks/sample_info")
  # load("data/7_24_mike/metabolomics/data_preparation/peaks/variable_info")
  # load("data/7_24_mike/metabolomics/data_preparation/peaks/expression_data")
  # metabolomicspeak_sample_info = sample_info
  # metabolomicspeak_variable_info = variable_info
  # metabolomicspeak_expression_data = expression_data
  # 
  # ###cortisol
  # metabolomicspeak_variable_info$mz
  # metabolomicspeak_variable_info$rt
  
  ###metabolomics
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  metabolomics_expression_data = expression_data
  
  ##remove QC and blank samples
  metabolomics_sample_info = 
    metabolomics_sample_info %>% 
    dplyr::filter(!is.na(accurate_time))
  
  metabolomics_expression_data = 
    metabolomics_expression_data[,metabolomics_sample_info$sample_id]
  
  # ###remove metabolites from other databases
  # table(metabolomics_variable_info$Database)
  # 
  # metabolomics_variable_info = 
  #   metabolomics_variable_info %>% 
  #   dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2", "msDatabase_hilic0.0.2",
  #                                 "msDatabase_rplc0.0.2", "nistDatabase0.0.2"))
  # 
  # metabolomics_expression_data =
  #   metabolomics_expression_data[metabolomics_variable_info$variable_id, ]
  
  ###remove outlier samples
  metabolomics_sample_info = 
    metabolomics_sample_info %>% 
    dplyr::filter(!sample_id %in% c("67", "91", "104", "74", "58", "47", "2")) 
  
  metabolomics_expression_data = 
    metabolomics_expression_data[,metabolomics_sample_info$sample_id]
  
  ###proteomics
  load("data/7_24_mike/proteomics/data_preparation/sample_info")
  load("data/7_24_mike/proteomics/data_preparation/variable_info")
  load("data/7_24_mike/proteomics/data_preparation/expression_data")
  proteomics_sample_info = sample_info
  proteomics_variable_info = variable_info
  proteomics_expression_data = expression_data
  
  ###total protein
  load("data/7_24_mike/total_protein/data_preparation/sample_info")
  load("data/7_24_mike/total_protein/data_preparation/variable_info")
  load("data/7_24_mike/total_protein/data_preparation/expression_data")
  total_protein_sample_info = sample_info
  total_protein_variable_info = variable_info
  total_protein_expression_data = expression_data
}


setwd("data/7_24_mike/summary_info")

{
  dim(proteomics_expression_data)
  dim(metabolomics_expression_data)
  dim(metabolic_panel_expression_data)
  dim(lipidomics_expression_data)
  dim(cytokine_expression_data)
  dim(cortisol_expression_data)
  dim(total_protein_expression_data)
  
  proteomics_sample_info$accurate_time
  metabolomics_sample_info$accurate_time
  metabolic_panel_sample_info$accurate_time
  lipidomics_sample_info$accurate_time
  cytokine_sample_info$accurate_time
  cortisol_sample_info$accurate_time
  total_protein_sample_info$accurate_time
  
  colnames(proteomics_expression_data) = as.character(proteomics_sample_info$accurate_time)
  colnames(metabolomics_expression_data) = as.character(metabolomics_sample_info$accurate_time)
  colnames(metabolic_panel_expression_data) = as.character(metabolic_panel_sample_info$accurate_time)
  colnames(lipidomics_expression_data) = as.character(lipidomics_sample_info$accurate_time)
  colnames(cytokine_expression_data) = as.character(cytokine_sample_info$accurate_time)
  colnames(cortisol_expression_data) = as.character(cortisol_sample_info$accurate_time)
  colnames(total_protein_expression_data) = as.character(total_protein_sample_info$accurate_time)
  
  intersect_sample_id =
    Reduce(f = union, 
           x = list(colnames(proteomics_expression_data),
                    colnames(metabolomics_expression_data),
                    colnames(metabolic_panel_expression_data),
                    colnames(lipidomics_expression_data),
                    colnames(cytokine_expression_data),
                    colnames(cortisol_expression_data),
                    colnames(total_protein_expression_data)
           )) %>% 
    sort()
  
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
  
  proteomics_expression_data1 = 
    cbind(proteomics_expression_data,
          setNames(
            lapply(setdiff(
              intersect_sample_id,
              colnames(proteomics_expression_data)
            ),
            function(x)
              x = NA),
            setdiff(intersect_sample_id,
                    colnames(proteomics_expression_data))
          )) %>% 
    dplyr::select(intersect_sample_id)
  
  metabolomics_expression_data1 = 
    cbind(metabolomics_expression_data,
          setNames(
            lapply(setdiff(
              intersect_sample_id,
              colnames(metabolomics_expression_data)
            ),
            function(x)
              x = NA),
            setdiff(intersect_sample_id,
                    colnames(metabolomics_expression_data))
          )) %>% 
    dplyr::select(intersect_sample_id)
  
  metabolic_panel_expression_data1 = 
    cbind(metabolic_panel_expression_data,
          setNames(
            lapply(setdiff(
              intersect_sample_id,
              colnames(metabolic_panel_expression_data)
            ),
            function(x)
              x = NA),
            setdiff(intersect_sample_id,
                    colnames(metabolic_panel_expression_data))
          )) %>% 
    dplyr::select(intersect_sample_id)
  
  lipidomics_expression_data1 = 
    cbind(lipidomics_expression_data,
          setNames(
            lapply(setdiff(
              intersect_sample_id,
              colnames(lipidomics_expression_data)
            ),
            function(x)
              x = NA),
            setdiff(intersect_sample_id,
                    colnames(lipidomics_expression_data))
          )) %>% 
    dplyr::select(intersect_sample_id)
  
  cytokine_expression_data1 = 
    cbind(cytokine_expression_data,
          setNames(
            lapply(setdiff(
              intersect_sample_id,
              colnames(cytokine_expression_data)
            ),
            function(x)
              x = NA),
            setdiff(intersect_sample_id,
                    colnames(cytokine_expression_data))
          )) %>% 
    dplyr::select(intersect_sample_id)
  
  cortisol_expression_data1 = 
    cbind(cortisol_expression_data,
          setNames(
            lapply(setdiff(
              intersect_sample_id,
              colnames(cortisol_expression_data)
            ),
            function(x)
              x = NA),
            setdiff(intersect_sample_id,
                    colnames(cortisol_expression_data))
          )) %>% 
    dplyr::select(intersect_sample_id)
  
  total_protein_expression_data1 = 
    cbind(total_protein_expression_data,
          setNames(
            lapply(setdiff(
              intersect_sample_id,
              colnames(total_protein_expression_data)
            ),
            function(x)
              x = NA),
            setdiff(intersect_sample_id,
                    colnames(total_protein_expression_data))
          )) %>% 
    dplyr::select(intersect_sample_id)
  
}

temp_data = 
  rbind(data.frame(proteomics_expression_data1, class = "proteomics", check.names = FALSE),
        data.frame(metabolomics_expression_data1, class = "metabolomics", check.names = FALSE),
        data.frame(metabolic_panel_expression_data1, class = "metabolic_panel", check.names = FALSE),
        data.frame(lipidomics_expression_data1, class = "lipidomics", check.names = FALSE),
        data.frame(cytokine_expression_data1, class = "cytokine", check.names = FALSE),
        data.frame(cortisol_expression_data1, class = "cortisol", check.names = FALSE),
        data.frame(total_protein_expression_data1, class = "total_protein", check.names = FALSE))

dim(temp_data)

library(ComplexHeatmap)

class = temp_data$class
temp_data = temp_data %>% 
  dplyr::select(-class)

library(lubridate)

temp_data = 
  temp_data %>% 
  apply(1, function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }) %>% 
  t() %>% 
  as.data.frame()

remove_idx = 
temp_data %>% 
  apply(1, function(x){
    all(is.na(x))
  }) %>% 
  which()
remove_idx
temp_data = temp_data[-remove_idx,]
class = class[-remove_idx]

range(temp_data, na.rm = TRUE)
temp_data[temp_data > 4] = 4
temp_data[temp_data < -4] = -4

dim(temp_data)

dim(sample_info)

# sample_info2 =
#   sample_info %>%
#   dplyr::mutate(sample_id = as.character(accurate_time))

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

table(class)
class[class == "cortisol" | class == "metabolic_panel" | class == "total_protein" |
        class == "cytokine"] =
  "Targeted assay"

library(circlize)

col_fun = colorRamp2(
  breaks = c(-3, 0, 3),
  colors =
    c("#366A9FFF", "white", "red"),
  transparency = 0
)

table(class)
omics_color = c(class_color["proteomics"],
                class_color["metabolomics"],
                class_color["cytokine"],
                class_color["lipidomics"]
                )

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
  # row_split = class,
  # row_title = rep("", 4),
  # row_title_rot = 0,
  na_col = "grey"
  # right_annotation = rowAnnotation(foo = anno_block(
  #   gp = gpar(fill = omics_color),
  #   labels = c("Lipidomics", "Targeted assay", "Proteomics", "Metabolomics"),
  #   labels_gp = gpar(col = "white", fontsize = 10)
  # ))
)

rownames(temp_data)[row_order(plot)]

plot

temp_data =
  temp_data[rownames(temp_data)[row_order(plot)],]

#####use ggplot2 to get the heatmap
colnames(temp_data)

plot1 = ggplotify::as.ggplot(plot)

plot1

ggsave(plot1, filename = "heatmap.dent.pdf", width = 7, height = 7)

time = lubridate::as_datetime(colnames(temp_data), tz = "America/Los_Angeles")

min_time = lubridate::as_datetime("2019-04-29 00:00:01", tz = "America/Los_Angeles")
max_time = lubridate::as_datetime("2019-05-07 23:59:00", tz = "America/Los_Angeles")

time[length(time)] = max_time

colnames(temp_data)

time1 = c(min_time, time[-length(time)])
time2 = c(time)

names(time1) = 
  names(time2) = 
  colnames(temp_data)

time = 
  data.frame(sample_name = colnames(temp_data),
             time1, 
             time2)

temp_data1 = 
  temp_data

temp_data1 = 
temp_data1 %>% 
  tibble::rowid_to_column(var = "variable_id") %>% 
  dplyr::mutate(variable_id = as.numeric(variable_id)) %>% 
  tidyr::pivot_longer(cols = -variable_id, names_to = "sample_name", values_to = "fill") %>% 
  dplyr::left_join(time, by = "sample_name")

range(temp_data1$fill, na.rm = TRUE)

plot = 
ggplot(temp_data1, aes(
  xmin = time1,
  xmax = time2,
  ymin = variable_id,
  ymax = variable_id + 1
)) +
  geom_rect(aes(fill = fill), colour = NA) +
  scale_fill_gradient2(low = "#366A9FFF", mid = "white", high = "red") +
  scale_x_datetime(
    breaks = date_breaks("12 hour"),
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
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

plot
ggsave(plot, filename = "molecular_heatmap.pdf", width = 15, height = 5)


table(class)
dim(lipido)

table(class)
dim(lipidomics_variable_info)
dim(metabolomics_variable_info)
dim(proteomics_variable_info)
dim(total_protein_variable_info)
dim(cortisol_variable_info)
dim(metabolic_panel_variable_info)
dim(cytokine_variable_info)

df = 
data.frame(value = as.numeric(table(class)),
           group = names(table(class)))

library(ggpubr)
plot = 
ggpie(
  data = df,
  x = "value",
  label = "group",
  fill = "group",
  color = "white",
  palette = omics_color[df$group]
)

ggsave(plot, file = "feature_number.pdf", width = 7, height = 7)


df = 
  data.frame(value = c(nrow(cytokine_variable_info), 
                       nrow(total_protein_variable_info), 
                       nrow(cortisol_variable_info), 
                       nrow(metabolic_panel_variable_info)),
             group = c("cytokine", "total_protein", "cortisol", "metabolic_panel"))

library(ggpubr)
plot = 
  ggpie(
    data = df,
    x = "value",
    label = "group",
    fill = "group",
    color = "white",
    palette = class_color[df$group]
  )
plot
ggsave(plot, file = "feature_number2.pdf", width = 7, height = 7)



