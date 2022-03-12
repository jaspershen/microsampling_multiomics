no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

{
  ####load data
  ###Step
  load("data/7_24_mike/steps/data_preparation/sample_info")
  load("data/7_24_mike/steps/data_preparation/variable_info")
  load("data/7_24_mike/steps/data_preparation/expression_data")
  steps_expression_data = expression_data
  steps_sample_info = sample_info
  steps_variable_info = variable_info
  
  ###lipidomics
  load("data/7_24_mike/lipidomics/data_preparation/sample_info")
  load("data/7_24_mike/lipidomics/data_preparation/variable_info")
  load("data/7_24_mike/lipidomics/data_preparation/expression_data")
  lipidomics_sample_info = sample_info
  lipidomics_variable_info = variable_info
  lipidomics_expression_data = expression_data
  
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

######steps vs lipidomics
dir.create("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_lipidomics")
setwd("data/7_24_mike/wearable_omics_correlation/steps_omics_correlation/steps_lipidomics")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between steps and lipidomicss
#global correlation

dir.create("lagged_correlation")

# lagged_cor = rep(NA, nrow(lipidomics_expression_data))
# global_cor = rep(NA, nrow(lipidomics_expression_data))
# 
# lagged_result = vector(mode = "list", length = nrow(lipidomics_expression_data))
# 
# for(i in 1:nrow(lipidomics_expression_data)){
#   cat(i, " ")
#   x = as.numeric(lipidomics_expression_data[i, ])
#   time1 = lipidomics_sample_info$accurate_time
#   y = as.numeric(steps_expression_data[1, ])
#   time2 = steps_sample_info$accurate_time
# 
#   result = lagged_correlation(
#     x = x,
#     y = y,
#     time1 = time1,
#     time2 = time2,
#     time_tol = 60/60,
#     step = 5/60
#   )
# 
#   lagged_cor[i] = result$max_cor
#   global_cor[i] = result$global_cor
#   lagged_result[[i]] = result
# }
# names(lagged_result) = rownames(lipidomics_expression_data)
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
#   lipidomics_variable_info$variable_id
# 
# cor_data =
#   data.frame(wearable = "Step",
#              lipidomics_variable_info,
#              global_cor = global_cor,
#              lagged_cor = lagged_cor,
#              shift_time = shift_time) %>% 
#   dplyr::filter(abs(lagged_cor) > 0.2)
# 
# cor_data$KEGG_ID = stringr::str_replace(cor_data$KEGG_ID, "\xa0", "")
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
# addWorksheet(wb, sheetName = "Step lipidomics global cor",
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
  dplyr::filter(lagged_cor > 0) %>% 
  tail(10)

neg_top_10 =
  cor_data %>% 
  dplyr::arrange(lagged_cor) %>% 
  dplyr::filter(lagged_cor < 0) %>% 
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
#                         internal_omics_color = class_color["lipidomics"],
#                         wearable_color = wearable_color["step"],
#                         internal_omics_name = temp$mol_name[i],
#                         warable_name = "Step",
#                         which = "max",
#                         x_limit = c(1,1000),
#                         non_matched_point_size = 0.1,
#                         wearable_point_size = 0.5,
#                         internal_omics_point_size = 2,
#                         integrated = FALSE)
# 
#   plot2 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["lipidomics"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
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
#                           internal_omics_color = class_color["lipidomics"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
#                           which = "max",
#                           x_limit = c(1,1000),
#                           non_matched_point_size = 0.1,
#                           wearable_point_size = 0.5,
#                           internal_omics_point_size = 2,
#                           integrated = TRUE)
# 
# 
#   plot4 =
#     lagged_alignment_plot(object = lagged_result[[temp$variable_id[i]]],
#                           day_night_df = day_night_df,
#                           internal_omics_color = class_color["lipidomics"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$mol_name[i],
#                           warable_name = "Step",
#                           which = "max",
#                           x_limit = c(1,30),
#                           non_matched_point_size = 3,
#                           wearable_point_size = 3,
#                           internal_omics_point_size = 3,
#                           integrated = TRUE)
# 
#   name = paste("Step vs",temp$mol_name[i])
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

plot =
  cor_data %>%
  dplyr::mutate(direction = case_when(
    shift_time > 0 ~ "After",
    shift_time < 0 ~ "Before",
    shift_time == 0 ~ "Synchronization"
  )) %>%
  ggplot(aes(global_cor, lagged_cor)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(fill = direction,
                 size = -log(lagged_cor_p_adjust, 10)),
             shape = 21) +
  labs(x = "Global Pearson correlation",
       y = "Lagged correlation") +
  scale_fill_manual(values = c(
    "After" = ggsci::pal_aaas()(n = 10)[1],
    "Before" = ggsci::pal_aaas()(n = 10)[2],
    "Synchronization" = "grey"
  )) +
  guides(size = guide_legend(title = "-log10(p.adjust)")) +
  scale_size_continuous(range = c(0.1,5)) +
  shadowtext::geom_shadowtext(aes(label = ifelse(lagged_cor_p_adjust < 0.05,
                                                 mol_name, NA)),
                              check_overlap = TRUE,
                              bg.colour='white',
                              color = "black") +
  base_theme

plot

# ggsave(plot, filename = "global_lagged_correlation.pdf", width = 9, height = 7)

####output the shift time vs lagged cor plot for the important lipids
important_lipid =
  cor_data %>%
  dplyr::filter(lagged_cor > quantile(cor_data$lagged_cor[cor_data$lagged_cor > 0], 0.75) |
                  lagged_cor < quantile(cor_data$lagged_cor[cor_data$lagged_cor < 0], 0.25))

dir.create("shift_time_vs_cor")

# for(i in 1:nrow(important_lipid)){
#   cat(i, "")
#   result = lagged_result[[important_lipid$variable_id[i]]]
#   
#   temp_data =
#     result[c("shift_time", "all_cor")] %>%
#     do.call(cbind, .) %>%
#     as.data.frame() %>%
#     dplyr::mutate(shift_time = stringr::str_replace(shift_time, "\\(", "")) %>%
#     dplyr::mutate(shift_time = stringr::str_replace(shift_time, "\\]", "")) %>%
#     dplyr::mutate(shift_time = stringr::str_split(shift_time, ",")) %>%
#     dplyr::mutate(all_cor = round(as.numeric(all_cor), 4))
#   
#   temp_data$shift_time =
#     temp_data$shift_time %>%
#     purrr::map(function(x){
#       mean(as.numeric(x))
#     }) %>%
#     unlist()
#   
#   plot =
#     temp_data %>%
#     ggplot(aes(x = shift_time, y = all_cor)) +
#     geom_vline(xintercept = temp_data$shift_time[result$which_max_idx],
#                color = "red") +
#     geom_hline(yintercept = 0) +
#     annotate(geom = "text",
#              x = temp_data$shift_time[result$which_max_idx],
#              y = result$max_cor,
#              label = round(result$max_cor, 4)) +
#     annotate(geom = "text",
#              x = temp_data$shift_time[result$which_global_idx],
#              y = result$global_cor,
#              label = round(result$global_cor, 4)) +
#     geom_point() +
#     geom_line(aes(group = 1)) +
#     base_theme +
#     labs(x = "Shift time (Omics - HR, min)",
#          y = "Pearsom correlation") +
#     theme()
#   
#   name =
#     lipidomics_variable_info$mol_name[match(names(lagged_result)[i], lipidomics_variable_info$variable_id)]
#   
#   ggsave(
#     plot,
#     file = file.path("shift_time_vs_cor", paste(name, ".pdf", sep = "")),
#     width = 8,
#     height = 7
#   )
# }


##Pathway analysis for the negative or positive correlation analysis
##because the correlation is too small, so here we just use the cutoff from
##0.75
##Pathway analysis for the negative or positive correlation analysis
library(clusterProfiler)
library(org.Hs.eg.db)

dir.create("pathway_enrichment")

important_lipid =
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

# save(important_lipid, file = "important_lipid")
load("important_lipid")

table(important_lipid$class1)

dir.create("lipidminion")

dir.create("lipidminion/positive")
dir.create("lipidminion/negative")

variable_info = 
  variable_info %>% 
  dplyr::mutate(Lipid_Name = case_when(
    is.na(Lipid_Name) ~ mol_name,
    !is.na(Lipid_Name) ~ Lipid_Name
  ))

variable_info$Lipid_Name = 
  variable_info$Lipid_Name %>% 
  stringr::str_replace_all("\\_", "\\/")

variable_info$Lipid_Name[grep("TAG", variable_info$Lipid_Name)] =
  variable_info$Lipid_Name[grep("TAG", variable_info$Lipid_Name)] %>% 
  purrr::map(function(x){
    # main = 
    stringr::str_extract(x, "TAG[0-9]{1,3}\\.[0-9]{1,2}") %>% 
      stringr::str_replace("TAG", "") %>% 
      stringr::str_split("\\.") %>% 
      `[[`(1) %>% 
      paste(collapse = ":") %>% 
      paste("TAG(",., ")",sep = "")
    
    # chain = 
    #   stringr::str_extract(x, "FA[0-9]{1,3}\\.[0-9]{1,2}") %>% 
    #   stringr::str_replace("FA", "") %>% 
    #   stringr::str_split("\\.") %>% 
    #   `[[`(1) %>% 
    #   paste(collapse = ":") %>% 
    #   paste("(FA ",.,")", sep = "")
    # paste(main, chain, sep = " ")
  }) %>% 
  unlist()

###TAG to TG and DAG to TG
variable_info$Lipid_Name = 
  variable_info$Lipid_Name %>% 
  stringr::str_replace_all("TAG", "TG") %>% 
  stringr::str_replace_all("DAG", "DG") 

important_lipid = 
  important_lipid %>% 
  dplyr::mutate(Lipid_Name = case_when(
    is.na(Lipid_Name) ~ mol_name,
    !is.na(Lipid_Name) ~ Lipid_Name
  ))

important_lipid$Lipid_Name = 
  important_lipid$Lipid_Name %>% 
  stringr::str_replace_all("\\_", "\\/")

important_lipid$Lipid_Name[grep("TAG", important_lipid$Lipid_Name)] =
  important_lipid$Lipid_Name[grep("TAG", important_lipid$Lipid_Name)] %>% 
  purrr::map(function(x){
    # main = 
    stringr::str_extract(x, "TAG[0-9]{1,3}\\.[0-9]{1,2}") %>% 
      stringr::str_replace("TAG", "") %>% 
      stringr::str_split("\\.") %>% 
      `[[`(1) %>% 
      paste(collapse = ":") %>% 
      paste("TAG(",., ")",sep = "")
    
    # chain = 
    #   stringr::str_extract(x, "FA[0-9]{1,3}\\.[0-9]{1,2}") %>% 
    #   stringr::str_replace("FA", "") %>% 
    #   stringr::str_split("\\.") %>% 
    #   `[[`(1) %>% 
    #   paste(collapse = ":") %>% 
    #   paste("(FA ",.,")", sep = "")
    # paste(main, chain, sep = " ")
  }) %>% 
  unlist()

important_lipid$Lipid_Name = 
  important_lipid$Lipid_Name %>% 
  stringr::str_replace_all("TAG", "TG") %>% 
  stringr::str_replace_all("DAG", "DG") 


positive_lipid = important_lipid %>% 
  dplyr::filter(class1 == "positive correlation") %>% 
  # dplyr::select(Lipid_Name) %>% 
  dplyr::filter(!is.na(Lipid_Name)) %>% 
  dplyr::rename(lipid = Lipid_Name) %>% 
  dplyr::distinct(lipid, .keep_all = TRUE) 

negative_lipid = important_lipid %>% 
  dplyr::filter(class1 == "negative correlation") %>% 
  # dplyr::select(Lipid_Name) %>% 
  dplyr::filter(!is.na(Lipid_Name)) %>% 
  dplyr::rename(lipid = Lipid_Name) %>% 
  dplyr::distinct(lipid, .keep_all = TRUE)


##because we have removed the carbon chain of the TG, so we need to remove the overlaping between
####remove the same lipid (TG) which are in both 
name = 
intersect(positive_lipid$lipid,
          negative_lipid$lipid) %>%
  data.frame(name = .) %>%
  dplyr::left_join(positive_lipid[, c("lipid", "lagged_cor")], by = c("name" = "lipid")) %>%
  dplyr::rename(cor_pos = lagged_cor) %>%
  dplyr::left_join(negative_lipid[, c("lipid", "lagged_cor")], by = c("name" = "lipid")) %>%
  dplyr::rename(cor_neg = lagged_cor)

positve_remove_name = name$name[which(name$cor_pos < abs(name$cor_neg))]  
negative_remove_name = name$name[which(name$cor_pos > abs(name$cor_neg))]  

positive_lipid = important_lipid %>% 
  dplyr::filter(class1 == "positive correlation") %>% 
  dplyr::select(Lipid_Name) %>% 
  dplyr::filter(!is.na(Lipid_Name)) %>% 
  dplyr::rename(lipid = Lipid_Name) %>% 
  dplyr::distinct(lipid) %>% 
  dplyr::filter(!lipid %in% positve_remove_name)

negative_lipid = important_lipid %>% 
  dplyr::filter(class1 == "negative correlation") %>% 
  dplyr::select(Lipid_Name) %>% 
  dplyr::filter(!is.na(Lipid_Name)) %>% 
  dplyr::rename(lipid = Lipid_Name) %>% 
  dplyr::distinct(lipid) %>% 
  dplyr::filter(!lipid %in% negative_remove_name)


# write.table(positive_lipid, file = "lipidminion/positive/positive_lipid.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(negative_lipid, file = "lipidminion/negative/negative_lipid.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

universe_lipid = 
  variable_info %>% 
  dplyr::select(Lipid_Name) %>% 
  dplyr::filter(!is.na(Lipid_Name)) %>% 
  dplyr::rename(lipid = Lipid_Name) %>% 
  dplyr::distinct(lipid)

# write.table(universe_lipid, file = "lipidminion/universe_lipid.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)




#####read the lipid minion results
pos_enriched = 
  readr::read_delim("lipidminion/positive/Fisher output table (6 pvals _ 0.05).txt", delim = " ")

pos_network_node = readr::read_delim("lipidminion/positive/Network_nodes.txt", delim = "\t")
pos_network_edge = readr::read_delim("lipidminion/positive/Network_edges.txt", delim = "\t")
pos_network_edge_attr = readr::read_delim("lipidminion/positive/Network_edge_attributes.txt", delim = "\t")

neg_enriched = 
  readr::read_delim("lipidminion/negative/Fisher output table (3 pvals _ 0.05).txt", delim = " ")

neg_network_node = readr::read_delim("lipidminion/negative/Network_nodes.txt", delim = "\t")
neg_network_edge = readr::read_delim("lipidminion/negative/Network_edges.txt", delim = "\t")
neg_network_edge_attr = readr::read_delim("lipidminion/negative/Network_edge_attributes.txt", delim = "\t")

unique(pos_network_node$title)

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
addWorksheet(wb, sheetName = "Positive correlation", gridLines = TRUE)
addWorksheet(wb, sheetName = "Negative correlation", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
writeDataTable(wb, sheet = 1, x = pos_enriched,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = neg_enriched,
               colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "lipidminion/lipid_enrichment.xlsx", overwrite = TRUE)

###network
###positive network
pos_network_node

library(igraph)
library(ggraph)
library(tidygraph)

pos_network_node

node = 
  pos_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>% 
  dplyr::mutate(class = case_when(
    shape == "#cccccc" ~ "lipid",
    TRUE ~ "class"
  ))

node$title[node$title == "Uncategorized"] = "Free facty acid"

edge = 
  pos_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

pos_netwrok = 
  tidygraph::tbl_graph(nodes = node, 
                       edges = edge)

plot = 
  ggraph(pos_netwrok,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "grey"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title
    ),
    color = "black",
    bg.color = "white",
    size = 4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "lipid" = unname(class_color["lipidomics"]),
    "class" = "red"
  )) +
  scale_size_manual(values = c(
    "lipid" = 8,
    "class" = 10
  )) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) 

plot

# ggsave(plot, filename = "lipidminion/positive/pos_network.pdf", width = 7, height = 7)


###negative network
neg_network_node

library(igraph)
library(ggraph)
library(tidygraph)

neg_network_node

node = 
  neg_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>% 
  dplyr::mutate(class = case_when(
    shape == "#cccccc" ~ "lipid",
    TRUE ~ "class"
  ))

edge = 
  neg_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

neg_netwrok = 
  tidygraph::tbl_graph(nodes = node, 
                       edges = edge)

plot = 
  ggraph(neg_netwrok,
         layout = 'stress',
         circular = FALSE) +
  geom_edge_link(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "grey"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = FALSE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title
    ),
    color = "black",
    bg.color = "white",
    size = 4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "lipid" = unname(class_color["lipidomics"]),
    "class" = "red"
  )) +
  scale_size_manual(values = c(
    "lipid" = 8,
    "class" = 10
  )) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) 

plot

# ggsave(plot, filename = "lipidminion/negative/neg_network.pdf", width = 7, height = 7)



