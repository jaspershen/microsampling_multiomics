no_function()

library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")
source("code/modified_dtw.R")
source("code/lagged_correlation.R")

{
  ####load data
  ###Step
  load("data/24_7_study/steps/data_preparation/sample_info")
  load("data/24_7_study/steps/data_preparation/variable_info")
  load("data/24_7_study/steps/data_preparation/expression_data")
  steps_expression_data = expression_data
  steps_sample_info = sample_info
  steps_variable_info = variable_info
  
  ###metabolomics
  load("data/24_7_study/metabolomics/data_preparation/metabolites/sample_info")
  load("data/24_7_study/metabolomics/data_preparation/metabolites/variable_info")
  load("data/24_7_study/metabolomics/data_preparation/metabolites/expression_data")
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  metabolomics_expression_data = expression_data
  
  plot(density(metabolomics_expression_data$`52`))
  plot(density(metabolomics_expression_data$`53`))
  
  ###log2
  metabolomics_expression_data = 
    log(metabolomics_expression_data + 1, 2)
  
  plot(density(metabolomics_expression_data$`52`))
  plot(density(metabolomics_expression_data$`53`))
  
  metabolomics_variable_info$Compound.name
  
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
######steps vs metabolomics
dir.create("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_metabolomics")
setwd("data/24_7_study/wearable_omics_correlation/steps_omics_correlation/steps_metabolomics")

#####---------------------------------------------------------------------------
#####---------------------------------------------------------------------------
###correlation between steps and metabolomicss
#global correlation

dir.create("lagged_correlation")

# lagged_result = vector(mode = "list", length = nrow(metabolomics_expression_data))
# 
# for(i in 1:nrow(metabolomics_expression_data)){
#   cat(i, " ")
#   x = as.numeric(metabolomics_expression_data[i, ])
#   time1 = metabolomics_sample_info$accurate_time
#   y = as.numeric(steps_expression_data[1, ])
#   time2 = steps_sample_info$accurate_time
#   
#   result = try(lagged_correlation(
#     x = x,
#     y = y,
#     time1 = time1,
#     time2 = time2,
#     time_tol = 60/60,
#     step = 5/60
#   ), silent = TRUE)
#   
#   if(class(result) == "try-error"){
#     lagged_result[[i]] = NA
#   }else{
#     lagged_result[[i]] = result  
#   }
#   
#   
# }
# names(lagged_result) = rownames(metabolomics_expression_data)
# save(lagged_result, file = "lagged_correlation/lagged_result")

load("lagged_correlation/lagged_result")

# lagged_cor = 
#   lagged_result %>% 
#   purrr::map(function(x){
#     if(is.na(x)){
#       return(0)
#     }
#     x = x$max_cor
#     if(is.null(x)){
#       x = 0
#     }
#     x
#   }) %>% 
#   unlist()
# 
# global_cor = 
#   lagged_result %>% 
#   purrr::map(function(x){
#     if(is.na(x)){
#       return(0)
#     }
#     x = x$global_cor
#     if(is.null(x)){
#       x = 0
#     }
#     x
#   }) %>% 
#   unlist()
# 
# shift_time = 
#   lagged_result %>% 
#   purrr::map(function(x){
#     if(is.na(x)){
#       return(0)
#     }
#     if(is.null(x)){
#       return(0)
#     }
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
#   metabolomics_variable_info$variable_id
# 
# cor_data =
#   data.frame(wearable = "Step",
#              metabolomics_variable_info,
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
# addWorksheet(wb, sheetName = "Step metabolomics global cor",
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
#                         internal_omics_color = class_color["metabolomics"],
#                         wearable_color = wearable_color["step"],
#                         internal_omics_name = temp$Compound.name[i],
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
#                           internal_omics_color = class_color["metabolomics"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$Compound.name[i],
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
#                           internal_omics_color = class_color["metabolomics"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$Compound.name[i],
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
#                           internal_omics_color = class_color["metabolomics"],
#                           wearable_color = wearable_color["step"],
#                           internal_omics_name = temp$Compound.name[i],
#                           warable_name = "Step",
#                           which = "max",
#                           x_limit = c(1,30),
#                           non_matched_point_size = 3,
#                           wearable_point_size = 3,
#                           internal_omics_point_size = 3,
#                           integrated = TRUE)
# 
#   name = paste("Step vs",temp$Compound.name[i])
#   name = name %>% 
#   stringr::str_replace("\\/", "_")
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
                                                 Compound.name, NA)),
                              check_overlap = TRUE,
                              bg.colour='white',
                              color = "black") +
  base_theme

plot

# ggsave(plot, filename = "global_lagged_correlation.pdf", width = 9, height = 7)


####output the shift time vs lagged cor plot for the important metabolites
important_metabolite =
  cor_data %>%
  dplyr::filter(lagged_cor > quantile(cor_data$lagged_cor[cor_data$lagged_cor > 0], 0.75) |
                  lagged_cor < quantile(cor_data$lagged_cor[cor_data$lagged_cor < 0], 0.25))

dir.create("shift_time_vs_cor")

# for(i in 1:nrow(important_metabolite)){
#   cat(i, "")
#   result = lagged_result[[important_metabolite$variable_id[i]]]
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
#     metabolomics_variable_info$Compound.name[match(names(lagged_result)[i], metabolomics_variable_info$variable_id)]
# 
#   name =
#     name %>%
#     stringr::str_replace("\\/", "_")
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

important_metabolite =
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

# save(important_metabolite, file = "important_metabolite")
load("important_metabolite")


important_metabolite_pos = 
  important_metabolite %>% 
  dplyr::filter(class1 == "positive correlation")

important_metabolite_neg = 
  important_metabolite %>% 
  dplyr::filter(class1 == "negative correlation")

###metpath package
important_metabolite_pos$KEGG.ID
important_metabolite_pos$HMDB.ID

########positive correlation
library(metPath)
library(tidyverse)

####KEGG pathway
data("kegg_hsa_pathway", package = "metPath")
kegg_hsa_pathway
get_pathway_class(kegg_hsa_pathway)

remain_idx =
  kegg_hsa_pathway@pathway_class %>%
  unlist() %>%
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

pathway_database =
  filter_pathway(object = kegg_hsa_pathway, remain_idx = remain_idx)

# kegg_result_pos =
#   enrich_kegg(query_id = unique(important_metabolite_pos$KEGG.ID[!is.na(important_metabolite_pos$KEGG.ID)]),
#               query_type = "compound",
#               id_type = "KEGG",
#               pathway_database = pathway_database,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
# 
# kegg_result_pos
# 
# save(kegg_result_pos, file = "kegg_result_pos")

load("kegg_result_pos")
data("hmdb_pathway", package = "metPath")
hmdb_pathway

remain_idx =
  hmdb_pathway@pathway_class %>%
  unlist() %>%
  stringr::str_detect("primary_pathway") %>%
  # `!`() %>%
  which()

hmdb_pathway =
  filter_pathway(object = hmdb_pathway, remain_idx = remain_idx)

# hmdb_result_pos =
#   enrich_hmdb(query_id = unique(important_metabolite_pos$HMDB.ID[!is.na(important_metabolite_pos$HMDB.ID)]),
#               query_type = "compound",
#               id_type = "HMDB",
#               pathway_database = hmdb_pathway,
#               only_primary_pathway = TRUE,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
# 
# hmdb_result_pos
# 
# save(hmdb_result_pos, file = "hmdb_result_pos")
load("hmdb_result_pos")
plot = 
  enrich_bar_plot(object = kegg_result_pos,
                  x_axis = "p_value",
                  top = 10)

plot
ggsave(plot, filename = "kegg_pos_pathway.pdf", width = 7, height = 7)
enrich_bar_plot(object = hmdb_result_pos,
                x_axis = "p_value",
                top = 10)

# kegg_result_neg =
#   enrich_kegg(query_id = unique(important_metabolite_neg$KEGG.ID[!is.na(important_metabolite_neg$KEGG.ID)]),
#               query_type = "compound",
#               id_type = "KEGG",
#               pathway_database = pathway_database,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
# save(kegg_result_neg, file = "kegg_result_neg")
load("kegg_result_neg")
kegg_result_neg

# hmdb_result_neg =
#   enrich_hmdb(query_id = unique(important_metabolite_neg$HMDB.ID[!is.na(important_metabolite_neg$HMDB.ID)]),
#               query_type = "compound",
#               id_type = "HMDB",
#               pathway_database = hmdb_pathway,
#               only_primary_pathway = TRUE,
#               p_cutoff = 0.05,
#               p_adjust_method = "BH",
#               threads = 3)
# 
# hmdb_result_neg
# 
# save(hmdb_result_neg, file = "hmdb_result_neg")
load("hmdb_result_neg")

plot = 
  enrich_bar_plot(object = kegg_result_neg,
                  top = 10,
                  x_axis = "p_value")
plot
ggsave(plot, filename = "kegg_neg_pathway.pdf", width = 7, height = 7)

enrich_bar_plot(object = hmdb_result_neg,
                top = 10,
                x_axis = "p_value")

###output the pathways
kegg_pos_pathway = 
  kegg_result_pos@result %>% 
  dplyr::arrange(p_value_adjust) %>% 
  dplyr::filter(p_value < 0.05)

kegg_neg_pathway = 
  kegg_result_neg@result %>% 
  dplyr::arrange(p_value_adjust) %>% 
  dplyr::filter(p_value < 0.05)

library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
addWorksheet(wb, sheetName = "Positive correlation", gridLines = TRUE)
addWorksheet(wb, sheetName = "Negative correlation", gridLines = TRUE)
freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) 
freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) 
writeDataTable(wb, sheet = 1, x = kegg_pos_pathway,
               colNames = TRUE, rowNames = FALSE)
writeDataTable(wb, sheet = 2, x = kegg_neg_pathway,
               colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "kegg_pathway_enrichment.xlsx", overwrite = TRUE)







