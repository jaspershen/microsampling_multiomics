##
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

###load data
##lipidomics
load("data/shake_study/lipidomics_data_analysis/data_preparation/expression_data")
load("data/shake_study/lipidomics_data_analysis/data_preparation/sample_info")
load("data/shake_study/lipidomics_data_analysis/data_preparation/variable_info")

load("data/shake_study/lipidomics_data_analysis/DEG/anova_marker_name")

expression_data <-
  expression_data[anova_marker_name, ] 

variable_info <-
  variable_info[match(anova_marker_name, variable_info$variable_id), ] %>% 
  dplyr::mutate(mol_class = "lipid")

lipid_info <- read.table("data/shake_study/lipidomics_data_analysis/DEG/Lipomat05.txt", header = TRUE, sep = "\t")

variable_info <-
  variable_info %>%
  dplyr::left_join(lipid_info, by = c("mol_name" = "Lipid_Name"))

sxtTools::setwd_project()
setwd("data/shake_study/3_omics/individual_scores/fat_score/")

variable_info$mol_name[!is.na(variable_info$Metabolite)] <- 
  variable_info$Metabolite[!is.na(variable_info$Metabolite)]
  
dim(expression_data)
dim(sample_info)
dim(variable_info)

###remove outliers
value <- 
  c("Lipid" = ggsci::pal_aaas()(10)[1],
    "Metabolite" = ggsci::pal_aaas()(10)[3],
    "Cytokine" = ggsci::pal_aaas()(10)[4]
  )

##remove some outliers
sample_info <-
  sample_info %>% 
  dplyr::filter(!subject_id %in% c("S18"))

expression_data <-
  expression_data %>% 
  dplyr::select(sample_info$sample_id)

####fat score
fat_variable_info = variable_info %>%
  dplyr::filter(stringr::str_detect(mol_name, "TAG")) 

fat_expression_data <-
  expression_data[match(fat_variable_info$variable_id, rownames(expression_data)),] %>% 
  `+`(1) %>% 
  log(2) %>% 
  apply(1, function(x){
    x / sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

load("../subject_col")

library(plyr)

plot = 
  fat_expression_data[names(tail(sort(apply(fat_expression_data, 1, mean)), 6)),] %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "sample_id",
                      values_to = "value") %>% 
  dplyr::left_join(variable_info, by = "variable_id") %>% 
  dplyr::left_join(sample_info, by = "sample_id")  %>% 
  dplyr::select(TP, mol_name, subject_id, sample_id, value) %>% 
  plyr::dlply(.variables = .(subject_id, mol_name)) %>% 
  purrr::map(.f = function(x){
    x$value = 
      x$value - x$value[x$TP == 0]
    x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  ggplot(aes(TP, value, group = subject_id)) +
  geom_hline(yintercept = 0) +
  # geom_point(aes(color = subject_id), show.legend = FALSE) +
  geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  geom_smooth(method = 'loess', aes(group = 1), se = FALSE, color = "black") +
  scale_color_manual(values = subject_col) +
  theme_bw() +
  labs(x = "Time point (min)", y = "Scaled log2 intensity") +
  scale_x_continuous(breaks = c(0, 30, 60, 120, 240), 
                     labels = c(0, 30, 60, 120, 240)) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  facet_wrap(facets = vars(mol_name), scales = "free_y")
plot
# ggsave(plot, filename = "fat_plot.pdf", width = 14, height = 7)

plot =
  fat_expression_data[names(head(sort(apply(fat_expression_data, 1, mean)), 6)),] %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "sample_id",
                      values_to = "value") %>% 
  dplyr::left_join(variable_info, by = "variable_id") %>% 
  dplyr::left_join(sample_info, by = "sample_id")  %>% 
  dplyr::select(TP, mol_name, subject_id, sample_id, value) %>% 
  plyr::dlply(.variables = .(subject_id, mol_name)) %>% 
  purrr::map(.f = function(x){
    x$value = 
      x$value - x$value[x$TP == 0]
    x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  ggplot(aes(TP, value, group = subject_id)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(color = subject_id), show.legend = FALSE) +
  geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  geom_area(aes(group = subject_id, fill = subject_id), alpha = 0.7, 
            show.legend = FALSE) +
  scale_color_manual(values = subject_col) +
  scale_fill_manual(values = subject_col) +
  theme_bw() +
  labs(x = "Time point (min)", y = "Scaled log2 intensity") +
  scale_x_continuous(breaks = c(0, 30, 60, 120, 240), 
                     labels = c(0, 30, 60, 120, 240)) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )  +
  facet_grid(vars(mol_name), vars(subject_id), scales = "free_y")
plot
# ggsave(plot, filename = "fat_plot2.pdf", width = 21, height = 10)

# 
# ###calculate fat score for each person
# fat_score =
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(temp_subject_id) {
#     cat(temp_subject_id, " ")
#     temp_data =
#       fat_expression_data %>%
#       tibble::rownames_to_column(var = "variable_id") %>%
#       tidyr::pivot_longer(
#         cols = -variable_id,
#         names_to = "sample_id",
#         values_to = "value"
#       ) %>%
#       dplyr::left_join(variable_info, by = "variable_id") %>%
#       dplyr::left_join(sample_info, by = "sample_id")  %>%
#       dplyr::select(TP, mol_name, subject_id, sample_id, value, variable_id) %>%
#       dplyr::filter(subject_id == temp_subject_id) %>%
#       dplyr::arrange(TP)
#     
#     if (nrow(temp_data) == 0) {
#       return(NULL)
#     }
#     
#     plot =
#       temp_data %>%
#       dplyr::filter(variable_id %in% names(tail(sort(
#         apply(fat_expression_data, 1, mean)
#       ), 6))) %>%
#       plyr::dlply(.variables = .(mol_name)) %>%
#       purrr::map(
#         .f = function(x) {
#           x =
#             x %>%
#             dplyr::mutate(value = value - x$value[x$TP == 0])
#           x
#         }
#       ) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       ggplot(aes(TP, value, group = subject_id)) +
#       geom_hline(yintercept = 0) +
#       geom_point(
#         aes(fill = subject_id),
#         shape = 21,
#         size = 3,
#         show.legend = FALSE
#       ) +
#       geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
#       geom_area(
#         aes(group = subject_id, fill = subject_id),
#         alpha = 0.5,
#         position = "stack",
#         show.legend = FALSE
#       ) +
#       scale_color_manual(values = subject_col) +
#       scale_fill_manual(values = subject_col) +
#       theme_bw() +
#       labs(x = "Time point (min)", y = "Scaled log2 intensity") +
#       scale_x_continuous(breaks = c(0, 30, 60, 120, 240),
#                          labels = c(0, 30, 60, 120, 240)) +
#       theme(
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 12),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 12),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA)
#       )  +
#       facet_wrap(facets = vars(mol_name))
#     
#     ggsave(plot,
#            file = file.path(paste(temp_subject_id, ".pdf", sep = "")),
#            width = 14,
#            height = 7)
#     
#     area =
#       temp_data %>%
#       plyr::dlply(.variables = .(mol_name)) %>%
#       purrr::map(function(x) {
#         x <-
#           x %>%
#           dplyr::arrange(TP)
#         trapz(x = x$TP,
#               y = x$value - x$value[x$TP == 0])
#       }) %>%
#       unlist()
#     area
#   })
# 
# names(fat_score) = unique(sample_info$subject_id)
# 
# fat_score =
#   fat_score %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# save(fat_score, file = "fat_score")
load("fat_score")

plot = 
  fat_score %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "mol_name",
                      values_to = "value") %>%
  dplyr::mutate(subject_id = factor(subject_id, 
                                    levels = stringr::str_sort(unique(subject_id), numeric = TRUE))) %>% 
  ggplot(aes(subject_id, value)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(fill = mol_name), shape = 21,
             size = 4, show.legend = FALSE) +
  geom_line(aes(group = mol_name, color = mol_name), show.legend = FALSE) +
  # ggsci::scale_color_jama() +
  # ggsci::scale_fill_jama() +
  theme_bw() +
  labs(x = "", y = "Area under curve") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) 
plot
# ggsave(plot, filename = "fat_score_plot.pdf", width = 10, height = 7)

fat_score %>% 
  apply(1, mean) %>% 
  sort()

