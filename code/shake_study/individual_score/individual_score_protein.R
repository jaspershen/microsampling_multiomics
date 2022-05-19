##
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

###load data
##metabolomics
load("data/shake_study/metabolomics_data_analysis/data_preparation/expression_data")
load("data/shake_study/metabolomics_data_analysis/data_preparation/sample_info")
load("data/shake_study/metabolomics_data_analysis/data_preparation/variable_info")

load("data/shake_study/metabolomics_data_analysis/metabolites/DEG/anova_marker_name")

expression_data1 <-
  expression_data[anova_marker_name, ]

variable_info1 <-
  variable_info[match(anova_marker_name, variable_info$variable_id), ] %>% 
  dplyr::mutate(mol_class = "metabolite")

masstools::setwd_project()
setwd("data/shake_study/3_omics/individual_scores/protein_score/")

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


####protein score
protein_variable_info <- variable_info %>%
  dplyr::filter(
    mol_name %in% c(
      "D-Alloisoleucine",
      "L-Alanine",
      "L-Isoleucine|L-Leucine",
      "L-Methionine",
      "L-Norvaline",
      "L-Phenylalanine",
      "L-Tryptophan",
      "L-Tyrosine",
      "Phenylalanine"
    )
  )

protein_expression_data <-
  expression_data[match(protein_variable_info$variable_id, rownames(expression_data)), ] %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    x / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

# load("../subject_col")

library(plyr)

plot =
  protein_expression_data %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "sample_id",
                      values_to = "value") %>%
  dplyr::left_join(variable_info, by = "variable_id") %>%
  dplyr::left_join(sample_info, by = "sample_id")  %>%
  dplyr::select(TP, mol_name, subject_id, sample_id, value) %>%
  plyr::dlply(.variables = .(subject_id, mol_name)) %>%
  purrr::map(
    .f = function(x) {
      x$value =
        x$value - x$value[x$TP == 0]
      x
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  ggplot(aes(TP, value, group = subject_id)) +
  geom_hline(yintercept = 0) +
  # geom_point(aes(color = subject_id), show.legend = FALSE) +
  geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  geom_smooth(method = 'loess',
              aes(group = 1),
              se = FALSE,
              color = "black") +
  facet_wrap(facets = vars(mol_name)) +
  scale_color_manual(values = subject_col) +
  theme_bw() +
  labs(x = "Time point (min)", y = "Scaled log2 intensity") +
  scale_x_continuous(breaks = c(0, 30, 60, 120, 240),
                     labels = c(0, 30, 60, 120, 240)) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot
# ggsave(plot,
#        filename = "protein_plot.pdf",
#        width = 7,
#        height = 7)

plot =
  protein_expression_data %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -variable_id,
                      names_to = "sample_id",
                      values_to = "value") %>%
  dplyr::left_join(variable_info, by = "variable_id") %>%
  dplyr::left_join(sample_info, by = "sample_id")  %>%
  dplyr::select(TP, mol_name, subject_id, sample_id, value) %>%
  plyr::dlply(.variables = .(subject_id, mol_name)) %>%
  purrr::map(
    .f = function(x) {
      x$value =
        x$value - x$value[x$TP == 0]
      x
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::mutate(subject_id = factor(subject_id, 
                                    levels = stringr::str_sort(unique(subject_id), numeric = TRUE))) %>% 
  ggplot(aes(TP, value, group = subject_id)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(color = subject_id), show.legend = FALSE) +
  geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  geom_area(aes(group = subject_id, fill = subject_id),
            alpha = 0.7,
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
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )  +
  facet_grid(vars(mol_name), vars(subject_id))
plot
# ggsave(plot, filename = "protein_plot2.pdf", width = 14, height = 7)


# protein_score =
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(temp_subject_id) {
#     cat(temp_subject_id, " ")
#     temp_data =
#       protein_expression_data %>%
#       tibble::rownames_to_column(var = "variable_id") %>%
#       tidyr::pivot_longer(
#         cols = -variable_id,
#         names_to = "sample_id",
#         values_to = "value"
#       ) %>%
#       dplyr::left_join(variable_info, by = "variable_id") %>%
#       dplyr::left_join(sample_info, by = "sample_id")  %>%
#       dplyr::select(TP, mol_name, subject_id, sample_id, value) %>%
#       dplyr::filter(subject_id == temp_subject_id) %>%
#       dplyr::arrange(TP)
# 
#     if (nrow(temp_data) == 0) {
#       return(NULL)
#     }
# 
#     plot =
#       temp_data %>%
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
#            width = 9,
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
# names(protein_score) = unique(sample_info$subject_id)
# 
# protein_score =
#   protein_score %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
# 
# save(protein_score, file = "protein_score")
load("protein_score")

plot =
  protein_score %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "mol_name",
                      values_to = "value") %>%
  dplyr::mutate(subject_id = factor(subject_id,
                                    levels = stringr::str_sort(unique(subject_id), numeric = TRUE))) %>%
  ggplot(aes(subject_id, value)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(fill = mol_name), shape = 21,
             size = 4) +
  geom_line(aes(group = mol_name, color = mol_name)) +
  ggsci::scale_color_uchicago() +
  ggsci::scale_fill_uchicago() +
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
# ggsave(plot, filename = "protein_score_plot.pdf", width = 10, height = 7)

protein_score %>%
  apply(1, mean) %>%
  sort()





protein_score1 = 
  protein_score %>% 
  apply(2, function(x){
    1 - ((x - min(x))/(max(x) - min(x)))
  }) %>% 
  as.data.frame()

temp_data <-
  protein_score1 %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "class",
                      values_to = 'value') %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(
    unique(sample_info$subject_id), numeric = TRUE
  )))

library(plyr)
rsd <- 
  temp_data %>% 
  group_by(subject_id) %>% 
  plyr::dlply(.variables = .(subject_id)) %>% 
  purrr::map(function(x){
    rsd = sd(x$value)*100/mean(x$value)
    data.frame(subject_id = unique(x$subject_id),
               rsd = rsd)
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

library(ggside)

plot <- 
  temp_data %>% 
  dplyr::left_join(rsd, by = "subject_id") %>% 
  ggplot(aes(subject_id, value)) +
  geom_boxplot(aes(color = subject_id),
               show.legend = FALSE) +
  geom_jitter(aes(fill = subject_id), size = 3,
              shape = 21,
              show.legend = FALSE) +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  base_theme +
  labs(x = "", y = "Score") +
  geom_text(aes(x = subject_id, 
                y = 1,label = paste0(round(rsd, 2), "%")),
            angle = 90)

plot  

ggsave(plot, filename = "protein_score_boxplot.pdf", width = 14, height = 3)







