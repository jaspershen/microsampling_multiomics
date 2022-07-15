##
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

###load data
##cytokine
load(
  "data/ensure_shake_study/cytokine_data_analysis/data_preparation/expression_data"
)
load("data/ensure_shake_study/cytokine_data_analysis/data_preparation/sample_info")
load(
  "data/ensure_shake_study/cytokine_data_analysis/data_preparation/variable_info"
)

load("data/ensure_shake_study/cytokine_data_analysis/DEG/anova_marker_name")

masstools::setwd_project()
dir.create("data/ensure_shake_study/3_omics/individual_scores/cytokine_score/")
setwd("data/ensure_shake_study/3_omics/individual_scores/cytokine_score/")

variable_info$mol_name[!is.na(variable_info$Metabolite)] <-
  variable_info$Metabolite[!is.na(variable_info$Metabolite)]

dim(expression_data)
dim(sample_info)
dim(variable_info)

###remove outliers
value <-
  c(
    "Lipid" = ggsci::pal_aaas()(10)[1],
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

###cytokine
cytokine_variable_info <- variable_info %>%
  dplyr::filter(subclass == "H41")

cytokine_expression_data <-
  expression_data[match(cytokine_variable_info$variable_id,
                        rownames(expression_data)),] %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    x / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

# load("../subject_col")

library(plyr)

plot <-
  cytokine_expression_data[names(tail(sort(apply(
    cytokine_expression_data, 1, mean
  )), 6)), ] %>%
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
# ggsave(plot, filename = "cytokine_plot.pdf", width = 14, height = 7)

plot <-
  cytokine_expression_data[names(tail(sort(apply(
    cytokine_expression_data, 1, mean
  )), 6)), ] %>%
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
                                    levels = stringr::str_sort(
                                      unique(sample_info$subject_id, numeric = TRUE)
                                    ))) %>%
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
    axis.text.x = element_text(
      angle = 45,
      size = 10,
      hjust = 1,
      vjust = 1
    ),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  )  +
  facet_grid(vars(mol_name), vars(subject_id))
plot

# ggsave(plot, filename = "cytokine_plot2.pdf", width = 14, height = 7)

# cytokine_score <-
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(temp_subject_id) {
#     cat(temp_subject_id, " ")
#     temp_data <-
#       cytokine_expression_data %>%
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
#     plot <-
#       temp_data %>%
#       dplyr::filter(variable_id %in% names(tail(sort(
#         apply(cytokine_expression_data, 1, mean)
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
#            width = 9,
#            height = 7)
#
#     area <-
#       temp_data %>%
#       plyr::dlply(.variables = .(mol_name)) %>%
#       purrr::map(function(x) {
#         x <-
#           x %>%
#           dplyr::arrange(TP)
#         pracma::trapz(x = x$TP,
#               y = x$value - x$value[x$TP == 0])
#       }) %>%
#       unlist()
#     area
#   })
#
# names(cytokine_score) = unique(sample_info$subject_id)
#
# cytokine_score =
#   cytokine_score %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# save(cytokine_score, file = "cytokine_score")
load("cytokine_score")

temp_data <-
  seq_len(nrow(cytokine_score)) %>%
  purrr::map(function(i) {
    x <- cytokine_score[i, , drop = TRUE] %>%
      unlist()
    mean_value <- x %>% as.numeric() %>% mean()
    data.frame(dist = rev(sort(abs(x - mean_value))),
               rank = 1:length(x)) %>%
      tibble::rownames_to_column(var = "name")
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::arrange(name)


library(plyr)

temp_data %>%
  plyr::dlply(.variables = .(name)) %>%
  lapply(function(x) {
    sd(x$rank) / mean(x$rank)
  }) %>%
  unlist() %>%
  data.frame(rsd = .) %>%
  tibble::rownames_to_column(var = "name")

temp_data %>%
  ggplot(aes(name, rank)) +
  geom_boxplot() +
  geom_jitter() +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

####FGFB, IL7, MDC, RANTES
plot <-
  cytokine_score %>%
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
  # ggsci::scale_color_uchicago() +
  # ggsci::scale_fill_uchicago() +
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
# ggsave(plot, filename = "cytokine_score_plot.pdf", width = 10, height = 7)

cytokine_score %>%
  apply(1, mean) %>%
  sort()

cytokine_score1 <-
  cytokine_score %>%
  apply(2, function(x) {
    ((x - min(x)) / (max(x) - min(x)))
  }) %>%
  as.data.frame()

temp_data <-
  cytokine_score1 %>%
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
  purrr::map(function(x) {
    rsd = sd(x$value) * 100 / mean(x$value)
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
  geom_jitter(
    aes(fill = subject_id),
    size = 3,
    shape = 21,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  base_theme +
  labs(x = "", y = "Score") +
  geom_text(aes(
    x = subject_id,
    y = 1,
    label = paste0(round(rsd, 2), "%")
  ),
  angle = 90)

plot

# ggsave(plot,
#        filename = "cytokine_score_boxplot.pdf",
#        width = 14,
#        height = 3)
