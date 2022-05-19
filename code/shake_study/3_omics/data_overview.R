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

metabolite_num = sum(!is.na(variable_info$Metabolite))

expression_data1 <-
  expression_data[anova_marker_name, ]

variable_info1 <-
  variable_info[match(anova_marker_name, variable_info$variable_id), ]

##lipidomics
load("data/shake_study/lipidomics_data_analysis/data_preparation/expression_data")
load("data/shake_study/lipidomics_data_analysis/data_preparation/sample_info")
load("data/shake_study/lipidomics_data_analysis/data_preparation/variable_info")

load("data/shake_study/lipidomics_data_analysis/DEG/anova_marker_name")

lipid_num = nrow(variable_info)

expression_data2 <-
  expression_data[anova_marker_name, ]

variable_info2 <-
  variable_info[match(anova_marker_name, variable_info$variable_id), ]

##cytokine
load("data/shake_study/cytokine_data_analysis/data_preparation/expression_data")
load("data/shake_study/cytokine_data_analysis/data_preparation/sample_info")
load("data/shake_study/cytokine_data_analysis/data_preparation/variable_info")

load("data/shake_study/cytokine_data_analysis/DEG/anova_marker_name")

cytokine_num = nrow(variable_info)

expression_data3 <-
  expression_data[anova_marker_name, ]

variable_info3 <-
  variable_info[match(anova_marker_name, variable_info$variable_id), ]

masstools::setwd_project()
setwd("data/shake_study/3_omics/data_overview")

value <- 
  c("Lipid" = unname(class_color["lipidomics"]),
    "Metabolite" = unname(class_color["metabolomics"]),
    "Cytokine" = unname(class_color["cytokine"])
  )

temp_data = 
  data.frame(class = c("Metabolite", "Lipid", "Cytokine"),
             number = c(metabolite_num, lipid_num, cytokine_num))

plot = 
temp_data %>% 
  dplyr::mutate(class = factor(class, levels = class)) %>% 
  ggplot(aes(class, number)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_bar(stat = "identity", aes(fill = class), show.legend = FALSE) +
  theme_classic() +
  scale_fill_manual(values = value) +
  labs(x = "", y = "Variable number") +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.ticks.x = element_blank())
plot
# ggsave(plot, file = "variable_number.pdf", width = 7, height = 7)

intersect_name <- 
  Reduce(f = intersect, x = list(colnames(expression_data1),
                                 colnames(expression_data2),
                                 colnames(expression_data3)))

expression_data <-
  rbind(expression_data1[, intersect_name],
        expression_data2[, intersect_name],
        expression_data3[, intersect_name])

variable_info <- 
  variable_info1 %>%
  dplyr::full_join(variable_info2, 
                   by = intersect(colnames(variable_info1), 
                                  colnames(variable_info2))) %>% 
  dplyr::full_join(variable_info3, 
                   by = intersect(colnames(.), 
                                  colnames(variable_info3)))

variable_info$mol_name[!is.na(variable_info$Metabolite)] <- 
  variable_info$Metabolite[!is.na(variable_info$Metabolite)]

sample_info <- 
  sample_info[match(colnames(expression_data), sample_info$sample_id),]

dim(expression_data)
dim(sample_info)
dim(variable_info)

###PCA analysis
##log transformation
subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

subject_data <-
  log(subject_data + 1, 2)

subject_data <- 
  apply(subject_data, 1, function(x){
    (x  - mean(x)) / sd(x)
  })

pca_object <-
  prcomp(x = subject_data, center = FALSE, scale. = FALSE)

library(ggfortify)

data <- cbind(pca_object$x[, 1:2],
              sample_info) %>%
  as.data.frame() %>%
  dplyr::mutate(TP = factor(TP, levels = stringr::str_sort(unique(TP), numeric = TRUE))) %>% 
  dplyr::mutate(subject_id = factor(subject_id, 
                                    stringr::str_sort(unique(sample_info$subject_id), numeric = TRUE)))

plot <-
  autoplot(
    object = pca_object,
    data = data,
    fill = "TP",
    frame.colour  = "TP",
    variance_percentage = TRUE,
    size = 4,
    shape = 21,
    alpha = 1,
    frame = FALSE,
    frame.type = 'norm',
  ) +
  geom_path(aes(PC1, PC2),
            arrow = grid::arrow(
              type = "closed",
              angle = 30,
              length = unit(0.1, "inches")
            )) +
  guides(fill = guide_legend(title = "Time point")) +
  ggrepel::geom_text_repel(aes(label = TP), size = 3) +
  scale_fill_manual(values = tp_color) +
  scale_color_manual(values = tp_color) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 12)
  ) +
  facet_wrap(facets = vars(subject_id), scales = "free")

plot

# ggsave(filename = "PCA_plot_individaul.pdf",
#        width = 14,
#        height = 10)


# library(ComplexHeatmap)
# 
# Heatmap(subject_data,
#         show_column_names = FALSE,
#         row_names_gp = gpar(fontsize = 5),
#         clustering_distance_rows = "euclidean",
#         clustering_method_rows = "ward.D",
#         clustering_distance_columns = "euclidean",
#         clustering_method_columns = "ward.D"
#         )

###tsne
tsne_object =
  Rtsne::Rtsne(
    X = subject_data,
    pca = TRUE,
    pca_center = FALSE,
    pca_scale = FALSE,
    perplexity = 30,
    theta = 0.0
  )

temp_data <- cbind(tsne_object$Y,
              sample_info) %>%
  as.data.frame() %>%
  dplyr::mutate(TP = factor(TP, levels = stringr::str_sort(unique(TP), numeric = TRUE)))

colnames(temp_data)[1:2] = c("x", "y")

library(ggforce)

temp_data2 <-
  temp_data %>% 
  group_by(subject_id) %>% 
  dplyr::summarise(x = mean(x),
            y = mean(y)) %>% 
  ungroup()

plot <-
  temp_data %>%
  ggplot(aes(x, y)) +
  geom_point(
    aes(fill = subject_id),
    shape = 21,
    size = 7,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  guides(fill = guide_legend(title = "Time point")) +
  ggrepel::geom_label_repel(aes(x = x,
                                y = y,
                                label = subject_id,
                                color = subject_id),
                            size = 5,
                            show.legend = FALSE,
                           data = temp_data2) +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  theme_bw() +
  labs(x = "t-SNE1", y = "t-SNE2") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13)
  ) 
# facet_wrap(facets = vars(TP), nrow = 1)

plot

# ggsave(plot, filename = "tsne.pdf", width = 7, height = 7)

temp_data$subject_id <-
  factor(temp_data$subject_id, 
         levels = stringr::str_sort(unique(temp_data$subject_id), numeric = TRUE))

plot <-
  temp_data %>%
  dplyr::arrange(subject_id, TP) %>% 
  ggplot(aes(x, y)) +
  geom_point(
    aes(fill = TP),
    shape = 21,
    size = 5,
    show.legend = FALSE
  ) +
  geom_path(aes(x,y), arrow = grid::arrow(type = "closed", 
                                          angle = 30, 
                                          length = unit(0.1, "inches"))) +
  guides(fill = guide_legend(title = "Time point")) +
  ggrepel::geom_text_repel(aes(label = TP),
                            size = 3,
                            show.legend = FALSE) +
  theme_bw() +
  scale_fill_manual(values = tp_color) +
  scale_x_continuous(expand = expansion(mult = c(0.2,0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2,0.2))) +
  labs(x = "t-SNE1", y = "t-SNE2") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 12)
  ) +
facet_wrap(facets = vars(subject_id), scales = "free")

plot

# ggsave(plot, filename = "tsne_individual.pdf", width = 14, height = 10)


# Show the objects in the 2D tsne representation
plot(tsne_out$Y, col = iris_unique$Species, asp = 1)


####
value <- 
  c("Lipid" = ggsci::pal_aaas()(10)[1],
    "Metabolite" = ggsci::pal_aaas()(10)[3],
    "Cytokine" = ggsci::pal_aaas()(10)[4]
  )

plot <- 
  expression_data %>% 
  dplyr::select(contains("T0")) %>% 
  `+`(1) %>% 
  log(2) %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "variable_id") %>% 
  tidyr::pivot_longer(cols = -variable_id, 
                      names_to = "subject_id", 
                      values_to = "value") %>% 
  dplyr::mutate(subject_id = stringr::str_replace_all(subject_id, "_T0", "")) %>% 
  dplyr::mutate(class = case_when(
    variable_id %in% variable_info1$variable_id ~ "Metabolite",
    variable_id %in% variable_info2$variable_id ~ "Lipid",
    variable_id %in% variable_info3$variable_id ~ "Cytokine"
  )) %>% 
  ggplot(aes(x = subject_id, value)) +
  geom_jitter(aes(subject_id, value, color = class, shape = class),
              size = 1.8, alpha = 0.7) +
  geom_boxplot(fill = "transparent", color = "red", outlier.shape = NA) +
  scale_color_manual(values = value) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_bw() +
  labs(x = "", y = "Scaled intensity") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
        )

plot

# ggsave(plot, filename = "t0_box_plot.pdf", width = 15, height = 4)