##
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

load("data/shake_study/lipidomics_data_analysis/data_preparation/expression_data")
load("data/shake_study/lipidomics_data_analysis/data_preparation/sample_info")
load("data/shake_study/lipidomics_data_analysis/data_preparation/variable_info")
load("data/shake_study/lipidomics_data_analysis/DEG/anova_marker_name")

load("data/shake_study/subject_info/subject_info")

sample_info <-
  sample_info %>% 
  dplyr::left_join(subject_info, by = "subject_id")

sxtTools::setwd_project()
setwd("data/shake_study/lipidomics_data_analysis/data_overview/")

variable_info <- 
  variable_info
  # dplyr::filter(variable_id %in% anova_marker_name)

expression_data <- 
  expression_data[match(variable_info$variable_id, rownames(expression_data)),]

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
  dplyr::mutate(TP = factor(TP, levels = stringr::str_sort(unique(TP), numeric = TRUE)))

plot <-
  autoplot(
    object = pca_object,
    data = data,
    fill = "TP",
    # frame.colour  = "TP",
    variance_percentage = TRUE,
    size = 4,
    shape = 21,
    alpha = 1,
    frame = TRUE,
    frame.type = 'norm',
  ) +
  geom_vline(xintercept = 0,
             linetype = 2,
             color = "black") +
  geom_hline(yintercept = 0,
             linetype = 2,
             color = "black") +
  # geom_point(size = 4, shape = 16, alpha = 0.7) +
  guides(fill = guide_legend(title = "Time point")) +
  ggrepel::geom_text_repel(aes(label = subject_id), size = 5) +
  ggsci::scale_fill_d3() +
  ggsci::scale_color_d3() +
  theme_bw() +
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
  ) +
  facet_wrap(facets = vars(TP), nrow = 1)

plot

ggsave(filename = "PCA_plot.pdf",
       width = 15,
       height = 4)
