##
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

###load data
load("data/ensure_shake_study/3_omics/individual_variation/cluster")
cluster <-
  data.frame(subject_id = names(cluster),
             cluster)
##metabolomics
load(
  "data/ensure_shake_study/metabolomics_data_analysis/data_preparation/expression_data"
)
load(
  "data/ensure_shake_study/metabolomics_data_analysis/data_preparation/sample_info"
)
load(
  "data/ensure_shake_study/metabolomics_data_analysis/data_preparation/variable_info"
)
load("data/ensure_shake_study/subject_info/subject_info")

load(
  "data/ensure_shake_study/metabolomics_data_analysis/metabolites/DEG/anova_marker_name"
)

expression_data1 <-
  expression_data[anova_marker_name,]

variable_info1 <-
  variable_info[match(anova_marker_name, variable_info$variable_id),] %>%
  dplyr::mutate(mol_class = "metabolite")

##lipidomics
load(
  "data/ensure_shake_study/lipidomics_data_analysis/data_preparation/expression_data"
)
load(
  "data/ensure_shake_study/lipidomics_data_analysis/data_preparation/sample_info"
)
load(
  "data/ensure_shake_study/lipidomics_data_analysis/data_preparation/variable_info"
)

load("data/ensure_shake_study/lipidomics_data_analysis/DEG/anova_marker_name")

expression_data2 <-
  expression_data[anova_marker_name,]

variable_info2 <-
  variable_info[match(anova_marker_name, variable_info$variable_id),] %>%
  dplyr::mutate(mol_class = "lipid")

lipid_info <-
  read.table(
    "data/ensure_shake_study/lipidomics_data_analysis/DEG/Lipomat05.txt",
    header = TRUE,
    sep = "\t"
  )

variable_info2 <-
  variable_info2 %>%
  dplyr::left_join(lipid_info, by = c("mol_name" = "Lipid_Name"))

##cytokine
load(
  "data/ensure_shake_study/cytokine_data_analysis/data_preparation/expression_data"
)
load("data/ensure_shake_study/cytokine_data_analysis/data_preparation/sample_info")
load(
  "data/ensure_shake_study/cytokine_data_analysis/data_preparation/variable_info"
)

load("data/ensure_shake_study/cytokine_data_analysis/DEG/anova_marker_name")

expression_data3 <-
  expression_data[anova_marker_name,]

variable_info3 <-
  variable_info[match(anova_marker_name, variable_info$variable_id),] %>%
  dplyr::mutate(mol_class = "cytokine")

masstools::setwd_project()
setwd("data/ensure_shake_study/3_omics/individual_scores")

intersect_name <-
  Reduce(f = intersect, x = list(
    colnames(expression_data1),
    colnames(expression_data2),
    colnames(expression_data3)
  ))
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
  sample_info[match(colnames(expression_data), sample_info$sample_id), ]

dim(expression_data)
dim(sample_info)
dim(variable_info)

sum(colnames(expression_data) == sample_info$sample_id)

sample_info =
  sample_info %>%
  dplyr::left_join(subject_info, by = "subject_id")

###remove outliers
###color
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

######load all individual score
load('carb_score/carb_score')
load('protein_score/protein_score')
load('fat_score/fat_score')
load('inslulin_secreation_score/inslulin_secreation_score')
load('inslulin_sensitivity_score/inslulin_sensitivity_score')

# load("subject_col")

intersect_name =
  Reduce(f = intersect,
         x = list(
           rownames(carb_score),
           rownames(fat_score),
           rownames(protein_score),
           rownames(inslulin_secreation_score),
           rownames(inslulin_sensitivity_score)
         ))

carb_score1 =
  carb_score[intersect_name, ] %>%
  apply(2, function(x) {
    1 - ((x - min(x)) / (max(x) - min(x)))
  }) %>%
  as.data.frame()

fat_score1 =
  fat_score[intersect_name, ] %>%
  apply(2, function(x) {
    1 - ((x - min(x)) / (max(x) - min(x)))
  }) %>%
  as.data.frame()

protein_score1 =
  protein_score[intersect_name, ] %>%
  apply(2, function(x) {
    1 - ((x - min(x)) / (max(x) - min(x)))
  }) %>%
  as.data.frame()

inslulin_secreation_score1 =
  inslulin_secreation_score[intersect_name, ] %>%
  apply(2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  }) %>%
  as.data.frame()

inslulin_sensitivity_score1 =
  inslulin_sensitivity_score[intersect_name, ] %>%
  apply(2, function(x) {
    1 - ((x - min(x)) / (max(x) - min(x)))
  }) %>%
  as.data.frame()

all_feature_score =
  rbind(
    carb_score1 %>%
      tibble::rownames_to_column(var = "subject_id") %>%
      tidyr::pivot_longer(-subject_id, names_to = "variable_id") %>%
      dplyr::mutate(class = "carb"),
    fat_score1 %>%
      tibble::rownames_to_column(var = "subject_id") %>%
      tidyr::pivot_longer(-subject_id, names_to = "variable_id") %>%
      dplyr::mutate(class = "fat"),
    protein_score1 %>%
      tibble::rownames_to_column(var = "subject_id") %>%
      tidyr::pivot_longer(-subject_id, names_to = "variable_id") %>%
      dplyr::mutate(class = "protein"),
    inslulin_secreation_score1 %>%
      tibble::rownames_to_column(var = "subject_id") %>%
      tidyr::pivot_longer(-subject_id, names_to = "variable_id") %>%
      dplyr::mutate(class = "inslulin_secreation"),
    inslulin_sensitivity_score1 %>%
      tibble::rownames_to_column(var = "subject_id") %>%
      tidyr::pivot_longer(-subject_id, names_to = "variable_id") %>%
      dplyr::mutate(class = "inslulin_sensitivity")
  )

#####how to combine area for each individual score
carb_variable_info <- variable_info %>%
  dplyr::filter(mol_name %in% c("D(-)-Fructose", "L-Lactic acid", "Pyruvic acid"))

carb_expression_data <-
  expression_data[match(carb_variable_info$variable_id, rownames(expression_data)),] %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    x / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

# carb_cor_matrix =
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(x) {
#     temp_sample_info =
#       sample_info %>%
#       dplyr::filter(subject_id == x)
#
#     if (nrow(temp_sample_info) == 0) {
#       return(NULL)
#     }
#
#     library(ComplexHeatmap)
#     temp_data =
#       carb_expression_data[, temp_sample_info$sample_id]
#
#     rownames(temp_data) = carb_variable_info$mol_name
#
#     library(corrplot)
#     library(circlize)
#     col_fun = circlize::colorRamp2(breaks = c(-1, 0, 1),
#                                    colors = c("blue", "white", "red"))
#     temp_data =
#       temp_data %>%
#       t() %>%
#       cor()
#
#     plot =
#       temp_data %>%
#       Heatmap(
#         col = col_fun,
#         rect_gp = gpar(col = "white"),
#         row_names_side = "left",
#         row_dend_side = "right",
#         clustering_method_columns = "ward.D",
#         clustering_method_rows = "ward.D",
#         clustering_distance_rows = "euclidean",
#         clustering_distance_columns = "euclidean",
#         border = TRUE,
#         column_title = paste("Mean correlation:", round(mean(temp_data[upper.tri(temp_data)]), 3))
#       ) %>%
#       ggplotify::as.ggplot()
#
#     ggsave(
#       plot,
#       filename = file.path("carb_score/correltion_plot/",
#                            paste(x, ".pdf", sep = "")),
#       width = 8,
#       height = 7
#     )
#     return(as.data.frame(temp_data))
#   })
#
# names(carb_cor_matrix) = unique(sample_info$subject_id)
# save(carb_cor_matrix, file = "carb_score/carb_cor_matrix")
load("carb_score/carb_cor_matrix")


###fat
#####how to combine area for each individual score
fat_variable_info <- variable_info %>%
  dplyr::filter(stringr::str_detect(mol_name, "TAG"))

fat_expression_data <-
  expression_data[match(fat_variable_info$variable_id, rownames(expression_data)),] %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    x / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

# fat_cor_matrix =
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(x) {
#     temp_sample_info =
#       sample_info %>%
#       dplyr::filter(subject_id == x)
#
#     if (nrow(temp_sample_info) == 0) {
#       return(NULL)
#     }
#
#     library(ComplexHeatmap)
#     temp_data =
#       fat_expression_data[, temp_sample_info$sample_id]
#
#     rownames(temp_data) = fat_variable_info$mol_name
#
#     library(corrplot)
#     library(circlize)
#     col_fun = circlize::colorRamp2(breaks = c(-1, 0, 1),
#                                    colors = c("blue", "white", "red"))
#     temp_data =
#       temp_data %>%
#       t() %>%
#       cor()
#
#     plot =
#       temp_data %>%
#       Heatmap(
#         col = col_fun,
#         # rect_gp = gpar(col = "white"),
#         row_names_side = "left",
#         row_dend_side = "right",
#         clustering_method_columns = "ward.D",
#         clustering_method_rows = "ward.D",
#         clustering_distance_rows = "euclidean",
#         clustering_distance_columns = "euclidean",
#         border = TRUE,
#         column_names_gp = gpar(fontsize = 4),
#         row_names_gp = gpar(fontsize = 4),
#         column_title = paste("Mean correlation:", round(mean(temp_data[upper.tri(temp_data)]), 3))
#       ) %>%
#       ggplotify::as.ggplot()
#
#     ggsave(
#       plot,
#       filename = file.path("fat_score/correltion_plot/",
#                            paste(x, ".pdf", sep = "")),
#       width = 8,
#       height = 7
#     )
#     return(as.data.frame(temp_data))
#   })
#
# names(fat_cor_matrix) = unique(sample_info$subject_id)
# save(fat_cor_matrix, file = "fat_score/fat_cor_matrix")
load("fat_score/fat_cor_matrix")

###protein
#####how to combine area for each individual score
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
  expression_data[match(protein_variable_info$variable_id,
                        rownames(expression_data)),] %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    x / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

# protein_cor_matrix =
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(x) {
#     temp_sample_info =
#       sample_info %>%
#       dplyr::filter(subject_id == x)
#
#     if (nrow(temp_sample_info) == 0) {
#       return(NULL)
#     }
#
#     library(ComplexHeatmap)
#     temp_data =
#       protein_expression_data[, temp_sample_info$sample_id]
#
#     rownames(temp_data) = protein_variable_info$mol_name
#
#     library(corrplot)
#     library(circlize)
#     col_fun = circlize::colorRamp2(breaks = c(-1, 0, 1),
#                                    colors = c("blue", "white", "red"))
#     temp_data =
#       temp_data %>%
#       t() %>%
#       cor()
#
#     plot =
#       temp_data %>%
#       Heatmap(
#         col = col_fun,
#         # rect_gp = gpar(col = "white"),
#         row_names_side = "left",
#         row_dend_side = "right",
#         clustering_method_columns = "ward.D",
#         clustering_method_rows = "ward.D",
#         clustering_distance_rows = "euclidean",
#         clustering_distance_columns = "euclidean",
#         border = TRUE,
#         column_names_gp = gpar(fontsize = 10),
#         row_names_gp = gpar(fontsize = 10),
#         column_title = paste("Mean correlation:", round(mean(temp_data[upper.tri(temp_data)]), 3))
#       ) %>%
#       ggplotify::as.ggplot()
#
#     ggsave(
#       plot,
#       filename = file.path("protein_score/correltion_plot/",
#                            paste(x, ".pdf", sep = "")),
#       width = 8,
#       height = 7
#     )
#     return(as.data.frame(temp_data))
#   })
#
# names(protein_cor_matrix) = unique(sample_info$subject_id)
# save(protein_cor_matrix, file = "protein_score/protein_cor_matrix")
load("protein_score/protein_cor_matrix")


###inslulin_secreation
#####how to combine area for each individual score
inslulin_secreation_variable_info <- variable_info %>%
  dplyr::filter(mol_name %in% c("CPEPTIDE", "INSULIN"))

inslulin_secreation_expression_data <-
  expression_data[match(inslulin_secreation_variable_info$variable_id,
                        rownames(expression_data)),] %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    x / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

# inslulin_secreation_cor_matrix =
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(x) {
#     temp_sample_info =
#       sample_info %>%
#       dplyr::filter(subject_id == x)
#
#     if (nrow(temp_sample_info) == 0) {
#       return(NULL)
#     }
#
#     library(ComplexHeatmap)
#     temp_data =
#       inslulin_secreation_expression_data[, temp_sample_info$sample_id]
#
#     rownames(temp_data) = inslulin_secreation_variable_info$mol_name
#
#     library(corrplot)
#     library(circlize)
#     col_fun = circlize::colorRamp2(breaks = c(-1, 0, 1),
#                                    colors = c("blue", "white", "red"))
#     temp_data =
#       temp_data %>%
#       t() %>%
#       cor()
#
#     plot =
#       temp_data %>%
#       Heatmap(
#         col = col_fun,
#         rect_gp = gpar(col = "white"),
#         row_names_side = "left",
#         row_dend_side = "right",
#         clustering_method_columns = "ward.D",
#         clustering_method_rows = "ward.D",
#         clustering_distance_rows = "euclidean",
#         clustering_distance_columns = "euclidean",
#         border = TRUE,
#         column_names_gp = gpar(fontsize = 10),
#         row_names_gp = gpar(fontsize = 10),
#         column_title = paste("Mean correlation:", round(mean(temp_data[upper.tri(temp_data)]), 3))
#       ) %>%
#       ggplotify::as.ggplot()
#
#     ggsave(
#       plot,
#       filename = file.path("inslulin_secreation_score/correltion_plot/",
#                            paste(x, ".pdf", sep = "")),
#       width = 8,
#       height = 7
#     )
#     return(as.data.frame(temp_data))
#   })
#
# names(inslulin_secreation_cor_matrix) = unique(sample_info$subject_id)
# save(inslulin_secreation_cor_matrix, file = "inslulin_secreation_score/inslulin_secreation_cor_matrix")
load("inslulin_secreation_score/inslulin_secreation_cor_matrix")


###inslulin_sensitivity
#####how to combine area for each individual score
inslulin_sensitivity_variable_info <- variable_info %>%
  dplyr::filter(stringr::str_detect(mol_name, "FFA"))

inslulin_sensitivity_expression_data <-
  expression_data[match(inslulin_sensitivity_variable_info$variable_id,
                        rownames(expression_data)),] %>%
  `+`(1) %>%
  log(2) %>%
  apply(1, function(x) {
    x / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

# inslulin_sensitivity_cor_matrix =
#   unique(sample_info$subject_id) %>%
#   purrr::map(function(x) {
#     temp_sample_info =
#       sample_info %>%
#       dplyr::filter(subject_id == x)
#
#     if (nrow(temp_sample_info) == 0) {
#       return(NULL)
#     }
#
#     library(ComplexHeatmap)
#     temp_data =
#       inslulin_sensitivity_expression_data[, temp_sample_info$sample_id]
#
#     rownames(temp_data) = inslulin_sensitivity_variable_info$mol_name
#
#     library(corrplot)
#     library(circlize)
#     col_fun = circlize::colorRamp2(breaks = c(-1, 0, 1),
#                                    colors = c("blue", "white", "red"))
#     temp_data =
#       temp_data %>%
#       t() %>%
#       cor()
#
#     plot =
#       temp_data %>%
#       Heatmap(
#         col = col_fun,
#         rect_gp = gpar(col = "white"),
#         row_names_side = "left",
#         row_dend_side = "right",
#         clustering_method_columns = "ward.D",
#         clustering_method_rows = "ward.D",
#         clustering_distance_rows = "euclidean",
#         clustering_distance_columns = "euclidean",
#         border = TRUE,
#         column_names_gp = gpar(fontsize = 10),
#         row_names_gp = gpar(fontsize = 10),
#         name = "cor",
#         column_title = paste("Mean correlation:", round(mean(temp_data[upper.tri(temp_data)]), 3))
#       ) %>%
#       ggplotify::as.ggplot()
#
#     ggsave(
#       plot,
#       filename = file.path("inslulin_sensitivity_score/correltion_plot/",
#                            paste(x, ".pdf", sep = "")),
#       width = 8,
#       height = 7
#     )
#     return(as.data.frame(temp_data))
#   })
#
# names(inslulin_sensitivity_cor_matrix) = unique(sample_info$subject_id)
# save(inslulin_sensitivity_cor_matrix, file = "inslulin_sensitivity_score/inslulin_sensitivity_cor_matrix")
load("inslulin_sensitivity_score/inslulin_sensitivity_cor_matrix")

temp_data =
  rbind(
    carb_cor_matrix %>%
      purrr::map(function(x) {
        mean(x[upper.tri(x)])
      }) %>%
      unlist(),
    fat_cor_matrix %>%
      purrr::map(function(x) {
        mean(x[upper.tri(x)])
      }) %>%
      unlist(),
    protein_cor_matrix %>%
      purrr::map(function(x) {
        mean(x[upper.tri(x)])
      }) %>%
      unlist(),
    inslulin_secreation_cor_matrix %>%
      purrr::map(function(x) {
        mean(x[upper.tri(x)])
      }) %>%
      unlist(),
    inslulin_sensitivity_cor_matrix %>%
      purrr::map(function(x) {
        mean(x[upper.tri(x)])
      }) %>%
      unlist()
  ) %>%
  as.data.frame()

rownames(temp_data)  = c("carb",
                         "fat",
                         "protein",
                         "inslulin_secreation",
                         "inslulin_sensitivity")

###combine score together
carb_score2 =
  carb_score1 %>%
  apply(1, median)

fat_score2 =
  fat_score1 %>%
  apply(1, median)

protein_score2 =
  protein_score1 %>%
  apply(1, median)

inslulin_secreation_score2 =
  inslulin_secreation_score1 %>%
  apply(1, median)

inslulin_sensitivity_score2 =
  inslulin_sensitivity_score1 %>%
  apply(1, median)

all_score =
  rbind(
    carb_score2,
    fat_score2,
    protein_score2,
    inslulin_secreation_score2,
    inslulin_sensitivity_score2
  ) %>%
  as.data.frame() %>%
  apply(1, function(x) {
    (x - min(x)) / (max(x) - min(x))
  }) %>%
  t() %>%
  as.data.frame()

rownames(all_score) =
  c("carb",
    "fat",
    "protein",
    "inslulin_secreation",
    "inslulin_sensitivity")

# save(all_score, file = "all_score")
load("all_score")

write.csv(t(all_score), "all_score.csv")

text_df =
  all_score %>%
  tibble::rownames_to_column(var = "class") %>%
  tidyr::pivot_longer(cols = -class,
                      names_to = "subject_id",
                      values_to = "value") %>%
  dplyr::mutate(subject_id = factor(subject_id,
                                    levels = stringr::str_sort(unique(subject_id), numeric = TRUE))) %>%
  dplyr::mutate(class = factor(
    class,
    levels = c(
      "carb",
      "fat",
      "protein",
      "inslulin_secreation",
      "inslulin_sensitivity"
    )
  ))

plot =
  all_feature_score %>%
  dplyr::mutate(subject_id = factor(subject_id,
                                    levels = stringr::str_sort(unique(subject_id), numeric = TRUE))) %>%
  dplyr::mutate(class = factor(
    class,
    levels = c(
      "carb",
      "fat",
      "protein",
      "inslulin_secreation",
      "inslulin_sensitivity"
    )
  )) %>%
  ggplot(aes(subject_id, value)) +
  geom_boxplot(aes(color = subject_id), show.legend = FALSE) +
  geom_jitter(
    shape = 21,
    aes(fill = subject_id),
    alpha = 0.6,
    show.legend = FALSE
  ) +
  geom_line(aes(subject_id, value, group = class),
            data = text_df,
            show.legend = FALSE) +
  geom_point(
    aes(subject_id, value, fill = subject_id),
    shape = 21,
    size = 4,
    data = text_df,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  theme_bw() +
  labs(x = "", y = "Scaled area under curve (AUC)") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  facet_grid(rows = vars(class))

plot

# ggsave(plot, filename = "all_score.pdf", width = 14, height = 10)


library(GGally)

# Check correlations (as scatterplots), distribution and print corrleation coefficient
# load("subject_col")
my_plot <-
  function(data,
           mapping,
           ...,
           point_fill) {
    ggplot(data = data, mapping = mapping) +
      geom_point(..., shape = 21, size = 3) +
      geom_smooth(se = FALSE,
                  color = "black",
                  method = "lm")
    # scale_fill_manual(values = subject_color)
  }

plot =
  ggpairs(
    as.data.frame(t(all_score)) %>% tibble::rownames_to_column(var = "subject_id"),
    columns = c(
      "carb",
      "fat",
      "protein",
      "inslulin_secreation",
      "inslulin_sensitivity"
    ),
    columnLabels = c(
      "Carbohydrate",
      "Fat",
      "Protein",
      "Insulin secreation",
      "Insulin sensitivity"
    ),
    lower = list(continuous = wrap(
      my_plot,
      mapping = aes(fill = subject_id),
      point_fill = subject_col
    )),
    upper = list(continuous = wrap("cor", method = "spearman"))
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

plot
# ggsave(plot, filename = "score_cor_matirx.pdf", width = 9, height = 7)

idx1 = "fat"
idx2 = "inslulin_sensitivity"

temp_data =
  all_score %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  dplyr::select(one_of(c("subject_id", idx1, idx2)))
# dplyr::filter(subject_id != "S13")

cor_value =
  cor.test(as.numeric(temp_data[, 2]),
           as.numeric(temp_data[, 3]),
           method = "spearman")

cor_value =
  paste(
    round(cor_value$estimate, 3),
    case_when(
      cor_value$p.value < 0.001 ~ "***",
      cor_value$p.value > 0.001 &
        cor_value$p.value <= 0.01 ~ "**",
      cor_value$p.value > 0.01 &
        cor_value$p.value <= 0.05 ~ "*",
      cor_value$p.value > 0.05 ~ "NS"
    )
  )


temp_data %>%
  ggplot() +
  geom_point(
    aes_string(idx1, idx2, fill = "subject_id"),
    size = 4,
    shape = 21,
    show.legend = FALSE
  ) +
  geom_smooth(
    aes_string(idx1, idx2),
    se = FALSE,
    method = "lm",
    color = "black"
  ) +
  ggrepel::geom_text_repel(aes_string(idx1, idx2, label = "subject_id")) +
  annotate(
    geom = "text",
    x = -Inf,
    y = -Inf,
    label = cor_value,
    hjust = -1,
    vjust = -1,
    size = 4
  ) +
  scale_fill_manual(values = subject_col) +
  theme_bw() +
  theme(panel.grid = element_blank())

# ggsave(
#   filename =
#     paste(idx1, "_", idx2, ".pdf", sep = ""),
#   width = 7,
#   height = 7
# )

all_score

library(ComplexHeatmap)
library(circlize)
col_fun = circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "red"))

bmi <-
  subject_info[match(colnames(all_score), subject_info$subject_id), ]$bmi
age <-
  subject_info[match(colnames(all_score), subject_info$subject_id), ]$age
sspg <-
  subject_info[match(colnames(all_score), subject_info$subject_id), ]$sspg
sex <-
  subject_info[match(colnames(all_score), subject_info$subject_id), ]$sex
ethnicity <-
  subject_info[match(colnames(all_score), subject_info$subject_id), ]$ethnicity

subject_info <-
  subject_info %>%
  dplyr::left_join(cluster, by = "subject_id") %>%
  dplyr::rename(group = cluster)

group <-
  cluster[match(colnames(all_score), cluster$subject_id), ]$cluster

library(ComplexHeatmap)

ha1 = HeatmapAnnotation(
  gp = gpar(col = "black"),
  sex = factor(sex, levels = c("F", "M")),
  ethnicity = factor(ethnicity, levels = c("A", "B", "C", "H")),
  # group = factor(as.character(group), levels = c("1", "2")),
  col = list(
    sex =
      c("M" = unname(sex_color["M"]),
        "F" = unname(sex_color["F"])),
    
    ethnicity = c(
      "A" =  unname(ethnicity_color["A"]),
      "B" = unname(ethnicity_color["B"]),
      "C" = unname(ethnicity_color["C"]),
      "H" = unname(ethnicity_color["H"])
    )
    # group = c(
    #   "1" =  ggsci::pal_d3()(10)[1],
    #   "2" = ggsci::pal_d3()(10)[2]
    # )
    
  )
)

ha2 =
  HeatmapAnnotation(
    bmi = anno_lines(
      bmi,
      add_points = TRUE,
      ylim = c(min(bmi, na.rm = TRUE), max(bmi, na.rm = TRUE)),
      height = unit(3, "cm"),
      size = unit(5, "mm"),
      pch = 21,
      extend = 0.1,
      pt_gp = gpar(fill = 3, fill = ggsci::pal_aaas()(n = 10)[4])
    ),
    age = anno_lines(
      age,
      add_points = TRUE,
      ylim = c(min(age, na.rm = TRUE), max(age, na.rm = TRUE)),
      height = unit(3, "cm"),
      pch = 21,
      size = unit(5, "mm"),
      extend = 0.1,
      pt_gp = gpar(cex = 3, fill = ggsci::pal_aaas()(n = 10)[5])
    ),
    sspg = anno_lines(
      sspg,
      add_points = TRUE,
      ylim = c(min(sspg, na.rm = TRUE), max(sspg, na.rm = TRUE)),
      height = unit(3, "cm"),
      pch = 21,
      size = unit(5, "mm"),
      extend = 0.1,
      pt_gp = gpar(cex = 3,
                   fill = ggsci::pal_aaas()(n = 10)[6])
    )
  )

###annotation
###top annotation
library(circlize)

plot =
  Heatmap(
    all_score,
    rect_gp = gpar(col = "white"),
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    col = col_fun,
    name = "Individual score",
    border = TRUE,
    top_annotation = ha1,
    bottom_annotation = ha2,
    column_names_rot = 45
  )

library(dendextend)

name1 = colnames(plot@matrix)[column_order(object = plot)]
dend1 = ComplexHeatmap::column_dend(object = plot)
dend1 = dendextend::as.ggdend(dend1)

plot = ggplotify::as.ggplot(plot)

plot

# ggsave(plot, filename = "individual_score_heatmap.pdf", width = 10, height = 7)


###all feature score
temp_data =
  all_feature_score %>%
  tidyr::pivot_wider(names_from = subject_id, values_from = value)

plot =
  Heatmap(
    temp_data[, -c(1, 2)],
    column_names_rot = 45,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "average",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "average",
    col = col_fun,
    name = "Individual score",
    row_split = factor(temp_data$class, levels = c(unique(temp_data$class))),
    border = TRUE,
    show_row_names = FALSE,
    row_title = NULL,
    right_annotation = rowAnnotation(foo = anno_block(
      gp = gpar(fill = rep("grey", 5)),
      labels = unique(temp_data$class),
      labels_gp = gpar(col = "white", fontsize = 10),
      labels_rot = 0
    ))
  )

name2 = colnames(plot@matrix)[column_order(object = plot)]
dend2 = ComplexHeatmap::column_dend(object = plot)
dend2 = dendextend::as.ggdend(dend2)

plot = ggplotify::as.ggplot(plot)

plot

# ggsave(plot, filename = "all_individual_score_heatmap.pdf", width = 10, height = 7)

cbind(name1, name2)

name1 =
  data.frame(
    x = 1:length(name1),
    subject_id = name1,
    class = "all_score",
    stringsAsFactors = FALSE
  )

name2 =
  data.frame(
    x = 1:length(name2),
    subject_id = name2,
    class = "combine_score",
    stringsAsFactors = FALSE
  )

temp_data =
  rbind(name1, name2)

text_data =
  cbind(
    name1 %>%
      dplyr::arrange(subject_id) %>%
      dplyr::rename(x1 = x, y1 = class),
    name2 %>%
      dplyr::arrange(subject_id) %>%
      dplyr::rename(x2 = x, y2 = class) %>%
      dplyr::select(-subject_id)
  )

plot =
  temp_data %>%
  dplyr::mutate(class = factor(class, levels = c("combine_score", "all_score"))) %>%
  ggplot(aes(x, class)) +
  # geom_segment(x = 1, y = "all_score", xend = 4, yend = "combine_score") +
  geom_segment(
    aes(
      x = x1,
      y = y1,
      xend = x2,
      yend = y2,
      color = subject_id,
    ),
    show.legend = FALSE,
    data = text_data
  ) +
  # geom_segment(aes(x, y, xend = x, yend = y)) +
  geom_point(
    aes(fill = subject_id),
    shape = 21,
    size = 10,
    show.legend = FALSE
  ) +
  geom_text(aes(label = subject_id), color = "white") +
  theme_void() +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  theme(plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"))

plot

library(patchwork)

plot1 = ggplot(dend1, labels = FALSE) +
  theme(plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"))
plot2 = ggplot(dend2, labels = FALSE) +
  scale_y_reverse() +
  theme(plot.margin = unit(x = c(0, 0, 0, 0), units = "mm"))

plot =
  plot1 + plot + plot2 + patchwork::plot_layout(nrow = 3, heights = c(1, 5, 1))
plot
# ggsave(plot, filename = "cluster_useing_all_and_combine_score.pdf", height = 7, width = 10)

####barplot to show the five score

plot =
  text_df %>%
  ggplot(aes(class, value)) +
  geom_boxplot(fill = "transparent") +
  geom_line(aes(group = subject_id, color = subject_id), show.legend = FALSE) +
  geom_point(
    aes(fill = subject_id),
    shape = 21,
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  labs(x = "", y = "Individual score") +
  # scale_x_discrete(expand = expansion(mult = c(0.05,0.05))) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
plot
# ggsave(plot, filename = "individual_score_plot.pdf", width = 10, height = 7)


####radar to show the plot for each person
# load("subject_col")
library(ggiraphExtra)

dir.create("radar_plot")
# for(idx in 1:ncol(all_score)){
#   cat(idx, " ")
#   plot =
#     all_score[, idx, drop = FALSE] %>%
#     t() %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column(var = "subject_id") %>%
#     dplyr::mutate(subject_id = factor(subject_id,
#                                       levels = stringr::str_sort(unique(subject_id), numeric = TRUE))) %>%
#     ggRadar(
#       aes(group = subject_id, facet = subject_id),
#       interactive = FALSE,
#       show.legend = FALSE,
#       size = 5,
#       shape = "21",
#       rescale = FALSE,
#       legend.position = "none"
#     ) +
#     theme_bw() +
#     scale_fill_manual(values = subject_col) +
#     scale_color_manual(values = subject_col) +
#     theme(panel.grid.minor = element_blank())
#
#   ggsave(plot, filename = file.path("radar_plot", paste(colnames(all_score)[idx], "_radar.pdf", sep = "")),
#          width = 8, height = 7)
# }


#####score vs characteristics
colnames(all_score)
plot(as.numeric(all_score[1, ]), as.numeric(all_score[4, ]))
cor.test(as.numeric(all_score[1, ]), as.numeric(all_score[4, ]), method = "spearman")

sspg = subject_info$sspg[match(colnames(all_score), subject_info$subject_id)]
age = subject_info$age[match(colnames(all_score), subject_info$subject_id)]
bmi = subject_info$bmi[match(colnames(all_score), subject_info$subject_id)]

temp_data1 =
  all_score %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "score_class",
                      values_to = "score")

temp_data2 =
  subject_info %>%
  dplyr::select(subject_id, sspg, age, bmi) %>%
  tidyr::pivot_longer(cols = -subject_id,
                      names_to = "info_class",
                      values_to = "value")

temp_data =
  temp_data1 %>%
  dplyr::left_join(temp_data2, by = "subject_id") %>%
  dplyr::filter(!is.na(info_class)) %>%
  dplyr::mutate(subject_id = factor(subject_id, levels = stringr::str_sort(unique(subject_id), numeric = TRUE)))

plot =
  temp_data %>%
  # dplyr::filter(subject_id != "S1") %>%
  ggplot(aes(score, value)) +
  geom_point(shape = 21, aes(fill = subject_id), size = 4) +
  geom_smooth(se = FALSE, method = "lm", color = "black") +
  scale_fill_manual(values = subject_col) +
  theme_bw() +
  labs(x = "Individual score", y = "Characteristics score") +
  facet_wrap(vars(info_class, score_class),
             scales = "free_y",
             nrow = 3) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12)
  )
plot
# ggsave(plot, filename = "score_vs_value.pdf", width = 10, height = 7)

purrr::map(
  .x = c(unique(temp_data$score_class)),
  .f = function(x) {
    purrr::map(
      .x = c(unique(temp_data$info_class)),
      .f = function(y) {
        temp =
          dplyr::filter(temp_data, score_class == x &
                          info_class == y)  %>%
          dplyr::filter(!is.na(value)) %>%
          dplyr::select(score, value)
        
        c(
          x,
          y,
          cor.test(temp$score, temp$value, method = "spearman")$estimate,
          cor.test(temp$score, temp$value, method = "spearman")$p.value
        )
      }
    ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
  }
) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::rename(x = V1, y = V2, p = V4) %>%
  dplyr::mutate(rho = as.numeric(rho), p = as.numeric(p)) %>%
  dplyr::arrange(y)


###remove outlier S1
purrr::map(
  .x = c(unique(temp_data$score_class)),
  .f = function(x) {
    purrr::map(
      .x = c(unique(temp_data$info_class)),
      .f = function(y) {
        temp =
          dplyr::filter(temp_data, score_class == x &
                          info_class == y)  %>%
          dplyr::filter(subject_id != "S1") %>%
          dplyr::filter(!is.na(value)) %>%
          dplyr::select(score, value)
        
        c(
          x,
          y,
          cor.test(temp$score, temp$value, method = "spearman")$estimate,
          cor.test(temp$score, temp$value, method = "spearman")$p.value
        )
      }
    ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
  }
) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::rename(x = V1, y = V2, p = V4) %>%
  dplyr::mutate(rho = as.numeric(rho), p = as.numeric(p)) %>%
  dplyr::arrange(y)



idx1 = "inslulin_sensitivity"
idx2 = "age"

plot =
  temp_data %>%
  dplyr::filter(subject_id != "S1") %>%
  dplyr::filter(score_class == idx1 & info_class == idx2) %>%
  ggplot(aes(value, score)) +
  geom_point(
    shape = 21,
    aes(fill = subject_id),
    size = 5,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(mapping = aes(label = subject_id)) +
  geom_smooth(se = FALSE, method = "lm", color = "black") +
  scale_fill_manual(values = subject_col) +
  theme_bw() +
  labs(
    y = paste(idx1, ":Individual score"),
    x = paste(idx2, ":Characteristics score")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    strip.text = element_text(size = 13)
  )
plot

# ggsave(plot, filename = paste(idx1, "_",idx2, ".pdf",sep = ""), width = 7, height = 7)




####PCA and t-sne analysis
library(Rtsne)
##PCA
all_score

temp_data <-
  all_score

temp_data <-
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  })

# #PCA analysis
# pca_object1 <- prcomp(x = temp_data, scale. = FALSE)
# save(pca_object1, file = "pca_object1")
load("pca_object1")

library(ggfortify)
x <- pca_object1$x

x <-
  x[, 1:2] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id")

library(ggforce)

library(wesanderson)

plot =
  ggplot(x, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(
    aes(fill = subject_id),
    shape = 21,
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = subject_col) +
  ggrepel::geom_text_repel(mapping = aes(label = subject_id)) +
  labs(
    x = paste("PC1 (", round(summary(pca_object1)$importance[2, 1] * 100, 2), "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object1)$importance[2, 2] * 100, 2), "%)", sep = "")
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "pca_plot1.pdf", width = 7, height = 7)


###all feature score
temp_data <-
  all_feature_score %>%
  dplyr::select(-class) %>%
  tidyr::pivot_wider(names_from = variable_id, values_from = "value") %>%
  tidyr::drop_na() %>%
  tibble::column_to_rownames(var = "subject_id") %>%
  t() %>%
  as.data.frame()

temp_data <-
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  })

#PCA analysis
# pca_object2 <- prcomp(x = temp_data, scale. = FALSE)
# save(pca_object2, file = "pca_object2")
load("pca_object2")

library(ggfortify)
x <- pca_object2$x

x <-
  x[, 1:2] %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id")

library(ggforce)

# child weight
library(wesanderson)

plot =
  ggplot(x, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(
    aes(fill = subject_id),
    shape = 21,
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = subject_col) +
  ggrepel::geom_text_repel(mapping = aes(label = subject_id)) +
  labs(
    x = paste("PC1 (", round(summary(pca_object2)$importance[2, 1] * 100, 2), "%)", sep = ""),
    y = paste("PC2 (", round(summary(pca_object2)$importance[2, 2] * 100, 2), "%)", sep = "")
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot, filename = "pca_plot2.pdf", width = 7, height = 7)



###radar plot
idx = c("S13", "S36", "S17", "S30", "S8")
plot =
  all_score[, idx, drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  dplyr::mutate(subject_id = factor(subject_id,
                                    levels = stringr::str_sort(unique(subject_id), numeric = TRUE))) %>%
  ggRadar(
    aes(group = subject_id),
    interactive = FALSE,
    show.legend = FALSE,
    size = 5,
    shape = "21",
    rescale = FALSE,
    legend.position = "none"
  ) +
  theme_bw() +
  scale_fill_manual(values = subject_col) +
  scale_color_manual(values = subject_col) +
  theme(panel.grid.minor = element_blank())

plot

# ggsave(
#   plot,
#   filename = file.path("radar_plot",
#                        paste(
#                          paste(idx, collapse = "_"),
#                          "_radar.pdf", sep = ""
#                        )),
#   width = 8,
#   height = 7
# )

mean_score =
  apply(all_score, 1, mean) %>%
  data.frame() %>%
  t() %>%
  as.data.frame()

mean_score$subject_id = paste(idx, collapse = ",")

mean_score =
  mean_score %>%
  separate_rows(subject_id, convert = TRUE) %>%
  dplyr::select(subject_id, everything()) %>%
  dplyr::mutate(class = "mean")


plot =
  all_score[, idx, drop = FALSE] %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "subject_id") %>%
  dplyr::mutate(class = "individual") %>%
  rbind(mean_score, .) %>%
  dplyr::mutate(subject_id = factor(subject_id,
                                    levels = name1$subject_id)) %>%
  ggRadar(
    aes(
      group = class,
      facet = subject_id,
      size = class,
      shape = class
    ),
    interactive = FALSE,
    show.legend = FALSE,
    size = 5,
    shape = 21,
    rescale = FALSE,
    legend.position = "none"
  ) +
  theme_bw() +
  # scale_size_manual(values = c(mean = 1, individual = 5)) +
  scale_fill_manual(values = c(
    "mean" = "grey",
    "individual" = ggsci::pal_aaas()(n = 10)[2]
  )) +
  scale_color_manual(values = c(
    "mean" = "grey",
    "individual" = ggsci::pal_aaas()(n = 10)[2]
  )) +
  theme(panel.grid.minor = element_blank())
plot

# ggsave(
#   plot,
#   filename = file.path("radar_plot",
#                        paste(paste(idx, collapse = "_"),
#                              "_radar2.pdf", sep = "")),
#   width = 10,
#   height = 7
# )
