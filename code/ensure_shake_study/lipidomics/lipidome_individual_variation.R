##
no_function()

masstools::setwd_project()
library(tidyverse)
rm(list = ls())

source("code/tools.R")

###load data
load(
  "data/ensure_shake_study/lipidomics_data_analysis/data_preparation/expression_data"
)
load(
  "data/ensure_shake_study/lipidomics_data_analysis/data_preparation/sample_info"
)
load(
  "data/ensure_shake_study/lipidomics_data_analysis/data_preparation/variable_info"
)

load("data/ensure_shake_study/subject_info/subject_info")

sample_info <-
  sample_info %>%
  dplyr::left_join(subject_info, by = "subject_id")

masstools::setwd_project()
setwd("data/ensure_shake_study/lipidomics_data_analysis/invidual_variation")

dim(expression_data)
dim(sample_info)
dim(variable_info)

dim(expression_data)
dim(sample_info)
dim(variable_info)

###calculate distance for each person
##for lipidomics, we use the euclidean distance
plot(density(as.numeric(expression_data[1, ])))

expression_data1 <- log(expression_data + 1, 2)
plot(density(as.numeric(expression_data1[1, ])))

distance_individual <-
  purrr::map(
    unique(sample_info$subject_id),
    .f = function(x) {
      temp_sample_info <- sample_info %>%
        dplyr::filter(subject_id == x) %>%
        dplyr::arrange(TP)
      
      temp_data <- expression_data1 %>%
        dplyr::select(temp_sample_info$sample_id)
      
      tp_list <-
        c(0, 30, 60, 120, 240)
      
      distance <-
        purrr::map(
          2:length(tp_list),
          .f = function(idx) {
            temp_tp <- tp_list[c(1, idx)]
            temp_idx <-
              which(temp_sample_info$TP %in% temp_tp)
            if (length(temp_idx) != 2) {
              return(NA)
            } else{
              dist(x = t(as.matrix(temp_data[, temp_idx])), method = "euclidean") %>%
                as.numeric()
            }
          }
        ) %>%
        unlist()
      
      distance <- c(0, distance)
      
    }
  )

names(distance_individual) <- unique(sample_info$subject_id)

distance_individual

# save(distance_individual, file = "distance_individual")
load("distance_individual")

tp_list <- c(0, 30, 60, 120, 240)

distance_individual_old <-
  distance_individual %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(distance_individual_old) <-
  c("0", "30", "60", "120", "240")

which(is.na(distance_individual_old), arr.ind = TRUE)

distance_individual_old[is.na(distance_individual_old)] <- 0

###consensus clustering
library(CancerSubtypes)

# result <-
#   ExecuteCC(
#     clusterNum = 2,
#     d = as.matrix(t(distance_individual_old)),
#     maxK = 6,
#     reps = 1000,
#     pItem = 0.8,
#     pFeature = 0.8,
#     title = "k_means_consensus",
#     clusterAlg = "km",
#     distance = "euclidean",
#     plot = "png",
#     writeTable = TRUE
#   )
#
# save(result, file = "result")
load("result")
idx <- 2
sil = silhouette_SimilarityMatrix(result$originalResult[[idx]]$consensusClass,
                                  result$originalResult[[idx]]$consensusMatrix)
sil_plot <-
  plot_silhouette(sil)

sil_plot

# ggsave(
#   sil_plot,
#   file = paste("sil_plot", idx, ".pdf", sep = ""),
#   width = 7,
#   height = 7
# )

name2 <-
  rownames(distance_individual_old)[result$originalResult[[idx]]$consensusTree$order]

distance_individual_old2 <- distance_individual_old[name2, ]
cluster <- result$originalResult[[idx]]$consensusClass

temp_sample_info <-
  subject_info[match(rownames(distance_individual_old2), subject_info$subject_id), ]

cbind(temp_sample_info$subject_id,
      rownames(distance_individual_old2))

temp_sample_info$subject_id <- rownames(distance_individual_old2)

cluster <-
  cluster[match(temp_sample_info$subject_id, names(cluster))]

names(cluster) == rownames(distance_individual_old2)

##reorder cluster
cluster2 <- cluster[cluster == 2]
cluster1 <- cluster[cluster == 1]
# cluster3 <- cluster[cluster == 3]
# cluster4 <- cluster[cluster == 4]
# cluster <- c(cluster3, cluster1, cluster2)
# # cluster <- rev(cluster)
# temp_sample_info <-
#   temp_sample_info[match(names(cluster), temp_sample_info$sample_id),]
#
# temp_subject_data2 <-
#   temp_subject_data2[,names(cluster)]
#
# names(cluster) == temp_sample_info$sample_id
# names(cluster) == colnames(temp_subject_data2)

###complex heatamp
library(ComplexHeatmap)

##bmi
bmi <- temp_sample_info$bmi

##age
age <- temp_sample_info$age

##sspg
sspg <- temp_sample_info$sspg

##sex
sex <- temp_sample_info$sex

##ethnicity
ethnicity <- temp_sample_info$ethnicity

library(ComplexHeatmap)

###annotation
###top annotation
library(circlize)
col_fun = colorRamp2(c(
  range(distance_individual_old)[1],
  range(distance_individual_old)[2] / 2,
  range(distance_individual_old)[2]
),
c("green", "white", "red"))

cluster_col_fun = colorRamp2(c(0, 2, 3), c("green", "white", "red"))

ha1 = rowAnnotation(
  bmi = bmi,
  age = age,
  sspg = sspg,
  sex = factor(sex, levels = c("F", "M")),
  ethnicity = factor(ethnicity, levels = c("A", "B", "C", "H")),
  col = list(
    bmi = circlize::colorRamp2(
      breaks = c(
        min(bmi, na.rm = TRUE),
        mean(bmi, na.rm = TRUE),
        max(bmi, na.rm = TRUE)
      ),
      colors = c(
        alpha(ggsci::pal_aaas()(10)[1], 1),
        "white",
        alpha(ggsci::pal_aaas()(10)[2], 1)
      )
    ),
    
    age = circlize::colorRamp2(
      breaks = c(
        min(age, na.rm = TRUE),
        mean(age, na.rm = TRUE),
        max(age, na.rm = TRUE)
      ),
      colors = c(
        alpha(ggsci::pal_aaas()(10)[1], 1),
        "white",
        alpha(ggsci::pal_aaas()(10)[2], 1)
      )
    ),
    
    sspg = circlize::colorRamp2(
      breaks = c(
        min(sspg, na.rm = TRUE),
        mean(sspg, na.rm = TRUE) ,
        max(sspg, na.rm = TRUE)
      ),
      colors = c(
        alpha(ggsci::pal_aaas()(10)[1], 1),
        "white",
        alpha(ggsci::pal_aaas()(10)[2], 1)
      )
    ),
    sex =
      c(
        "M" = ggsci::pal_aaas()(n = 10)[2],
        "F" = ggsci::pal_aaas()(n = 10)[10]
      ),
    
    ethnicity = c(
      "A" = ggsci::pal_aaas()(n = 10)[3],
      "B" = ggsci::pal_aaas()(n = 10)[2],
      "C" = ggsci::pal_d3()(10)[2],
      "H" = "snow1"
    )
  )
)

ha2 = rowAnnotation(cluster = as.character(cluster),
                    
                    col = list(
                      cluster = c(
                        "1" = ggsci::pal_d3()(10)[1],
                        "2" = ggsci::pal_d3()(10)[2],
                        "3" = ggsci::pal_d3()(10)[3],
                        "4" = ggsci::pal_d3()(10)[4],
                        "5" = ggsci::pal_d3()(10)[5]
                      )
                    ))


temp_data <- distance_individual_old2

plot <-
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    col = col_fun,
    left_annotation = ha1,
    right_annotation = ha2,
    column_names_rot = 0,
    rect_gp = gpar(col = "black"),
    border = TRUE,
    name = "Euclidean Distance",
    clustering_method_rows = "ward.D",
    clustering_distance_rows = "euclidean"
  )

plot

library(ggplotify)

plot <- as.ggplot(plot)

# ggsave(plot,
#        filename = "heatmap.pdf",
#        width = 10,
#        height = 7)

cluster_group <-
  cluster %>%
  data.frame(cluster = .) %>%
  tibble::rownames_to_column(var = "subject_id")

distance_individual <-
  distance_individual %>%
  do.call(cbind, .) %>%
  data.frame(time = c("0",
                      "30",
                      "60",
                      "120",
                      "240")) %>%
  tidyr::pivot_longer(cols = -time,
                      names_to = "subject_id",
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = c("0",
                                               "30",
                                               "60",
                                               "120",
                                               "240"))) %>%
  dplyr::left_join(cluster_group, by = "subject_id")

cluster_color = c(
  "1" = ggsci::pal_d3()(10)[1],
  "2" = ggsci::pal_d3()(10)[2],
  "3" = ggsci::pal_d3()(10)[3],
  "4" = ggsci::pal_d3()(10)[4],
  "5" = ggsci::pal_d3()(10)[5]
)

plot <-
  distance_individual %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  ggplot(aes(time, value, group = cluster)) +
  geom_point(aes(fill = cluster, group = subject_id), shape = 21) +
  geom_line(aes(color = cluster, group = subject_id)) +
  scale_fill_manual(values = cluster_color) +
  scale_color_manual(values = cluster_color) +
  theme_bw() +
  labs(x = "Time", y = "Euclidean Distance") +
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
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 15),
    panel.grid.minor = element_blank()
  )

plot

# ggsave(plot,
#        filename = "lipidomics_inter_vist_distance.pdf",
#        width = 9,
#        height = 7)

####cluster 1 3 should be one group, and cluster 2 and 4 should be one group
cluster[cluster == 1 | cluster == 3] <- "group1"
cluster[cluster == 2 | cluster == 4] <- "group2"

cluster <-
  cluster %>%
  data.frame(group = .) %>%
  tibble::rownames_to_column(var = "subject_id")

subject_info <-
  subject_info %>%
  dplyr::left_join(cluster, by = "subject_id")

####bmi
library(ggpubr)

p <- ggboxplot(
  subject_info %>%
    dplyr::filter(!is.na(group)),
  x = "group",
  y = "bmi",
  color = "group",
  palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
) +
  theme_bw() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  ) +
  labs(x = "", y = "BMI")
p

my_comparisons <- list(c("group1", "group2"))
p <-
  p +
  stat_compare_means(comparisons = my_comparisons)# Add pairwise comparisons p-value

p

# ggsave(p,
#        filename = "bmi_compare.pdf",
#        width = 7,
#        height = 7)

####sspg
library(ggpubr)

p <- ggboxplot(
  subject_info %>%
    dplyr::filter(!is.na(group)),
  x = "group",
  y = "sspg",
  color = "group",
  palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
) +
  theme_bw() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  ) +
  labs(x = "", y = "sspg")
p

my_comparisons <- list(c("group1", "group2"))
p <-
  p +
  stat_compare_means(comparisons = my_comparisons)# Add pairwise comparisons p-value

p

# ggsave(p,
#        filename = "sspg_compare.pdf",
#        width = 7,
#        height = 7)


####age
library(ggpubr)

p <- ggboxplot(
  subject_info %>%
    dplyr::filter(!is.na(group)),
  x = "group",
  y = "age",
  color = "group",
  palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
) +
  theme_bw() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank()
  ) +
  labs(x = "", y = "age")
p

my_comparisons <- list(c("group1", "group2"))
p <-
  p +
  stat_compare_means(comparisons = my_comparisons)# Add pairwise comparisons p-value

p

# ggsave(p,
#        filename = "age_compare.pdf",
#        width = 7,
#        height = 7)

####sex
temp_data <-
  subject_info %>%
  dplyr::filter(!is.na(group))

table(temp_data$group, temp_data$sex) %>%
  chisq.test()

table(temp_data$group, temp_data$ethnicity) %>%
  chisq.test()

# ###find the difference expression metabolites in different groups
# library(OmicsLonDA)
# library(SummarizedExperiment)
#
# temp_data <-
#   expression_data1 %>%
#   apply(1, function(x){
#     (x - mean(x)) / sd(x)
#   }) %>%
#   t()
#
# temp_sample_info <-
#   sample_info %>%
#   dplyr::left_join(cluster, by = "subject_id") %>%
#   dplyr::select(sample_id, subject_id, TP, group) %>%
#   tibble::column_to_rownames(var = "sample_id") %>%
#   dplyr::select(Subject = subject_id, Group = group, Time = TP)
#
# se_ome_matrix = as.matrix(temp_data)
# se_metadata = DataFrame(temp_sample_info)
# omicslonda_se_object = SummarizedExperiment(assays=list(se_ome_matrix),
#                                             colData = se_metadata)
#
# # omicslonda_se_object_adjusted = adjustBaseline(se_object = omicslonda_se_object)
#
# assay(omicslonda_se_object_adjusted)[1:5, 1:5]
#
# important_feature <- rep(NA, nrow(variable_info))
# for(idx in 17:nrow(variable_info)){
#   cat(idx, " ")
#   omicslonda_test_object = omicslonda_se_object[idx,]
#   visualizeFeature(
#     se_object = omicslonda_test_object,
#     text = stringr::str_replace_all(variable_info$Metabolite[idx],"/", '_'),
#     unit = "mins",
#     ylabel = "Normalized intensity",
#     col = c("#1F77B4FF", "#D62728FF"),
#     prefix = "OmicsLonDA"
#   )
#
#   points = seq(0, 240, by = 1)
#
#   res = omicslonda(
#     se_object = omicslonda_test_object,
#     n.perm = 100,
#     fit.method = "ssgaussian",
#     points = points,
#     text = stringr::str_replace_all(variable_info$Metabolite[idx],"/", '_'),
#     parall = TRUE,
#     pvalue.threshold = 0.05,
#     adjust.method = "BH",
#     time.unit = "mins",
#     ylabel = "Normalized Count",
#     col = c("#1F77B4FF", "#D62728FF"),
#     prefix = "OmicsLonDA"
#   )
#
#   if(length(res$start) == 0 &
#      length(res$end) == 0){
#     next()
#   }else{
#     important_feature[idx] <- idx
#     visualizeFeatureSpline(
#       se_object = omicslonda_test_object,
#       omicslonda_object = res,
#       fit.method = "ssgaussian",
#       text = stringr::str_replace_all(variable_info$Metabolite[idx],"/", '_'),
#       unit = "mins",
#       ylabel = "Normalized Count",
#       col = c("#1F77B4FF", "#D62728FF"),
#       prefix = "OmicsLonDA"
#     )
#   }
# }

###simple method
temp_data <-
  expression_data1

temp_sample_info <-
  sample_info %>%
  dplyr::left_join(cluster, by = "subject_id") %>%
  dplyr::select(sample_id, subject_id, TP, group) %>%
  tibble::column_to_rownames(var = "sample_id") %>%
  dplyr::select(Subject = subject_id,
                Group = group,
                Time = TP)

p_fc <-
  purrr::map(
    as.data.frame(t(temp_data)),
    .f = function(x) {
      temp_p_fc <-
        purrr::map(
          c(0, 30, 60, 120, 240),
          .f = function(y) {
            test <-
              wilcox.test(x[temp_sample_info$Group == "group1" &
                              temp_sample_info$Time == y],
                          x[temp_sample_info$Group == "group2" &
                              temp_sample_info$Time == y])
            
            c(broom::tidy(test)$p.value,
              mean(x[temp_sample_info$Group == "group2" &
                       temp_sample_info$Time == y]) /
                mean(x[temp_sample_info$Group == "group1" &
                         temp_sample_info$Time == y]))
          }
        ) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      temp_p_fc
    }
  )

p_fc

idx <-
  purrr::map(p_fc, function(x) {
    any(x$V1 < 0.05)
  }) %>%
  unlist() %>%
  which()

p_fc2 <- p_fc

names(p_fc2)

###output all the data
temp_data <-
  expression_data1
# apply(1, function(x){
#   (x - mean(x)) / sd(x)
# }) %>%
# t()

temp_sample_info <-
  sample_info %>%
  dplyr::left_join(cluster, by = "subject_id") %>%
  dplyr::select(sample_id, subject_id, TP, group) %>%
  dplyr::select(sample_id,
                Subject = subject_id,
                Group = group,
                Time = TP)


dir.create("important_features")

group_color = c("group1" = ggsci::pal_d3()(10)[1],
                "group2" = ggsci::pal_d3()(10)[2])

library(plyr)

# purrr::walk(
#   names(p_fc2),
#   .f = function(x) {
#     cat(x, " ")
#     name <- variable_info$mol_name[variable_info$variable_id == x]
#     plot <-
#       data.frame(value = as.numeric(temp_data[x, ]),
#                  sample_id = colnames(temp_data)) %>%
#       dplyr::left_join(temp_sample_info, by = "sample_id") %>%
#       dplyr::mutate(Time = factor(as.character(Time),
#                                   levels = c("0", "30", "60", '120', "240"))) %>%
#       plyr::dlply(.variables = .(Subject)) %>%
#       purrr::map(
#         .f = function(y) {
#           y %>%
#             dplyr::mutate(value = value - y$value[y$Time == 0])
#         }
#       ) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       ggplot(aes(x = Time, y = value, group = Subject)) +
#       geom_line(aes(group = Subject, col = Group)) +
#       geom_point(aes(fill = Group, group = Subject), shape = 21) +
#       scale_fill_manual(values = group_color) +
#       scale_color_manual(values = group_color) +
#       theme_bw() +
#       labs(x = "Time", y = "Scaled intensity", title = name) +
#       theme(
#         axis.title = element_text(size = 13),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 13),
#         legend.text = element_text(size = 12),
#         legend.position = c(0, 1),
#         legend.justification = c(0, 1),
#         panel.background = element_rect(fill = "transparent", color = NA),
#         plot.background = element_rect(fill = "transparent", color = NA),
#         legend.background = element_rect(fill = "transparent", color = NA),
#         strip.background = element_rect(fill = "#0099B47F"),
#         strip.text = element_text(color = "white", size = 15),
#         panel.grid.minor = element_blank()
#       )
#     
#     name <- stringr::str_replace_all(name, "/", "_")
#     ggsave(
#       plot,
#       filename = file.path("important_features",
#                            paste(name, ".pdf", sep = "")),
#       width = 9,
#       height = 7
#     )
#   }
# )
