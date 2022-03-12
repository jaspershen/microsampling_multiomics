no_source()

library(tidyverse)
library(data.table)
library(tidymass)

tinyTools::setwd_project()

###load data
load("data/shake_study/metabolome_data_analysis/data_preparation/peaks/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/peaks/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/peaks/variable_info")

load("data/shake_study/metabolome_data_analysis/metabolites/DEG/anova_marker_name")

###load plasma and milk match table
load("data/shake_study/milk_component_in_body/plasma_milk_match_table")
load("data/shake_study/milk_component_in_body/idx")

annotation_table =
  readr::read_csv("data/shake_study/milk_component_in_body/annotation_table_manual.csv")

plasma_milk_match_table = 
plasma_milk_match_table %>% 
  dplyr::left_join(annotation_table[,c("plasma_variable_id", "Check")], 
                   by = c("plasma_variable_id")) %>% 
  dplyr::filter(in_baseline == "No")

###remove the duplicated annotations for milk annotation
library(plyr)
plasma_milk_match_table = 
plasma_milk_match_table %>% 
  mutate(class = case_when(
    is.na(milk_annotation) ~ "no",
    !is.na(milk_annotation) ~ milk_annotation
  )) %>% 
  plyr::dlply(.variables = .(class)) %>% 
  purrr::map(function(x){
    if(x$class[1] == 'no'){
      return(x)
    }
    x %>% 
      dplyr::filter(rt.error == min(rt.error)) %>% 
      dplyr::filter(mz.error == min(mz.error))
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()
  
variable_info = 
variable_info %>% 
  dplyr::left_join(plasma_milk_match_table[,c("milk_variable_id", "plasma_variable_id", "milk_annotation", "Check")],
                   by = c("variable_id" = "plasma_variable_id"))

###remove the annotations which are not different
variable_info = 
variable_info %>% 
  dplyr::filter(is.na(Check) | Check == "Same") %>% 
  dplyr::filter(!is.na(milk_variable_id))

expression_data = expression_data[variable_info$variable_id,]

###ANOVA test
###find marker which are change according to pregnancy
expression_data[5,] %>% 
  as.numeric() %>% 
  density() %>%
  plot()

subject_data <- 
  expression_data %>% 
  dplyr::select(one_of(sample_info$sample_id))

##log transformation
subject_data <-
  log(subject_data + 1, 2)

###ANOVA analysis
library(tidyverse)
library(ggpubr)
library(rstatix)

# anova_p <-
#   purrr::map(as.data.frame(t(subject_data)), .f = function(x){
#     temp_data <-
#       data.frame(
#         subject_id = sample_info$subject_id,
#         tp = factor(
#           sample_info$TP,
#           levels = stringr::str_sort(unique(sample_info$TP),
#                                      numeric = TRUE)),
#         x = x,
#         stringsAsFactors = FALSE
#       )
# 
#     library(plyr)
#     intersect_name <-
#       temp_data %>%
#       plyr::dlply(.variables = .(tp)) %>%
#       purrr::map(.f = function(x){
#         x$subject_id
#       }) %>%
#       Reduce(intersect, .) %>%
#       unique()
# 
#     temp_data <-
#       temp_data %>%
#       dplyr::filter(subject_id %in% intersect_name)
# 
#     temp_data$subject_id <-
#       factor(temp_data$subject_id, levels = unique(temp_data$subject_id))
# 
#     res.aov <-
#       anova_test(
#         data = temp_data,
#         dv = x,
#         wid = subject_id,
#         within = tp
#       )
# 
#     p <-
#       as.data.frame(get_anova_table(res.aov))$p
# 
#     # pairwise comparisons
#     pwc <- temp_data %>%
#       pairwise_t_test(x ~ tp,
#                       paired = TRUE,
#                       p.adjust.method = "fdr")
#     pwc
# 
#     p2 <-
#       pwc[,c(2,3,8,9)]
# 
#     name <- paste(pwc$group2, pwc$group1, sep = "_")
# 
#     p2 <- matrix(pwc$p, nrow = 1)
# 
#     colnames(p2) <- name
# 
#     p <-
#       data.frame(p, p2, stringsAsFactors = FALSE, check.names = FALSE)
# 
#     p
# 
#   })
# 
# anova_p <- anova_p %>% do.call(rbind, .) %>% as.data.frame()
# save(anova_p, file = "data/shake_study/milk_component_in_body/fuzzy_c_means/anova_p")
load("data/shake_study/milk_component_in_body/fuzzy_c_means/anova_p")

anova_marker_name <-
  rownames(anova_p)[which(anova_p$p < 0.05)]

# expression_data = 
#   expression_data[anova_marker_name,]
# 
# variable_info =
#   variable_info[match(anova_marker_name,variable_info$variable_id),]

sum(!is.na(variable_info$milk_annotation))

####metabolite cluster information
cluster1_info = 
  readr::read_csv("data/shake_study/3_omics/k_means_clustering/cluster_3/cluster3_metabolite_manual.csv")

cluster2_info = 
readr::read_csv("data/shake_study/3_omics/k_means_clustering/cluster_1/cluster1_metabolite_manual.csv")

cluster3_info = 
  readr::read_csv("data/shake_study/3_omics/k_means_clustering/cluster_2/cluster2_metabolite_manual.csv")

tinyTools::setwd_project()
setwd("data/shake_study/milk_component_in_body/fuzzy_c_means")

dim(expression_data)
dim(sample_info)
dim(variable_info)

intersect(variable_info$variable_id, cluster1_info$variable_id)
intersect(variable_info$variable_id, cluster2_info$variable_id)
intersect(variable_info$variable_id, cluster3_info$variable_id)

###cluster
temp_data <-
  log(expression_data + 1, 2)

sample_info$TP

tp <- c(0, 30, 60, 120, 240)

temp_data_mean <-
  purrr::map(tp, function(x) {
    temp_idx <- which(sample_info$TP == x)
    temp_data[, temp_idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame() %>% 
  apply(1, function(x)
    (x - mean(x)) / sd(x)) %>%
  t() %>%
  as.data.frame()

colnames(temp_data_mean) <- paste("Time", tp, sep = "_")

library(ComplexHeatmap)
Heatmap(temp_data_mean, show_row_names = FALSE, cluster_columns = FALSE)

#######fuzzy k means
library(Mfuzz)
library(e1071)

sample_info$sample_id == colnames(expression_data)
sample_info$TP

temp_data <-
  log(expression_data + 1, 2)

sample_info$TP

tp <- c(0, 30, 60, 120, 240)

temp_data_mean <-
  purrr::map(tp, function(x) {
    temp_idx <- which(sample_info$TP == x)
    temp_data[, temp_idx] %>%
      apply(1, mean)
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(temp_data_mean) <- paste("Time", tp, sep = "_")

library(Mfuzz)
#first get the time point data together:
# bind that to the data frame
##scale
temp_data <-
  temp_data_mean %>% 
  apply(1, function(x)
    (x - mean(x)) / sd(x)) %>%
  t() %>%
  as.data.frame()

time <- c(0, 30, 60, 120, 240)

temp_data <- rbind(time, temp_data)

row.names(temp_data)[1] <- "time"

# write.table(
#   temp_data,
#   file = "temp_data.txt",
#   sep = '\t',
#   quote = FALSE,
#   col.names = NA
# )

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

plot <-
  Dmin(
    data.s,
    m = m1,
    crange = seq(2, 22, 1),
    repeats = 3,
    visu = TRUE
  )

plot <-
  plot %>%
  data.frame(distance = plot,
             k = seq(2, 22, 1)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21, size = 4, fill = "black") +
  # geom_smooth() +
  geom_segment(aes(
    x = k,
    y = 0,
    xend = k,
    yend = distance
  )) +
  theme_bw() +
  theme(
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    panel.grid = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(x = "Cluster number",
       y = "Min. centroid distance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

plot

# ggsave(plot,
#        filename = "distance_k_number.pdf",
#        width = 7,
#        height = 7)

cluster = 2

c <- mfuzz(data.s, c = cluster, m = m1)

mfuzz.plot(
  eset = data.s,
  # min.mem = 0.6,
  cl = c,
  mfrow = c(1, 1),
  time.labels = time,
  new.window = FALSE
)

# save(c, file = "c")
load("c")
#
#
# ####any two clusters with correlation > 0.8 should be considered as one
library(corrplot)
layout(1)
center <- c$centers
rownames(center) <- paste("Cluster", rownames(center), sep = ' ')

corrplot::corrplot(
  corr = cor(t(center)),
  type = "full",
  diag = TRUE,
  order = "hclust",
  hclust.method = "ward.D",
  # addrect = 5,
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

library(ComplexHeatmap)

###
cluster_color <-
  ggsci::pal_jama()(n = 7)[1:4]

names(cluster_color) <- as.character(1:4)

centers <- c$centers

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

# xlsx::write.xlsx(cluster_info,
#                  "cluster_info.xlsx",
#                  row.names = FALSE)

####plot for each cluster
idx <- 1

value <- 
  c("Lipid" = ggsci::pal_aaas()(10)[1],
    "Metabolite" = ggsci::pal_aaas()(10)[3],
    "Cytokine" = ggsci::pal_aaas()(10)[4]
  )

# for (idx in 1:4) {
#   cat(idx, " ")
#   cluster_data <-
#     cluster_info %>%
#     dplyr::filter(cluster == idx) %>%
#     dplyr::select(1, 1 + idx, cluster)
# 
#   colnames(cluster_data)[2] <- c("membership")
# 
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.5)
# 
#   if(nrow(cluster_data) == 0){
#     next()
#   }
#   
#   path <- paste("cluster", idx, sep = "_")
#   dir.create(path)
# 
#   xlsx::write.xlsx(cluster_data,
#                    file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
#                    row.names = FALSE)
# 
#   temp_center <-
#     centers[idx, , drop = TRUE] %>%
#     data.frame(time = names(.),
#                value = .,
#                stringsAsFactors = FALSE) %>%
#     dplyr::mutate(time = factor(time, levels = time))
# 
#   plot <-
#     temp_data[cluster_data$variable_id,] %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     dplyr::mutate(
#       class = "Metabolite"
#     ) %>%
#     dplyr::left_join(variable_info[,c("variable_id", "milk_annotation")], by = "variable_id") %>% 
#     dplyr::mutate(class2 = case_when(
#       is.na(milk_annotation) ~ "NO",
#       !is.na(milk_annotation) ~ "Yes"
#     )) %>% 
#     dplyr::select(-milk_annotation) %>% 
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership, class, class2),
#       names_to = "time",
#       values_to = "value"
#     ) %>%
#     dplyr::mutate(time = factor(time, levels = unique(time))) %>%
#     ggplot(aes(time, value, group = variable_id)) +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(color = class2,
#                   alpha = class2),
#               show.legend = FALSE) +
#     scale_color_manual(values = c("Yes" = ggsci::pal_aaas()(10)[3], 
#                                   "NO" = "grey")) +
#     scale_alpha_manual(values = c("Yes" = 1, "NO" = 1)) +
#     theme_bw() +
#     theme(
#       legend.position = c(0, 1),
#       legend.justification = c(0, 1),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size = 13),
#       axis.text = element_text(size = 12),
#       axis.text.x = element_text(
#         angle = 45,
#         hjust = 1,
#         vjust = 1,
#         size = 12
#       ),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background = element_rect(fill = "transparent", color = NA),
#       legend.background = element_rect(fill = "transparent", color = NA)
#     ) +
#     labs(
#       x = "",
#       y = "Z-score",
#       title = paste("Cluster ",
#                     idx,
#                     " (",
#                     nrow(cluster_data),
#                     " features)",
#                     sep = "")
#     ) +
#     geom_line(
#       mapping = aes(time, value, group = 1),
#       data = temp_center,
#       size = 2
#     )
# 
#   plot
# 
#   ggsave(
#     plot,
#     filename = file.path(path, paste("cluster", idx, ".pdf", sep = "")),
#     width = 8,
#     height = 7
#   )
# }

dim(cluster_data)

###annotation for each cluster
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx")
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx")
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx")

nrow(cluster1) +
  nrow(cluster2) +
  nrow(cluster3)

cluster1 %>%
  dplyr::left_join(plasma_milk_match_table[, c("plasma_variable_id", "in_baseline")],
                   by = c("variable_id" = "plasma_variable_id")) %>% 
  pull(in_baseline) %>% 
  table()

cluster2 %>%
  dplyr::left_join(plasma_milk_match_table[, c("plasma_variable_id", "in_baseline")],
                   by = c("variable_id" = "plasma_variable_id")) %>% 
  pull(in_baseline) %>% 
  table()

cluster3 %>%
  dplyr::left_join(plasma_milk_match_table[, c("plasma_variable_id", "in_baseline")],
                   by = c("variable_id" = "plasma_variable_id")) %>% 
  pull(in_baseline) %>% 
  table()

##metabolites that not in clusters
non_metabolite <-
  setdiff(
    cluster_info$variable_id,
    c(
      cluster1$variable_id,
      cluster2$variable_id,
      cluster3$variable_id
    )
  )

plot <-
  temp_data[non_metabolite, ] %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(
    cols = -c(variable_id),
    names_to = "time",
    values_to = "value"
  ) %>%
  dplyr::mutate(time = factor(time, levels = unique(time))) %>%
  ggplot(aes(time, value, group = variable_id)) +
  geom_line(alpha = 1) +
  # scale_color_gradientn(colours = c(
  #   RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)]
  # )) +
  theme_bw() +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12
    ),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(
    x = "",
    y = "Z-score",
    title = paste(
      "Unclustered metabolic peaks (",
      length(non_metabolite),
      " metabolic peaks)",
      sep = ""
    )
  )

plot

# ggsave(plot,
#        filename = "non_cluster_metabolites.pdf",
#        width = 8,
#        height = 7)


dim(cluster1)
dim(cluster2)
dim(cluster3)


variable_info %>% 
  dplyr::filter(variable_id %in% cluster1$variable_id) %>% 
  dplyr::filter(!is.na(milk_annotation)) %>% 
  dplyr::pull(milk_annotation) %>% 
  unique()

variable_info %>% 
  dplyr::filter(variable_id %in% cluster2$variable_id) %>% 
  dplyr::filter(!is.na(milk_annotation)) %>% 
  dplyr::pull(milk_annotation) %>% 
  unique()

variable_info %>% 
  dplyr::filter(variable_id %in% cluster3$variable_id) %>% 
  dplyr::filter(!is.na(milk_annotation)) %>% 
  dplyr::pull(milk_annotation) %>% 
  unique()









# ###heatmap for all the samples
# library(ComplexHeatmap)
# 
# temp_data <-
#   temp_data_mean[c(cluster2$variable_id,
#                    cluster1$variable_id,
#                    cluster3$variable_id),]
# 
# temp_data <-
#   apply(temp_data, 1, function(x) {
#     (x - mean(x)) / sd(x)
#   }) %>%
#   t() %>%
#   as.data.frame()
# 
# rownames(temp_data) <-
#   variable_info$mol_name[match(rownames(temp_data), variable_info$variable_id)]
# 
# pal <-
#   wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")
# 
# range(temp_data)
# library(circlize)
# col_fun = colorRamp2(c(-2, 0, 2), c(pal[1], "white", pal[100]))
# 
# cluster <-
#   c(rep(2, nrow(cluster2)),
#     rep(1, nrow(cluster1)),
#     rep(3, nrow(cluster3)))
# 
# cluster_color <-
#   c(ggsci::pal_d3()(n = 10)[1:3])
# 
# names(cluster_color) <- c(as.character(1:3))
# 
# text_color <- 
#   cluster_color[cluster]
# 
# ha <-
#   rowAnnotation(Cluster = cluster,
#                 col = list(Cluster = cluster_color))
# 
# plot <-
#   Heatmap(
#     temp_data,
#     cluster_columns = FALSE,
#     cluster_rows = FALSE,
#     show_row_names = FALSE,
#     border = TRUE,
#     row_names_side = "left",
#     row_names_gp = gpar(fontsize = 7, 
#                         col = text_color),
#     col = col_fun,
#     right_annotation = ha,
#     name = "Z score",
#     column_names_rot = 0
#     # row_split = factor(as.character(cluster), levels = c(6,3,8,1,4,5,2,9,7))
#   )
# 
# plot <- ggplotify::as.ggplot(plot)
# 
# plot
# 
# # ggsave(
# #   plot = plot,
# #   filename = "heatmap_for_all_cluster.pdf",
# #   width = 7,
# #   height = 16
# # )
# 
# 
# ###functional annotation for different cluster
# sxtTools::setwd_project()
# load("data/shake_study/metabolome_data_analysis/metabolites/DEG/subject_data_mean")
# load("data/shake_study/metabolome_data_analysis/metabolites/DEG/subject_data_sd")
# load("data/shake_study/metabolome_data_analysis/metabolites/DEG/subject_data_sem")
# load("data/shake_study/metabolome_data_analysis/metabolites/DEG/fc_p_value")
# metabolome_subject_data_mean <-
#   subject_data_mean
# metabolome_subject_data_sd <-
#   subject_data_sd
# metabolome_subject_data_sem <-
#   subject_data_sem
# metabolome_fc_p <- fc_p_value
# 
# load("data/shake_study/lipidomics_data_analysis/DEG/subject_data_mean")
# load("data/shake_study/lipidomics_data_analysis/DEG/subject_data_sd")
# load("data/shake_study/lipidomics_data_analysis/DEG/subject_data_sem")
# load("data/shake_study/lipidomics_data_analysis/DEG/fc_p_value")
# lipidomics_subject_data_mean <-
#   subject_data_mean
# lipidomics_subject_data_sd <-
#   subject_data_sd
# lipidomics_subject_data_sem <-
#   subject_data_sem
# lipidomics_fc_p <- fc_p_value
# 
# load("data/shake_study/cytokine_data_analysis/DEG/subject_data_mean")
# load("data/shake_study/cytokine_data_analysis/DEG/subject_data_sd")
# load("data/shake_study/cytokine_data_analysis/DEG/subject_data_sem")
# load("data/shake_study/cytokine_data_analysis/DEG/fc_p_value")
# cytokine_subject_data_mean <-
#   subject_data_mean
# cytokine_subject_data_sd <-
#   subject_data_sd
# cytokine_subject_data_sem <-
#   subject_data_sem
# cytokine_fc_p <- fc_p_value
# 
# setwd("data/shake_study/3_omics/k_means_clustering/")
# 
# ##metabolomics
# cluster1_metabolite <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster1$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info1$variable_id)
# 
# cluster2_metabolite <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster2$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info1$variable_id)
# 
# cluster3_metabolite <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster3$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info1$variable_id)
# 
# non_cluster_metabolite <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% non_metabolite) %>% 
#   dplyr::filter(variable_id %in% variable_info1$variable_id)
# 
# cluster1_metabolite$Metabolite
# cluster2_metabolite$Metabolite
# cluster3_metabolite$Metabolite
# non_cluster_metabolite$Metabolite
# 
# ####get the KEGG ID and HMDB ID
# cluster1_metabolite$Metabolite
# cluster1_metabolite$HMDB
# getwd()
# # write.csv(cluster1_metabolite, "cluster_1/cluster1_metabolite.csv", row.names = FALSE)
# # write.csv(cluster2_metabolite, "cluster_2/cluster2_metabolite.csv", row.names = FALSE)
# # write.csv(cluster3_metabolite, "cluster_3/cluster3_metabolite.csv", row.names = FALSE)
# 
# ###cluster2
# cluster2_path <- 
#   readr::read_csv("cluster_2/pathway_enrichment/pathway_results.csv")
# 
# plot <-
#   cluster2_path %>%
#   dplyr::mutate(fdr = -log(p.adjust(`Raw p`, method = "fdr"), 10)) %>%
#   dplyr::mutate(class = case_when(fdr > 1.30103 & Hits >= 1 ~ "Yes",
#                                   TRUE ~ "No")) %>%
#   ggplot(aes(Impact, fdr)) +
#   geom_hline(yintercept = 1.30103) +
#   geom_point(aes(size = Hits, fill = class), shape = 21) +
#   scale_fill_manual(values = c(
#     "Yes" = ggsci::pal_aaas()(n=10)[2],
#     "No" = ggsci::pal_aaas()(n=10)[9]
#   )) +
#   scale_size_continuous(range = c(10,20)) +
#   ggrepel::geom_text_repel(aes(label = ifelse(class == "Yes", X1, NA))) +
#   labs(x = "Pathway Impact", y = "-log10(FDR adjusted P value)") +
#   theme_bw() +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 13)
#   )
# plot
# # ggsave(plot, filename = "cluster_2/cluster2_pathway.pdf", width = 8, height = 7)
# 
# 
# ###cluster3
# cluster3_path <- 
#   readr::read_csv("cluster_3/pathway_enrichment/pathway_results.csv")
# 
# plot <-
#   cluster3_path %>%
#   dplyr::mutate(fdr = -log(p.adjust(`Raw p`, method = "fdr"), 10)) %>%
#   dplyr::mutate(class = case_when(fdr > 1.30103 & Hits > 1 ~ "Yes",
#                                   TRUE ~ "No")) %>%
#   ggplot(aes(Impact, fdr)) +
#   geom_hline(yintercept = 1.30103) +
#   geom_point(aes(size = Hits, fill = class), shape = 21) +
#   scale_fill_manual(values = c(
#     "Yes" = ggsci::pal_aaas()(n=10)[2],
#     "No" = ggsci::pal_aaas()(n=10)[9]
#   )) +
#   scale_size_continuous(range = c(5,15)) +
#   ggrepel::geom_text_repel(aes(label = ifelse(class == "Yes", X1, NA))) +
#   labs(x = "Pathway Impact", y = "-log10(FDR adjusted P value)") +
#   theme_bw() +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 13)
#   )
# plot
# # ggsave(plot, filename = "cluster_3/cluster3_pathway.pdf", width = 8, height = 7)
# 
# 
# ###functional annotation for different cluster
# #lipidomics
# cluster1_lipid <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster1$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info2$variable_id)
# 
# cluster2_lipid <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster2$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info2$variable_id)
# 
# cluster3_lipid <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster3$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info2$variable_id)
# 
# non_cluster_lipid <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% non_metabolite) %>% 
#   dplyr::filter(variable_id %in% variable_info2$variable_id)
# 
# cluster1_lipid$mol_name
# cluster2_lipid$mol_name
# cluster3_lipid$mol_name
# non_cluster_lipid$mol_name
# 
# 
# ###lipid class changes
# 
# 
# # getwd()
# # write.csv(cluster1_lipid, "cluster_1/cluster1_lipid.csv", row.names = FALSE)
# # write.csv(cluster2_lipid, "cluster_2/cluster2_lipid.csv", row.names = FALSE)
# # write.csv(cluster3_lipid, "cluster_3/cluster3_lipid.csv", row.names = FALSE)
# 
# 
# 
# 
# ###functional annotation for different cluster
# color <- ggsci::pal_futurama()(n=10)[1:5]
# names(color) <- c("SM", "CE", "DAG", "TAG", "FFA")
# 
# temp1 <-
#   cluster1_lipid %>%
#   dplyr::pull(subclass) %>%
#   table() %>%
#   as.data.frame()
# 
# colnames(temp1)[1] <- "class"
# 
# temp1 <- 
#   temp1 %>%
#   dplyr::arrange(Freq) %>% 
#   dplyr::mutate(class = factor(class, levels = class)) %>% 
#   dplyr::mutate(Freq = Freq * 100/sum(Freq)) %>% 
#   mutate(prop = round(Freq / sum(Freq) * 100, 2)) %>%
#   mutate(ypos = cumsum(prop) - 0.5 * prop)
# 
# library(scales)
# 
# plot1 <- 
#   temp1 %>%
#   ggplot(aes(x = "", y = Freq, fill = class)) +
#   geom_bar(width = 1,
#            stat = "identity",
#            aes(fill = class),
#            color = "white") +
#   scale_fill_manual(values = color) +
#   geom_text(aes(y = ypos, label = prop), color = "white", size=6) +
#   theme_void() +
#   coord_polar("y", start = 0) 
# 
# plot1  
# 
# # ggsave(plot1, filename = "cluster_1/cluster1_annotation.pdf", width = 7, height = 7)
# 
# temp2 <-
#   cluster2_lipid %>% 
#   dplyr::pull(subclass) %>%
#   table() %>%
#   as.data.frame()
# 
# colnames(temp2)[1] <- "class"
# 
# temp2 <- 
#   temp2 %>%
#   dplyr::arrange(Freq) %>% 
#   dplyr::mutate(class = factor(class, levels = class)) %>% 
#   dplyr::mutate(Freq = Freq * 100/sum(Freq)) %>% 
#   mutate(prop = round(Freq / sum(Freq) * 100, 2)) %>%
#   mutate(ypos = cumsum(prop) - 0.5 * prop)
# 
# library(scales)
# 
# plot2 <- 
#   temp2 %>%
#   ggplot(aes(x = "", y = prop, fill = class)) +
#   geom_bar(width = 1,
#            stat = "identity",
#            aes(fill = class),
#            color = "white") +
#   scale_fill_manual(values = color) +
#   coord_polar("y", start = 0) +
#   geom_text(aes(y = ypos, label = prop), color = "white", size=6) +
#   theme_void() 
# 
# plot2
# 
# # ggsave(plot2, filename = "cluster_2/cluster2_annotation.pdf", width = 7, height = 7)
# 
# temp3 <-
#   cluster3_lipid %>%
#   dplyr::pull(subclass) %>%
#   table() %>%
#   as.data.frame()
# 
# colnames(temp3)[1] <- "class"
# 
# temp3 <- 
#   temp3 %>%
#   dplyr::arrange(Freq) %>% 
#   dplyr::mutate(class = factor(class, levels = class)) %>% 
#   dplyr::mutate(Freq = Freq * 100/sum(Freq)) %>% 
#   mutate(prop = round(Freq / sum(Freq) * 100, 2)) %>%
#   mutate(ypos = cumsum(prop) - 0.5 * prop)
# 
# library(scales)
# 
# plot3 <- 
#   temp3 %>%
#   ggplot(aes(x = "", y = prop, fill = class)) +
#   geom_bar(width = 1,
#            stat = "identity",
#            aes(fill = class),
#            color = "white") +
#   scale_fill_manual(values = color) +
#   coord_polar("y", start = 0) +
#   geom_text(aes(y = ypos, label = prop), color = "white", size=6) +
#   theme_void() 
# 
# plot3
# 
# # ggsave(plot3, filename = "cluster3_annotation.pdf", width = 7, height = 7)
# 
# ###functional annotation for different cluster
# #cytokine
# cluster1_cytokine <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster1$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info3$variable_id)
# 
# cluster2_cytokine <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster2$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info3$variable_id)
# 
# cluster3_cytokine <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% cluster3$variable_id) %>% 
#   dplyr::filter(variable_id %in% variable_info3$variable_id)
# 
# non_cluster_cytokine <- 
#   variable_info %>% 
#   dplyr::filter(variable_id %in% non_metabolite) %>% 
#   dplyr::filter(variable_id %in% variable_info3$variable_id)
# 
# cluster1_cytokine$mol_name
# cluster2_cytokine$mol_name
# cluster3_cytokine$mol_name
# non_cluster_cytokine$mol_name
# 
# getwd()
# # write.csv(cluster1_cytokine, "cluster_1/cluster1_cytokine.csv", row.names = FALSE)
# # write.csv(cluster2_cytokine, "cluster_2/cluster2_cytokine.csv", row.names = FALSE)
# # write.csv(cluster3_cytokine, "cluster_3/cluster3_cytokine.csv", row.names = FALSE)
# 
# idx <- 
#   match(cluster3_cytokine$variable_id, 
#         rownames(cytokine_subject_data_mean))
# 
# mean_value <- 
#   cytokine_subject_data_mean[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "mean")
# 
# sem_value <- 
#   cytokine_subject_data_sem[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "sem")
# 
# plot <- 
#   mean_value %>%
#   dplyr::left_join(sem_value, by = c("variable_id", "TP")) %>%
#   dplyr::mutate(TP = factor(as.character(TP),
#                             levels = c("0", "30", "60", "120", "240"))) %>% 
#   dplyr::left_join(variable_info3, by = "variable_id") %>% 
#   ggplot(aes(TP, mean, group = mol_name)) +
#   geom_line(aes(color = mol_name)) +
#   geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem,
#                     color = mol_name), width = 0) +
#   theme_bw() +
#   ggsci::scale_color_lancet() +
#   labs(x = "Time point (min)", y = "Mean(log2 intensity)") +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) 
# 
# plot
# 
# # ggsave(plot, filename = "cluster_3/cluster3_cytokine_plot.pdf", width = 9, height = 7)
# 
# 
# ##amino acid
# amino_acid_variable_id <- 
#   cluster3_metabolite$variable_id[c(1, 2, 3, 10, 12, 17, 19, 33, 39)]
# 
# idx <- 
#   match(amino_acid_variable_id, 
#         rownames(metabolome_subject_data_mean))
# 
# mean_value <- 
#   metabolome_subject_data_mean[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "mean")
# 
# sem_value <- 
#   metabolome_subject_data_sem[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "sem")
# 
# plot <- 
#   mean_value %>%
#   dplyr::left_join(sem_value, by = c("variable_id", "TP")) %>%
#   dplyr::mutate(TP = factor(as.character(TP),
#                             levels = c("0", "30", "60", "120", "240"))) %>% 
#   dplyr::left_join(variable_info1, by = "variable_id") %>% 
#   ggplot(aes(TP, mean, group = Metabolite)) +
#   geom_line(aes(color = Metabolite)) +
#   geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem,
#                     color = Metabolite), width = 0) +
#   theme_bw() +
#   ggsci::scale_color_jco() +
#   labs(x = "Time point (min)", y = "Mean(log2 intensity)") +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) 
# 
# plot
# # ggsave(plot, filename = "cluster_3/amino_acid_plot.pdf", width = 9, height = 7)
# 
# 
# 
# 
# 
# ##amino acid
# amino_acid_variable_id <- 
#   cluster3_metabolite$variable_id[c(1, 2, 3, 10, 12, 17, 19, 33, 39)]
# 
# idx <- 
#   match(amino_acid_variable_id, 
#         rownames(metabolome_subject_data_mean))
# 
# mean_value <- 
#   metabolome_subject_data_mean[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "mean")
# 
# sem_value <- 
#   metabolome_subject_data_sem[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "sem")
# 
# plot <- 
#   mean_value %>%
#   dplyr::left_join(sem_value, by = c("variable_id", "TP")) %>%
#   dplyr::mutate(TP = factor(as.character(TP),
#                             levels = c("0", "30", "60", "120", "240"))) %>% 
#   dplyr::left_join(variable_info1, by = "variable_id") %>% 
#   ggplot(aes(TP, mean, group = Metabolite)) +
#   geom_line(aes(color = Metabolite)) +
#   geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem,
#                     color = Metabolite), width = 0) +
#   theme_bw() +
#   ggsci::scale_color_jco() +
#   labs(x = "Time point (min)", y = "Mean(log2 intensity)") +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) 
# 
# plot
# # ggsave(plot, filename = "cluster_3/amino_acid_plot.pdf", width = 9, height = 7)
# 
# 
# #Carbohydrates
# carbohydrates_variable_id <- 
#   c(cluster3_metabolite$variable_id[c(15, 20)], "4.83_87.0088m/z")
# 
# idx <- 
#   match(carbohydrates_variable_id, 
#         rownames(metabolome_subject_data_mean))
# 
# mean_value <- 
#   metabolome_subject_data_mean[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "mean")
# 
# sem_value <- 
#   metabolome_subject_data_sem[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "sem")
# 
# plot <- 
#   mean_value %>%
#   dplyr::left_join(sem_value, by = c("variable_id", "TP")) %>%
#   dplyr::mutate(TP = factor(as.character(TP),
#                             levels = c("0", "30", "60", "120", "240"))) %>% 
#   dplyr::left_join(variable_info1, by = "variable_id") %>% 
#   ggplot(aes(TP, mean, group = Metabolite)) +
#   geom_line(aes(color = Metabolite)) +
#   geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem,
#                     color = Metabolite), width = 0) +
#   theme_bw() +
#   ggsci::scale_color_jco() +
#   labs(x = "Time point (min)", y = "Mean(log2 intensity)") +
#   theme(
#     legend.position = c(0, 1),
#     legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) 
# 
# plot
# # ggsave(plot, filename = "cluster_3/carbohydrates_plot.pdf", width = 7, height = 7)
# 
# 
# 
# 
# 
# 
# 
# #acylcarnitine
# acylcarnitine_variable_id <- 
#   cluster2_metabolite$variable_id[c(2, 3, 5, 6, 14:17,20)]
# 
# idx <- 
#   match(acylcarnitine_variable_id, 
#         rownames(metabolome_subject_data_mean))
# 
# mean_value <- 
#   metabolome_subject_data_mean[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "mean")
# 
# sem_value <- 
#   metabolome_subject_data_sem[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "sem")
# 
# plot <- 
#   mean_value %>%
#   dplyr::left_join(sem_value, by = c("variable_id", "TP")) %>%
#   dplyr::mutate(TP = factor(as.character(TP),
#                             levels = c("0", "30", "60", "120", "240"))) %>% 
#   dplyr::left_join(variable_info1, by = "variable_id") %>% 
#   ggplot(aes(TP, mean, group = Metabolite)) +
#   geom_line(aes(color = Metabolite)) +
#   geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem,
#                     color = Metabolite), width = 0) +
#   theme_bw() +
#   ggsci::scale_color_jco() +
#   labs(x = "Time point (min)", y = "Mean(log2 intensity)") +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) 
# 
# plot
# ggsave(plot, filename = "cluster_2/acylcarnitine_plot.pdf", width = 9, height = 7)
# 
# 
# 
# 
# 
# 
# #cytokine
# cytokine_variable_id3 <- 
#   cluster3_cytokine$variable_id
# cytokine_variable_id2 <- 
#   cluster2_cytokine$variable_id
# 
# idx <- 
#   match(c(cytokine_variable_id3, cytokine_variable_id2), 
#         rownames(cytokine_subject_data_mean))
# 
# mean_value <- 
#   cytokine_subject_data_mean[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "mean")
# 
# sem_value <- 
#   cytokine_subject_data_sem[idx, ] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "variable_id") %>%
#   tidyr::pivot_longer(cols = -variable_id,
#                       names_to = "TP",
#                       values_to = "sem")
# 
# plot <- 
#   mean_value %>%
#   dplyr::left_join(sem_value, by = c("variable_id", "TP")) %>%
#   dplyr::mutate(TP = factor(as.character(TP),
#                             levels = c("0", "30", "60", "120", "240"))) %>% 
#   dplyr::left_join(variable_info3, by = "variable_id") %>% 
#   dplyr::mutate(cluster = case_when(
#     variable_id %in% cytokine_variable_id3 ~ "Cluster 3",
#     variable_id %in% cytokine_variable_id2 ~ "Cluster 2"
#   )) %>% 
#   ggplot(aes(TP, mean, group = mol_name)) +
#   geom_line(aes(color = mol_name)) +
#   geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem,
#                     color = mol_name), width = 0) +
#   theme_bw() +
#   ggsci::scale_color_locuszoom() +
#   labs(x = "Time point (min)", y = "Mean(log2 intensity)") +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) +
#   facet_wrap(facets = vars(cluster), scales = "free_y")
# 
# plot
# ggsave(plot, filename = "cytokine_plot.pdf", width = 12, height = 7)
# 
# 
# ###lipid class chanfges
# ##PC, PE, DAG TAG FFA
# pc_idx <- cluster1_metabolite$variable_id[c(5,11)]
# pe_idx <- cluster1_metabolite$variable_id[c(7,8,9,12,13,19)]
# dag_idx <- cluster1_lipid %>% 
#   dplyr::filter(stringr::str_detect(mol_name, "DAG")) %>% 
#   dplyr::pull(variable_id)
# 
# tag_idx <- cluster1_lipid %>% 
#   dplyr::filter(stringr::str_detect(mol_name, "TAG")) %>% 
#   dplyr::pull(variable_id)  
# 
# ffa_idx <- cluster2_lipid %>% 
#   dplyr::filter(stringr::str_detect(mol_name, "FFA")) %>% 
#   dplyr::pull(variable_id)  
# 
# pc_data <- 
#   metabolome_fc_p %>% 
#   purrr::map(.f = function(x){
#     x %>% 
#       dplyr::filter(variable_id %in% pc_idx) %>% 
#       dplyr::select(p_value, fc) %>% 
#       apply(2, mean)
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame() %>% 
#   dplyr::mutate(fc = log(fc, 2)) %>% 
#   dplyr::mutate(p_value = -log(p_value, 10)) %>% 
#   dplyr::add_row(p_value = 0, fc = 0, .before = 1) %>% 
#   dplyr::mutate(TP = c(0, 30, 60, 120, 240)) %>% 
#   dplyr::mutate(class = "PC")
# 
# pe_data <- 
#   metabolome_fc_p %>% 
#   purrr::map(.f = function(x){
#     x %>% 
#       dplyr::filter(variable_id %in% pe_idx) %>% 
#       dplyr::select(p_value, fc) %>% 
#       apply(2, mean)
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame() %>% 
#   dplyr::mutate(fc = log(fc, 2)) %>% 
#   dplyr::mutate(p_value = -log(p_value, 10)) %>% 
#   dplyr::add_row(p_value = 0, fc = 0, .before = 1) %>% 
#   dplyr::mutate(TP = c(0, 30, 60, 120, 240)) %>% 
#   dplyr::mutate(class = "PE")
# 
# dag_data <- 
#   lipidomics_fc_p %>% 
#   purrr::map(.f = function(x){
#     x %>% 
#       dplyr::filter(variable_id %in% dag_idx) %>% 
#       dplyr::select(p_value, fc) %>% 
#       apply(2, mean)
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame() %>% 
#   dplyr::mutate(fc = log(fc, 2)) %>% 
#   dplyr::mutate(p_value = -log(p_value, 10)) %>% 
#   dplyr::add_row(p_value = 0, fc = 0, .before = 1) %>% 
#   dplyr::mutate(TP = c(0, 30, 60, 120, 240)) %>% 
#   dplyr::mutate(class = "DAG")
# 
# tag_data <- 
#   lipidomics_fc_p %>% 
#   purrr::map(.f = function(x){
#     x %>% 
#       dplyr::filter(variable_id %in% tag_idx) %>% 
#       dplyr::select(p_value, fc) %>% 
#       apply(2, mean)
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame() %>% 
#   dplyr::mutate(fc = log(fc, 2)) %>% 
#   dplyr::mutate(p_value = -log(p_value, 10)) %>% 
#   dplyr::add_row(p_value = 0, fc = 0, .before = 1) %>% 
#   dplyr::mutate(TP = c(0, 30, 60, 120, 240)) %>% 
#   dplyr::mutate(class = "TAG")
# 
# ffa_data <- 
#   lipidomics_fc_p %>% 
#   purrr::map(.f = function(x){
#     x %>% 
#       dplyr::filter(variable_id %in% ffa_idx) %>% 
#       dplyr::select(p_value, fc) %>% 
#       apply(2, mean)
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame() %>% 
#   dplyr::mutate(fc = log(fc, 2)) %>% 
#   dplyr::mutate(p_value = -log(p_value, 10)) %>% 
#   dplyr::add_row(p_value = 0, fc = 0, .before = 1) %>% 
#   dplyr::mutate(TP = c(0, 30, 60, 120, 240)) %>% 
#   dplyr::mutate(class = "FFA")
# 
# temp_data =
#   rbind(pc_data,
#         pe_data,
#         dag_data,
#         tag_data,
#         ffa_data)
# plot <- 
#   temp_data %>% 
#   dplyr::mutate(TP = factor(as.character(TP), levels = c("0", "30", "60", "120", "240"))) %>% 
#   dplyr::mutate(class = factor(class, levels = c("PC", "PE", 'DAG', "TAG", "FFA"))) %>% 
#   ggplot(aes(TP, class)) +
#   geom_point(aes(size = p_value, 
#                  fill = fc),
#              shape = 21) +
#   scale_fill_gradient2(low = ggsci::pal_aaas()(n=10)[1], 
#                        high = ggsci::pal_aaas()(n=10)[2], 
#                        mid = "white") +
#   scale_size_continuous(range = c(1, 20)) +
#   theme_bw() +
#   labs(x = "Time point (min)", y = "") +
#   theme(
#     # legend.position = c(0, 1),
#     # legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA)
#   ) 
# 
# 
# # ggsave(plot, filename = "lipid_change.pdf", width = 9, height = 5)
# 
# 
# ###lipid changes according to carbon and un number
# sxtTools::setwd_project()
# lipid_info <- read.table("data/shake_study/lipidomics_data_analysis/DEG/Lipomat05.txt", header = TRUE, sep = "\t")
# load("data/shake_study/lipidomics_data_analysis/data_preparation/expression_data")
# lipidomics_expression_data <-
#   expression_data
# setwd("data/shake_study/3_omics/k_means_clustering/")
# 
# tag_variable_info <-
#   variable_info2[match(cluster1_lipid$variable_id, variable_info2$variable_id), ] %>%
#   dplyr::filter(subclass == 'TAG') %>%
#   dplyr::left_join(lipid_info, by = c("mol_name" = "Lipid_Name"))
# 
# 
# plot(tag_variable_info$Total_Carb,
#      tag_variable_info$Total_Unsat)
# 
# temp_data1 <- 
#   tag_variable_info[,c("Total_Carb", "Total_Unsat")] %>% 
#   dplyr::mutate(variable_id = paste(Total_Carb, Total_Unsat, sep = "_")) %>% 
#   dplyr::group_by(variable_id) %>% 
#   dplyr::summarise(n = n()) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(Total_Carb = stringr::str_split(variable_id, "_") %>% 
#                   lapply(function(x)x[1])
#                 %>% unlist()) %>% 
#   dplyr::mutate(Total_Unsat = stringr::str_split(variable_id, "_") %>% 
#                   lapply(function(x)x[2])
#                 %>% unlist()) 
# 
# plot1 <- 
#   temp_data1 %>% 
#   ggplot(aes(x = Total_Carb, y = Total_Unsat)) +
#   geom_point(aes(size = n), 
#              shape = 21,
#              fill = ggsci::pal_aaas()(n=10)[9]) +
#   scale_size_continuous(range = c(5, 10)) +
#   theme_bw() +
#   labs(x = "Number of carbons", y = "Number of unsaturations") +
#   theme(
#     legend.position = c(0, 1),
#     legend.justification = c(0, 1),
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 13),
#     axis.text = element_text(size = 12),
#     panel.background = element_rect(fill = "transparent", color = NA),
#     plot.background = element_rect(fill = "transparent", color = NA),
#     legend.background = element_rect(fill = "transparent", color = NA),
#     plot.margin = unit(x = c(0, 0, 0, 0),
#                        units = "cm")
#   ) 
# 
# temp_data2 <- 
#   tag_variable_info$Total_Carb %>% 
#   unique() %>% 
#   purrr::map(function(x){
#     temp_index <- 
#       tag_variable_info %>% 
#       dplyr::filter(Total_Carb == x) %>% 
#       dplyr::pull(variable_id)
#     
#     temp_data <- 
#       lipidomics_expression_data[temp_index,] %>% 
#       `+`(1) %>% 
#       log(2) %>% 
#       apply(1, function(x){
#         x / sd(x)
#       }) %>% 
#       apply(1, mean)
#     
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# rownames(temp_data2) <- 
#   tag_variable_info$Total_Carb %>% 
#   unique()
# 
# library(plyr)
# 
# 
# ##mean value
# temp_data2_mean <- 
#   temp_data2 %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   purrr::map(.f = function(x){
#     data.frame(id =  colnames(temp_data2), x) %>% 
#       dplyr::mutate(TP = stringr::str_split(id, "_") %>% 
#                       lapply(function(z)z[2]) %>% 
#                       unlist()) %>% 
#       dplyr::mutate(subject_id = stringr::str_split(id, "_") %>% 
#                       lapply(function(z)z[1]) %>% 
#                       unlist()) %>% 
#       plyr::dlply(.variables = .(TP)) %>% 
#       purrr::map(function(y){
#         y[,"x"] %>% 
#           mean()
#       }) %>% 
#       unlist()
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# plot2 <- 
#   temp_data2_mean %>% 
#   apply(1, function(x){
#     (x - mean(x)) / sd(x)
#   }) %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "carbon_number") %>% 
#   tidyr::pivot_longer(cols = -carbon_number, names_to = "TP", values_to = "value") %>% 
#   dplyr::mutate(TP = factor(TP, levels = stringr::str_sort(unique(TP), numeric = TRUE))) %>% 
#   ggplot(aes(x = carbon_number, y = TP)) +
#   geom_tile(aes(fill = value), color = "white", show.legend = FALSE) +
#   scale_fill_gradient2(low = ggsci::pal_aaas()(n=10)[1],
#                        high = ggsci::pal_aaas()(n=10)[2],
#                        mid = "white") +
#   scale_x_discrete(expand = expansion(mult = c(0,0))) +
#   scale_y_discrete(expand = expansion(mult = c(0,0))) +
#   theme_bw() +
#   labs(y = "", x = "") +
#   theme(
#     panel.grid = element_blank(),
#     legend.position = "top",
#     axis.ticks = element_blank(),
#     axis.text.x = element_blank(),
#     plot.margin = unit(x = c(0, 0, 0, 0),
#                        units = "cm")
#   )
# 
# temp_data3 <- 
#   tag_variable_info$Total_Unsat %>% 
#   unique() %>% 
#   purrr::map(function(x){
#     temp_index <- 
#       tag_variable_info %>% 
#       dplyr::filter(Total_Unsat == x) %>% 
#       dplyr::pull(variable_id)
#     
#     temp_data <- 
#       lipidomics_expression_data[temp_index,] %>% 
#       `+`(1) %>% 
#       log(2) %>% 
#       apply(1, function(x){
#         x / sd(x)
#       }) %>% 
#       apply(1, mean)
#     
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# rownames(temp_data3) <- 
#   tag_variable_info$Total_Unsat %>% 
#   unique()
# 
# library(plyr)
# 
# 
# ##mean value
# temp_data3_mean <- 
#   temp_data3 %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   purrr::map(.f = function(x){
#     data.frame(id =  colnames(temp_data3), x) %>% 
#       dplyr::mutate(TP = stringr::str_split(id, "_") %>% 
#                       lapply(function(z)z[2]) %>% 
#                       unlist()) %>% 
#       dplyr::mutate(subject_id = stringr::str_split(id, "_") %>% 
#                       lapply(function(z)z[1]) %>% 
#                       unlist()) %>% 
#       plyr::dlply(.variables = .(TP)) %>% 
#       purrr::map(function(y){
#         y[,"x"] %>% 
#           mean()
#       }) %>% 
#       unlist()
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# plot3 <- 
#   temp_data3_mean %>% 
#   apply(1, function(x){
#     (x - mean(x)) / sd(x)
#   }) %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   tibble::rownames_to_column(var = "un_number") %>% 
#   tidyr::pivot_longer(cols = -un_number, names_to = "TP", values_to = "value") %>% 
#   dplyr::mutate(TP = factor(TP, levels = stringr::str_sort(unique(TP), numeric = TRUE))) %>% 
#   ggplot(aes(y = un_number, x = TP)) +
#   geom_tile(aes(fill = value), color = "white") +
#   scale_fill_gradient2(low = ggsci::pal_aaas()(n=10)[1],
#                        high = ggsci::pal_aaas()(n=10)[2],
#                        mid = "white") +
#   scale_x_discrete(expand = expansion(mult = c(0,0))) +
#   scale_y_discrete(expand = expansion(mult = c(0,0))) +
#   theme_bw() +
#   labs(y = "", x = "") +
#   theme(
#     panel.grid = element_blank(),
#     legend.position = "top",
#     axis.ticks = element_blank(),
#     axis.text.y = element_blank(),
#     plot.margin = unit(x = c(0, 0, 0, 0),
#                        units = "cm")
#   )
# 
# library(patchwork)
# plot2 + plot1 + patchwork::plot_layout(ncol = 1, heights = c(1,2))
# plot1 + plot3 + patchwork::plot_layout(widths = c(2,1))
# 
# 
# 
# 
# plot <-
#   {
#     plot2 + plot2 + plot_layout(ncol = 2, widths = c(2, 1))
#   } -
#   {
#     plot1 + plot3 + plot_layout(ncol = 2, widths = c(2, 1))
#   } +
#   plot_layout(ncol = 1, heights = c(1, 2))
# 
# plot
# 
# 
# ggsave(plot, filename = "lipid_carbon_number.pdf", width = 8, height = 7)