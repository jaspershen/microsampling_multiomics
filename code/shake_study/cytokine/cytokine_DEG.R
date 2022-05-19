#to avoid source
no_exist_function()

masstools::setwd_project()
rm(list = ls())
library(tidyverse)
source("code/tools.R")

##load data
load("data/shake_study/cytokine_data_analysis/data_preparation/expression_data")
load("data/shake_study/cytokine_data_analysis/data_preparation/sample_info")
load("data/shake_study/cytokine_data_analysis/data_preparation/variable_info")

load("data/shake_study/subject_info/subject_info")

sample_info <-
  sample_info %>% 
  dplyr::left_join(subject_info, by = "subject_id")

masstools::setwd_project()
setwd("data/shake_study/cytokine_data_analysis/DEG")

dim(expression_data)

sample_info <- 
  sample_info %>% 
  dplyr::arrange(sample_id, TP)

expression_data <- 
  expression_data %>% 
  dplyr::select(sample_info$sample_id)

## for each person, just combine the samples in the same TP range
###combine different samples in one time together
library(plyr)

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
# 
# anova_p <-
# purrr::map(as.data.frame(t(subject_data)), .f = function(x){
#   temp_data <-
#     data.frame(
#       subject_id = sample_info$subject_id,
#       tp = factor(
#         sample_info$TP,
#                   levels = stringr::str_sort(unique(sample_info$TP),
#                                              numeric = TRUE)),
#       x = x,
#       stringsAsFactors = FALSE
#     )
# 
#   intersect_name <-
#   temp_data %>%
#     plyr::dlply(.variables = .(tp)) %>%
#     purrr::map(.f = function(x){
#       x$subject_id
#     }) %>%
#     Reduce(intersect, .) %>%
#     unique()
# 
#   temp_data <-
#     temp_data %>%
#     dplyr::filter(subject_id %in% intersect_name)
# 
#   temp_data$subject_id <-
#     factor(temp_data$subject_id, levels = unique(temp_data$subject_id))
# 
#   res.aov <-
#     anova_test(
#       data = temp_data,
#       dv = x,
#       wid = subject_id,
#       within = tp
#     )
# 
#   p <-
#     as.data.frame(get_anova_table(res.aov))$p
# 
#   # pairwise comparisons
#   pwc <- temp_data %>%
#     pairwise_t_test(x ~ tp,
#                     paired = TRUE,
#                     p.adjust.method = "fdr")
#   pwc
# 
#   p2 <-
#   pwc[,c(2,3,8,9)]
# 
#   name <- paste(pwc$group2, pwc$group1, sep = "_")
# 
#   p2 <- matrix(pwc$p, nrow = 1)
# 
#   colnames(p2) <- name
# 
#   p <-
#     data.frame(p, p2, stringsAsFactors = FALSE, check.names = FALSE)
# 
#   p
# 
#   })
# 
# anova_p <- anova_p %>% do.call(rbind, .) %>% as.data.frame()
# 
# save(anova_p, file = "anova_p")
load("anova_p")

anova_marker_name <-
  rownames(anova_p)[which(anova_p$p < 0.05)]

# save(anova_marker_name, file = "anova_marker_name")

load("anova_marker_name")


# ###permutation test
# permutation_marker_number <- 
#   purrr::map(1:100, function(i){
#     cat(i, "")
#     anova_p <-
#       purrr::map(as.data.frame(t(subject_data)), .f = function(x){
#         temp_data <-
#           data.frame(
#             subject_id = sample_info$subject_id,
#             tp = factor(
#               sample_info$TP,
#               levels = stringr::str_sort(unique(sample_info$TP),
#                                          numeric = TRUE)),
#             x = x,
#             stringsAsFactors = FALSE
#           )
#         
#         intersect_name <-
#           temp_data %>%
#           plyr::dlply(.variables = .(tp)) %>%
#           purrr::map(.f = function(x){
#             x$subject_id
#           }) %>%
#           Reduce(intersect, .) %>%
#           unique()
#         
#         temp_data <-
#           temp_data %>%
#           dplyr::filter(subject_id %in% intersect_name)
#         
#         temp_data$subject_id <-
#           factor(temp_data$subject_id, levels = unique(temp_data$subject_id))
#         
#         temp_data$x = sample(temp_data$x)
#         
#         res.aov <-
#           anova_test(
#             data = temp_data,
#             dv = x,
#             wid = subject_id,
#             within = tp
#           )
#         
#         p <-
#           as.data.frame(get_anova_table(res.aov))$p
#         
#         p <-
#           data.frame(p, stringsAsFactors = FALSE, check.names = FALSE)
#         p
#       }) %>% 
#       do.call(rbind, .) %>% 
#       as.data.frame()
#     
#     sum(anova_p$p < 0.05)
#   }) %>% 
#   unlist()
# 
# save(permutation_marker_number, file = "permutation_marker_number")
load(permutation_marker_number, file = "permutation_marker_number")
sum(permutation_marker_number > 7)/100

###calculate the FC 30/0
idx0 <- which(stringr::str_detect(colnames(subject_data), "T0"))
idx30 <- which(stringr::str_detect(colnames(subject_data), "T30"))

fc <- 
  subject_data %>% 
  apply(1,function(x){
    x <- as.numeric(x)
    mean(x[idx30])/mean(x[idx0])
  })

names(fc) == rownames(anova_p)

##volcano plot
plot <- 
  volcano_plot(
    fc = fc,
    p_value = anova_p$p,
    p.cutoff = 0.05,
    fc.cutoff = 1,
    alpha = 1,
    point.size = 5,
    text = TRUE, 
    variable_id = names(fc)
  )

# plot <- 
#   plot +
#   scale_x_continuous(limits = c(-0.3, 0.3)) +
#   scale_y_continuous(limits = c(0, 8))

plot

# ggsave(plot, filename = "volcano_plot.pdf", width = 5, height = 7)


subject_data <- 
  apply(subject_data, 1, function(x){
    (x) / sd(x)
  })

subject_data2 <- 
  subject_data %>% 
  data.frame(., TP = sample_info$TP, 
             subject_id = sample_info$subject_id,
             stringsAsFactors = FALSE, check.names = FALSE) %>% 
  mutate(TP = factor(TP,levels = sample_info$TP %>% 
                       unique() %>% 
                       stringr::str_sort(numeric = TRUE))) %>% 
  plyr::dlply(.variables = .(TP))

# temp_data <-
#   expression_data %>%
#   dplyr::select(one_of(sample_info$sample_id))
# 
# ##log transformation
# temp_data <-
#   log(temp_data + 1, 2) %>%
#   t()
#   # apply(1, function(x) {
#   #   (x - mean(x)) / sd(x)
#   # })
# 
# temp_data2 <-
#   temp_data %>%
#   data.frame(., TP = sample_info$TP,
#              subject_id = sample_info$subject_id,
#              stringsAsFactors = FALSE,
#              check.names = FALSE) %>%
#   mutate(TP = factor(TP, levels = sample_info$TP %>%
#                        unique() %>%
#                        stringr::str_sort(numeric = TRUE))) %>%
#   plyr::dlply(.variables = .(TP))
# 
# temp_data2 <-
# temp_data2 %>%
#   purrr::map(.f = function(x){
#     x <-
#       x %>%
#       tibble::column_to_rownames(var = "subject_id") %>%
#       dplyr::select(-c(TP))
#     baseline <- temp_data2[[1]] %>%
#       tibble::column_to_rownames(var = "subject_id") %>%
#       dplyr::select(-c(TP))
#     x - baseline
#   })
# 
# subject_data_mean <-
#   lapply(temp_data2, function(x) {
#     apply(x, 2, function(y){
#       mean(y)
#     })
#   }) %>%
#   do.call(cbind, .)
# 
# subject_data_sd <-
#   lapply(temp_data2, function(x) {
#     apply(x, 2, function(y){
#       sd(y)
#     })
#   }) %>%
#   do.call(cbind, .)
# 
# subject_data_sem <-
#   lapply(temp_data2, function(x) {
#     apply(x, 2, function(y) {
#       sd(y) / sqrt(nrow(x))
#     })
#   }) %>%
#   do.call(cbind, .)
# 
# save(subject_data_mean, file = "subject_data_mean")
# save(subject_data_sd, file = "subject_data_sd")
# save(subject_data_sem, file = "subject_data_sem")

load("subject_data_mean")
load("subject_data_sd")
load("subject_data_sem")

subject_data2 <-
  lapply(subject_data2, function(x) {
    x <-
      x %>%
      tibble::column_to_rownames(var = "subject_id") %>% 
      dplyr::select(-TP)
  })

# ##find all the peaks in different time points
# fc_p_value <-
#   pbapply::pblapply(subject_data2[-1], function(x) {
#     y <- subject_data2[[1]]
#     intersect_name <- intersect(rownames(x), rownames(y))
#     p_value <- lapply(1:ncol(x), function(idx) {
#       wilcox.test(x[intersect_name, idx],
#                   subject_data2[[1]][intersect_name, idx],
#                   paired = TRUE)$p.value
#     }) %>%
#       unlist() %>%
#       p.adjust(method = "fdr")
# 
#     fc <- lapply(1:ncol(x), function(idx) {
#       mean(x[, idx]) / mean(subject_data2[[1]][, idx])
#     }) %>%
#       unlist()
# 
#     fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
# 
#     data.frame(p_value,
#                fc,
#                variable_id = variable_info$variable_id,
#                stringsAsFactors = FALSE)
#   })
# 
# save(fc_p_value, file = "fc_p_value")

load("fc_p_value")

names(fc_p_value)

dir.create("marker_in_different_points")

##find markers for each time points
marker_each_point <- 
  lapply(fc_p_value, function(x){
    idx1 <- which(x$p_value < 0.05 & x$fc > 1)
    idx2 <- which(x$p_value < 0.05 & x$fc < 1)
    
    gene1 <- 
      try(
        data.frame(variable_id = variable_info$variable_id[idx1],
                   x[idx1,],
                   class = "increase",
                   stringsAsFactors = FALSE
        ),silent = TRUE 
      )
    
    if(class(gene1) == "try-error"){
      gene1 <- NULL
    }
    
    gene2 <- 
      try(
        data.frame(variable_id = variable_info$variable_id[idx2],
                   x[idx2,],
                   class = "decrease",
                   stringsAsFactors = FALSE
        ),silent = TRUE
      )
    
    if(class(gene2) == "try-error"){
      gene2 <- NULL
    }
    
    rbind(gene1, gene2)
  })

marker_each_point[[1]]

names(marker_each_point)

# save(marker_each_point, file = "marker_each_point")
load("marker_each_point")

#####a sankey 
marker_each_point %>% 
  lapply(nrow) %>% 
  unlist()

all_marker_name <- 
  lapply(marker_each_point, function(x){
    x$variable_id
  }) %>% 
  unlist() %>% 
  unique()

length(all_marker_name)

library(ggalluvial)

temp_data <- 
  lapply(marker_each_point, function(x){
    if(is.null(x)){
      return(NULL)
    }
    x <- 
      data.frame(variable_id = all_marker_name,
                 stringsAsFactors = FALSE) %>% 
      left_join(x, by = "variable_id") %>% 
      dplyr::select(variable_id, class)
    
    x$class[is.na(x$class)] <- "no"
    x$freq <- 1
    x
    
  })

temp_data <-
  purrr::map2(.x = temp_data, .y = names(temp_data), .f = function(x,y){
    if(is.null(x)){
      return(NULL)
    }
    data.frame(x, point = y, stringsAsFactors = FALSE)
  })

temp_data <- 
  do.call(rbind, temp_data)

temp_data$point <- 
  factor(temp_data$point, levels = unique(temp_data$point))

RColorBrewer::display.brewer.all()

plot1 <- 
  ggplot(temp_data,
         aes(x = point, 
             y = freq,
             stratum = class, 
             alluvium = variable_id,
             fill = class, 
             label = class)) +
  scale_x_discrete(expand = c(.1, .1)) +
  ggalluvial::geom_flow() +
  labs(x = "", y = "") +
  scale_fill_manual(values = c(
    "increase" = ggsci::pal_d3()(10)[2],
    "decrease" = ggsci::pal_d3()(10)[1],
    "no" = "azure2"
  )) +
  ggalluvial::geom_stratum(alpha = 1, color = "black") +
  # geom_text(stat = "stratum", size = 3) +
  theme_bw() +
  theme(legend.position = "top", 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot1

# ggsave(
#   plot1,
#   file = file.path("marker_in_different_points", "gene_sankey_light.pdf"),
#   width = 14,
#   height = 7,
#   bg = "transparent"
# )
# 
# ggsave(
#   plot1,
#   file = file.path("marker_in_different_points", "gene_sankey_light.png"),
#   width = 14,
#   height = 7,
#   bg = "transparent"
# )

length(all_marker_name)

# save(all_marker_name, file = "all_marker_name")

length(all_marker_name)


####cytokine information
temp_data <- 
  purrr::map(subject_data2, .f = function(x){
    x[,anova_marker_name]
  })

class <- 
  variable_info$subclass[match(anova_marker_name, variable_info$variable_id)]

temp_data <-
  temp_data %>% 
  purrr::map(function(x){
    # x <-
      t(x) %>% 
      as.data.frame() %>% 
      data.frame(class, .) %>% 
      plyr::dlply(.variables = .(class)) %>% 
      purrr::map(function(y){
        y %>% 
          dplyr::select(-class) %>% 
          colMeans()
      }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
  })

temp_data_p_fc <- 
  purrr::map(temp_data[-1], function(x){
    y <- temp_data[[1]]
    intersect_name <- intersect(colnames(x), colnames(y))
    p_value <- lapply(1:nrow(x), function(idx) {
      wilcox.test(as.numeric(x[idx, intersect_name]),
                  as.numeric(temp_data[[1]][idx, intersect_name]),
                  paired = TRUE)$p.value
    }) %>%
      unlist() %>%
      p.adjust(method = "fdr")
    
    fc <- lapply(1:nrow(x), function(idx) {
      mean(as.numeric(x[idx,])) / mean(as.numeric(y[idx,]))
    }) %>%
      unlist()
    
    fc[is.infinite(fc)] <- max(fc[!is.infinite(fc)])
    
    data.frame(variable_id = rownames(x),
               p_value,
               fc,
               stringsAsFactors = FALSE)
  })

temp_data <- 
  purrr::map2(
    .x = names(temp_data_p_fc),
    .y = temp_data_p_fc,
    .f = function(x, y) {
      data.frame(TP = x, y, stringsAsFactors = FALSE)
    }
  ) %>%
  do.call(rbind, .) %>%
  as.data.frame()

plot <- 
temp_data %>% 
  dplyr::mutate(TP = factor(x = TP, levels = unique(TP))) %>% 
  dplyr::mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.1 & p_value > 0.001 ~ "**",
    p_value < 0.05 & p_value > 0.01 ~ "*",
    p_value > 0.05 ~ "N.S."
  )) %>% 
  ggplot(aes(TP, variable_id)) +
    geom_point(aes(size = -log(p_value, 10), 
                   fill = log(fc, 2)), shape = 21) +
  geom_text(aes(label = sig)) +
  scale_size_continuous(range = c(1, 25)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_color_manual(values = c("Yes" = "blue", "No" = "green")) +
  theme_bw() +
  labs(x = "TP", y = "") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        panel.grid.minor = element_blank())

plot

# ggsave(plot, filename = "class_change.pdf", width = 10, height = 7)
# ggsave(plot, filename = "class_change.png", width = 10, height = 7)

