#to avoid source
no_exist_function()

sxtTools::setwd_project()
rm(list = ls())
library(tidyverse)

##load data
load("data/shake_study/metabolome_data_analysis/data_preparation/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/variable_info")

load("data/shake_study/subject_info/subject_info")

sample_info <- 
sample_info %>% 
  dplyr::left_join(subject_info, by = "subject_id")

##only remain metabolites
variable_info <-
  variable_info %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level != 3)

expression_data <- 
  expression_data[variable_info$variable_id,]

sxtTools::setwd_project()
setwd("data/shake_study/metabolome_data_analysis/DEG")

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

subject_data <-
  t(subject_data)
  # apply(subject_data, 1, function(x){
  #   (x) / sd(x)
  # })

subject_data2 <- 
  subject_data %>% 
  data.frame(., TP = sample_info$TP, 
             subject_id = sample_info$subject_id,
             stringsAsFactors = FALSE, check.names = FALSE) %>% 
  mutate(TP = factor(TP, levels = sample_info$TP %>% 
                       unique() %>% 
                       stringr::str_sort(numeric = TRUE))) %>% 
  plyr::dlply(.variables = .(TP))

subject_data_mean <-
  lapply(subject_data2, function(x) {
    x <-
      x %>%
      dplyr::select(-c(TP, subject_id))
    apply(x, 2, function(y){
      y <- 
        (y - mean(y))/sd(y)
      mean(y)
    })
  }) %>%
  do.call(cbind, .)

subject_data_sd <-
  lapply(subject_data2, function(x) {
    x <-
      x %>%
      dplyr::select(-c(TP, subject_id))
    apply(x, 2, function(y){
      y <- 
        (y - mean(y))/sd(y)
      sd(y)
    })
  }) %>%
  do.call(cbind, .)

subject_data_sem <-
  lapply(subject_data2, function(x) {
    x <-
      x %>%
      dplyr::select(-c(TP, subject_id))
    apply(x, 2, function(y) {
      y <- (y - mean(y))/sd(y)
      sd(y) / sqrt(nrow(x))
    })
  }) %>%
  do.call(cbind, .)

save(subject_data_mean, file = "subject_data_mean")
save(subject_data_sd, file = "subject_data_sd")
save(subject_data_sem, file = "subject_data_sem")

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
#       wilcox.test(x[, idx],
#                   subject_data2[[1]][, idx],
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
#     data.frame(variable_id = variable_info$variable_id,
#                p_value,
#                fc,
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

save(marker_each_point, file = "marker_each_point")
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

ggsave(
  plot1,
  file = file.path("marker_in_different_points", "gene_sankey_light.pdf"),
  width = 14,
  height = 7,
  bg = "transparent"
)

ggsave(
  plot1,
  file = file.path("marker_in_different_points", "gene_sankey_light.png"),
  width = 14,
  height = 7,
  bg = "transparent"
)

length(all_marker_name)

save(all_marker_name, file = "all_marker_name")

length(all_marker_name)




















