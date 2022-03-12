##
no_function()

sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

###load data
load("data/shake_study/cytokine_data_analysis/data_preparation/expression_data")
load("data/shake_study/cytokine_data_analysis/data_preparation/sample_info")
load("data/shake_study/cytokine_data_analysis/data_preparation/variable_info")

load("data/shake_study/cytokine_data_analysis/DEG/anova_marker_name")

sxtTools::setwd_project()
setwd("data/shake_study/cytokine_data_analysis/k_means_clustering")

expression_data <- 
  expression_data[anova_marker_name,]

dim(expression_data)
dim(sample_info)
dim(variable_info)

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
  temp_data[,temp_idx] %>% 
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
  apply(1, function(x) (x - mean(x))/sd(x)) %>% 
  t() %>% 
  as.data.frame()

time <- c(0, 30, 60, 120, 240)

temp_data <- rbind(time, temp_data)

row.names(temp_data)[1] <- "time"

write.table(
  temp_data,
  file = "temp_data.txt",
  sep = '\t',
  quote = FALSE,
  col.names = NA
)

#read it back in as an expression set
data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

plot <-
Dmin(data.s, m=m1, crange=seq(2,6,1), repeats=3, visu=TRUE)

plot <-
plot %>%
  data.frame(distance = plot,
             k = seq(2,6,1)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21, size = 4, fill = "black") +
  geom_smooth() +
  geom_segment(aes(x = k, y = 0, xend = k, yend = distance)) +
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
  labs(
    x = "Cluster number",
    y = "Min. centroid distance"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))

plot

ggsave(plot, filename = "distance_k_number.pdf", width = 7, height = 7)

cluster = 2

c <- mfuzz(data.s, c = cluster, m = m1)

mfuzz.plot(eset = data.s,
           # min.mem = 0.6,
           cl = c,
           mfrow=c(4,4),
           time.labels = time,
           new.window = FALSE)

save(c, file = "c")
load("c")

####any two clusters with correlation > 0.8 should be considered as one
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
  col = colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(n=100),
  number.cex = .7, 
  addCoef.col = "black"
)

###cor_plot_cluster

center %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "cluster") %>% 
  tidyr::pivot_longer(cols = -cluster, names_to = "point", values_to = "value") %>% 
  # dplyr::filter(cluster %in% c("Cluster 2", "Cluster 4")) %>% 
  dplyr::mutate(point = factor(point, levels = unique(point))) %>% 
  ggplot(aes(point, value, group = cluster)) +
  geom_point(aes(color = cluster)) +
  geom_line(aes(group = cluster, color = cluster)) 
# geom_smooth(aes(group = cluster), se = FALSE)

library(ComplexHeatmap)

plot <- 
  Heatmap(
    matrix = center,
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D",
    name = "Z-score",
    border = TRUE,
    rect_gp = gpar(col= "white")
  )


plot <- ggplotify::as.ggplot(plot)
plot
ggsave(plot, filename = "cluster_heatmap.pdf", width = 9, height = 7)

###
cluster_color <- 
  colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = "Spectral"))(8)

names(cluster_color) <- as.character(1:8)

plot <- 
  center %>%
  as.data.frame() %>%
  tibble::rowid_to_column(var = "cluster") %>%
  tidyr::pivot_longer(cols = -cluster,
                      names_to = "time",
                      values_to = "value") %>%
  dplyr::mutate(time = factor(time, levels = unique(time))) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  ggplot(aes(time, value)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(group = cluster, color = cluster), show.legend = FALSE) +
  # geom_point(aes(group = cluster, fill = cluster), show.legend = FALSE, shape = 21, size = 3) +
  # geom_smooth(aes(color = cluster, group = cluster), se = FALSE) +
  scale_color_manual(values = cluster_color) +
  scale_fill_manual(values = cluster_color) +
  # facet_grid(rows = vars(class)) +
  theme_bw() +
  labs(x = "", y = "Z-score") +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
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
  )

plot

ggsave(plot, filename = "cluster_center_plot.pdf", width = 7, height = 4)

centers <- c$centers

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

xlsx::write.xlsx(cluster_info,
                 "cluster_info.xlsx",
                 row.names = FALSE)

####plot for each cluster
idx <- 3
for(idx in 1:2) {
  cat(idx, " ")
  cluster_data <-
    cluster_info %>%
    dplyr::filter(cluster == idx) %>%
    dplyr::select(1, 1 + idx, cluster)
  
  colnames(cluster_data)[2] <- c("membership")
  
  cluster_data <-
    cluster_data %>%
    dplyr::filter(membership > 0.5)
  path <- paste("cluster", idx, sep = "_")
  dir.create(path)
  
  xlsx::write.xlsx(cluster_data,
                   file = file.path(path, paste("cluster", idx, ".xlsx", sep = "")),
                   row.names = FALSE)
  
  temp_center <-
    centers[idx, , drop = TRUE] %>%
    data.frame(time = names(.),
               value = .,
               stringsAsFactors = FALSE) %>%
    dplyr::mutate(time = factor(time, levels = time))
  
  plot <-
    temp_data[cluster_data$variable_id, ] %>%
    data.frame(
      membership = cluster_data$membership,
      .,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) %>%
    tibble::rownames_to_column(var = "variable_id") %>%
    tidyr::pivot_longer(
      cols = -c(variable_id, membership),
      names_to = "time",
      values_to = "value"
    ) %>%
    dplyr::mutate(time = factor(time, levels = unique(time))) %>%
    ggplot(aes(time, value, group = variable_id)) +
    geom_line(aes(color = membership)) +
    scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
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
        "Cluster ",
        idx,
        " (",
        nrow(cluster_data),
        " metabolic peaks)",
        sep = ""
      )
    ) +
    geom_line(
      mapping = aes(time, value, group = 1),
      data = temp_center,
      size = 2
    )
  
  plot

  ggsave(plot, filename = file.path(path, paste("cluster",idx, ".pdf", sep = "")),
         width = 8, height = 7)

}






###annotation for each cluster
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx") 
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx") 

nrow(cluster1) +
  nrow(cluster2) 

##metabolites that not in clusters
non_metabolite <- 
  setdiff(cluster_info$variable_id,
          c(cluster1$variable_id,
            cluster2$variable_id))

plot <-
  temp_data[non_metabolite,] %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  tidyr::pivot_longer(cols = -c(variable_id),
                      names_to = "time", 
                      values_to = "value") %>%
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
    title = paste("Unclustered metabolic peaks (", length(non_metabolite), " metabolic peaks)", sep = "")
  ) 
  
plot

ggsave(plot, filename = "non_cluster_metabolites.pdf",
       width = 8, height = 7)

dim(cluster1)
dim(cluster2)

###heatmap for all the samples
library(ComplexHeatmap)
temp_data <-
  temp_data_mean[c(
    cluster1$variable_id,
    cluster2$variable_id
  ), ]

temp_data <- 
  apply(temp_data, 1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

rownames(temp_data) <-
  variable_info$mol_name[match(rownames(temp_data), variable_info$variable_id)]

range(temp_data)
library(circlize)
pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

range(temp_data)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c(pal[1], "white", pal[100]))

cluster <-
  c(
    rep(1, nrow(cluster1)),
    rep(2, nrow(cluster2))
  )

cluster_color <-
  c(ggsci::pal_d3()(n = 10)[1:2])

names(cluster_color) <- c(as.character(1:2))

text_color <- 
  cluster_color[cluster]

ha <-
  rowAnnotation(
    Cluster = cluster,
    col = list(
      Cluster = cluster_color
    )
    # annotation_name_side = c("left")
  )

plot <-
  Heatmap(
    temp_data,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    border = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10, 
                        col = text_color),
    col = col_fun,
    right_annotation = ha,
    name = "Z score",
    column_names_rot = 0, rect_gp = gpar(col = "white")
    # row_split = factor(as.character(cluster), levels = c(6,3,8,1,4,5,2,9,7))
  )

plot <- ggplotify::as.ggplot(plot)

plot

ggsave(plot = plot, filename = "heatmap_for_all_cluster.pdf", width = 7, height = 7)
ggsave(plot = plot, filename = "heatmap_for_all_cluster.png", width = 7, height = 7)

###functional annotation for different cluster
cluster1 %>% 
  dplyr::left_join(variable_info, by = "variable_id")

cluster2 %>% 
  dplyr::left_join(variable_info, by = "variable_id")
