no_function()

library(tidyverse)
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")

####load data all the omics data
{
  ###
  load(
    here::here(
      "data/24_7_study/combine_omics/data_preparation/new_expression_data"
    )
  )
  load(here::here(
    "data/24_7_study/combine_omics/data_preparation/new_sample_info"
  ))
  load(
    here::here(
      "data/24_7_study/combine_omics/data_preparation/new_variable_info"
    )
  )
  expression_data = new_expression_data
  sample_info = new_sample_info
  variable_info = new_variable_info
}

####load consistence score for each molecules
load(
  here::here(
    "data/24_7_study/circadian_analysis/all_omics/day_consistence/consistence_score"
  )
)

consistence_score$consistence_score1 <-
  consistence_score$consistence_score * 0.5 + consistence_score$new_consistence_score * 0.5

variable_info =
  variable_info %>%
  dplyr::left_join(consistence_score, by = "variable_id")

library(MetaCycle)

setwd("data/24_7_study/circadian_analysis/all_omics_loess_data")

######manual check
load(here::here("data/24_7_study/summary_info/day_night_df"))

####read metadata
metdata =
  readr::read_csv(here::here(
    "data/24_7_study/raw_data_from_box/sample_registration.csv"
  ))

milk_time =
  metdata %>%
  dplyr::filter(!is.na(food)) %>%
  dplyr::filter(stringr::str_detect(food, "milk|Milk")) %>%
  pull(date_time)

grep("Salicylic", variable_info$mol_name, value = TRUE)
grep("Salicylic", variable_info$mol_name, value = FALSE)

time_plot(
  x = as.numeric(expression_data[1233,]),
  time = sample_info$accurate_time,
  day_night_df = day_night_df,
  add_point = TRUE,
  x_name = "metabolite",
  y_axis_name = "Scaled intensity",
) +
  geom_vline(xintercept = lubridate::as_datetime(milk_time, tz = "PDT"))

plot =
  plot +
  geom_vline(xintercept = lubridate::as_datetime(coffee_time, tz = "PDT"))

###scale for omics
dim(expression_data)

####here we should scale for each day
expression_data =
  unique(sample_info$day) %>%
  purrr::map(function(day) {
    temp_expression_data =
      expression_data[, which(sample_info$day == day)]
    
    temp_expression_data =
      temp_expression_data %>%
      apply(1, function(x) {
        (x - mean(x)) / sd(x)
      }) %>%
      t() %>%
      as.data.frame()
    
    temp_expression_data
    
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

expression_data =
  expression_data[, sample_info$sample_id]

###reference https://abego.cn/2019/05/31/the-rule-of-gene-expression-in-the-day-and-nigth/
###https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html

sample_info$time

data.frame(
  time = sample_info$time,
  value = as.numeric(expression_data["lipid_135",]),
  day = as.character(sample_info$day)
) %>%
  ggplot(aes(time, value)) +
  geom_line(aes(group = day, color = day)) +
  # facet_wrap(facets = vars(day), ncol = 1) +
  geom_point(aes(color = day))

##temp_data is the expression data, one row is one variable, and one column is
##one sample (time point)
temp_data =
  expression_data %>%
  tibble::rownames_to_column(var = "variable_id")

sample_info$day

###as,numeric transfer time to seconds of this day.
time_point =
  as.numeric(sample_info$time) / (60 * 60)

plot(time_point)

colnames(temp_data)[-1] = time_point

# write.csv(temp_data, "temp_data.csv", row.names = FALSE)
# meta2d(
#     infile = "temp_data.csv",
#     filestyle = "csv",
#     outdir = ".",
#     timepoints = time_point,
#     # cycMethod = "LS",
#     outputFile = TRUE,
#     outRawData = FALSE
#   )

####read result
result1 = readr::read_csv("LSresult_temp_data.csv")
result2 = readr::read_csv("meta2d_temp_data.csv")

plot =
  result1 %>%
  dplyr::left_join(variable_info, by = c("CycID" = "variable_id")) %>%
  ggplot(aes(-log(BH.Q, 10), consistence_score1)) +
  geom_point() +
  geom_vline(xintercept = 0.69897, color = "red") +
  geom_hline(yintercept = 0, color = "red") +
  base_theme +
  labs(y = "Consistence score") +
  ggrepel::geom_label_repel(aes(label = ifelse(consistence_score1 > 0.6, mol_name, NA)))

plot

variable_info %>%
  dplyr::filter(consistence_score1 > 0.4)

# ggsave(plot, filename = "circadian_p_vs_consistence_score.pdf", width = 9, height = 7)

result1 %>%
  dplyr::filter(BH.Q < 0.001)

result1 %>%
  dplyr::left_join(variable_info, by = c("CycID" = "variable_id")) %>%
  dplyr::filter(consistence_score1 < 0 & BH.Q < 0.2)

result1$BH.Q

####so the molecules with consistence score < 0 will not be considered for circiadian
####analysis.
cutoff = as.numeric(quantile(variable_info$new_consistence_score, probs = 0.75))

new_result =
  result1 %>%
  dplyr::left_join(variable_info[, c("variable_id", "new_consistence_score", "consistence_score1")],
                   by = c("CycID" = "variable_id")) %>%
  dplyr::filter(new_consistence_score > cutoff)

new_result$BH.Q = p.adjust(new_result$p, method = "BH")

idx = which(new_result$BH.Q < 0.05)
idx

# save(new_result, file = "new_result")
load("new_result")

result <-
  new_result[idx, ]

# write.csv(result, "circidian_marker.csv", row.names = FALSE)

load(here::here("data/24_7_study/summary_info/day_night_df"))

day_night_df =
  data.frame(
    start = hms::hms(
      seconds = 0,
      minutes = 0,
      hours = 6
    ),
    end = hms::hms(
      seconds = 0,
      minutes = 0,
      hours = 18
    )
  )

#####output plot
dir.create("plot")

# purrr::walk(idx, function(temp_idx) {
#     cat(temp_idx, " ")
#     p = new_result$p[temp_idx]
#     bh = new_result$BH.Q[temp_idx]
#     mol_id = new_result$CycID[temp_idx]
#     mol_name = variable_info$mol_name[match(mol_id, variable_info$variable_id)]
#
#     value = as.numeric(expression_data[mol_id, ])
#
#     temp_data =
#       data.frame(
#         time = sample_info$time,
#         value = value,
#         day = as.character(sample_info$day)
#       )
#
#     plot =
#       ggplot() +
#       geom_rect(
#         mapping = aes(
#           xmin = start,
#           xmax = end,
#           ymin = -Inf,
#           ymax = Inf
#         ),
#         fill = "lightyellow",
#         data = day_night_df,
#         show.legend = FALSE
#       ) +
#       geom_point(
#         aes(x = time,
#             y = value,
#             color = day),
#         alpha = 1,
#         size = 2,
#         data = temp_data
#       ) +
#       geom_smooth(
#         aes(x = time,
#             y = value),
#         color = "black",
#         size = 2,
#         se = FALSE,
#         data = temp_data
#       ) +
#       ggsci::scale_color_lancet() +
#       labs(
#         y = "Continuous glucose monitoring",
#         x = "",
#         title = paste(mol_name, ",p:", round(p, 5), ",BH", round(bh, 5))
#       ) +
#       base_theme
#
#     if (bh < 0.2) {
#       plot =
#         plot +
#         stat_smooth(
#           geom = 'line',
#           aes(x = time,
#               y = value,
#               color = day),
#           alpha = 0.3,
#           se = FALSE,
#           data = temp_data
#         )
#     } else{
#       plot =
#         plot +
#         stat_smooth(
#           geom = 'line',
#           aes(x = time,
#               y = value,
#               color = day),
#           alpha = 0.3,
#           se = FALSE,
#           data = temp_data
#         )
#     }
#
#     name =
#       paste(mol_id, "pdf", sep = ".")
#     name = stringr::str_replace(name, "\\/", "_")
#     ggsave(
#       plot,
#       filename = file.path("plot", name),
#       width = 10,
#       height = 7
#     )
#   })

library(tictoc)
# tic()
# temp_data_loess_day =
# purrr::map(idx, function(temp_idx) {
#   cat(temp_idx, " ")
#   p = new_result$p[temp_idx]
#   bh = new_result$BH.Q[temp_idx]
#   mol_id = new_result$CycID[temp_idx]
#   mol_name = variable_info$mol_name[match(mol_id, variable_info$variable_id)]
#
#   value = as.numeric(expression_data[mol_id, ])
#
#   temp_data =
#     data.frame(
#       time = sample_info$time,
#       value = value,
#       day = as.character(sample_info$day)
#     )
#
#   temp_data$new_time = as.numeric(temp_data$time)/(60*60)
#   temp_data =
#     temp_data %>%
#     dplyr::arrange(new_time)
#   optimization_result =
#   optimize_loess_span(x = temp_data$new_time,
#                       y = temp_data$value,
#                       span_range = c(0.4,0.5,0.6,0.7,0.8))
#   span =
#   optimization_result[[1]]$span[which.min(optimization_result[[1]]$rmse)]
#
#   loess_rg =
#   loess(formula = value ~ new_time, span = span, data = temp_data)
#
#   predicted_value =
#     predict(object = loess_rg,
#             newdata = data.frame(new_time = unique(temp_data$new_time)))
#
#   data.frame(variable_id = mol_id,
#     time = unique(temp_data$time),
#              value = unname(predicted_value))
#
# })
#   # dplyr::bind_cols()
# toc()
#
# save(temp_data_loess_day, file = "temp_data_loess_day")
load("temp_data_loess_day")

temp_data = temp_data_loess_day

####heatmap
temp_data =
  temp_data %>%
  dplyr::bind_rows() %>%
  tidyr::pivot_wider(names_from = "time", values_from = "value") %>%
  tibble::column_to_rownames(var = "variable_id")

rownames(temp_data) = variable_info$mol_name[match(rownames(temp_data), variable_info$variable_id)]

sample_info =
  sample_info %>%
  dplyr::arrange(time)

colnames(temp_data) ==
  unique(as.character(sample_info$time))

temp_data =
  temp_data %>%
  apply(1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>%
  t()

library(ComplexHeatmap)
library(circlize)

range(temp_data)

col_fun = colorRamp2(
  breaks = c(-3, 0, 3),
  colors =
    c("#366A9FFF", "white", "red"),
  transparency = 0
)

p = -log(new_result$BH.Q[idx], 10)
bar_col = ifelse(p > 1.30103, "red", "grey")

text_cor =
  class_color[variable_info$data_type[idx]]

# sample_number = as.numeric(table(sample_info$time))

library(factoextra)
fviz_nbclust(temp_data,
             FUN = hcut,
             method = "silhouette",
             k.max = 30)
fviz_nbclust(temp_data,
             FUN = hcut,
             method = "wss",
             k.max = 30)

log_p =
  -log(new_result$BH.Q[idx], 10)
log_p[is.infinite(log_p)] = max(log_p[!is.infinite(log_p)])

plot =
  Heatmap(
    as.matrix(temp_data),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    name = "z-score",
    border = TRUE,
    col = col_fun,
    column_names_rot = 45,
    clustering_method_rows = "ward.D",
    clustering_distance_rows = "euclidean",
    # rect_gp = gpar(col = "white"),
    row_names_gp = gpar(cex = 0.8, col = text_cor),
    column_names_gp = gpar(cex = 0.8),
    # row_km = 4,
    show_row_names = FALSE,
    row_names_side = "left"
    # top_annotation =
    #   columnAnnotation("Sample number" = anno_barplot(
    #   x = sample_number,
    #   gp = gpar(col = "black",
    #             fill = ggsci::pal_lancet()(n=10)[4])
    # ))
  ) +
  rowAnnotation("-log(10, BH)" = anno_barplot(x = log_p,
                                              gp = gpar(col = "black",
                                                        fill = "black")))

plot

row_orders = ComplexHeatmap::row_order(plot)

library(ggplotify)

# pdf(file = "heatmap.pdf", width = 10, height = 15)
# plot
# dev.off()


#####then we need to fuzzy-c means clustering

time <- c(1:ncol(temp_data))

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
library(Mfuzz)
data <- table2eset(filename = "temp_data.txt")

data.s <- standardise(data)
m1 <- mestimate(data.s)
m1

plot <-
  Dmin(
    data.s,
    m = m1,
    crange = seq(2, 40, 1),
    repeats = 3,
    visu = TRUE
  )

plot <-
  plot %>%
  data.frame(distance = plot,
             k = seq(2, 40, 1)) %>%
  ggplot(aes(k, distance)) +
  geom_point(shape = 21, size = 4, fill = "black") +
  geom_smooth() +
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

# ggsave(plot, filename = "distance_k_number.pdf", width = 7, height = 7)

clust = 5

# c <- mfuzz(data.s, c = clust, m = m1)
#
# mfuzz.plot(eset = data.s,
#            min.mem = 0.8,
#            cl = c,
#            mfrow=c(2,3),
#            time.labels = time,
#            new.window = FALSE)
#
# names(c$cluster) <- rownames(temp_data)[-1]
# rownames(c$membership) <- rownames(temp_data)[-1]
# save(c, file = "c")
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
  col = colorRampPalette(colors = rev(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  ))(n = 100),
  number.cex = .7,
  addCoef.col = "black"
)

acore <- acore(data.s, c, min.acore = 0)
acore

centers <- c$centers
names(c$cluster) == rownames(c$membership)

cluster_info <-
  data.frame(
    variable_id = names(c$cluster),
    c$membership,
    cluster = c$cluster,
    stringsAsFactors = FALSE
  ) %>%
  arrange(cluster)

# openxlsx::write.xlsx(x = cluster_info,
#                      file = "cluster_info.xlsx",
#                      asTable = TRUE, overwrite = TRUE)

cluster_info <- readxl::read_xlsx("cluster_info.xlsx")
#####output the expression data of different clusters

######plot for each cluster
# for (cluster_idx in 1:clust) {
#   cat(cluster_idx, " ")
#   dir.create(paste("cluster", cluster_idx, sep = "_"))
#   cluster_data <-
#     cluster_info %>%
#     dplyr::filter(cluster_idx == cluster_idx) %>%
#     dplyr::select(1, 1 + cluster_idx)
# 
#   colnames(cluster_data) <- c("variable_id", "membership")
# 
#   cluster_data <-
#     cluster_data %>%
#     dplyr::filter(membership > 0.9)
# 
#   openxlsx::write.xlsx(x = cluster_data,
#                        file = file.path(
#                          paste("cluster", cluster_idx, sep = "_"),
#                          paste("cluster", cluster_idx, ".xlsx", sep = "")
#                        ),
#                        asTable = TRUE, overwrite = TRUE)
# 
# ###cluster plot
# 
#   temp =
#     temp_data[cluster_data$variable_id, ] %>%
#     data.frame(
#       membership = cluster_data$membership,
#       .,
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     ) %>%
#     tibble::rownames_to_column(var = "variable_id") %>%
#     tidyr::pivot_longer(
#       cols = -c(variable_id, membership),
#       names_to = "sample_id",
#       values_to = "value"
#     ) %>%
#     dplyr::left_join(sample_info[, c("time", "time")] %>%
#                        dplyr::mutate(time = as.character(time)) %>%
#                        dplyr::rename(accurate_time = time.1),
#                      by = c("sample_id" = "time"))
# 
#   temp =
#   temp %>%
#     dplyr::left_join(variable_info[,c("mol_name", "data_type")],
#                      by = c("variable_id" =
#                               "mol_name"))
# 
#   title =
#   temp %>% dplyr::distinct(variable_id, .keep_all = TRUE) %>%
#     pull(data_type) %>%
#     table()
#   title =
#   paste(paste(names(title), as.numeric(title)), collapse = ";")
# 
#   plot <-
#     ggplot() +
#     geom_rect(
#       mapping = aes(
#         xmin = start,
#         xmax = end,
#         ymin = -Inf,
#         ymax = Inf
#       ),
#       fill = "lightyellow",
#       data = day_night_df,
#       show.legend = FALSE
#     ) +
#     geom_hline(yintercept = 0) +
#     geom_line(aes(accurate_time, value,
#                   group = variable_id,
#                   color = data_type),
#               data = temp) +
#     scale_color_manual(values = class_color) +
#     # scale_color_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "OrRd")[c(1:7)])) +
#     theme_bw() +
#     theme(
#       panel.grid = element_blank(),
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
#       title = title
#     )
# 
#   plot
# 
#   ggsave(
#     plot,
#     filename = file.path(paste("cluster", cluster_idx, sep = "_"),
#                          paste("cluster", cluster_idx, ".pdf", sep = "")),
#     width = 10,
#     height = 7
#   )
# }

## (4) feature number
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx")
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx")
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx")
cluster4 <- readxl::read_xlsx("cluster_4/cluster4.xlsx")
cluster5 <- readxl::read_xlsx("cluster_5/cluster5.xlsx")

###only remain the 2,3,4
cluster = list(cluster2,
               cluster3,
               cluster4)

cluster_number =
  c(nrow(cluster2),
    nrow(cluster3),
    nrow(cluster4))

total_number = sum(cluster_number)

track_heights =
  cluster_number * 0.6 / total_number

library(circlize)

time_data = matrix(nrow = 1, ncol = ncol(temp_data))
colnames(time_data) = colnames(temp_data)
# colnames(time_data) = sort(c(colnames(temp_data),
#                         c("22:30:00", "23:00:00", "23:30:00", "24:00:00",
#                           "1:00:00", "0:30:00")))

time_data[which(colnames(time_data) == "06:00:00"):which(colnames(time_data) == "17:30:00")] = 1
time_data[is.na(time_data)] = 0

col_fun1 = colorRamp2(c(0, 1), c("grey", "lightyellow"))

colnames(time_data) =
  stringr::str_replace(string = colnames(time_data), pattern = "\\:00$", "")

# circos.clear()
# circos.par(
#   start.degree = 90-11.25,
#   clock.wise = TRUE,
#   gap.after = 45
# )
# circos.par(gap.after = 0)
# circos.heatmap(
#   t(as.matrix(time_data)),
#   rownames.side = "outside",
#   col = col_fun1,
#   dend.side = "inside",
#   track.height = 0.02,
#   cluster = FALSE,
#   bg.border = "black",
#   cell.border = "black"
# )
#
# for(i in 1:length(cluster_number)){
#   cat(i, " ")
#   circos.heatmap(
#     t(as.matrix(temp_data)[cluster[[i]]$variable_id, ]),
#     # rownames.side = "outside",
#     col = col_fun,
#     dend.side = "inside",
#     track.height = track_heights[i],
#     cluster = FALSE,
#     bg.border = "black"
#   )
#
#   circos.track(
#     track.index = get.current.track.index(),
#     panel.fun = function(x, y) {
#       if (CELL_META$sector.numeric.index == 1) {
#         # the last sector
#         circos.rect(xleft = CELL_META$cell.xlim[2] + convert_x(1, "mm"),
#                     ybottom = 0,
#                     xright =  CELL_META$cell.xlim[2] + convert_x(5, "mm"),
#                     ytop = cluster_number[i],
#                     col = "grey",border = "black"
#         )
#         circos.text(x = CELL_META$cell.xlim[2] + convert_x(3, "mm"),
#                     y = cluster_number[i]/2,
#                     col = "white",
#                     labels = c("module 2", "module 3", "module 4")[i],
#           cex = 0.5,
#           facing = "clockwise"
#         )
#       }
#     },
#     bg.border = NA
#   )
# }
#
# lgd = Legend(title = "mat1", col_fun = col_fun)
#
# grid.draw(lgd)

#######function for different cluster
cluster2 %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "mol_name")) %>%
  pull(data_type) %>%
  table()

cluster3 %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "mol_name")) %>%
  pull(data_type) %>%
  table()

cluster4 %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "mol_name")) %>%
  pull(data_type) %>%
  table()


# ######organize plot in plot folder
# 1:length(list(cluster1, cluster2, cluster3, cluster4, cluster5)) %>%
#   purrr::map(function(idx){
#     cat(idx, " ")
#     temp_cluster = list(cluster1, cluster2, cluster3, cluster4, cluster5)[[idx]]
#     name =
#     temp_cluster %>%
#       dplyr::left_join(variable_info[, c("variable_id", "mol_name")],
#                        by = c("variable_id" = "mol_name")) %>%
#       pull(variable_id.y) %>%
#       paste("pdf", sep = ".") %>%
#       stringr::str_replace("\\/", "_")
#
#     dir.create(file.path("plot", paste("cluster", idx, sep = "_")))
#     file.copy(file.path("plot", name),
#               to = file.path("plot", paste("cluster", idx, sep = "_")),
#               overwrite = TRUE, recursive = TRUE)
#   })

dim(cluster2)

load(
  here::here(
    "data/24_7_study/combine_omics/data_preparation/lipidomics_variable_info"
  )
)

cluster2_new =
  cluster2 %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "mol_name")) %>%
  dplyr::filter(data_type == "lipidomics") %>%
  dplyr::left_join(lipidomics_variable_info,
                   by = c("variable_id.y" = "variable_id"))

cluster3_new =
  cluster3 %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "mol_name")) %>%
  dplyr::filter(data_type == "lipidomics") %>%
  dplyr::left_join(lipidomics_variable_info,
                   by = c("variable_id.y" = "variable_id"))

cluster4_new =
  cluster4 %>%
  dplyr::left_join(variable_info, by = c("variable_id" = "mol_name")) %>%
  dplyr::filter(data_type == "lipidomics") %>%
  dplyr::left_join(lipidomics_variable_info,
                   by = c("variable_id.y" = "variable_id"))

cluster2_new$subclass %>% table()
cluster3_new$subclass %>% table()
cluster4_new$subclass %>% table()

####
####mosaic plot to show the lipid distributation
library(ggmosaic)

plot =
  rbind(
    data.frame(cluster2_new[, c("subclass")], cluster = '2'),
    data.frame(cluster3_new[, c("subclass")], cluster = '3'),
    data.frame(cluster4_new[, c("subclass")], cluster = '4')
  ) %>%
  dplyr::mutate(subclass = as.character(subclass)) %>%
  ggplot() +
  geom_mosaic(aes(x = product(subclass, cluster), fill = subclass),
              offset = 0.01) +
  theme_mosaic() +
  scale_fill_manual(values = lipid_class_color) +
  labs(x = "", y = "") +
  theme(
    panel.border = element_rect(color = "black",
                                fill = "transparent"),
    legend.position = "bottom"
  )

plot

# ggsave(plot, filename = "module2_3_4_lipid.pdf", width = 5, height = 12)

table(cluster2_new$subclass) * 100 / nrow(cluster2_new)
table(cluster3_new$subclass) * 100 / nrow(cluster3_new)
table(cluster4_new$subclass) * 100 / nrow(cluster4_new)

####
cluster2_new %>%
  dplyr::filter(stringr::str_detect(mol_name, "LPC")) %>%
  dplyr::filter(new_consistence_score > 0.5) %>%
  pull(variable_id.y)

######component for each cluster
###change cluster name
###1 - 4
###2 - 1
###3 - 2
###4 - 3
###5 - 5

plot =
  rbind(
    data.frame(cluster1 %>%
                 dplyr::left_join(
                   variable_info[, c("mol_name", "data_type")],
                   by = c("variable_id" = "mol_name")
                 ),
               cluster = '4'),
    data.frame(cluster2 %>%
                 dplyr::left_join(
                   variable_info[, c("mol_name", "data_type")],
                   by = c("variable_id" = "mol_name")
                 ),
               cluster = '1'),
    data.frame(cluster3 %>%
                 dplyr::left_join(
                   variable_info[, c("mol_name", "data_type")],
                   by = c("variable_id" = "mol_name")
                 ),
               cluster = '2'),
    data.frame(cluster4 %>%
                 dplyr::left_join(
                   variable_info[, c("mol_name", "data_type")],
                   by = c("variable_id" = "mol_name")
                 ),
               cluster = '3'),
    data.frame(cluster5 %>%
                 dplyr::left_join(
                   variable_info[, c("mol_name", "data_type")],
                   by = c("variable_id" = "mol_name")
                 ),
               cluster = '5')
  ) %>%
  ggplot() +
  geom_mosaic(aes(x = product(data_type, cluster), fill = data_type),
              offset = 0.01) +
  theme_mosaic() +
  scale_fill_manual(values = class_color) +
  labs(x = "", y = "") +
  theme(
    panel.border = element_rect(color = "black",
                                fill = "transparent"),
    legend.position = "right"
  )

plot

# ggsave(plot, filename = "module1-5 component.pdf", width = 5, height = 8)

#####lipid minion analysis
dir.create("lipidminion")

table(important_lipid$class1)

dir.create("lipidminion/module2")
dir.create("lipidminion/module3")
dir.create("lipidminion/module4")

lipidomics_variable_info =
  lipidomics_variable_info %>%
  dplyr::mutate(Lipid_Name = case_when(is.na(Lipid_Name) ~ mol_name,
                                       !is.na(Lipid_Name) ~ Lipid_Name))

lipidomics_variable_info$Lipid_Name =
  lipidomics_variable_info$Lipid_Name %>%
  stringr::str_replace_all("\\_", "\\/")

lipidomics_variable_info$Lipid_Name[grep("TAG", lipidomics_variable_info$Lipid_Name)] =
  lipidomics_variable_info$Lipid_Name[grep("TAG", lipidomics_variable_info$Lipid_Name)] %>%
  purrr::map(function(x) {
    # main =
    stringr::str_extract(x, "TAG[0-9]{1,3}\\.[0-9]{1,2}") %>%
      stringr::str_replace("TAG", "") %>%
      stringr::str_split("\\.") %>%
      `[[`(1) %>%
      paste(collapse = ":") %>%
      paste("TAG(", ., ")", sep = "")
  }) %>%
  unlist()

###TAG to TG and DAG to TG
lipidomics_variable_info$Lipid_Name =
  lipidomics_variable_info$Lipid_Name %>%
  stringr::str_replace_all("TAG", "TG") %>%
  stringr::str_replace_all("DAG", "DG")


cluster2_lipid =
  lipidomics_variable_info$Lipid_Name[match(cluster2_new$variable_id.y,
                                            lipidomics_variable_info$variable_id)] %>%
  data.frame(lipid = .)

cluster3_lipid =
  lipidomics_variable_info$Lipid_Name[match(cluster3_new$variable_id.y,
                                            lipidomics_variable_info$variable_id)] %>%
  data.frame(lipid = .)

cluster4_lipid =
  lipidomics_variable_info$Lipid_Name[match(cluster4_new$variable_id.y,
                                            lipidomics_variable_info$variable_id)] %>%
  data.frame(lipid = .)

# write.table(
#   cluster2_lipid,
#   file = "lipidminion/module2/positive_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )
#
# write.table(
#   cluster3_lipid,
#   file = "lipidminion/module3/positive_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )
#
# write.table(
#   cluster4_lipid,
#   file = "lipidminion/module4/positive_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )

universe_lipid =
  lipidomics_variable_info %>%
  dplyr::select(Lipid_Name) %>%
  dplyr::filter(!is.na(Lipid_Name)) %>%
  dplyr::rename(lipid = Lipid_Name) %>%
  dplyr::distinct(lipid)

# write.table(
#   universe_lipid,
#   file = "lipidminion/universe_lipid.txt",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote = FALSE
# )


#####read the lipid minion results
#####module 2
module2_result = readr::read_delim("lipidminion/module2/Fisher output table (10 pvals _ 0.05).txt",
                                   delim = " ")

module2_network_node =
  readr::read_delim("lipidminion/module2/Network_nodes.txt", delim = "\t")
module2_network_edge =
  readr::read_delim("lipidminion/module2/Network_edges.txt", delim = "\t")
module2_network_edge_attr =
  readr::read_delim("lipidminion/module2/Network_edge_attributes.txt", delim = "\t")

unique(module2_network_node$title)

###network
###module 2
library(igraph)
library(ggraph)
library(tidygraph)

module2_network_node

node =
  module2_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class"))

edge =
  module2_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

module2_netwrok =
  tidygraph::tbl_graph(nodes = node,
                       edges = edge)

plot =
  ggraph(module2_netwrok,
         layout = 'kk',
         circular = FALSE) +
  geom_edge_link(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "black"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 8,
                               "class" = 15)) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(plot, filename = "lipidminion/module2/module2_network.pdf", width = 7, height = 7)


#####read the lipid minion results
##module3
module3_result = readr::read_delim("lipidminion/module3/Fisher output table (4 pvals _ 0.05).txt",
                                   delim = " ")

module3_network_node =
  readr::read_delim("lipidminion/module3/Network_nodes.txt", delim = "\t")
module3_network_edge =
  readr::read_delim("lipidminion/module3/Network_edges.txt", delim = "\t")
module3_network_edge_attr =
  readr::read_delim("lipidminion/module3/Network_edge_attributes.txt", delim = "\t")

unique(module3_network_node$title)

###network
###module3
library(igraph)
library(ggraph)
library(tidygraph)

module3_network_node

node =
  module3_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class"))

edge =
  module3_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

module3_netwrok =
  tidygraph::tbl_graph(nodes = node,
                       edges = edge)

plot =
  ggraph(module3_netwrok,
         layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "black"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 8,
                               "class" = 15)) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(plot, filename = "lipidminion/module3/module3_network.pdf", width = 7, height = 7)


#####read the lipid minion results
##module4
module4_result = readr::read_delim("lipidminion/module4/Fisher output table (4 pvals _ 0.05).txt",
                                   delim = " ")

module4_network_node =
  readr::read_delim("lipidminion/module4/Network_nodes.txt", delim = "\t")
module4_network_edge =
  readr::read_delim("lipidminion/module4/Network_edges.txt", delim = "\t")
module4_network_edge_attr =
  readr::read_delim("lipidminion/module4/Network_edge_attributes.txt", delim = "\t")

unique(module4_network_node$title)

###network
###module4
library(igraph)
library(ggraph)
library(tidygraph)

module4_network_node

node =
  module4_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class"))

edge =
  module4_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

module4_netwrok =
  tidygraph::tbl_graph(nodes = node,
                       edges = edge)

plot =
  ggraph(module4_netwrok,
         layout = 'fr',
         circular = FALSE) +
  geom_edge_link(
    strength = 1,
    alpha = 1,
    show.legend = FALSE,
    color = "black"
  ) +
  geom_node_point(
    aes(color = class, size = class),
    shape = 16,
    alpha = 1,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 8,
                               "class" = 15)) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

# ggsave(plot, filename = "lipidminion/module4/module4_network.pdf", width = 7, height = 7)




#####combine three network together
module2_node =
  module2_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class")) %>%
  dplyr::select(node, title, class)

module2_edge =
  module2_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

module2_edge$from = module2_node$title[match(module2_edge$from, module2_node$node)]
module2_edge$to = module2_node$title[match(module2_edge$to, module2_node$node)]

module2_node$node = module2_node$title

module3_node =
  module3_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class")) %>%
  dplyr::select(node, title, class)

module3_edge =
  module3_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

module3_edge$from = module3_node$title[match(module3_edge$from, module3_node$node)]
module3_edge$to = module3_node$title[match(module3_edge$to, module3_node$node)]

module3_node$node = module3_node$title

module4_node =
  module4_network_node %>% dplyr::select(-id) %>% dplyr::rename(node = label) %>%
  dplyr::mutate(class = case_when(shape == "#cccccc" ~ "lipid",
                                  TRUE ~ "class")) %>%
  dplyr::select(node, title, class)

module4_edge =
  module4_network_edge %>% dplyr::select(to, color) %>% dplyr::rename(from = to, to = color)

module4_edge$from = module4_node$title[match(module4_edge$from, module4_node$node)]
module4_edge$to = module4_node$title[match(module4_edge$to, module4_node$node)]

module4_node$node = module4_node$title



module2_node$module = 1
module3_node$module = 1
module4_node$module = 1

node_data =
  module2_node %>%
  dplyr::full_join(module3_node, by = c("node", "title", "class")) %>%
  dplyr::full_join(module4_node, by = c("node", "title", "class")) %>%
  dplyr::rename(module2 = module.x,
                module3 = module.y,
                module4 = module) %>%
  dplyr::mutate(module2 = case_when(is.na(module2) ~ 0,
                                    TRUE ~ module2)) %>%
  dplyr::mutate(module3 = case_when(is.na(module3) ~ 0,
                                    TRUE ~ module3))  %>%
  dplyr::mutate(module4 = case_when(is.na(module4) ~ 0,
                                    TRUE ~ module4)) %>%
  dplyr::mutate(radius = case_when(class == "class" ~ 0.5,
                                   TRUE ~ 0.2))

edge_data =
  rbind(module2_edge,
        module3_edge,
        module4_edge) %>%
  dplyr::distinct(.keep_all = TRUE)

network =
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data)


library(scatterpie)
library(igraph)
library(graphlayouts)
xy =
  ggraph::create_layout(network, layout = "igraph", algorithm = 'fr')
V(network)$x <- xy$x
V(network)$y <- xy$y
library(scatterpie)

plot =
  ggraph(
    network,
    layout = 'manual',
    x = V(network)$x,
    y = V(network)$y,
    circular = FALSE
  ) +
  geom_edge_link(alpha = 1,
                 show.legend = FALSE,
                 color = "black") +
  geom_scatterpie(
    aes(x = x, y = y, r = radius),
    data = as_data_frame(network, "vertices"),
    cols = c("module2", "module3", "module4")
  ) +
  scale_fill_manual(
    values = c(
      "module2" = ggsci::pal_lancet()(n = 7)[2],
      "module3" = ggsci::pal_lancet()(n = 7)[4],
      "module4" = ggsci::pal_lancet()(n = 7)[5]
    )
  ) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_fixed()

plot

######seperate it
subnetwork =
  igraph::cluster_edge_betweenness(graph = network)


####subnetwork1
node_data1 =
  node_data[subnetwork$membership == 1, ]

edge_data1 =
  edge_data %>%
  dplyr::filter(from %in% node_data1$node &
                  to %in% node_data1$node)

network1 =
  tidygraph::tbl_graph(nodes = node_data1,
                       edges = edge_data1)
xy =
  ggraph::create_layout(network1, layout = "igraph", algorithm = 'kk')
V(network1)$x <- xy$x
V(network1)$y <- xy$y
library(scatterpie)

plot =
  ggraph(
    network1,
    layout = 'manual',
    x = V(network1)$x,
    y = V(network1)$y,
    circular = FALSE
  ) +
  geom_edge_link(alpha = 1,
                 show.legend = FALSE,
                 color = "black") +
  geom_scatterpie(
    aes(x = x, y = y, r = radius),
    data = as_data_frame(network1, "vertices"),
    cols = c("module2", "module3", "module4")
  ) +
  scale_fill_manual(
    values = c(
      "module2" = ggsci::pal_lancet()(n = 7)[2],
      "module3" = ggsci::pal_lancet()(n = 7)[4],
      "module4" = ggsci::pal_lancet()(n = 7)[5]
    )
  ) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_fixed()

plot

# ggsave(plot, filename = "lipidminion/subnetwork1_2.pdf", width = 7, height = 7)


####subnetwork2
node_data2 =
  node_data[subnetwork$membership == 2, ]

edge_data2 =
  edge_data %>%
  dplyr::filter(from %in% node_data2$node &
                  to %in% node_data2$node)

network2 =
  tidygraph::tbl_graph(nodes = node_data2,
                       edges = edge_data2)
xy =
  ggraph::create_layout(network2, layout = "igraph", algorithm = 'fr')
V(network2)$x <- xy$x
V(network2)$y <- xy$y
library(scatterpie)

plot =
  ggraph(
    network2,
    layout = 'manual',
    x = V(network2)$x,
    y = V(network2)$y,
    circular = FALSE
  ) +
  geom_edge_link(alpha = 1,
                 show.legend = FALSE,
                 color = "black") +
  geom_scatterpie(
    aes(x = x, y = y, r = radius),
    data = as_data_frame(network2, "vertices"),
    cols = c("module2", "module3", "module4")
  ) +
  scale_fill_manual(
    values = c(
      "module2" = ggsci::pal_lancet()(n = 7)[2],
      "module3" = ggsci::pal_lancet()(n = 7)[4],
      "module4" = ggsci::pal_lancet()(n = 7)[5]
    )
  ) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_fixed()

plot
# ggsave(plot, filename = "lipidminion/subnetwork2_2.pdf", width = 7, height = 7)










####subnetwork3
node_data3 =
  node_data[subnetwork$membership == 3, ]

edge_data3 =
  edge_data %>%
  dplyr::filter(from %in% node_data3$node &
                  to %in% node_data3$node)

network3 =
  tidygraph::tbl_graph(nodes = node_data3,
                       edges = edge_data3)
xy =
  ggraph::create_layout(network3, layout = "igraph", algorithm = 'fr')
V(network3)$x <- xy$x
V(network3)$y <- xy$y
library(scatterpie)

plot =
  ggraph(
    network3,
    layout = 'manual',
    x = V(network3)$x,
    y = V(network3)$y,
    circular = FALSE
  ) +
  geom_edge_link(alpha = 1,
                 show.legend = FALSE,
                 color = "black") +
  geom_scatterpie(
    aes(x = x, y = y, r = radius),
    data = as_data_frame(network3, "vertices"),
    cols = c("module2", "module3", "module4")
  ) +
  scale_fill_manual(
    values = c(
      "module2" = ggsci::pal_lancet()(n = 7)[2],
      "module3" = ggsci::pal_lancet()(n = 7)[4],
      "module4" = ggsci::pal_lancet()(n = 7)[5]
    )
  ) +
  ggnewscale::new_scale(new_aes = "size") +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = title,
      size = class
    ),
    color = "black",
    bg.color = "white",
    # size = 4,
    check_overlap = FALSE,
    show.legend = FALSE
  ) +
  scale_size_manual(values = c("lipid" = 3,
                               "class" = 5)) +
  scale_color_manual(values = c("lipid" = unname(class_color["lipidomics"]),
                                "class" = "red")) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_fixed()

plot

# ggsave(plot, filename = "lipidminion/subnetwork3-2.pdf", width = 7, height = 7)


temp_data <-
  new_result %>%
  dplyr::left_join(variable_info[, c("variable_id", "mol_name")],
                   by = c("CycID" = "variable_id")) %>%
  dplyr::left_join(cluster_info[, c("variable_id", "cluster")],
                   by = c("mol_name" = "variable_id")) %>%
  dplyr::filter(!is.na(cluster)) %>%
  dplyr::filter(cluster %in% c(2, 3, 4)) %>%
  dplyr::filter(mol_name %in% c(
    cluster2$variable_id,
    cluster3$variable_id,
    cluster4$variable_id
  ))

table(temp_data$cluster)

temp_data %>%
  dplyr::filter(cluster == 2) %>%
  pull(Period) %>%
  mean()

temp_data %>%
  dplyr::filter(cluster == 2) %>%
  pull(PhaseShift) %>%
  mean()

temp_data %>%
  dplyr::filter(cluster == 3) %>%
  pull(Period) %>%
  mean()

temp_data %>%
  dplyr::filter(cluster == 3) %>%
  pull(PhaseShift) %>%
  mean()

temp_data %>%
  dplyr::filter(cluster == 4) %>%
  pull(Period) %>%
  mean()

temp_data %>%
  dplyr::filter(cluster == 4) %>%
  pull(PhaseShift) %>%
  mean()
