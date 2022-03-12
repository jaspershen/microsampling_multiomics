no_function()
library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

####we need the lipidomics module information
lipidomics_variable_info =
  readxl::read_xlsx("data/7_24_mike/lipidomics/k_means_clustering/final_info_manual.xlsx")

load("data/7_24_mike/lipidomics/k_means_clustering/new_variable_info")
load("data/7_24_mike/lipidomics/k_means_clustering/new_expression_data")

lipidomics_variable_info =
  new_variable_info

####load all omics data loess data
load(
  here::here(
    "data/7_24_mike/combine_omics/data_preparation/new_expression_data"
  )
)
load(here::here(
  "data/7_24_mike/combine_omics/data_preparation/new_sample_info"
))
load(here::here(
  "data/7_24_mike/combine_omics/data_preparation/new_variable_info"
))

variable_info1 =
  new_variable_info %>%
  dplyr::filter(data_type == "lipidomics")

expression_data1 =
  new_expression_data[variable_info1$variable_id, ]

variable_info2 =
  new_variable_info %>%
  dplyr::filter(data_type != "lipidomics")

expression_data2 =
  new_expression_data[variable_info2$variable_id, ]

expression_data1 =
  lipidomics_variable_info$old_variable_id %>%
  purrr::map(function(x) {
    x = stringr::str_split(x, ";")[[1]]
    expression_data1[x, , drop = FALSE] %>%
      colMeans()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

rownames(expression_data1) = lipidomics_variable_info$variable_id

variable_info1 =
  data.frame(lipidomics_variable_info[, c("variable_id", "mol_name")],
             data_type = "lipidomics")

variable_info =
  rbind(variable_info1,
        variable_info2)

expression_data =
  rbind(expression_data1,
        expression_data2)

sample_info =
  new_sample_info

rownames(expression_data) == variable_info$variable_id
colnames(expression_data) == sample_info$sample_id

####this is for the day night time
load("data/7_24_mike/summary_info/day_night_df")

day_night_df =
  day_night_df %>%
  dplyr::mutate(
    start_time = as.POSIXct(hms::as_hms(start)),
    end_time = as.POSIXct(hms::as_hms(end)),
    week = format(day, "%a")
  ) %>%
  dplyr::mutate(week = paste(week,
                             lubridate::month(day),
                             lubridate::day(day),
                             sep = "-")) %>%
  dplyr::mutate(week = factor(week, unique(week)))

######
setwd("data/7_24_mike/inter_omics_correlation/inter_all_omics_loess_data")

lagged_cor =
  readxl::read_xlsx("lagged_correlation/lagged_cor.xlsx")

table(paste(lagged_cor$from_data_type, lagged_cor$to_data_type, sep = "_"))

####
lagged_cor1 =
  lagged_cor %>%
  dplyr::filter(shift_time != "(-15,15]")

lagged_cor2 =
  lagged_cor %>%
  dplyr::filter(shift_time == "(-15,15]")

####construct the network for lagged correlation
##############total network
temp_data =
  lagged_cor1 %>%
  dplyr::filter(p_adjust < 0.05)

dim(temp_data)

###there are too many edges,
###so we need to make filtering

temp_data =
  temp_data %>%
  dplyr::filter(cor > quantile(temp_data$cor[temp_data$cor > 0], 0.75) |
                  cor < quantile(temp_data$cor[temp_data$cor < 0], 0.25))

dim(temp_data)

shift_time2 =
  purrr::map(stringr::str_split(temp_data$shift_time, ","), function(x) {
    mean(as.numeric(stringr::str_replace(x, "\\(|\\]", "")))
  }) %>%
  unlist()

temp_data$shift_time2 =
  shift_time2

table(temp_data$shift_time2)

plot(temp_data$cor)

####example network
###Network shows top 300
###significant correlations (FDR P < 0.05) between each pair of measurement
###types
edge_data =
  temp_data %>%
  dplyr::select(from,
                to,
                shift_time,
                shift_time2,
                cor,
                p_adjust,
                dplyr::everything()) %>%
  dplyr::rename(time = shift_time2) %>%
  dplyr::mutate(shift = case_when(time > 0 ~ "1",
                                  time < 0 ~ "-1")) %>%
  dplyr::arrange(desc(abs(cor)))

edge_data$edge_class =
  apply(edge_data, 1, function(x) {
    paste(sort(as.character(x[c(11, 13)])), collapse = "_")
  })

edge_data %>%
  dplyr::count(edge_class)

library(plyr)

edge_data =
  edge_data %>%
  plyr::dlply(.variables = .(edge_class)) %>%
  purrr::map(function(x) {
    head(x, 50)
  }) %>%
  dplyr::bind_rows()

dim(edge_data)

edge_data$p_adjust[edge_data$p_adjust == 0] =
  min(edge_data$p_adjust[edge_data$p_adjust != 0])

node_data =
  data.frame(node = unique(c(edge_data$from, edge_data$to))) %>%
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))

node_data =
  node_data %>%
  dplyr::mutate(true_name = case_when(stringr::str_detect(node, "module") ~ node,
                                      TRUE ~ mol_name))

dim(node_data)
dim(edge_data)

library(igraph)
library(tidygraph)

total_lagged_graph <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  tidygraph::mutate(Degree = centrality_degree(mode = 'all'))

library(ggraph)

####just add the node and edge labels based on if they are high-degree nodes or
###high correlations

node_data =
  igraph::vertex_attr(total_lagged_graph) %>%
  dplyr::bind_cols() %>%
  as.data.frame()

label_node1 =
  node_data %>%
  plyr::dlply(.variables = .(data_type)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::arrange(desc(Degree)) %>%
      head(1)
  }) %>%
  dplyr::bind_rows() %>%
  pull(node)

label_node2 =
  edge_data %>%
  plyr::dlply(.variables = .(edge_class)) %>%
  purrr::map(function(x) {
    x =
      x %>%
      dplyr::arrange(desc(abs(cor))) %>%
      head(1) %>%
      dplyr::select(from, to)
    unique(c(x$from, x$to))
  }) %>%
  unlist() %>%
  unname() %>%
  unique()

label_node = unique(c(label_node1, label_node2))

plot <-
  ggraph(total_lagged_graph,
         layout = 'fr') +
  geom_edge_link(
    aes(
      color = cor,
      strength = 0.1,
      width = -log(p_adjust, 10)
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_node_point(
    aes(fill = data_type,
        size = Degree),
    shape = 21,
    alpha = 0.9,
    show.legend = TRUE
  ) +
  # geom_node_text(
  #   aes(
  #     x = x * 1.05,
  #     y = y * 1.05,
  #     label = ifelse(node %in% label_node, true_name, NA)
  #   ),
  #   color = "black",
  #   size = 2,
  #   show.legend = FALSE,
  #   check_overlap = TRUE
  # ) +
shadowtext::geom_shadowtext(
  aes(x = x * 1.05,
      y = y * 1.05,
    label =
        ifelse(node %in% label_node, true_name, NA)), 
  check_overlap = TRUE, size = 3,
  bg.colour='white',
  colour = "black") +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color)) +
  scale_color_manual(values = c(class_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
                           override.aes = list(linetype = "none")),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = NA,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    fill = guide_legend(
      title = "Class",
      override.aes = list(
        shape = 21,
        size = 3,
        alpha = 1
      )
    )
  ) +
  scale_edge_width_continuous(range = c(0.05, 0.7)) +
  ggraph::scale_edge_color_gradient2(
    low = alpha("#3B4992FF", 0.7),
    mid = "white",
    high = alpha("#EE0000FF", 0.7)
  ) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

node_data =
  igraph::vertex_attr(graph = total_lagged_graph) %>%
  dplyr::bind_cols() %>%
  as.data.frame()

ggsave(
  plot,
  filename = (
    "example_lagged_network.pdf"
  ),
  width = 8.3,
  height = 7
)

# ###output the network data
# library(openxlsx)
# wb <- createWorkbook()
# modifyBaseFont(wb, fontSize = 12, fontName = "Time New Roma")
# addWorksheet(wb, sheetName = "Node data",
#              gridLines = TRUE)
# addWorksheet(wb, sheetName = "Edge data",
#              gridLines = TRUE)
# freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE)
# freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE)
# writeDataTable(wb, sheet = 1, x = node_data,
#                colNames = TRUE, rowNames = FALSE)
# writeDataTable(wb, sheet = 2, x = edge_data,
#                colNames = TRUE, rowNames = FALSE)
# saveWorkbook(wb, "example_network_data.xlsx", overwrite = TRUE)



####property of network
plot =
  node_data %>%
  dplyr::arrange(Degree) %>%
  dplyr::filter(Degree > 10 | data_type != "lipidomics") %>%
  dplyr::mutate(mol_name = factor(mol_name, levels = mol_name)) %>%
  ggplot(aes(Degree, mol_name)) +
  geom_segment(aes(
    x = 0,
    xend = Degree,
    y = mol_name,
    yend = mol_name,
    color = data_type
  ),
  show.legend = FALSE) +
  geom_point(aes(size = Degree,
                 color = data_type),
             shape = 16,
             show.legend = FALSE) +
  scale_color_manual(values = class_color) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  base_theme +
  geom_text(aes(x = Degree, mol_name, label = mol_name), hjust = -0.1) +
  labs(y = "", x = "Degree") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ggforce::facet_col(vars(data_type),
                     scales = "free_y", space = "free")
plot

ggsave(plot,
       filename = "lagged_cor_network/degree.pdf",
       width = 7,
       height = 7)

###cor distributation
plot(edge_data$time, edge_data$cor)
plot1 =
  edge_data %>%
  dplyr::filter(cor > 0) %>%
  ggplot(aes(time, cor)) +
  # geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(
    aes(size = -log(p_adjust, 10),
        fill = cor),
    alpha = 0.5,
    shape = 21,
    show.legend = FALSE
  ) +
  base_theme +
  scale_fill_gradient2(low = "#366A9FFF", mid = "white", high = "red") +
  labs(x = "",
       y = "Lagged correlation") +
  scale_x_continuous(
    breaks = sort(unique(edge_data$time)),
    labels = sort(unique(edge_data$time)),
    limits = c(-180, 250)
  ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot2 =
  edge_data %>%
  dplyr::filter(cor < 0) %>%
  ggplot(aes(time, cor)) +
  # geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(
    aes(size = -log(p_adjust, 10),
        fill = cor),
    alpha = 0.5,
    shape = 21,
    show.legend = FALSE
  ) +
  base_theme +
  scale_fill_gradient2(low = "#366A9FFF", mid = "white", high = "red") +
  labs(x = "Shift time (min)",
       y = "Lagged correlation") +
  scale_x_continuous(
    breaks = sort(unique(edge_data$time)),
    labels = sort(unique(edge_data$time)),
    limits = c(-180, 250)
  )

library(patchwork)

plot =
  plot1 + plot2 + patchwork::plot_layout(ncol = 1)
plot
# ggsave(plot,
#        filename = "lagged_cor_network/shift_time_vs_cor.pdf",
#        width = 7,
#        height = 7)


table(edge_data$time[edge_data$cor > 0])
table(edge_data$time[edge_data$cor < 0])

######

####output protein and hormone with their lipids correlation
node_data %>%
  dplyr::filter(data_type != "lipidomics")

###GIP
edge_data %>%
  dplyr::filter(from_mol_name == "GIP")

edge_data %>%
  dplyr::filter(to_mol_name == "GIP")

temp_data =
  edge_data %>%
  dplyr::filter(to_mol_name == "GIP" |
                  to_mol_name == "P02649_APOE" |
                  to_mol_name == "P06727_APOA4")

plot =
  temp_data %>%
  dplyr::filter(score > 0.5) %>%
  ggplot(aes(time, from_mol_name)) +
  geom_segment(aes(
    x = 0,
    xend = time,
    y = from_mol_name,
    yend = from_mol_name,
    size = score
  ),
  alpha = 0.5) +
  scale_size_continuous(range = c(0.01, 1)) +
  ggnewscale::new_scale(new_aes = "size") +
  geom_vline(xintercept = 0) +
  geom_point(aes(size = -log(p_adjust, 10),
                 fill = cor), shape = 21) +
  scale_fill_gradient2(low = "#366A9FFF", mid = "white", high = "red") +
  scale_size_continuous(range = c(0.5, 9)) +
  geom_text(
    aes(
      x = 0,
      y = from_mol_name,
      color = cor,
      label = ifelse(time > 0, round(cor, 2), NA)
    ),
    hjust = 1,
    size = 4
  ) +
  geom_text(
    aes(
      0,
      from_mol_name,
      color = cor,
      label = ifelse(time < 0, round(cor, 2), NA)
    ),
    hjust = -0.5,
    size = 4
  ) +
  scale_color_gradient2(low = "#366A9FFF", mid = "white", high = "red") +
  base_theme +
  labs(y = "", x = "Shift time (GIP - molecules, min)") +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  # scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  facet_grid(
    rows = vars(to_mol_name),
    scales = "free_y",
    shrink = TRUE,
    as.table = TRUE,
    space = "free"
  )

plot

# ggsave(plot, filename = "lagged_cor_network/protein_hormone_lipid.pdf", width = 7, height = 10)


###top degree lipid
# important_lipid =
# node_data %>%
#   dplyr::filter(Degree > 20 & data_type == "lipidomics") %>%
#   dplyr::pull(node)

temp_data =
  edge_data %>%
  # dplyr::filter(from %in% important_lipid | to %in% important_lipid) %>%
  dplyr::filter(from_data_type == "lipidomics" &
                  to_data_type == "lipidomics")

temp_edge_data =
  temp_data

temp_node_data =
  data.frame(node = unique(c(temp_edge_data$from, temp_edge_data$to))) %>%
  dplyr::left_join(variable_info, by = c("node" = "variable_id"))

temp_edge_data =
  temp_edge_data %>%
  dplyr::filter(from %in% temp_node_data$node &
                  to %in% temp_node_data$node)

lipid_lagged_graph <-
  tidygraph::tbl_graph(nodes = temp_node_data,
                       edges = temp_edge_data,
                       directed = FALSE) %>%
  tidygraph::mutate(Degree = centrality_degree(mode = 'all')) %>%
  tidygraph::mutate(class = ifelse(Degree > 20, "yes", "no"))

g <- lipid_lagged_graph
library(igraph)
library(ggraph)
V(g)$type <- rep(TRUE, nrow(temp_node_data))

coords <-
  ggraph::create_layout(g, layout = "bipartite") %>%
  dplyr::select(node, class, mol_name, x, y)

coords$y[coords$class == "yes"] = 0.2
coords$y[coords$class == "no"] = 1

table(coords$class)

coords$x[coords$y == 0.2]
coords$mol_name[coords$y == 0.2]

coords$x[coords$y == 0.2] = c(1, 15, 25, 35, 45)

coords <-
  coords %>%
  dplyr::select(x, y) %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

my_graph <-
  create_layout(
    graph = g,
    layout = "manual",
    x = coords$x,
    y = coords$y
  )

library(ggraph)
plot <-
ggraph(my_graph,
       layout = 'bipartite') +
  geom_edge_diagonal(
    aes(
      # label = ifelse(cor > 0.8, round(cor, 2), ""),
      color = cor,
      strength = 0.1,
      width = -log(p_adjust, 10)
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  shadowtext::geom_shadowtext(
    aes(
      x = x,
      y = y,
      label = ifelse(class == "yes", mol_name, NA),
      color = data_type
    ),
    bg.color = "white",
    size = 3,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  geom_node_point(
    aes(fill = data_type,
        size = Degree),
    shape = 21,
    alpha = 0.5,
    show.legend = TRUE
  ) +
  geom_node_text(
    aes(
      x = x * 1.05,
      y = y * 1.05,
      label = ifelse(class == "no", mol_name, NA),
      hjust = 'outward',
      angle = -((-node_angle(x, y) + 90) %% 180) + 90
    ),
    color = "black",
    size = 2,
    show.legend = FALSE,
    check_overlap = FALSE
  ) +
  scale_size_continuous(range = c(3, 10)) +
  scale_fill_manual(values = c(class_color)) +
  scale_color_manual(values = c(class_color)) +
  guides(
    linetype = "none",
    color = guide_colorbar(title = "Lagged correlation",
                           override.aes = list(linetype = "none")),
    size = guide_legend(
      title = "Degree",
      override.aes = list(
        linetype = NA,
        fill = "transparent",
        shape = 21,
        color = "black"
      )
    ),
    fill = guide_legend(
      title = "Class",
      override.aes = list(
        shape = 21,
        size = 3,
        alpha = 1
      )
    )
  ) +
  scale_edge_width_continuous(range = c(0.1, 1)) +
  ggraph::scale_edge_color_gradient2(
    low = alpha("#3B4992FF", 0.7),
    mid = "white",
    high = alpha("#EE0000FF", 0.7)
  ) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

extrafont::loadfonts()

# ggsave(
#   plot,
#   filename = (
#     "lagged_cor_network/total_lipid_lagged_network.pdf"
#   ),
#   width = 8.3,
#   height = 7
# )

