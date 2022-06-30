no_function()

masstools::setwd_project()
library(tidyverse)

rm(list=ls())

source("code/tools.R")

load("data/ensure_shake_study/subject_info/subject_info")

load("data/ensure_shake_study/lipidomics_data_analysis/data_preparation/sample_info")
load("data/ensure_shake_study/lipidomics_data_analysis/data_preparation/variable_info")

lipidomics_sample_info <-
  sample_info

lipidomics_variable_info <-
  variable_info

load("data/ensure_shake_study/metabolomics_data_analysis/data_preparation/sample_info")
load("data/ensure_shake_study/metabolomics_data_analysis/data_preparation/variable_info")
metabolomics_sample_info <-
  sample_info
metabolomics_variable_info <-
  variable_info

load("data/ensure_shake_study/cytokine_data_analysis/data_preparation/sample_info")
load("data/ensure_shake_study/cytokine_data_analysis/data_preparation/variable_info")
cytokine_sample_info <-
  sample_info
cytokine_variable_info <-
  variable_info

all_subject_id <-
c(
  unique(cytokine_sample_info$subject_id),
  unique(metabolomics_sample_info$subject_id),
  unique(lipidomics_sample_info$subject_id)
) %>% 
  unique()

setdiff(subject_info$subject_id, all_subject_id)

setdiff(all_subject_id, subject_info$subject_id)

setwd("data/ensure_shake_study/study_information")

###omics sample overlap
dim(metabolomics_sample_info)
dim(lipidomics_sample_info)
dim(cytokine_sample_info)

library(ComplexUpset)

####the upset plot to show the overlap of different omics data
####we have six omics data
library(ComplexUpset)

dim(metabolomics_sample_info)
dim(cytokine_sample_info)
dim(lipidomics_sample_info)

temp_data1 =
  data.frame(metabolomics_sample_info[, c("subject_id", "sample_id")],
             class = "Metabolomics") %>%
  dplyr::select(sample_id, class) %>% 
  dplyr::distinct(sample_id, .keep_all = TRUE)

temp_data2 =
  data.frame(cytokine_sample_info[, c("subject_id", "sample_id")],
             class = "Cytokine") %>%
  dplyr::select(sample_id, class) %>% 
  dplyr::distinct(sample_id, .keep_all = TRUE)

temp_data3 =
  data.frame(lipidomics_sample_info[, c("subject_id", "sample_id")],
             class = "Lipidomics") %>%
  dplyr::select(sample_id, class) %>% 
  dplyr::distinct(sample_id, .keep_all = TRUE)

temp_data = 
  rbind(temp_data1,
        temp_data2,
        temp_data3) %>% 
  tidyr::pivot_wider(names_from = "class", values_from = "class") %>% 
  dplyr::mutate(Metabolomics = case_when(
    is.na(Metabolomics) ~ FALSE,
    TRUE ~ TRUE
  ),
  Cytokine = case_when(
    is.na(Cytokine) ~ FALSE,
    TRUE ~ TRUE
  ),
  Lipidomics = case_when(
    is.na(Lipidomics) ~ FALSE,
    TRUE ~ TRUE
  )
  )

library(ComplexUpset)

plot = 
  upset(
    data = temp_data,
    name = "",
    intersect = colnames(temp_data)[-1],
    # set_sizes = TRUE,
    min_degree = 1,
    # group_by='sets',
    # min_size = 3,
    sort_sets = FALSE,
    sort_intersections = "ascending",
    sort_intersections_by=c('degree', "cardinality"),
    stripes = alpha(unname(shake_omics_color[colnames(temp_data)[-1]]), 0.8),
    base_annotations = list(
      'Intersection size' = intersection_size(counts = TRUE,
                                              mode = "intersect",
                                              text = list(angle = 0),
                                              mapping = aes(fill = 'bars_color')) +
        scale_fill_manual(values = c('bars_color' = 'black'), 
                          guide = 'none') +
        theme_classic() +
        labs(x = "", y = "Sample number") +
        scale_y_continuous(expand = expansion(mult = c(0,0))) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    ),
    matrix = (
      intersection_matrix(
        geom=geom_point(
          shape=16,
          size = 13
        ),
        segment=geom_segment(
          linetype = 1
        ),
        outline_color=list(
          active='white',
          inactive='white'
        ) 
      )
      # scale_color_manual(
      #   values=c('TRUE' = 'black', 'FALSE' = 'grey')
      #   # labels=c('TRUE'='yes', 'FALSE'='no'),
      #   # breaks=c('TRUE', 'FALSE'),
      #   # name=''
      # )
    )
    # queries=list(
    #   upset_query(
    #     intersect=c('RNA'),
    #     color=omics_color["RNA"],
    #     fill=omics_color["RNA"],
    #     only_components=c('intersections_matrix', 'Intersection size')
    #   )
    # )
  ) 

plot
# ggsave(plot, filename = "upset_plot_of_all_samples.pdf", width = 9, height = 7)

lipidomics_sample_info %>% 
  ggplot(aes(TP, subject_id)) +
  geom_point()


library(scatterpie)

temp_data1 <- 
metabolomics_sample_info %>% 
  dplyr::select(subject_id, TP) %>% 
  dplyr::mutate(Metabolomics = 1)

temp_data2 <- 
  lipidomics_sample_info %>% 
  dplyr::select(subject_id, TP) %>% 
  dplyr::mutate(Lipidomics = 1)

temp_data3 <- 
  cytokine_sample_info %>% 
  dplyr::select(subject_id, TP) %>% 
  dplyr::mutate(Cytokine = 1)

temp_data <- 
temp_data1 %>% 
  dplyr::full_join(temp_data2, by = c("subject_id", "TP")) %>% 
  dplyr::full_join(temp_data3, by = c("subject_id", "TP"))

temp_data[is.na(temp_data)] <- 0

temp_data <- 
temp_data %>% 
  rowwise %>% 
  dplyr::mutate(None = sum(Metabolomics, Lipidomics, Cytokine)) %>% 
  dplyr::mutate(None = 3 - None)

new_subject_id <- 
  data.frame(subject_id = unique(temp_data$subject_id)) %>% 
  dplyr::arrange(as.numeric(stringr::str_replace(subject_id, "S", ""))) %>% 
  dplyr::mutate(new_subject_id = 1:28)

temp_data <- 
temp_data %>% 
  dplyr::left_join(new_subject_id, by = c("subject_id"))

temp_data$x <- temp_data$new_subject_id * 10

plot <- 
ggplot() +
  geom_scatterpie(
    mapping = aes(x = TP, y = x),
    data = temp_data,
    cols = factor(c("Metabolomics", "Lipidomics", "Cytokine", "None")),
    alpha = 0.8
  ) +
  coord_fixed() +
  scale_y_continuous(breaks = c(1:28)*10, labels = new_subject_id$subject_id,
                     expand = expansion(mult = c(0.01,0.01))) +
  scale_x_continuous(breaks = c(0, 30, 60, 120, 240), 
                     labels = c(0, 30, 60, 120, 240)) +
  labs(x = "Time (min)", y = "") +
  scale_fill_manual(values = c(shake_omics_color, "None" = "grey")) +
  base_theme +
  theme(legend.position = "none")

plot

plot1 <- 
temp_data %>% 
  dplyr::select(-c(TP, x)) %>% 
  tidyr::pivot_longer(cols = -c(subject_id,new_subject_id), 
                      names_to = "class", values_to = "number") %>% 
  dplyr::filter(number != 0 & class != "None") %>% 
  ggplot(aes(y = new_subject_id)) +
  geom_bar(aes(fill = class), show.legend = FALSE) +
  scale_fill_manual(values = c(shake_omics_color)) +
  scale_y_continuous(breaks = 1:28, labels = new_subject_id$subject_id,
                     expand = expansion(mult = c(0.01,0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = c(0, 5, 10, 15),
                     labels = c(0, 5, 10, 15)) +
  labs(x = "Sample number", y = "") + 
  base_theme +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot2 <- 
  temp_data %>% 
  dplyr::select(-c(subject_id, new_subject_id, x)) %>% 
  tidyr::pivot_longer(cols = -c(TP), 
                      names_to = "class", values_to = "number") %>% 
  dplyr::filter(number != 0 & class != "None") %>% 
  ggplot(aes(x = TP)) +
  geom_bar(aes(fill = class), show.legend = FALSE) +
  scale_fill_manual(values = c(shake_omics_color)) +
  scale_x_continuous(breaks = c(0, 30, 60, 120, 240),
                     labels = c(0, 30, 60, 120, 240)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Sample number", x = "") + 
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot2

# ggsave(plot, filename = "omics_plot.pdf", width = 7, height = 7)
# ggsave(plot1, filename = "omics_plot1.pdf", width = 3, height = 7)
# ggsave(plot2, filename = "omics_plot2.pdf", width = 7, height = 2)

#####subject_info
masstools::setwd_project()
load("data/ensure_shake_study/subject_info/subject_info")
setwd("data/ensure_shake_study/study_information/")

dim(subject_info)

df <- 
  data.frame(subject_info, x = 1, y = 1, 
             factors = factor(subject_info$subject_id, 
                              levels = stringr::str_sort(subject_info$subject_id, numeric = TRUE)))

library(circlize)
circos.par(
  "track.height" = 0.05,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0.5, 20), 90),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = df$factors, 
                  x = df$x,
                  xlim = c(0.5,1.5))

##sspg
temp_value <- df$sspg
circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = range(temp_value, na.rm = TRUE),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.22,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    #plot country labels
    
    circos.text(
      x = mean(xlim),
      y = 300,
      labels = name,
      facing = "bending.outside",
      niceFacing = TRUE, 
      cex = 0.8,
      adj = aa
    )
    
    circos.yaxis(
      side = "left",
      at = c(min(temp_value, na.rm = TRUE),
             round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE))/2,2),
             max(temp_value, na.rm = TRUE)),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = temp_value[i],
      col = ggsci::pal_aaas()(n = 10)[6],
      bg.border = "black"
    )
  }
)

###age
temp_value <- df$age
circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = range(temp_value, na.rm = TRUE),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.22,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = c(min(temp_value, na.rm = TRUE),
             round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE))/2,2),
             max(temp_value, na.rm = TRUE)),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = temp_value[i],
      col = ggsci::pal_aaas()(n = 10)[5],
      bg.border = "black"
    )
  }
)


###bmi
temp_value <- df$bmi
circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = range(temp_value, na.rm = TRUE),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.22,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = c(min(temp_value, na.rm = TRUE),
             round((min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE))/2,2),
             max(temp_value, na.rm = TRUE)),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = temp_value[i],
      col = ggsci::pal_aaas()(n = 10)[4],
      bg.border = "black"
    )
  }
)

####sex
temp_col <- df$sex
temp_col[temp_col == "M"] <- sex_color["M"]
temp_col[temp_col == "F"] <- sex_color["F"]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_col[i],
      bg.border = "black"
    )
  }
)

#####ethnicity
temp_col <- df$ethnicity
temp_col[temp_col == "A"] <- ethnicity_color["A"]
temp_col[temp_col == "B"] <- ethnicity_color["B"]
temp_col[temp_col == "C"] <- ethnicity_color["C"]
temp_col[temp_col == "H"] <- ethnicity_color["H"]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_col[i],
      bg.border = "black"
    )
  }
)

circos.clear()


####sspg
plot <- 
df %>%
  ggplot(aes(y = "SSPG", x = sspg)) +
  geom_boxplot(color = ggsci::pal_aaas()(n = 10)[6]) +
  geom_jitter(size = 5, color = ggsci::pal_aaas()(n = 10)[6]) +
  labs(y = "", x = "") +
  base_theme
plot
ggsave(plot, filename = "sspg_plot.pdf", width = 6, height = 2)


####age
plot <- 
  df %>%
  ggplot(aes(y = "SSPG", x = age)) +
  geom_boxplot(color = ggsci::pal_aaas()(n = 10)[5]) +
  geom_jitter(size = 5, color = ggsci::pal_aaas()(n = 10)[5]) +
  labs(y = "", x = "") +
  base_theme
plot
ggsave(plot, filename = "age_plot.pdf", width = 6, height = 2)

####BMI
plot <- 
  df %>%
  ggplot(aes(y = "SSPG", x = bmi)) +
  geom_boxplot(color = ggsci::pal_aaas()(n = 10)[4]) +
  geom_jitter(size = 5,
              color = ggsci::pal_aaas()(n = 10)[4]) +
  labs(y = "", x = "") +
  base_theme
plot
ggsave(plot, filename = "bmi_plot.pdf", width = 6, height = 2)


####sex
plot <- 
  df %>%
  ggplot(aes(y = "SSPG")) +
  geom_bar(stat = "count", aes(fill = sex),
           show.legend = FALSE) +
  labs(y = "", x = "") +
  scale_fill_manual(values = sex_color) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  base_theme
plot
ggsave(plot, filename = "sex_plot.pdf", width = 6, height = 1.5)


####ethnicity
plot <- 
  df %>%
  ggplot(aes(y = "SSPG")) +
  geom_bar(stat = "count", aes(fill = ethnicity),
           show.legend = FALSE) +
  labs(y = "", x = "") +
  scale_fill_manual(values = ethnicity_color) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  base_theme
plot
ggsave(plot, filename = "ethnicity_plot.pdf", width = 6, height = 1.5)



