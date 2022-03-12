
no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

{
  ###wearbale
  
  ##sleep
  load("data/7_24_mike/sleep/data_preparation/sample_info")
  load("data/7_24_mike/sleep/data_preparation/variable_info")
  load("data/7_24_mike/sleep/data_preparation/expression_data")
  sleep_sample_info = sample_info
  sleep_variable_info = variable_info
  sleep_expression_data = expression_data
  
  ###cgm
  load("data/7_24_mike/cgm/data_preparation/sample_info")
  load("data/7_24_mike/cgm/data_preparation/variable_info")
  load("data/7_24_mike/cgm/data_preparation/expression_data")
  cgm_expression_data = expression_data
  cgm_sample_info = sample_info
  cgm_variable_info = variable_info
  
  ###hr
  load("data/7_24_mike/hr/data_preparation/sample_info")
  load("data/7_24_mike/hr/data_preparation/variable_info")
  load("data/7_24_mike/hr/data_preparation/expression_data")
  hr_sample_info = sample_info
  hr_variable_info = variable_info
  hr_expression_data = expression_data
  
  ###steps
  load("data/7_24_mike/steps/data_preparation/sample_info")
  load("data/7_24_mike/steps/data_preparation/variable_info")
  load("data/7_24_mike/steps/data_preparation/expression_data")
  
  steps_sample_info = sample_info
  steps_variable_info = variable_info
  steps_expression_data = expression_data  
}



###internal omics data
{
  ###cortisol
  load("data/7_24_mike/cortisol/data_preparation/sample_info")
  load("data/7_24_mike/cortisol/data_preparation/variable_info")
  load("data/7_24_mike/cortisol/data_preparation/expression_data")
  cortisol_expression_data = expression_data
  cortisol_sample_info = sample_info
  cortisol_variable_info = variable_info
  
  ###cytokine
  load("data/7_24_mike/cytokine/data_preparation/sample_info")
  load("data/7_24_mike/cytokine/data_preparation/variable_info")
  load("data/7_24_mike/cytokine/data_preparation/expression_data")
  cytokine_sample_info = sample_info
  cytokine_variable_info = variable_info
  cytokine_expression_data = expression_data
  
  ###food_log
  load("data/7_24_mike/food_log/data_preparation/sample_info")
  load("data/7_24_mike/food_log/data_preparation/variable_info")
  load("data/7_24_mike/food_log/data_preparation/expression_data")
  food_sample_info = sample_info
  food_variable_info = variable_info
  food_expression_data = expression_data
  
  dim(food_expression_data)
  
  ###lipidomics
  load("data/7_24_mike/lipidomics/data_preparation/sample_info")
  load("data/7_24_mike/lipidomics/data_preparation/variable_info")
  load("data/7_24_mike/lipidomics/data_preparation/expression_data")
  lipidomics_sample_info = sample_info
  lipidomics_variable_info = variable_info
  lipidomics_expression_datao = expression_data
  
  ###metabolic_panel
  load("data/7_24_mike/metabolic_panel/data_preparation/sample_info")
  load("data/7_24_mike/metabolic_panel/data_preparation/variable_info")
  load("data/7_24_mike/metabolic_panel/data_preparation/expression_data")
  metabolic_panel_sample_info = sample_info
  metabolic_panel_variable_info = variable_info
  metabolic_panel_expression_data = expression_data
  
  ###metabolomics
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  metabolomics_expression_data = expression_data
  
  load("data/7_24_mike/metabolomics/data_preparation/peaks/sample_info")
  load("data/7_24_mike/metabolomics/data_preparation/peaks/variable_info")
  load("data/7_24_mike/metabolomics/data_preparation/peaks/expression_data")
  metabolomicspeak_sample_info = sample_info
  metabolomicspeak_variable_info = variable_info
  metabolomicspeak_expression_data = expression_data
  
  ###cortisol
  metabolomicspeak_variable_info$mz
  metabolomicspeak_variable_info$rt
  
  ###metabolomics
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/sample_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/variable_info")
  load("data/7_24_mike/metabolomics/data_preparation/metabolites/expression_data")
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
  metabolomics_expression_data = expression_data
  
  ##remove QC and blank samples
  metabolomics_sample_info = 
    metabolomics_sample_info %>% 
    dplyr::filter(!is.na(accurate_time))
  
  metabolomics_expression_data = 
    metabolomics_expression_data[,metabolomics_sample_info$sample_id]
  
  ###proteomics
  load("data/7_24_mike/proteomics/data_preparation/sample_info")
  load("data/7_24_mike/proteomics/data_preparation/variable_info")
  load("data/7_24_mike/proteomics/data_preparation/expression_data")
  proteomics_sample_info = sample_info
  proteomics_variable_info = variable_info
  proteomics_expression_data = expression_data
  
  ###total protein
  load("data/7_24_mike/total_protein/data_preparation/sample_info")
  load("data/7_24_mike/total_protein/data_preparation/variable_info")
  load("data/7_24_mike/total_protein/data_preparation/expression_data")
  total_protein_sample_info = sample_info
  total_protein_variable_info = variable_info
  total_protein_expression_data = expression_data
  
}


###set work directory
setwd("data/7_24_mike/summary_info")

grep("Cortisol$",metabolomicspeak_variable_info$Compound.name, value = F)
grep("Corti",metabolomicspeak_variable_info$Compound.name, value = F)

x = cortisol_expression_data[1,,drop = FALSE]
y = metabolomicspeak_expression_data[9544,,drop = FALSE]

plot(
  as.numeric(x[1, intersect(colnames(x), colnames(y))]), 
  as.numeric(y[1, intersect(colnames(x), colnames(y))])  
)

##wearable information
##sleep
temp_data_sleep = 
  data.frame(accurate_time = sleep_sample_info$accurate_time,
             day = as.character(sleep_sample_info$day),
             time = sleep_sample_info$time, 
             hour = sleep_sample_info$hour,
             second = sleep_sample_info$seconds,
             value = as.character(sleep_expression_data[1,]))
##CGM
temp_data_cgm = 
  data.frame(accurate_time = cgm_sample_info$accurate_time,
             day = as.character(cgm_sample_info$day),
             time = cgm_sample_info$time, 
             value = as.numeric(cgm_expression_data[1,]))

##HR
temp_data_hr = 
  data.frame(accurate_time = hr_sample_info$accurate_time,
             day = as.character(hr_sample_info$day),
             time = hr_sample_info$time, 
             value = as.numeric(hr_expression_data[1,]))

##step
temp_data_step = 
  data.frame(accurate_time = steps_sample_info$accurate_time,
             day = as.character(steps_sample_info$day),
             time = steps_sample_info$time, 
             value = as.numeric(steps_expression_data[1,]))

##food
temp_data_food = 
  tinyTools::convert2long(expression_data = food_expression_data, 
                          sample_info = food_sample_info, 
                          variable_info = food_variable_info) %>% 
  dplyr::filter(value != 0)


head(temp_data_sleep,1)
head(temp_data_cgm,1)
head(temp_data_food,1)
head(temp_data_hr,1)
head(temp_data_step,1)

tail(temp_data_sleep,1)



####
all_accurate_time =
  unique(
    c(
      temp_data_sleep$accurate_time,
      temp_data_cgm$accurate_time,
      temp_data_hr$accurate_time,
      temp_data_step$accurate_time,
      temp_data_food$accurate_time
    )
  )

day_night_df =
  data.frame(time = all_accurate_time) %>%
  dplyr::mutate(day = lubridate::date(time),
                hour = lubridate::hour(time)) %>%
  dplyr::arrange(day, hour) %>%
  dplyr::group_by(day, hour) %>%
  dplyr::filter(time == min(time)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(start = time[hour == 6], end = time[hour == 18]) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(class = "day")

# save(all_accurate_time, file = "all_accurate_time")
# save(day_night_df, file = "day_night_df")
load("all_accurate_time")
load("day_night_df")

# data.frame(sleep_sample_info, value = as.character(sleep_expression_data[1,])) %>%
sleep_color =
  c("asleep" = alpha("#E41A1C", alpha = 0.3),
    "light" = alpha("#E41A1C", alpha = 0.5),
    "deep" = alpha("#E41A1C", alpha = 0.7),
    "rem" = alpha("#E41A1C", alpha = 1),
    "awake" = alpha("#377EB8", alpha = 0.4),
    "wake" = alpha("#377EB8", alpha = 0.7),
    "restless" = alpha("#377EB8", alpha = 1)
  )

temp_data_sleep$value = 
  factor(temp_data_sleep$value, 
         levels = c("asleep", "light", "deep", "rem",
                    "awake", "wake", "restless"))

library(scales)
plot_sleep = 
  ggplot() +
  geom_rect(
    aes(
      xmin = accurate_time,
      xmax = accurate_time + second,
      ymin = 0,
      ymax = 1,
      fill = value
    ),
    data = temp_data_sleep
  ) +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = 0,
      ymax = 1
    ),
    fill = "lightyellow",
    data = day_night_df,
    # alpha = 0.5, 
    show.legend = FALSE
  ) +
  labs(y = "Sleep", x = "") +
  scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    limits = c(min(all_accurate_time),
               max(all_accurate_time)),
    timezone = "America/Los_Angeles"
  ) +
  scale_fill_manual(values = sleep_color) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 90),
    legend.position = "top",
    legend.box = "horizontal",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = alpha("grey", 0.2)),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1)) 

plot_sleep

plot_cgm =
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_night_df,
    # alpha = 0.5, 
    show.legend = FALSE
  ) +
  geom_line(aes(x = accurate_time,
                y = value,
                group = 1),
            data = temp_data_cgm,
            show.legend = FALSE) +
  labs(y = "CGM", x = "") +
  scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    limits = c(min(all_accurate_time),
               max(all_accurate_time)),
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0))) +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 

plot_cgm

plot_cgm2 =
  ggplot() +
  # geom_rect(
  #   mapping = aes(
  #     xmin = start,
  #     xmax = end,
  #     ymin = -Inf,
  #     ymax = Inf
  #   ),
  #   fill = "lightyellow",
  #   data = day_night_df,
  #   # alpha = 0.5, 
  #   show.legend = FALSE
  # ) +
  geom_line(aes(x = time,
                y = value,
                color = day,
                group = 1),
            data = temp_data_cgm %>% dplyr::mutate(time = as.POSIXct(time)),
            show.legend = FALSE) +
  ggsci::scale_color_lancet() +
  labs(y = "CGM", x = "") +
  scale_x_datetime(
    breaks = date_breaks("6 hour"),
    date_labels = "%H:%M",
    limits = c(min(as.POSIXct(temp_data_cgm$time)),
               max(as.POSIXct(temp_data_cgm$time))),
    timezone = "America/Los_Angeles"
  ) +
  base_theme +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
  facet_grid(rows = vars(day), scales = "free")
  
plot_cgm2

library(patchwork)

plot_sleep + plot_cgm + plot_layout(ncol = 1)

plot_hr =
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_night_df,
    # alpha = 0.5, 
    show.legend = FALSE
  ) +
  geom_line(aes(x = accurate_time,
                y = value,
                group = 1),
            data = temp_data_hr[seq(1, nrow(temp_data_hr), 3),],
            show.legend = FALSE) +
  labs(y = "HR", x = "") +
  scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    limits = c(min(all_accurate_time),
               max(all_accurate_time)),
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0))) +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 

plot_hr

plot_step =
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_night_df,
    # alpha = 0.5, 
    show.legend = FALSE
  ) +
  geom_line(aes(x = accurate_time,
                y = value,
                group = 1),
            data = temp_data_step,
            show.legend = FALSE) +
  labs(y = "Step", x = "") +
  scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    limits = c(min(all_accurate_time),
               max(all_accurate_time)),
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0))) +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 

plot_step

library(patchwork)

food_color =
  c("Carbs_g" = ggsci::pal_jama()(n=10)[1], 
    "Fat_g" = ggsci::pal_jama()(n=10)[2],  
    "Protein_g" = ggsci::pal_jama()(n=10)[3], 
    "Alcohol_g"  = ggsci::pal_jama()(n=10)[4])

plot_food = 
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_night_df,
    # alpha = 0.5, 
    show.legend = FALSE
  ) +
  geom_point(aes(x = accurate_time,
                 y = value, 
                 fill = mol_name),
             size = 3,
             data = temp_data_food %>%
               dplyr::filter(mol_name %in% 
                               c("Carbs_g", "Fat_g", "Protein_g", "Alcohol_g")),
             show.legend = TRUE,
             shape = 22) +
  scale_fill_manual(values = food_color) +
  labs(y = "Food", x = "") +
  scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    limits = c(min(all_accurate_time),
               max(all_accurate_time)),
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1))) +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

plot_food

plot_wearable =
  plot_sleep + plot_cgm + plot_hr + plot_step + plot_food +
  patchwork::plot_layout(ncol = 1,
                         heights = c(0.5, 1, 1, 1, 1))
plot_wearable
# ggsave(plot_wearable, filename = "plot_wearbale.pdf", width = 15, height = 6)

######internal omics data
temp_data = 
  rbind(
    data.frame(time = total_protein_sample_info[, c("accurate_time")],
               class = "total_protein"),
    data.frame(time = cortisol_sample_info[, c("accurate_time")],
               class = "cortisol"),
    data.frame(time = cytokine_sample_info[, c("accurate_time")],
               class = "cytokine"),
    data.frame(time = lipidomics_sample_info[, c("accurate_time")],
               class = "lipidomics"),
    data.frame(time = metabolic_panel_sample_info[, c("accurate_time")],
               class = "metabolic_panel"),
    data.frame(time = metabolomics_sample_info[, c("accurate_time")],
               class = "metabolomics"),
    data.frame(time = proteomics_sample_info[, c("accurate_time")],
               class = "proteomics")
  ) %>% 
  as.data.frame() %>% 
  dplyr::filter(!is.na(time)) %>%
  dplyr::mutate(time = as.character(time)) %>% 
  dplyr::mutate(class = factor(class, levels = unique(class)))

plot1 = 
  temp_data %>% 
  ggplot(aes(time, class)) +
  geom_tile(aes(fill = class), color = "black", show.legend = FALSE) +
  scale_fill_manual(values = class_color) +
  labs(x = "Time", y = "") +
  scale_y_discrete(expand = expansion(mult = c(0,0))) +
  base_theme +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot1

# ggsave(plot1, file = "internal_omics.pdf", width = 15, height = 2)


plot2 = 
  temp_data %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(class = factor(class, levels = class)) %>% 
  ggplot(aes(x = n, y = class)) +
  geom_bar(stat = "identity", aes(fill = class), show.legend = FALSE) +
  labs(x = "Sample number", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0,0))) +
  scale_fill_manual(values = class_color) +
  geom_text(aes(x = n + 2, y = class, label = n)) +
  geom_text(aes(x = n + 2, y = class, label = n)) +
  base_theme +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

plot2

library(patchwork)

plot = 
  plot1 + plot2 + patchwork::plot_layout(ncol = 2, widths = c(3,1))
plot
# ggsave(plot, filename = "omics_sample.pdf", width = 14, height = 7)


plot2 = 
  temp_data %>% 
  dplyr::group_by(class) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(class = factor(class, levels = class)) %>% 
  ggplot(aes(x = n, y = class)) +
  geom_bar(stat = "identity", aes(fill = class), show.legend = FALSE) +
  labs(x = "Sample number", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  # scale_y_discrete(expand = expansion(mult = c(0,0))) +
  scale_fill_manual(values = class_color) +
  geom_text(aes(x = n + 2, y = class, label = n)) +
  geom_text(aes(x = n + 2, y = class, label = n)) +
  base_theme
# theme(axis.text.y = element_blank(),
#       axis.ticks.y = element_blank())

plot2

# ggsave(plot2, filename = "omics_sample_number.pdf", width = 7, height = 7)

dim(proteomics_expression_data)

date = 
  as.character(unique(lubridate::as_date(temp_data$time)))

plot = 
  data.frame(date) %>% 
  ggplot(aes(x = date, y = 1)) +
  geom_tile(aes(fill = date),
            show.legend = FALSE) +
  ggsci::scale_fill_uchicago() +
  labs(x = "", y = "") +
  scale_y_continuous(expand = expansion(mult = c(8,0))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_polar()

plot
# ggsave(plot = plot, file = "24_7_graph.pdf", width = 7, height = 7)


dim()


######sample collection overiview
temp_data = 
  rbind(
    data.frame(time = total_protein_sample_info[, c("accurate_time")],
               class = "total_protein"),
    data.frame(time = cortisol_sample_info[, c("accurate_time")],
               class = "cortisol"),
    data.frame(time = cytokine_sample_info[, c("accurate_time")],
               class = "cytokine"),
    data.frame(time = lipidomics_sample_info[, c("accurate_time")],
               class = "lipidomics"),
    data.frame(time = metabolic_panel_sample_info[, c("accurate_time")],
               class = "metabolic_panel"),
    data.frame(time = metabolomics_sample_info[, c("accurate_time")],
               class = "metabolomics"),
    data.frame(time = proteomics_sample_info[, c("accurate_time")],
               class = "proteomics")
  )

library(lubridate)
library(hms)
accurate_time = unique(temp_data$time)

time = 
  ymd_hms(accurate_time) %>%
  hms::as_hms() %>% 
  as.numeric()
# as.POSIXct()

time = time/60/60

temp_data = 
  data.frame(accurate_time, time) %>% 
  dplyr::mutate(day = lubridate::date(accurate_time),
                hour = lubridate::hour(accurate_time)) %>% 
  dplyr::arrange(day, time) %>% 
  dplyr::mutate(image = "blood.png")

segment_data =
  lubridate::ymd_hms(c("2019-04-29 00:00:00",
                       "2019-04-29 23:59:59",
                       "2019-04-30 23:59:59",
                       "2019-05-01 23:59:59",
                       "2019-05-02 23:59:59",
                       "2019-05-03 23:59:59",
                       "2019-05-04 23:59:59",
                       "2019-05-04 23:59:59",
                       "2019-05-05 23:59:59",
                       "2019-05-06 23:59:59",
                       "2019-05-07 23:59:59"), 
                     tz="America/Los_Angeles") %>% 
  data.frame(time = .)

library(ggimage)

plot = 
  ggplot(data = temp_data) +
  geom_rect(
    mapping = aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_night_df,
    show.legend = FALSE
  ) +
  geom_image(
    aes(x = accurate_time,
        y = time,
        image = image),
    asp = 2,
    size = 0.02,
    data = temp_data
  ) +
  geom_vline(aes(xintercept = time), data = segment_data, linetype = 2) +
  labs(y = "Time (houe)", x = "") +
  scale_x_datetime(
    breaks = date_breaks("12 hour"),
    date_labels = "%a %H:%M",
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(breaks = seq(0, 24, 2), 
                     labels = seq(0, 24, 2)
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 10
    ),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = alpha("grey", 0.2)),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0,
      unit = "pt"
    )
  ) 

plot
# ggsave(plot, filename = "sample_collection.pdf", width = 14, height = 7)


plot(temp_data$time, temp_data$day)
plot = 
  temp_data %>%
  dplyr::mutate(day = as.POSIXct(day)) %>%
  ggplot(aes(time, day)) +
  # annotate(geom = "rect", 
  #          xmin = 6, xmax = 18, 
  #          ymin = as.Date("2019-04-21"),
  #          ymax = as.Date("2019-05-10"),
  #          fill = "lightyellow") +
  geom_image(
    aes(x = time,
        y = day,
        image = image),
    asp = 2,
    size = 0.02,
    data = temp_data
  ) +
  
  scale_x_continuous(breaks = seq(0, 24, 2),
                     labels = seq(0, 24, 2)) +
  scale_y_date(date_breaks = "1 day",
               date_labels = "%A") +
  ggrepel::geom_text_repel(aes(label = round(time,2)), 
                           data = temp_data, ) +
  labs(x = "Time (hour of one day)", y = "") +
  base_theme 

plot

# ggsave(plot, filename = "sample_collection2.pdf", width = 14, height = 7)

#######sampling frequency
library(plyr)

day_class_info =
  day_night_df %>%
  dplyr::mutate(class = 1:9)

class = 
  temp_data$accurate_time %>%
  purrr::map(function(x) {
    class = which(x > day_class_info$start &
                    x < day_class_info$end)
    if(length(class) == 0){
      class = 0
    }
    class
  }) %>% 
  unlist()

temp_data$class = class  

library(gghalves)

temp_data2 =
  temp_data %>% 
  dplyr::filter(class != 0) %>% 
  plyr::dlply(.variables = .(class)) %>% 
  purrr::map(function(x){
    x %>% 
      dplyr::arrange(time) %>% 
      dplyr::pull(time) %>% 
      diff()
  }) %>% 
  unlist() %>% 
  data.frame(class = "yes", time = .) 

plot = 
  temp_data2 %>% 
  ggplot(aes(x = class, y = time)) +
  labs(x = "", y = "Sampling interval (hour)") +
  geom_half_violin() + 
  geom_half_boxplot(fill = "transparent", 
                    color = ggsci::pal_aaas()(n=10)[2],
                    outlier.shape = NA) + 
  geom_dotplot(binaxis = "y", method="histodot", stackdir="up") +
  base_theme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot

# ggsave(plot, filename = "sampling_frequency.pdf", width = 7, height = 7)


quantile(temp_data2$time)


#####frequency for wearable data
##CGM
temp = 
  cgm_sample_info %>%
  dplyr::arrange(accurate_time) %>%
  dplyr::pull(accurate_time) %>%
  diff() %>%
  as.numeric() %>%
  data.frame(time = .) %>% 
  dplyr::mutate(time = round(time, 3)) %>% 
  # dplyr::mutate(time = as.character(time)) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(time) %>% 
  dplyr::mutate(time = factor(as.character(time), levels = as.character(time))) 


library(ggbreak)

plot = 
  temp %>%
  ggplot() +
  geom_bar(aes(y = time, x = n), stat = "identity") +
  # scale_x_break(c(10, 2300)) +
  geom_text(aes(x = n + 0.4, y = time, label = n), size = 4) +
  base_theme +
  labs(y = "Sampling interval (minute)", x = "Number")
plot
ggsave(plot, filename = "cgm_sampling_frequency.pdf", width = 7, height = 7)

##HR
temp = 
  hr_sample_info %>%
  dplyr::arrange(accurate_time) %>%
  dplyr::pull(accurate_time) %>%
  diff() %>%
  as.numeric() %>%
  data.frame(time = .) %>% 
  dplyr::mutate(time = round(time, 3)) %>% 
  # dplyr::mutate(time = as.character(time)) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(time) %>% 
  dplyr::filter(n > 2) %>% 
  dplyr::mutate(time = factor(as.character(time), levels = as.character(time))) 

library(ggbreak)

plot = 
  temp %>%
  ggplot() +
  geom_bar(aes(y = time, x = n), stat = "identity") +
  # scale_x_break(c(10, 2300)) +
  geom_text(aes(x = n + 0.4, y = time, label = n), size = 4) +
  base_theme +
  labs(y = "Sampling interval (second)", x = "Number")
plot
ggsave(plot, filename = "hr_sampling_frequency.pdf", width = 7, height = 7)




##steps
temp = 
  steps_sample_info %>%
  dplyr::arrange(accurate_time) %>%
  dplyr::pull(accurate_time) %>%
  diff() %>%
  as.numeric() %>%
  data.frame(time = .) %>% 
  dplyr::mutate(time = round(time, 3)) %>% 
  # dplyr::mutate(time = as.character(time)) %>% 
  dplyr::group_by(time) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::arrange(time) %>% 
  # dplyr::filter(n > 2) %>% 
  dplyr::mutate(time = factor(as.character(time), levels = as.character(time))) 

library(ggbreak)

plot = 
  temp %>%
  ggplot() +
  geom_bar(aes(y = time, x = n), stat = "identity") +
  # scale_x_break(c(10, 2300)) +
  geom_text(aes(x = n + 0.4, y = time, label = n), size = 4) +
  base_theme +
  labs(y = "Sampling interval (minutes)", x = "Number")
plot
ggsave(plot, filename = "step_sampling_frequency.pdf", width = 7, height = 7)


####how many feature we get for each omics data
dim(proteomics_variable_info)
dim(metabolomics_variable_info)
dim(metabolic_panel_variable_info)
dim(lipidomics_variable_info)
dim(cytokine_variable_info)  
dim(cortisol_variable_info)
dim(total_protein_variable_info)

temp_data = 
  data.frame(
    number = c(nrow(total_protein_variable_info),
               nrow(cortisol_variable_info),
               nrow(cytokine_variable_info),
               nrow(lipidomics_variable_info),
               nrow(metabolic_panel_variable_info),
               nrow(metabolomics_variable_info),
               nrow(proteomics_variable_info)),
    class = c("total_protein", "cortisol", "cytokine", "lipidomics",
              "metabolic_panel", "metabolomics", "proteomics")) %>% 
  dplyr::mutate(class = factor(class, levels = class))


plot = 
  temp_data %>% 
  ggplot(aes(x = number, y = class)) +
  geom_segment(aes(x = 0, y = class, 
                   xend = number, yend = class,
                   color = class), show.legend = FALSE) +
  geom_point(aes(color = class), size = 4, show.legend = FALSE) +
  scale_color_manual(values = class_color) +
  base_theme +
  geom_text(aes(x = number, y = class, label = number), nudge_y = -0.2) +
  labs(x = "Feature number", y = "")

ggsave(plot, filename = "omics_feature_number.pdf", width = 7, height = 7)

