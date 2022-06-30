##total_proteins
no_function()

library(tidyverse)

###total_protein
masstools::setwd_project()
rm(list = ls())
source("code/tools.R")

load("data/24_7_study/total_protein/data_preparation/expression_data")
load("data/24_7_study/total_protein/data_preparation/sample_info")
load("data/24_7_study/total_protein/data_preparation/variable_info")

load("data/24_7_study/summary_info/day_night_df")
load("data/24_7_study/summary_info/all_accurate_time")

setwd("data/24_7_study/total_protein/data_overview")

expression_data[1,] %>% 
  as.numeric() %>% 
  density() %>% 
  plot()

expression_data =
  expression_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

day_night_df =
  day_night_df %>%
  dplyr::mutate(
    start_time = as.POSIXct(hms::as_hms(start)),
    end_time = as.POSIXct(hms::as_hms(end)),
    week = format(day, "%a")
  ) %>% 
  dplyr::mutate(week = paste(
    week,
    lubridate::month(day),
    lubridate::day(day),
    sep = "-"
  )) %>% 
  dplyr::mutate(week = factor(week, unique(week)))

##total_protein
temp_data_total_protein =
  data.frame(accurate_time = sample_info$accurate_time,
             day = as.character(lubridate::date(sample_info$accurate_time)),
             hour = lubridate::hour(sample_info$accurate_time),
             time = hms::as_hms(sample_info$accurate_time), 
             value = as.numeric(expression_data[1,])) %>% 
  dplyr::mutate(
    time = as.POSIXct(time),
    week = format(accurate_time, "%a")
  ) %>% 
  dplyr::mutate(week = paste(
    week,
    lubridate::month(day),
    lubridate::day(day),
    sep = "-"
  )) %>% 
  dplyr::mutate(week = factor(week, unique(week)))

library(plyr)
temp =
  temp_data_total_protein %>% plyr::dlply(.variables = .(day))

temp %>% 
  lapply(function(x){
    as.character(range(x$accurate_time))
  }) %>% 
  do.call(rbind, .)

temp %>% 
  lapply(function(x){
    as.character(range(x$time))
  }) %>% 
  do.call(rbind, .)

library(scales)

plot_total_protein1 =
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
            data = temp_data_total_protein,
            show.legend = FALSE) +
  geom_point(aes(x = accurate_time,
                y = value,
                fill = week),
             shape = 21,
             size = 3,
            data = temp_data_total_protein,
            show.legend = FALSE) +
  scale_fill_manual(values = week_color) +
  labs(y = "Value", x = "") +
  scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    limits = c(min(all_accurate_time),
               max(all_accurate_time)),
    timezone = "America/Los_Angeles"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0))) +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.line.x = element_blank(),
        # axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) 

plot_total_protein1

plot_total_protein2 =
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start_time,
      xmax = end_time,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_night_df,
    show.legend = FALSE
  ) +
  geom_line(aes(x = time,
                y = value,
                group = 1),
            data = temp_data_total_protein,
            show.legend = FALSE) +
  geom_point(aes(x = time,
                 y = value,
                 fill = week),
             shape = 21,
             size = 3,
             data = temp_data_total_protein,
             show.legend = FALSE) +
  labs(y = "Value", x = "") +
  scale_fill_manual(values = week_color) +
  scale_x_datetime(
    breaks = scales::date_breaks("2 hour"),
    date_labels = "%H:%M",
    # limits = c(as.POSIXct("1970-01-01 00:01:01 UTC"),
    #            as.POSIXct("1970-01-01 23:59:00 UTC")),
    expand = expansion(mult = c(0,0))
    # timezone = "America/Los_Angeles"
  ) +
  # scale_y_continuous(expand = expansion(mult = c(0,0))) +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.line.x = element_blank(),
        # axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
  facet_grid(rows = vars(week), scales = "free_y")

plot_total_protein2

plot_total_protein3 =
  ggplot() +
  geom_rect(
    mapping = aes(
      xmin = start_time,
      xmax = end_time,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow",
    data = day_night_df %>% 
      dplyr::filter(day == "2019-05-01"),
    show.legend = FALSE
  ) +
  geom_line(aes(x = time,
                y = value,
                group = week,
                color = week),
            data = temp_data_total_protein,
            show.legend = TRUE) +
  geom_point(aes(x = time,
                 y = value,
                 fill = week),
             shape = 21,
             size = 3,
             data = temp_data_total_protein,
             show.legend = TRUE) +
  labs(y = "Value", x = "") +
  scale_color_manual(values = week_color) +
  scale_fill_manual(values = week_color) +
  scale_x_datetime(
    breaks = scales::date_breaks("2 hour"),
    date_labels = "%H:%M",
    expand = expansion(mult = c(0,0))
  ) +
  # scale_y_continuous(expand = expansion(mult = c(0,0))) +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.line.x = element_blank(),
        legend.position = "top",
        panel.grid = element_blank(),
        panel.background = element_rect(fill = alpha("grey", 0.2)),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
  guides(color = guide_legend(nrow = 1))

plot_total_protein3

# ggsave(plot_total_protein1, filename = "plot_total_protein1.pdf", width = 14, height = 3)
# ggsave(plot_total_protein2, filename = "plot_total_protein2.pdf", width = 14, height = 7)
# ggsave(plot_total_protein3, filename = "plot_total_protein3.pdf", width = 14, height = 7)


