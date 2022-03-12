##steps
no_function()

library(tidyverse)

###sleep
sxtTools::setwd_project()
rm(list = ls())
source("R/tools.R")

load("data/7_24_mike/sleep/data_preparation/expression_data")
load("data/7_24_mike/sleep/data_preparation/sample_info")
load("data/7_24_mike/sleep/data_preparation/variable_info")

load("data/7_24_mike/summary_info/day_night_df")
load("data/7_24_mike/summary_info/all_accurate_time")

setwd("data/7_24_mike/sleep/data_overview")

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

##sleep
temp_data_sleep =
  data.frame(
    accurate_time = sample_info$accurate_time,
    day = as.character(sample_info$day),
    time = sample_info$time,
    hour = sample_info$hour,
    second = sample_info$seconds,
    value = as.character(expression_data[1,])
  ) %>%
  dplyr::mutate(
    time = as.POSIXct(time),
    value =
      factor(
        value,
        levels = c("asleep", "light", "deep", "rem",
                   "awake", "wake", "restless")
      ),
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
  temp_data_sleep %>% plyr::dlply(.variables = .(day))

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
plot_sleep1 = 
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
    show.legend = FALSE
  ) +
  labs(y = "Sleep", x = "") +
  scale_x_datetime(
    breaks = scales::date_breaks("4 hour"),
    date_labels = "%a %H:%M",
    # limits = c(min(all_accurate_time),
    #            max(all_accurate_time)),
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
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    panel.background = element_rect(fill = alpha("grey", 0.2)),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1)) 

plot_sleep1

plot_sleep2 = 
  ggplot() +
  geom_rect(
    aes(
      xmin = time, xmax = time + second, ymin = 0, ymax = 1, fill = value
    ),
    data = temp_data_sleep
  ) +
  geom_rect(
    mapping = aes(
      xmin = start_time,
      xmax = end_time,
      ymin = 0,
      ymax = 1
    ),
    fill = "lightyellow",
    data = day_night_df,
    show.legend = FALSE
  ) +
  labs(y = "", x = "") +
  scale_x_datetime(
    breaks = scales::date_breaks("2 hour"),
    date_labels = "%H:%M",
    # limits = c(as.POSIXct("1970-01-01 00:01:01 UTC"),
    #            as.POSIXct("1970-01-01 23:59:00 UTC")),
    expand = expansion(mult = c(0,0))
    # timezone = "America/Los_Angeles"
  ) +
  scale_fill_manual(values = sleep_color) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
    panel.background = element_rect(fill = alpha("grey", 0.2)),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  facet_grid(rows = vars(week), scales = "free_y")

plot_sleep2

# ggsave(plot_sleep2, filename = "plot_sleep2.pdf", width = 14, height = 7)


