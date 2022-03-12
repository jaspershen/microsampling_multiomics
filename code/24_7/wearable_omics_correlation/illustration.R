no_function()

library(tidyverse)
sxtTools::setwd_project()
rm(list = ls())

source("R/tools.R")
source("R/modified_dtw.R")
source("R/lagged_correlation.R")

setwd("Figures/other_figures")

set.seed(seed = 100)

x = sin(seq( 0,  2*pi,  0.05))
x = x + rnorm(n = length(x), 
              mean = mean(x)/5, 
              sd = sd(x)/5)
y = x + rnorm(n = length(x), 
              mean = mean(x)/5, 
              sd = sd(x)/5)

temp_data = 
purrr::map(
  .x = seq(-40, 40, 1),
  .f = function(i) {
    data1 = data.frame(value = x, index = 1:length(x),
                       class = "x", shift_time = i)
    data2 = data.frame(value = y, index = 1:length(x),
                       class = "y",
                      shift_time = i)
    
    data2$index = data2$index + i
  
    rbind(data1,
          data2)
    
  }
) %>% do.call(rbind, .) %>% 
  as.data.frame()

temp_data %>% 
  dplyr::filter(shift_time == -40) %>% 
  ggplot(aes(index, value)) +
  geom_line(aes(group = class, color = class))

library(ggplot2)
library(gganimate)

library(gapminder)
head(gapminder)

p <- 
  temp_data %>% 
  ggplot(
  aes(x = index, y = value, colour = class)
) +
  geom_line(aes(group = class, 
                color = class),
            show.legend = FALSE) +
  labs(x = "", 
       y = "") +
  ggsci::scale_color_aaas() +
  scale_x_continuous(limits = c(-40,167)) +
  base_theme +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
p

p = 
p + transition_time(shift_time) +
  labs(title = "Shift time: {frame_time}")

animate(p,height = 3, width = 7, units = "in", res = 150)

anim_save(filename = "line_shift.gif")
  


temp_data = 
  purrr::map(
    .x = seq(-40, 40, 1),
    .f = function(i) {
      data1 = data.frame(value = x, index = 1:length(x),
                         class = "x", shift_time = i)
      data2 = data.frame(value = y, index = 1:length(x),
                         class = "y",
                         shift_time = i)
      
      data2$index = data2$index + i
      
      data = 
      data1[,c("index", "value")] %>% 
        dplyr::left_join(data2[,c("index", "value")], by = "index") %>% 
        dplyr::filter(!is.na(value.x) & !is.na(value.y))
      c(shift_time = i, cor = cor(data$value.x, data$value.y))
    }
  ) %>% do.call(rbind, .) %>% 
  as.data.frame() %>% 
  dplyr::mutate(class = "1")

library(ggplot2)
library(gganimate)

p <- 
  temp_data %>% 
  ggplot(
    aes(x = shift_time, y = cor, colour = class)
  ) +
  geom_line(aes(group = class, 
                color = class),
            show.legend = FALSE) +
  labs(x = "Shift time", 
       y = "Correlation") +
  ggsci::scale_color_jama() +
  # scale_x_continuous(limits = c(-10,31)) +
  base_theme +
  theme(panel.grid = element_blank())
p

p = 
  p + transition_reveal(along = shift_time)

animate(p,height = 3, width = 7, units = "in", res = 150)

anim_save(filename = "cor_shift.gif")

