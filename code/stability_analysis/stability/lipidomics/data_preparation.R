##
no_function()

masstools::setwd_project()
setwd("data/")
library(tidyverse)
rm(list = ls())

masstools::setwd_project()
setwd("data/stability_analysis/danny/lipidomics/")

ls()

data <- readr::read_csv("stability.lipids.csv")
data <- data[, -1]

data <-
  data %>%
  tibble::column_to_rownames(var = "CollapsedSampleName") %>%
  t() %>%
  as.data.frame()

plot(as.numeric(data[3,]))

as.numeric(data[2,]) %>%
  density() %>%
  plot()

data <-
  data %>%
  purrr::map(function(x) {
    2 ^ x
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

plot(density(as.numeric(data[3,])))

write.csv(data, "stability_lipids_xiaotao.csv",
          row.names = FALSE)
