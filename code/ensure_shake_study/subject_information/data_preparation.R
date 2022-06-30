no_function()
masstools::setwd_project()
library(tidyverse)
rm(list = ls())

setwd("data/ensure_shake_study/subject_info/")
subject_info1 <- readr::read_csv("Shake_participants_Meta_rk.csv")
subject_info2 <- readr::read_csv("Shake_participants_Meta.csv")

subject_info <- subject_info1

subject_info <-
  subject_info %>%
  dplyr::select(
    subject_id = PID,
    sex = Sex,
    ethnicity = Ethnicity,
    sspg = SSPG,
    sspg_date = SSPG_date,
    age = AgeYears,
    weight = Weight_KG,
    bmi = Bmi
  )

subject_info3 <- readr::read_csv("Shake_participants.csv")

colnames(subject_info)

colnames(subject_info3)

save(subject_info, file = "subject_info")
