no_source()

##RPLC
##positive mode
masstools::setwd_project()
setwd("data/ensure_shake/metabolomics/RPLC/POS/")
library(tidyverse)
library(data.table)
library(tidymass)
#
# peak_table = readr::read_csv("210414_Ryan_Shake_RPLC_pos.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "msDatabase_rplc0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# result_rp_pos25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = grep("mgf", dir("NCE25/"), value = TRUE),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE25/"
# )
# 
# save(result_rp_pos25, file = "NCE25/result_rp_pos25")
# 
# result_rp_pos50 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = grep("mgf", dir("NCE50/"), value = TRUE),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE50/"
# )
# 
# save(result_rp_pos50, file = "NCE50/result_rp_pos50")
# 
# load("NCE25/result_rp_pos25")
# load("NCE50/result_rp_pos50")
# 
# annotation_table =
#   metID::get_identification_table_all(result_rp_pos25,
#                                       result_rp_pos50,
#                                       candidate.num = 1)
# head(annotation_table)
# 
# plot =
#   summary_annotation_result(object = annotation_table)
# 
# plot
# 
# ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)
# 
# peak_table = readr::read_csv("210414_Ryan_Shake_RPLC_pos.csv")
# 
# annotation_table$name == peak_table$Compound
# 
# annotation_table =
#   cbind(annotation_table, peak_table)
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)
# 
# ###RPLC
# ###negative mode
# masstools::setwd_project()
# setwd("data/ensure_shake/metabolomics/RPLC/NEG/")
# library(tidyverse)
# library(data.table)
# library(tidymass)
# 
# peak_table = readr::read_csv("210414_Ryan_Shake_RPLC_neg.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "msDatabase_rplc0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "rp",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 10,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# 
# result_rp_neg25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE25/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE25/"
# )
# 
# save(result_rp_neg25, file = "NCE25/result_rp_neg25")
# 
# 
# result_rp_neg50 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE50/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE50/"
# )
# 
# save(result_rp_neg50, file = "NCE50/result_rp_neg50")
# 
# load("NCE25/result_rp_neg25")
# load("NCE50/result_rp_neg50")
# 
# annotation_table =
#   metID::get_identification_table_all(result_rp_neg25,
#                                       result_rp_neg50,
#                                       candidate.num = 1)
# head(annotation_table)
# 
# plot =
#   summary_annotation_result(object = annotation_table)
# plot
# ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)
# 
# peak_table = readr::read_csv("210414_Ryan_Shake_RPLC_neg.csv")
# 
# annotation_table$name == peak_table$Compound
# 
# annotation_table =
#   cbind(annotation_table, peak_table)
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)
# 
# 
# 
# 
# ###HILIC
# ###positive mode
# masstools::setwd_project()
# setwd("data/ensure_shake/metabolomics/HILIC/POS/")
# library(tidyverse)
# library(data.table)
# library(metID)
# 
# peak_table = readr::read_csv("210414_Ryan_Shake_HILIC_pos.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_hilic0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "positive",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# result_hilic_pos25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE25/", pattern = "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE25/"
# )
# 
# save(result_hilic_pos25, file = "NCE25/result_hilic_pos25")
# 
# result_hilic_pos35 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir(path = "NCE35/", pattern = "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE35/"
# )
# 
# save(result_hilic_pos35, file = "NCE35/result_hilic_pos35")
# 
# load("NCE25/result_hilic_pos25")
# load("NCE35/result_hilic_pos35")
# 
# 
# annotation_table =
#   metID::get_identification_table_all(result_hilic_pos25,
#                                       result_hilic_pos35,
#                                       candidate.num = 1)
# head(annotation_table)
# 
# plot =
#   summary_annotation_result(object = annotation_table)
# plot
# ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)
# 
# peak_table = readr::read_csv("210414_Ryan_Shake_HILIC_pos.csv")
# 
# annotation_table$name == peak_table$Compound
# 
# annotation_table =
#   cbind(annotation_table, peak_table)
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)
# # 
# # 
# # load("NCE25/msDatabase_hilic0.0.2")
# # load("NCE25/metlinDatabase0.0.2")
# # 
# # which_has_identification(object = result_hilic_pos25[[4]])
# # example_plot = 
# # ms2plot(
# #   object = result_hilic_pos25[[4]],
# #   database = metlinDatabase0.0.2,
# #   which.peak = which_has_identification(object = result_hilic_pos25[[5]])$MS1.peak.name[4]
# # )
# # example_plot
# # ggsave(example_plot, filename = "example_plot.pdf", width = 10, height = 7)
# 
# 
# ##HILIC
# ##negative mode
# masstools::setwd_project()
# setwd("data/ensure_shake/metabolomics/HILIC/NEG/")
# library(tidyverse)
# library(data.table)
# library(metID)
# 
# peak_table = readr::read_csv("210414_Ryan_Shake_HILIC_neg.csv")
# 
# peak_table =
#   peak_table %>%
#   dplyr::select(name = Compound, mz = `m/z`, rt = `Retention time (min)`) %>%
#   dplyr::mutate(rt = rt * 60)
# 
# write.csv(peak_table, "peak_table.csv", row.names = FALSE)
# 
# ###level 1
# param1 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 30,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "msDatabase_hilic0.0.2"
#   )
# 
# ###level 2
# param2 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbDatabase0.0.2"
#   )
# 
# param3 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "massbankDatabase0.0.2"
#   )
# 
# param4 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "metlinDatabase0.0.2"
#   )
# 
# param5 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "monaDatabase0.0.2"
#   )
# 
# param6 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "nistDatabase0.0.2"
#   )
# 
# param7 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "orbitrapDatabase0.0.1"
#   )
# 
# ##level 3
# param8 <-
#   identify_metabolites_params(
#     ms1.match.ppm = 25,
#     rt.match.tol = 1000000,
#     polarity = "negative",
#     ce = "all",
#     column = "hilic",
#     total.score.tol = 0.5,
#     candidate.num = 3,
#     threads = 8,
#     database = "hmdbMS1Database0.0.1"
#   )
# 
# 
# 
# result_hilic_neg25 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE25/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE25/"
# )
# 
# save(result_hilic_neg25, file = "NCE25/result_hilic_neg25")
# 
# 
# result_hilic_neg35 <- identify_metabolite_all(
#   ms1.data = "peak_table.csv",
#   ms2.data = dir("NCE35/", "mgf"),
#   parameter.list = c(param1, param2, param3, param4, param5, param6, param7),
#   path = "NCE35/"
# )
# 
# save(result_hilic_neg35, file = "NCE35/result_hilic_neg35")
# 
# load("NCE25/result_hilic_neg25")
# load("NCE35/result_hilic_neg35")
# 
# 
# annotation_table =
#   metID::get_identification_table_all(result_hilic_neg25,
#                                       result_hilic_neg35,
#                                       candidate.num = 1)
# head(annotation_table)
# 
# plot =
#   summary_annotation_result(object = annotation_table)
# plot
# ggsave(plot, filename = "annotation_table_summary.pdf", width = 10, height = 7)
# 
# peak_table = readr::read_csv("210209_CarT_Plasma_HILIC_neg.csv")
# 
# annotation_table$name == peak_table$Compound
# 
# annotation_table =
#   cbind(annotation_table, peak_table)
# 
# write.csv(annotation_table, file = "annotation_table.csv", row.names = FALSE)





###combine them together
masstools::setwd_project()
setwd("data/ensure_shake/metabolomics/")
rplc_pos_data = readr::read_csv("RPLC/POS/annotation_table.csv")
rplc_neg_data = readr::read_csv("RPLC/NEG/annotation_table.csv")
hilic_pos_data = readr::read_csv("hilic/POS/annotation_table.csv")
hilic_neg_data = readr::read_csv("hilic/NEG/annotation_table.csv")

##rename variable
rplc_pos_data$name = paste("rplc_pos",rplc_pos_data$name, sep = "_")
rplc_neg_data$name = paste("rplc_neg",rplc_neg_data$name, sep = "_")
hilic_pos_data$name = paste("hilic_pos",hilic_pos_data$name, sep = "_")
hilic_neg_data$name = paste("hilic_neg",hilic_neg_data$name, sep = "_")

all_name = unique(c(
  colnames(rplc_pos_data),
  colnames(rplc_neg_data),
  colnames(hilic_pos_data),
  colnames(hilic_neg_data)
))

temp_data <-
  data.frame(
    class = c("rplc_pos", "rplc_neg", "hilic_pos", "hilic_neg"),
    level1 = c(
      sum(rplc_pos_data$Level == 1, na.rm = TRUE),
      sum(rplc_neg_data$Level == 1, na.rm = TRUE),
      sum(hilic_pos_data$Level == 1, na.rm = TRUE),
      sum(hilic_neg_data$Level == 1, na.rm = TRUE)
    ),
    level2 = c(
      sum(rplc_pos_data$Level == 2, na.rm = TRUE),
      sum(rplc_neg_data$Level == 2, na.rm = TRUE),
      sum(hilic_pos_data$Level == 2, na.rm = TRUE),
      sum(hilic_neg_data$Level == 2, na.rm = TRUE)
    )
  )

plot = 
  temp_data %>% 
  tidyr::pivot_longer(cols = -class, 
                      names_to = "level", values_to = "number") %>% 
  dplyr::mutate(class = factor(class, levels = unique(class))) %>% 
  ggplot(aes(class, y = number, fill = level, label = number)) +
  geom_bar(aes(fill = level), stat = "identity") +
  scale_fill_manual(values = c("level1" = ggsci::pal_aaas()(n=10)[2],
                               "level2" = ggsci::pal_aaas()(n=10)[1])) +
  geom_text(size = 5, position = position_stack(vjust = 0.5), color = "white") +
  labs(x = "", y = "Metabolite number") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13))

plot

ggsave(plot, filename = "annotation_info.pdf", width = 7, height = 7)

unique(c(
  rplc_pos_data %>%
    dplyr::filter(!is.na(Level)) %>%
    dplyr::filter(Level != 3) %>%
    pull(Compound.name),
  rplc_neg_data %>%
    dplyr::filter(!is.na(Level)) %>%
    dplyr::filter(Level != 3) %>%
    pull(Compound.name),
  hilic_pos_data %>%
    dplyr::filter(!is.na(Level)) %>%
    dplyr::filter(Level != 3) %>%
    pull(Compound.name),
  hilic_neg_data %>%
    dplyr::filter(!is.na(Level)) %>%
    dplyr::filter(Level != 3) %>%
    pull(Compound.name)
)) %>% 
  length()


###combine them together
##rename
colnames(hilic_pos_data) = 
  colnames(hilic_pos_data) %>%
  stringr::str_replace(pattern = "\\_[0-9]{10,15}", "")

colnames(hilic_pos_data)[31:32] = c("5F", "50F")

colnames(hilic_neg_data) = 
  colnames(hilic_neg_data) %>%
  stringr::str_replace(pattern = "\\_[0-9]{10,15}", "")

colnames(hilic_neg_data)[31:32] = c("5F", "50F")

colnames(rplc_neg_data)[31:32] = c("5F", "50F")

colnames(rplc_pos_data)[31:32] = c("5F", "50F")

###RPLC
rplc_data = 
  rbind(rplc_pos_data, rplc_neg_data)

hilic_data = 
  rbind(hilic_pos_data, hilic_neg_data)

dim(rplc_data)
dim(hilic_data)

expression_data = 
  rbind(rplc_data,
        hilic_data)

variable_info = 
  expression_data %>% 
  dplyr::select(c(name:`Minimum CV%`, `Accepted Compound ID`:`Compound Link`))

variable_info =   
  variable_info %>% 
  dplyr::rename(variable_id = name)

expression_data = 
  expression_data %>% 
  dplyr::select(-c(`Accepted Compound ID`:`Compound Link`, name:`Minimum CV%`))

sample_info = data.frame(sample_id = colnames(expression_data))

sample_info = 
  sample_info %>% 
  dplyr::mutate(class = case_when(
    stringr::str_detect(sample_id, "QC") ~ "QC",
    stringr::str_detect(sample_id, "dln") ~ "QC_DL",
    stringr::str_detect(sample_id, "blank") ~ "Blank",
    TRUE ~ "Subject"
  ))

colnames(expression_data) == sample_info$sample_id
rownames(expression_data) <- variable_info$variable_id

variable_info

dir.create("data_preparation")
dir.create("data_preparation/peak")
save(expression_data, file = "data_preparation/peak/expression_data")
save(sample_info, file = "data_preparation/peak/sample_info")
save(variable_info, file = "data_preparation/peak/variable_info")

write.csv(expression_data, file = "data_preparation/peak/expression_data.csv", row.names = FALSE)
write.csv(sample_info, file = "data_preparation/peak/sample_info.csv", row.names = FALSE)
write.csv(variable_info, file = "data_preparation/peak/variable_info.csv", row.names = FALSE)

####only remain metabolite
variable_info = 
  variable_info %>% 
  dplyr::filter(!is.na(Compound.name)) %>% 
  dplyr::filter(Level != 3)

library(plyr)

variable_info = 
  variable_info %>% 
  plyr::dlply(.variables = .(Compound.name)) %>% 
  purrr::map(function(x){
    x = 
      x %>% 
      dplyr::filter(Level == min(Level))
    
    if(nrow(x) == 1){
      return(x)
    }
    
    x = 
      x %>% 
      dplyr::filter(SS == max(SS))
    x
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

rownames(expression_data)

expression_data = 
  expression_data[variable_info$variable_id,]

rownames(expression_data) = variable_info$variable_id

dir.create("data_preparation/metabolite")
save(expression_data, file = "data_preparation/metabolite/expression_data")
save(sample_info, file = "data_preparation/metabolite/sample_info")
save(variable_info, file = "data_preparation/metabolite/variable_info")

write.csv(expression_data, file = "data_preparation/metabolite/expression_data.csv", row.names = FALSE)
write.csv(sample_info, file = "data_preparation/metabolite/sample_info.csv", row.names = FALSE)
write.csv(variable_info, file = "data_preparation/metabolite/variable_info.csv", row.names = FALSE)

dim(variable_info)
table(variable_info$Level)





