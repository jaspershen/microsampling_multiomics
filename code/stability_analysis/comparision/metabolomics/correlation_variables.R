no_function()

masstools::setwd_project()
rm(list = ls())
library(tidyverse)
library(massdataset)

load("data/stability_analysis/metabolomics/data_preparation/metabolomics_data")

dir.create("data/stability_analysis/metabolomics/correlation_variables",
           recursive = TRUE)
setwd("data/stability_analysis/metabolomics/correlation_variables")

metabolomics_data <-
  metabolomics_data %>%
  mutate_median_intensity() %>%
  activate_mass_dataset(what = "variable_info") %>%
  rename(median_intensity_total = median_intensity) %>%
  mutate_median_intensity(
    according_to_samples = metabolomics_data %>%
      activate_mass_dataset(what = "sample_info") %>%
      filter(class == "M") %>% pull(sample_id)
  ) %>%
  activate_mass_dataset(what = "variable_info") %>%
  rename(median_intensity_m = median_intensity) %>%
  mutate_median_intensity(
    according_to_samples = metabolomics_data %>%
      activate_mass_dataset(what = "sample_info") %>%
      filter(class == "P") %>% pull(sample_id)
  ) %>%
  activate_mass_dataset(what = "variable_info") %>%
  rename(median_intensity_p = median_intensity) %>%
  dplyr::filter(median_intensity_total > 0) %>%
  dplyr::filter(!is.na(MSI.level)) %>%
  dplyr::filter(MSI.level > 0)

ratio <-
  extract_variable_info((metabolomics_data))$median_intensity_p / extract_variable_info((metabolomics_data))$median_intensity_m

plot(log(ratio, 2))

metabolomics_data <-
  metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  mutate(ratio = ratio)

metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  ggplot(aes(log(median_intensity_total, 10), log(ratio, 2))) +
  geom_point() +
  theme_bw()

cor.test(
  log(metabolomics_data@variable_info$median_intensity_m + 1, 10),
  log(metabolomics_data@variable_info$median_intensity_p + 1, 10),
  method = "spearman"
)

lm_reg <-
  lm(
    log(metabolomics_data@variable_info$median_intensity_m + 1, 10) ~ log(metabolomics_data@variable_info$median_intensity_p + 1, 10)
  )

summary(lm_reg)

m_data <-
  metabolomics_data %>%
  activate_mass_dataset(what = "expression_data") %>%
  select(contains("M")) %>%
  `+`(1) %>%
  log(10)

p_data <-
  metabolomics_data %>%
  activate_mass_dataset(what = "expression_data") %>%
  select(contains("P")) %>%
  `+`(1) %>%
  log(10)

colnames(m_data)
colnames(p_data)

plot(unlist(m_data[1, , drop = TRUE]),
     unlist(p_data[1, , drop = TRUE]))

cor.test(unlist(m_data[1, , drop = TRUE]),
         unlist(p_data[1, , drop = TRUE]), method = "spearman")

cor_data <-
  purrr::map(1:nrow(m_data), function(idx) {
    cat(idx, " ")
    result <-
      cor.test(unlist(m_data[idx, , drop = TRUE]),
               unlist(p_data[idx, , drop = TRUE]), method = "spearman")
    data.frame(cor = unname(result$estimate),
               p = unname(result$p.value))
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

metabolomics_data <-
  metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor = cor_data$cor,
                p = cor_data$p)

metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  # dplyr::filter(p < 0.05) %>%
  ggplot(aes(variable_id, cor)) +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

dim(m_data)

temp_data <-
  metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(p < 0.05)

plot <-
  temp_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(cor > 0) %>%
  ggplot(aes(log(median_intensity_m + 1, 10),
             log(median_intensity_p + 1, 10))) +
  geom_point(aes(fill = cor,
                 size = -log(p, 10)),
             shape = 21) +
  theme_bw() +
  guides(color = guide_colorbar(title = "log10(intensity)")) +
  labs(x = "Microsampling", y = "Plasma") +
  theme(panel.grid.minor = element_blank()) +
  scale_fill_gradient(low = alpha(ggsci::pal_aaas()(n = 10)[2], 0.1),
                      high = ggsci::pal_aaas()(n = 10)[2]) +
  scale_size_continuous(range = c(0.5, 3)) +
  geom_smooth(method = "lm", se = TRUE)
# ggforce:: facet_zoom(xlim = c(0,0.5), ylim = c(0, 1))

plot

cor.test(
  log(temp_data@variable_info$median_intensity_m + 1, 10),
  log(temp_data@variable_info$median_intensity_p + 1, 10),
  method = "spearman"
)

# ggsave(plot, filename = "variable_cor_plot.pdf", width = 4.5, height = 3.5)

median(temp_data@variable_info$cor)
mean(temp_data@variable_info$cor)
#
#
# cor_data_sample <-
#   purrr::map(1:ncol(m_data), function(idx) {
#     cat(idx, " ")
#     result <-
#       cor.test(m_data[, idx, drop = TRUE],
#                p_data[, idx, drop = TRUE],
#                method = "spearman")
#     data.frame(cor = unname(result$estimate),
#                p = unname(result$p.value))
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# m_data
#
# temp_data <-
#   data.frame(median_intensity_m = apply(log(m_data + 1, 10), 2, median),
#              median_intensity_p = apply(log(p_data + 1, 10), 2, median),
#              cor_data_sample)
#
# plot <-
#   temp_data %>%
#   ggplot(aes(log(median_intensity_m + 1, 10),
#              log(median_intensity_p + 1, 10))) +
#   geom_point(aes(fill = cor,
#                  size = -log(p, 10)),
#              shape = 21) +
#   theme_bw() +
#   labs(x = "Microsampling", y = "Plasma") +
#   theme(panel.grid.minor = element_blank()) +
#   scale_fill_gradient(low = alpha(ggsci::pal_aaas()(n=10)[2], 0.1),
#                       high = ggsci::pal_aaas()(n=10)[2]) +
#   scale_size_continuous(range = c(0.5, 3)) +
#   geom_smooth(method = "lm", se = TRUE)
#
#
# plot
#
# median(temp_data$cor)
#
# cor.test(temp_data$median_intensity_m, temp_data$median_intensity_p, method = "spearman")
#
#
# cor.test(-log(metabolomics_data$M1 + 1, 10), -log(metabolomics_data$P1 + 1, 10))
#
# temp_data <-
#   metabolomics_data %>%
#   `+`(1) %>%
#   log(10) %>%
#   pivot_longer() %>%
#   dplyr::left_join(extract_sample_info(metabolomics_data)[,c("sample_id", "subject_id", "class")],
#                    by = "sample_id") %>%
#   dplyr::select(-sample_id) %>%
#   pivot_wider(names_from = class, values_from = value)
#
#
# median(cor_data_sample$cor)
#
# all_cor_sample <-
# metabolomics_data %>%
#   `+`(1) %>%
#   log(10) %>%
# cor_mass_dataset(margin = "sample",
#                  method = "spearman",
#                  data_type = "longer",
#                  p_adjust_method = 'BH')
#
# all_cor_sample %>%
#   dplyr::filter(from == "P1" | to == "P1") %>%
#   dplyr::filter(stringr::str_extract(from, "M|P") != stringr::str_extract(to, "M|P"))
#
#
# temp_data %>%
#   dplyr::mutate(subject_id = factor(subject_id,
#                                    levels = stringr::str_sort(unique(
#                                      temp_data$subject_id
#                                    ),
#                                    numeric = TRUE))) %>%
# ggplot(aes(M, P)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth(method = "lm") +
#   theme_bw() +
#   facet_wrap(facets = vars(subject_id)) +
#   theme(panel.grid.minor = element_blank()) +
#   labs(x = "Microsampling", y = "Plasma")



# ######PCA to show if the person can be clustered
# library(massstat)
# pca_object <-
#   metabolomics_data %>%
#   `+`(1) %>%
#   log(10) %>%
#   scale() %>%
#   run_pca()
#
# pca_score_plot(
#   object = metabolomics_data,
#   pca_object = pca_object,
#   color_by = "subject_id",
#   frame = FALSE
# ) +
#   ggrepel::geom_text_repel(aes(label = sample_id))


sum(temp_data@variable_info$cor < 0.5)
sum(temp_data@variable_info$cor > 0.5)

library(metpath)




query_id <- 
  query_id[!is.na(query_id)]

query_id <- query_id[query_id != ""]
   
query_id <- 
  query_id %>% 
  stringr::str_split(pattern = "\\|") %>% 
  unlist() %>% 
  unique()


data("kegg_hsa_pathway", package = "metpath")
kegg_hsa_pathway
pathway_class = 
  metpath::pathway_class(kegg_hsa_pathway)
remain_idx =
  pathway_class %>%
  unlist() %>%
  stringr::str_detect("Disease") %>%
  `!`() %>%
  which()

remain_idx

pathway_database =
  kegg_hsa_pathway[remain_idx]

result = 
  enrich_kegg(query_id = query_id, 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = pathway_database, 
              p_cutoff = 0.05, 
              p_adjust_method = "BH", 
              threads = 3)

result@result$p_value_adjust

plot <- 
enrich_scatter_plot(object = result)


ggsave(plot, file = "pathway_enrichment.pdf", width = 7, height = 7)
