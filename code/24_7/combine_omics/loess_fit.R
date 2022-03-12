no_function()

###
tinytools::setwd_project()
rm(list=ls())
source(here::here("R/tools.R"))
setwd("data/7_24_mike/combine_omics/data_preparation")


load("expression_data")
load("sample_info")
load("variable_info")

load("proteomics_variable_info")
load("metabolomics_variable_info")
load("metabolic_panel_variable_info")
load("lipidomics_variable_info")
load("cortisol_variable_info")
load("total_protein_variable_info")
load("cytokine_variable_info")


####fore each molecule, we just use the loess method to impute the space between the first sample and the last sample
####for each day

dim(expression_data)
library(plyr)
# new_expression_data = 
#   purrr::map(variable_info$variable_id, function(temp_variable_id) {
#     cat(temp_variable_id, " ")
#     temp_data =
#       data.frame(value = as.numeric(expression_data[temp_variable_id, ]),
#                  sample_info) %>%
#       plyr::dlply(.variables = .(day)) %>%
#       purrr::map(function(x) {
#         if(nrow(x) < 8){
#           return(NULL)
#         }
#         # cat(as.character(x$day[1]), " ")
#         x$new_time = as.numeric(x$time) / (60 * 60)
#         optimize_span =
#           optimize_loess_span(
#             x = x$new_time,
#             y = x$value,
#             span_range = c(0.3, 0.4, 0.5, 0.6)
#           )
#         
#         span =
#           optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
#         
#         value = x$value
#         new_time = x$new_time
#         
#         ls_reg =
#           loess(value ~ new_time,
#                 span = span)
#         
#         prediction_value =
#           predict(ls_reg,
#                   newdata = data.frame(new_time = seq(floor(min(
#                     x$new_time
#                   )) + 0.5,
#                   floor(max(
#                     x$new_time
#                   )) + 0.5,
#                   by = 0.5)))
#         
#         prediction = 
#           data.frame(new_time = seq(floor(min(
#             x$new_time
#           )) + 0.5,
#           floor(max(
#             x$new_time
#           )) + 0.5,
#           by = 0.5),
#           prediction = prediction_value) %>% 
#           dplyr::filter(!is.na(prediction))
#         
#         rbind(
#           data.frame(
#             variable_id = temp_variable_id,
#             time = prediction$new_time,
#             value = prediction$prediction,
#             type = "Predicted",
#             span = span,
#             day = x$day[1]
#           ),
#           data.frame(
#             variable_id = temp_variable_id,
#             time = x$new_time,
#             value = x$value,
#             type = "Real",
#             span = span,
#             day = x$day[1]
#           )
#         )
#       })
#     
#     plot_data = 
#       do.call(rbind, temp_data) %>% 
#       as.data.frame() 
#     
#     plot_data$accurate_time =
#       plot_data$time * 60 * 60 + 
#       lubridate::as_datetime(paste(plot_data$day, "00:00:00"), tz = "America/Los_Angeles")
#     
#     #####plot
#     plot =
#       plot_data %>%
#       ggplot(aes(accurate_time, value)) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[1]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[1])
#       ) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[2]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[2])
#       ) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[3]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[3])
#       ) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[4]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[4])
#       ) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[5]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[5])
#       ) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[6]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[6])
#       ) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[7]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[7])
#       ) +
#       geom_smooth(
#         se = FALSE,
#         span = plot_data$span[plot_data$day == unique(plot_data$day)[8]][1],
#         color = ggsci::pal_lancet()(n = 9)[1],
#         data = plot_data %>% dplyr::filter(type == "Real" & day == unique(plot_data$day)[8])
#       ) +
#       geom_point(aes(color = type), size = 2) +
#       ggsci::scale_color_d3() +
#       facet_grid(cols = vars(day), scales = "free_x") +
#       base_theme +
#       scale_x_datetime(
#         breaks = scales::date_breaks("6 hour"),
#         date_labels = "%H:%M",
#         timezone = "America/Los_Angeles"
#       ) 
#     
#     name = 
#       paste(plot_data$variable_id[1], "pdf", sep = ".") %>% 
#       stringr::str_replace("\\/", "_")
#     
#     ggsave(plot, filename = file.path("smooth_plot", name), width = 21, height = 5)
#     
#     purrr::map2(names(temp_data), .y = temp_data, function(x, y) {
#       if(is.null(y)){
#         return(y)
#       }
#       y = 
#         y %>% 
#         dplyr::filter(type == "Predicted") %>% 
#         dplyr::select(variable_id:value)
#       data.frame(day = x,  y)
#     }) %>% 
#       do.call(rbind, .) %>% 
#       as.data.frame()
#   }) %>% 
#   do.call(rbind, .) %>% 
#   as.data.frame()
# 
# dim(new_expression_data)
# head(new_expression_data)
# 
# new_expression_data$sample_id = 
#   paste(new_expression_data$day, new_expression_data$time, sep = "_")
# 
# new_sample_info = 
#   new_expression_data %>% 
#   dplyr::select(sample_id, day, time) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# new_expression_data = 
#   new_expression_data %>% 
#   dplyr::select(sample_id, variable_id, value) %>% 
#   tidyr::pivot_wider(names_from = sample_id, values_from = value) %>% 
#   column_to_rownames(var = "variable_id")
# 
# dim(new_sample_info)
# dim(new_expression_data)
# sum(is.na(new_expression_data))
# colnames(new_expression_data) == new_sample_info$sample_id
# 
# dim(new_expression_data)
# dim(variable_info)
# 
# rownames(new_expression_data) == variable_info$variable_id
# 
# new_variable_info = variable_info
# 
# time = 
#   as.POSIXct(new_sample_info$time * 60 * 60, 
#              tz = "PST",
#              origin = "2019-04-29") %>% 
#   hms::as_hms()
# 
# accurate_time = 
#   time +   
#   lubridate::as_datetime(paste(new_sample_info$day, "00:00:00"), tz = "America/Los_Angeles")
# 
# 
# new_sample_info$time = time
# new_sample_info$accurate_time = accurate_time
# 
# save(new_expression_data, file = "new_expression_data")
# save(new_sample_info, file = "new_sample_info")
# save(new_variable_info, file = "new_variable_info")


load('new_expression_data')
load("new_sample_info")
load("new_variable_info")

new_variable_info

dim(new_sample_info)
dim(sample_info)

sum(variable_info$variable_id == new_variable_info$variable_id)
dim(new_variable_info)

matched_idx = 
sample_info$accurate_time %>% 
  purrr::map(function(time1){
    diff_time_min = as.numeric(abs(
      difftime(
        time1 = time1,
        time2 = new_sample_info$accurate_time,
        units = "min"
      )
    ))
    idx = which(diff_time_min < 30)
    if(length(idx) > 1){
      idx = idx[which.min(diff_time_min[idx])][1]
    }
    if(length(idx) == 0){
      return(NA)
    }
    idx
  }) %>% 
  unlist()

matched_idx = 
  data.frame(idx1 = 1:nrow(sample_info),
             idx2 = matched_idx,
             time1 = sample_info$accurate_time,
             time2 = new_sample_info$accurate_time[matched_idx],
             day1 = sample_info$day,
             day2 = new_sample_info$day[matched_idx]
             ) %>% 
  dplyr::filter(!is.na(idx2))


cor_data = 
purrr::map(
  1:nrow(new_expression_data),
  .f = function(idx) {
    cat(idx, " ")
    temp_data =
      data.frame(value1 = as.numeric(expression_data[idx, matched_idx$idx1]),
                 value2 = as.numeric(new_expression_data[idx, matched_idx$idx2]),
                 day = matched_idx$day1)
    
    test =
      cor.test(temp_data$value1, temp_data$value2, method = "spearman")
    
    result = 
    temp_data %>% 
      plyr::dlply(.variables = .(day)) %>% 
      purrr::map(function(x){
        test =
          cor.test(x$value1, x$value2, method = "spearman")
        
        data.frame(
          variable_id = rownames(expression_data)[idx],
          cor = unname(test$estimate),
          p = unname(test$p.value),
          day = x$day[1]
        )
      }) %>% 
      dplyr::bind_rows()
    
    test =
      cor.test(temp_data$value1, temp_data$value2, method = "spearman")
    
    result2 = 
    data.frame(
      variable_id = rownames(expression_data)[idx],
      cor = unname(test$estimate),
      p = unname(test$p.value),
      day = "all"
    )
    
    result$day = as.character(result$day)
    rbind(result, result2)
  }
) %>% 
  dplyr::bind_rows()

plot = 
cor_data %>%
  # dplyr::filter(variable_id %in% c(paste("lipid", 1:10, sep = "_"))) %>%
  dplyr::group_by(variable_id) %>%
  dplyr::summarise(
    median = median(cor),
    upper = quantile(cor, 0.75),
    lower = quantile(cor, 0.25)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(new_variable_info, by = "variable_id") %>%
  dplyr::arrange(data_type, median) %>%
  dplyr::mutate(variable_id = factor(variable_id, levels = variable_id)) %>%
  ggplot() +
  geom_pointrange(aes(
    x = variable_id,
    y = median,
    ymin = lower,
    ymax = upper,
    color = data_type
  ), 
  fatten = 0.1,
  size = 0.5,
  alpha = 0.5) +
  scale_color_manual(values = class_color) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(x = "Variable ID", y = "Spearman correlation range")
plot
ggsave(plot, filename = "real_loess_comparison_each_feature.pdf", 
       height = 5, width = 16)








  






