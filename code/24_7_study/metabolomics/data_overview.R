##
no_function()

masstools::setwd_project()
source("code/tools.R")
setwd("data/24_7_study/metabolomics/data_preparation/")

load("metabolites/expression_data")
load("metabolites/variable_info")
load("metabolites/sample_info")

library(tidyverse)
library(data.table)

##remove blank samples
table(sample_info$class)
sample_info = 
sample_info %>% 
  dplyr::filter(class != "Blank")

expression_data = 
  expression_data[,sample_info$sample_id]

# ###remove metabolites from other databases
# table(variable_info$Database)
# 
# variable_info = 
#   variable_info %>% 
#   dplyr::filter(Database %in% c("hmdbDatabase0.0.2", "metlinDatabase0.0.2", "msDatabase_hilic0.0.2",
#                                 "msDatabase_rplc0.0.2", "nistDatabase0.0.2"))
# 
# expression_data =
#   expression_data[variable_info$variable_id, ]
# 
# dim(expression_data)

###remove zero > 50%
remain_idx = 
expression_data %>% 
  apply(1, function(x){
    sum(x == 0)/ncol(expression_data)
  }) %>% 
  "<"(0.5) %>% 
  which()

variable_info = 
  variable_info[remain_idx,]

expression_data = 
  expression_data[remain_idx,]

#######PCA to show the data quality
temp_data = 
expression_data %>% 
  apply(1, function(x){
    (x - mean(x))/sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

sum(is.na(temp_data))

pca_object = prcomp(x = t(as.matrix(temp_data)),
                    center = FALSE,
                    scale. = FALSE)

x = pca_object$x

x = x[,1:3] %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample_id") %>% 
  dplyr::left_join(sample_info, by = "sample_id")

library(ggfortify)

sample_info$day = as.character(sample_info$day)

plot = 
autoplot(object = pca_object,
         shape = "class",
         data = sample_info,
         fill = "class",
         size = 5,
         alpha = 0.5,
         frame = TRUE,
         frame.type = 'norm') +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  scale_shape_manual(values = c("QC" = 24,
                                "Subject" = 21)) +
  ggrepel::geom_text_repel(aes(PC1, PC2, 
                               label = ifelse(PC1 > 0.2, as.character(sample_id), NA))) +
  ggrepel::geom_text_repel(aes(PC1, PC2, 
                               label = ifelse(PC2 < -0.25, as.character(sample_id), NA))) +
  ggsci::scale_fill_jama() +
  base_theme 
plot
ggsave(plot, filename = "pca.pdf", width = 9, height = 7)









