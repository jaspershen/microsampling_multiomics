##
no_function()


sxtTools::setwd_project()
library(tidyverse)
rm(list = ls())

###load data
load("data/shake_study/metabolome_data_analysis/data_preparation/expression_data")
load("data/shake_study/metabolome_data_analysis/data_preparation/sample_info")
load("data/shake_study/metabolome_data_analysis/data_preparation/variable_info")

dim(expression_data)
dim(sample_info)
dim(variable_info)

setwd(dir = "data/shake_study/metabolome_data_analysis/k_means_clustering/")


###annotation for each cluster
cluster1 <- readxl::read_xlsx("cluster_1/cluster1.xlsx") 
cluster2 <- readxl::read_xlsx("cluster_2/cluster2.xlsx") 
cluster3 <- readxl::read_xlsx("cluster_3/cluster3.xlsx") 
cluster4 <- readxl::read_xlsx("cluster_4/cluster4.xlsx") 
cluster5 <- readxl::read_xlsx("cluster_5/cluster5.xlsx") 
cluster6 <- readxl::read_xlsx("cluster_6/cluster6.xlsx") 
cluster7 <- readxl::read_xlsx("cluster_7/cluster7.xlsx") 

###functional annotation for different cluster
cluster1 %>% 
  dplyr::left_join(variable_info, by = "variable_id") %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Metabolite)

cluster2 %>% 
  dplyr::left_join(variable_info, by = "variable_id") %>% 
  dplyr::filter(!is.na(Level)) %>% 
  dplyr::filter(Level == 1 | Level == 2) %>% 
  dplyr::pull(Metabolite)

cluster1_info <- read.table("PIUMet/cluster1_piumet.txt", sep = "\t")
cluster2_info <- read.table("PIUMet/cluster2_piumet.txt", sep = "\t")
cluster3_info <- read.table("PIUMet/cluster3_piumet.txt", sep = "\t")
cluster4_info <- read.table("PIUMet/cluster4_piumet.txt", sep = "\t")
cluster5_info <- read.table("PIUMet/cluster5_piumet.txt", sep = "\t")
cluster6_info <- read.table("PIUMet/cluster6_piumet.txt", sep = "\t")
cluster7_info <- read.table("PIUMet/cluster7_piumet.txt", sep = "\t")

colnames(cluster1_info) <- 
  colnames(cluster2_info) <- 
  colnames(cluster3_info) <- 
  colnames(cluster4_info) <- 
  colnames(cluster5_info) <- 
  colnames(cluster6_info) <- 
  colnames(cluster7_info) <- 
  c("mz", "polarity", "fdr")


cluster1 %>%
  dplyr::left_join(variable_info, by = "variable_id") %>%
  dplyr::filter(!is.na(Level)) %>%
  dplyr::filter(Level == 1 | Level == 2) %>%
  dplyr::pull(Metabolite)




######cluster1
marker <- 
  cluster1 %>% 
  dplyr::left_join(variable_info, by = "variable_id") %>% 
  dplyr::select("variable_id", "mz", "rt", "Mode") %>% 
  dplyr::rename(name = variable_id) %>% 
  dplyr::mutate(name = case_when(
    stringr::str_detect(Mode, "p") ~ paste(name, "_POS", sep = ""),
    stringr::str_detect(Mode, "n") ~ paste(name, "_NEG", sep = "")
  )) %>% 
  dplyr::select(-Mode)


readPIUMet(
  path = "PIUMet/cluster1_piumet_output/",
  marker = marker,
  text = FALSE,
  layout = "kk",
  size_range = c(2, 8)
)

###cluster1
load("PIUMet/cluster1_piumet_output/Result/node_data")

# hmdb_id <-
#   node_data %>%
#   dplyr::filter(!is.null(HMDB_ID)) %>%
#   dplyr::filter(HMDB_ID != " " & HMDB_ID != "NULL") %>%
#   dplyr::pull(HMDB_ID) %>%
#   unique()
# 
# kegg_id <-
#   lapply(hmdb_id, function(x) {
#     metflow2::transID(
#       query = x,
#       from = "Human Metabolome Database",
#       to = "KEGG",
#       top = 1
#     )
#   }) %>%
#   do.call(rbind, .)
# 
# kegg_id
# 
# load("PIUMet/hmdbMS1Database0.0.1")
# 
# hmdb_data <- hmdbMS1Database0.0.1@spectra.info
# 
# kegg_id2 <-
#   hmdb_data$KEGG.ID[match(kegg_id$`Human Metabolome Database`, hmdb_data$HMDB.ID)]
# 
# kegg_id <- data.frame(kegg_id, kegg_id2, stringsAsFactors = FALSE)
# 
# colnames(kegg_id) <- c("HMDB_ID", "KEGG_ID1", "KEGG_ID2")
# 
# KEGG_ID <-
#   apply(kegg_id, 1, function(x) {
#     x <- as.character(x)
#     if (is.na(x[2]) & is.na(x[3])) {
#       return(NA)
#     }
#     
#     if (!is.na(x[2])) {
#       return(x[2])
#     }
#     
#     if (!is.na(x[3])) {
#       return(x[3])
#     }
#     
#     
#   })
# 
# kegg_id <-
#   kegg_id %>%
#   dplyr::mutate(KEGG_ID = KEGG_ID) %>%
#   dplyr::select(-c(KEGG_ID1, KEGG_ID2))
# 
# kegg_id
# 
# node_data <-
#   node_data %>%
#   dplyr::left_join(kegg_id, by = c("HMDB_ID"))
# 
# annotation_result <- node_data

# save(annotation_result, file = "PIUMet/cluster1_piumet_output/Result/annotation_result")
load("PIUMet/cluster1_piumet_output/Result/annotation_result")

######pathway enrichment
######
######
######
##cluster1
load("PIUMet/cluster1_piumet_output/Result/annotation_result")
load("PIUMet/cluster1_piumet_output/Result/edge_data")

edge_data <-
  edge_data %>%
  dplyr::filter(stringr::str_detect(from, "POS") |
                  stringr::str_detect(from, "NEG") |
                  stringr::str_detect(to, "POS") |
                  stringr::str_detect(to, "NEG")
  ) %>%
  dplyr::select(from, to) %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(.f = function(x){
    if(stringr::str_detect(x[1], "POS") | stringr::str_detect(x[1], "NEG")){
      return(x)
    }else{
      return(rev(x))
    }
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  dplyr::arrange(V1) %>%
  dplyr::distinct()

colnames(edge_data) <- c('peak', "metabolite")

annotation_result <-
  annotation_result %>%
  dplyr::filter(node_class == "Metabolite") %>%
  dplyr::select(node, HMDB_ID, KEGG_ID, super.class)

annotation_result_cluster1 <-
  edge_data %>%
  dplyr::left_join(annotation_result, by = c("metabolite" = "node"))

# save(annotation_result_cluster1,
#      file = "PIUMet/cluster1_piumet_output/Result/annotation_result_cluster1")
load("PIUMet/cluster1_piumet_output/Result/annotation_result_cluster1")


########functional annotation
load("PIUMet/hsa_pathway")

##cluster1
load("PIUMet/cluster1_piumet_output/Result/annotation_result_cluster1")

kegg_id <-
  annotation_result_cluster1$KEGG_ID

kegg_id <-
  kegg_id[!is.na(kegg_id)]

enrichment_kegg_cluster1 <-
  enrich_pathway(id = kegg_id,
                 pathway_database = hsa_pathway)

# save(enrichment_kegg_cluster1, file = "PIUMet/cluster1_piumet_output/Result/enrichment_kegg_cluster1")
load("PIUMet/cluster1_piumet_output/Result/enrichment_kegg_cluster1")



enrichment_kegg_cluster1 <-
  enrichment_kegg_cluster1 %>%
  dplyr::arrange(p.value.fdr)


intersect(enrichment_kegg_cluster1$Pathway.name, 
          enrichment_kegg_cluster2$Pathway.name)

result_cluster <-
  list(
    enrichment_kegg_cluster1,
    enrichment_kegg_cluster2,
    enrichment_kegg_cluster3,
    enrichment_kegg_cluster4,
    enrichment_kegg_cluster5,
    enrichment_kegg_cluster6,
    enrichment_kegg_cluster7
  )

unique_name <- 
  unique(c(result_cluster[[1]]$Pathway.name,
           result_cluster[[2]]$Pathway.name,
           result_cluster[[3]]$Pathway.name,
           result_cluster[[4]]$Pathway.name,
           result_cluster[[5]]$Pathway.name,
           result_cluster[[6]]$Pathway.name,
           result_cluster[[7]]$Pathway.name)) %>% 
  sort()



for(x in unique_name){
  cat(x, "\n")
  idx1 <- which(result_cluster[[1]]$Pathway.name == x)
  idx2 <- which(result_cluster[[2]]$Pathway.name == x)
  idx3 <- which(result_cluster[[3]]$Pathway.name == x)
  idx4 <- which(result_cluster[[4]]$Pathway.name == x)
  idx5 <- which(result_cluster[[5]]$Pathway.name == x)
  idx6 <- which(result_cluster[[6]]$Pathway.name == x)
  idx7 <- which(result_cluster[[7]]$Pathway.name == x)
  
  temp_data <-
    rbind(
      data.frame(result_cluster[[1]][idx1, c(1:7)], class = rep("no", length(idx1)),cluster = rep(1, length(idx1))),
      data.frame(result_cluster[[2]][idx2, c(1:7)], class = rep("increase", length(idx2)),cluster = rep(2, length(idx2))),
      data.frame(result_cluster[[3]][idx3, c(1:7)], class = rep("increase", length(idx3)),cluster = rep(3, length(idx3))),
      data.frame(result_cluster[[4]][idx4, c(1:7)], class = rep("decrease", length(idx4)),cluster = rep(4, length(idx4))),
      data.frame(result_cluster[[5]][idx5, c(1:7)], class = rep("no", length(idx5)),cluster = rep(5, length(idx5))),
      data.frame(result_cluster[[6]][idx6, c(1:7)], class = rep("no", length(idx6)),cluster = rep(6, length(idx6))),
      data.frame(result_cluster[[7]][idx7, c(1:7)], class = rep("decrease", length(idx7)),cluster = rep(7, length(idx7)))
    ) %>%
    dplyr::arrange(p.value.bh) %>% 
    dplyr::filter(class != "no")
  
  if(nrow(temp_data) == 1){
    next()
  }
  
  if(all(c("decrease", "increase") %in% unique(temp_data$class))){
    temp_data <- 
      temp_data %>% 
      dplyr::filter(class != temp_data$class[1])
    
    ##remove wrong pathways in other cluster
    for(i in 1:nrow(temp_data)){
      # cat(i, "")
      temp_cluster <- temp_data$cluster[i]
      temp_id <- temp_data$Pathway.ID[i]
      result_cluster[[temp_cluster]] <-
        result_cluster[[temp_cluster]] %>%
        dplyr::filter(Pathway.ID != temp_id)
    }
  }
}







library(openxlsx)
wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roman")

purrr::walk(
  .x = 1:length(result_cluster),
  .f = function(i) {
    addWorksheet(wb, sheetName = paste("Cluster", i, sep = ""))
    freezePane(wb = wb, sheet = i, firstRow = TRUE, firstCol = FALSE) 
    writeDataTable(wb, sheet = i, x = result_cluster[[i]],
                   colNames = TRUE, rowNames = FALSE)
  }
)

saveWorkbook(wb, "cluster_annotation.xlsx", overwrite = TRUE)





save(enrichment_kegg_cluster1, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster1")
save(enrichment_kegg_cluster2, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster2")
save(enrichment_kegg_cluster4, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster3")
save(enrichment_kegg_cluster5, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster4")
save(enrichment_kegg_cluster6, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster6")
save(enrichment_kegg_cluster8, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster8")
save(enrichment_kegg_cluster9, file = "PIUMet/piumet_output_cluster1/Result/enrichment_kegg_cluster9")


enrichment_kegg_cluster1 %>%
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) 


plot <- 
  enrichment_kegg_cluster1 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail(10) %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/cluster1_piumet_output/Result/cluster1_annotation.pdf",
       width = 7,
       height = 7)




plot <- 
  enrichment_kegg_cluster2 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail(10) %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster2/Result/cluster2_annotation.pdf",
       width = 7,
       height = 7)



plot <- 
  enrichment_kegg_cluster4 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail(10) %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number", nrow = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster4/Result/cluster4_annotation.pdf",
       width = 7,
       height = 7)


plot <- 
  enrichment_kegg_cluster5 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail() %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number")) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster5/Result/cluster5_annotation.pdf",
       width = 7,
       height = 7)


plot <- 
  enrichment_kegg_cluster6 %>% 
  dplyr::mutate(Pathway.name = stringr::str_replace(Pathway.name, " - Homo sapiens \\(human\\)", "")) %>% 
  dplyr::mutate(class = case_when(
    p.value.fdr < 0.05 ~ "< 0.05",
    TRUE ~ "> 0.05"
  )) %>% 
  dplyr::mutate(fdr = -log(p.value.fdr, 10)) %>% 
  dplyr::arrange(fdr) %>% 
  dplyr::mutate(Pathway.ID = factor(Pathway.ID, levels = Pathway.ID)) %>% 
  tail() %>% 
  ggplot(aes(x = fdr, y = Pathway.ID)) +
  geom_point(aes(size = Overlap, fill = class), shape = 21) +
  geom_segment(aes(x = 0, xend = fdr, y = Pathway.ID, yend = Pathway.ID)) +
  geom_text(aes(x = fdr, y = Pathway.ID, label = Pathway.name)) +
  scale_size_continuous(range = c(2,10)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(FDR adjusted P-values)", y = "") + 
  scale_fill_manual(values = c(
    "< 0.05" = ggsci::pal_aaas()(n=10)[2],
    "> 0.05" = "black"
  )) +
  guides(fill = guide_legend(title = "FDR adjusted P value", 
                             override.aes = list(size = 7)),
         size = guide_legend(title = "Metabolite number", nrow = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)
  )

plot

ggsave(plot,
       filename = "PIUMet/piumet_output_cluster6/Result/cluster6_annotation.pdf",
       width = 7,
       height = 7)

