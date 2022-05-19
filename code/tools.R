###----------------------------------------------------------------------------
readPIUMet <-
  function(path = ".",
           marker,
           text = TRUE,
           layout = "kk",
           size_range = c(2, 5),
           width_range = c(0.2, 0.5),
           label.name) {
    library(ggraph)
    annotation_result <-
      read.table(
        file.path(
          path,
          "peaks_putative_metabolites_w10.0_b2.0_mu0.0005_R1.txt"
        ),
        sep = "\t",
        header = TRUE
      )
    
    if (nrow(annotation_result) == 0) {
      cat("No result.\n")
      return(NULL)
    }
    
    annotation_result <-
      annotation_result %>%
      dplyr::mutate(mz =
                      stringr::str_replace(mz.peak, "m/z=", "") %>%
                      as.numeric() %>%
                      round(4)) %>%
      dplyr::mutate(mz2 = as.character(mz))
    
    
    marker <-
      marker %>%
      dplyr::mutate(polarity = case_when(
        stringr::str_detect(name, "POS") ~ "positive",
        stringr::str_detect(name, "NEG") ~ "negative"
      )) %>%
      dplyr::mutate(mz2 = case_when(
        polarity == "positive" ~ as.character(round(mz, 4) - 1),
        polarity == "negative" ~ as.character(round(mz, 4) + 1)
      ))
    
    annotation_result <-
      annotation_result %>%
      dplyr::left_join(marker, by = "mz2") %>%
      dplyr::select(-c(mz.y)) %>%
      dplyr::rename(mz = mz.x)
    
    edge_attr <-
      read.table(
        file.path(path, "result_edge_frequency_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(edge = V1)
    
    edge_data <-
      read.table(
        file.path(path, "result_union_net_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(from = V1, to = V2) %>%
      dplyr::mutate(edge = paste(from, "(pp)", to, sep = " ")) %>%
      dplyr::left_join(edge_attr, by = "edge")
    
    node_data <-
      read.table(
        file.path(path, "result_node_frequency_w10.0_b2.0_mu0.0005_R1.txt"),
        sep = "\t",
        header = FALSE
      ) %>%
      dplyr::rename(node = V1,
                    node_class = V3,
                    HMDB_ID = V4) %>%
      dplyr::left_join(annotation_result[, c("name", "Metabolite.Name", "super.class")],
                       by = c("node" = "Metabolite.Name"))
    
    node <-
      node_data$node %>%
      stringr::str_replace("m/z=", "") %>%
      as.numeric() %>%
      round(4) %>%
      as.character()
    
    node <- marker$name[match(node, marker$mz2)]
    
    rename <-
      data.frame(name1 = node_data$node[!is.na(node)],
                 name2 = node[!is.na(node)])
    
    node_data$node <-
      sapply(node_data$node, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    edge_data$from <-
      sapply(edge_data$from, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    edge_data$to <-
      sapply(edge_data$to, function(x) {
        temp_idx <- match(x, rename$name1)
        if (is.na(temp_idx)) {
          return(x)
        } else{
          rename$name2[temp_idx]
        }
      }) %>%
      unname()
    
    
    edge_data$edge <-
      paste(edge_data$from, "(pp)", edge_data$to, sep = " ")
    
    node_data <-
      node_data %>%
      dplyr::select(-name) %>%
      dplyr::distinct()
    
    node_data$node_class[grep("Metabolite", node_data$node_class)] <-
      "Metabolite"
    
    graph <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(Degree = tidygraph::centrality_degree(mode = 'all'))
    
    fill <-
      c(
        "m/z Peak" = ggsci::pal_d3()(10)[2],
        "Metabolite" = ggsci::pal_aaas()(10)[3],
        "Protein" = "grey"
        # "Protein" = "tomato"
      )
    
    col <-
      c(
        "m/z Peak" = ggsci::pal_d3()(10)[2],
        "Metabolite" = ggsci::pal_aaas()(10)[3],
        "Protein" = "grey"
      )
    
    shape = c("m/z Peak" = 21,
              "Metabolite" = 22,
              # "Metabolite_others" = 22,
              "Protein" = 24)
    
    require(ggraph)
    if (text) {
      plot <-
        ggraph(graph,
               layout = layout) +
        geom_edge_link(
          aes(edge_width = V3),
          alpha = 1,
          color = "black",
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(
            size = Degree,
            fill = node_class,
            shape = node_class
          ),
          alpha = 1,
          show.legend = TRUE
        ) +
        scale_shape_manual(values = shape) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggraph::geom_node_text(
          aes(label = ifelse(node %in% label.name, node, NA)),
          color = "black",
          repel = TRUE,
          size = 3
        ) +
        ggraph::scale_edge_width(range = width_range) +
        scale_size_continuous(range = size_range) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = col) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
    } else{
      plot <-
        ggraph(graph,
               layout = layout) +
        geom_edge_link(
          aes(edge_width = V3),
          alpha = 1,
          color = "black",
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(
            size = Degree,
            fill = node_class,
            shape = node_class
          ),
          alpha = 1,
          show.legend = TRUE
        ) +
        scale_shape_manual(values = shape) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggraph::scale_edge_width(range = width_range) +
        scale_size_continuous(range = size_range) +
        scale_fill_manual(values = fill) +
        scale_color_manual(values = col) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "right",
          legend.background = element_rect(fill = "transparent", color = NA)
        )
    }
    
    
    
    output_path <- file.path(path, "Result")
    dir.create(output_path)
    save(edge_data, file = file.path(output_path, "edge_data"))
    save(node_data, file = file.path(output_path, "node_data"))
    save(graph, file = file.path(output_path, "graph"))
    save(annotation_result, file = file.path(output_path, "annotation_result"))
    
    ggsave(
      plot,
      filename = file.path(output_path, "graph_plog.pdf"),
      width = 7,
      height = 7
    )
    
    plot
  }






setGeneric(
  name = "enrich_pathway",
  def = function(id,
                 pathway_database,
                 method = c("hypergeometric", "fisher")) {
    id <- unique(id)
    method <- match.arg(method)
    
    all_id <-
      pathway_database %>%
      unlist() %>%
      unique() %>%
      unname() %>%
      as.character()
    
    ##remove the ID which is not in the all_id
    id <- id[which(is.element(id, all_id))]
    
    if (length(id) == 0) {
      return(NULL)
    }
    
    sig_id <- as.character(id)
    
    num_all <- length(all_id)
    num_sig <- length(sig_id)
    
    # --------------------------------------------------------------
    #Generating label matrix for detected metabolites
    # --------------------------------------------------------------
    
    all_matrix <- setlabel(all_id, pathway_database)
    
    # delete metabolite set
    all_matrix2 <-
      all_matrix[, colSums(all_matrix) != 0, drop = FALSE]
    
    # error handling
    if (ncol(all_matrix2) < 2) {
      # stop function
      return(NULL)
    }
    
    # --------------------------------------------------------------
    #Generating label matrix for significant metabolites
    # --------------------------------------------------------------
    
    sig_matrix <- setlabel(sig_id, pathway_database)
    sig_matrix <-
      sig_matrix[, colSums(all_matrix) != 0, drop = FALSE]
    
    # -------------------------------
    #Calculating  ORA
    # -------------------------------
    
    p_value <- rep(NA, ncol(all_matrix2))
    
    #for each pathway
    
    for (i in 1:ncol(all_matrix2)) {
      # ------------------------------------
      #Generating 2~2 table
      # -------------------------------------
      a1 <- sum(sig_matrix[, i])# significant and including pathway
      a2 <-
        sum(all_matrix2[, i]) - sum(sig_matrix[, i])# not significant and including pathway
      a3 <-
        length(sig_id) - a1# significant and not including pathway
      a4 <-
        (length(all_id) - length(sig_id)) - a2# not significant and not including pathway
      
      tab <- t(matrix(c(a1, a2, a3, a4), 2))
      
      if (method == "hypergeometric") {
        p_value[i] <-
          phyper(
            q = a1 - 1,
            m = sum(all_matrix2[, i]),
            n = num_all - sum(all_matrix2[, i]),
            k = num_sig,
            lower.tail = FALSE
          )
      } else{
        # ----------------------------------
        # Fisher's exact test
        # ----------------------------------
        check <- tryCatch({
          resfish <- fisher.test(tab, alternative = "greater")
        }, error = function(e) {
          NA
        })
        
        if (class(check) != "htest") {
          p_value[i] <- 1
        } else{
          resfish <- fisher.test(tab, alternative = "greater")
          p_value[i] <- resfish$p.value
        }
      }
      
    }
    
    # -----------------------
    #q-value
    # -----------------------
    p_value_bh <- p.adjust(p_value, method = "BH")
    p_value_fdr <- p.adjust(p_value, method = "fdr")
    
    # ----------------------
    #Result
    # ----------------------
    PQ <- data.frame(p_value, p_value_bh, p_value_fdr,
                     stringsAsFactors = FALSE)
    rownames(PQ) <- colnames(sig_matrix)
    colnames(PQ) <- c("p.value", "p.value.bh", "p.value.fdr")
    
    PQ <-
      PQ %>%
      rownames_to_column(., "Pathway.name.ID")
    
    ##calculate the impact of pathway
    info <- lapply(pathway_database, function(module) {
      overlap.number <- length(intersect(module, id))
      pathway.number <- length(module)
      c(pathway.number, overlap.number)
    })
    
    info <- do.call(rbind, info)
    colnames(info) <- c("Pathway.length", "Overlap")
    
    info <-
      info %>%
      as_tibble() %>%
      dplyr::mutate(Overlap.frac = Overlap / Pathway.length)
    
    rownames(info) <-
      names(pathway_database)
    
    info <-
      info %>%
      rownames_to_column(var = 'Pathway.name.ID')
    
    info <-
      info %>%
      dplyr::left_join(PQ, info, by = 'Pathway.name.ID')
    
    info <-
      info %>%
      dplyr::mutate(Pathway.name = unlist(lapply(stringr::str_split(info$Pathway.name.ID, ";"), function(x)
        x[1])),
        Pathway.ID = unlist(lapply(stringr::str_split(info$Pathway.name.ID, ";"), function(x)
          x[2]))) %>%
      dplyr::select(., Pathway.name, Pathway.ID, everything(), -Pathway.name.ID) %>%
      dplyr::filter(., Overlap != 0) %>%
      dplyr::arrange(., p.value)
    
    info <- info
    
  }
)



#' @title setlabel
#' @description setlabel
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @param M_ID M_ID.
#' @param M M
#' @return  The MSE analysis result.

setlabel <- function (M_ID, M) {
  temp <- lapply(M_ID, function(x) {
    unlist(lapply(M, function(y) {
      sum(is.element(x, y))
    }))
  })
  temp <- do.call(rbind, temp)
  rownames(temp) <- M_ID
  temp <- temp
}



volcano_plot <- function(fc,
                         p_value,
                         p.cutoff = 0.05,
                         fc.cutoff = 2,
                         theme = c("light", "dark"),
                         alpha = 0.7,
                         text = FALSE,
                         variable_id = NULL,
                         point.size = 1) {
  theme <- match.arg(theme)
  temp_data <- data.frame(
    fc = log(fc, 2),
    p_value = -log(p_value, 10),
    variable_id,
    stringsAsFactors = FALSE
  )
  temp_data <-
    temp_data %>%
    dplyr::mutate(
      class = case_when(
        p_value > -log(p.cutoff, 10) & fc > log(fc.cutoff, 2) ~ "Increase",
        p_value > -log(p.cutoff, 10) &
          fc < log(1 / fc.cutoff, 2) ~ "Decrease",
        TRUE ~ "No"
      )
    )
  
  plot <-
    temp_data %>%
    ggplot(aes(fc, p_value)) +
    geom_vline(xintercept = 0,
               color = "black",
               linetype = 2) +
    geom_hline(
      yintercept = -log(p.cutoff, 10),
      color = "#FB8072",
      linetype = 2
    ) +
    geom_point(
      shape = 16,
      aes(color = class),
      show.legend = FALSE,
      alpha = alpha,
      size = point.size
    ) +
    scale_color_manual(
      values = c(
        "Increase" = ggsci::pal_aaas()(10)[2],
        "Decrease" = ggsci::pal_aaas()(10)[1],
        "No" = "black"
      )
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 13),
      panel.grid.minor = element_blank()
    ) +
    xlab(label = expression(log[2](Fold ~ change))) +
    ylab(label = expression(log[10](FDR ~ adjusted ~ P ~ value)))
  
  if (theme == "light") {
    plot <- plot
  } else{
    plot <- plot + ggdark::dark_theme_classic() +
      theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.background = element_rect(fill = "#0099B47F"),
        strip.text = element_text(color = "white", size = 13)
      )
  }
  
  if (text) {
    plot <-
      plot +
      ggrepel::geom_text_repel(mapping = aes(fc, p_value, label = variable_id),
                               data = temp_data[temp_data$class == "Yes", ])
  }
  
  plot
  
}


plot_silhouette <- function(sil,
                            color = "red") {
  temp_data <-
    data.frame(
      cluster = sil[, 1],
      neighbor = sil[, 2],
      sil_width = sil[, 3],
      stringsAsFactors = FALSE
    )
  temp_data <-
    temp_data %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::arrange(desc(cluster), sil_width) %>%
    dplyr::mutate(index = 1:nrow(temp_data))
  
  plot <-
    temp_data %>%
    ggplot() +
    geom_bar(
      aes(
        y = sil_width,
        x = index,
        fill = cluster,
        color = cluster
      ),
      stat = "identity",
      show.legend = FALSE
    ) +
    geom_hline(yintercept = 0) +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    theme_classic() +
    labs(
      y = paste(
        "Silhousette width",
        "\nAverage silhousettle widht:",
        round(mean(temp_data$sil_width), 2)
      ),
      x = ""
    ) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    coord_flip()
  
  plot <-
    plot +
    ggplot2::annotate(
      geom = "text",
      y = 0,
      x = max(temp_data$index),
      label = paste("n =", nrow(temp_data)),
      color = "black",
      hjust = -0.5,
      vjust = -1,
      size = 4
    )
  
  cluster_num <- as.numeric(max(temp_data$cluster))
  title <-
    paste(cluster_num, "clusters Cj", "\n", "j: nj | avei<Cj Si")
  
  plot <-
    plot +
    ggplot2::annotate(
      geom = "text",
      y = max(temp_data$sil_width),
      x = max(temp_data$index),
      label = title,
      color = "black",
      hjust = 0,
      vjust = 0,
      size = 4
    )
  
  class <- temp_data$cluster %>%
    unique() %>%
    sort() %>%
    rev()
  
  cluster_num <-
    temp_data %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(cluster)) %>%
    dplyr::pull(n)
  
  cluster_num <-
    sapply(1:length(cluster_num), function(x) {
      if (x == 1) {
        cluster_num[x] / 2
      } else{
        tail(cumsum(cluster_num[1:(x - 1)]), 1)  +  cluster_num[x] / 2
      }
    })
  
  for (i in 1:length(class)) {
    label <-
      paste(class[i],
            ":",
            sum(temp_data$cluster == class[i]),
            "|",
            round(mean(temp_data$sil_width[temp_data$cluster == class[i]]), 2))
    plot <-
      plot +
      ggplot2::annotate(
        geom = "text",
        y = max(temp_data$sil_width),
        x = cluster_num[i],
        label = label,
        color = "black",
        hjust = 0,
        vjust = 0,
        size = 4
      )
  }
  
  plot
}



base_theme =
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"),
    strip.text = element_text(size = 12)
  )





lipid_class_color =
  c("CE" = ggsci::pal_npg()(n=10)[1],
    "CER" = ggsci::pal_npg()(n=10)[2],
    "LPC" = ggsci::pal_npg()(n=10)[3],
    "PC" = ggsci::pal_npg()(n=10)[4],
    "PE" = ggsci::pal_npg()(n=10)[5],
    "PI" = ggsci::pal_npg()(n=10)[6],
    "SM" = ggsci::pal_npg()(n=10)[7],
    "DAG" = ggsci::pal_npg()(n=10)[8],
    "LPE" = ggsci::pal_npg()(n=10)[9],
    "TAG" = ggsci::pal_npg()(n=10)[10])



######
class_color =
  c(
    "lipidomics" = ggsci::pal_aaas()(10)[1],
    "metabolomics" = ggsci::pal_aaas()(10)[3],
    "cytokine" = ggsci::pal_aaas()(10)[4],
    "total_protein" = ggsci::pal_aaas()(10)[5],
    "cortisol" = ggsci::pal_aaas()(10)[6],
    "metabolic_panel" = ggsci::pal_aaas()(10)[7],
    "proteomics" = ggsci::pal_aaas()(10)[8]
  )

######
shake_omics_color =
  c(
    "Lipidomics" = ggsci::pal_aaas()(10)[1],
    "Metabolomics" = ggsci::pal_aaas()(10)[3],
    "Cytokine" = ggsci::pal_aaas()(10)[4]
  )


wearable_color =
  c(
    "sleep" = ggsci::pal_d3()(n = 10)[2],
    "cgm"  = ggsci::pal_d3()(n = 10)[6],
    "hr" = ggsci::pal_d3()(n = 10)[7],
    "step"  = ggsci::pal_d3()(n = 10)[9],
    "food" = ggsci::pal_d3()(n = 10)[10]
  )


sleep_color =
  c(
    "asleep" = alpha("#E41A1C", alpha = 0.3),
    "light" = alpha("#E41A1C", alpha = 0.5),
    "deep" = alpha("#E41A1C", alpha = 0.7),
    "rem" = alpha("#E41A1C", alpha = 1),
    "awake" = alpha("#377EB8", alpha = 0.4),
    "wake" = alpha("#377EB8", alpha = 0.7),
    "restless" = alpha("#377EB8", alpha = 1)
  )

week_color =
  c(
    "Mon-4-29" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1],
    "Tue-4-30" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[2],
    "Wed-5-1" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[3],
    "Thu-5-2" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[4],
    "Fri-5-3" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[5],
    "Sat-5-4" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[6],
    "Sun-5-5" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[7],
    "Mon-5-6" = RColorBrewer::brewer.pal(n = 9, name = "Set1")[8],
    "Tue-5-7" = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[4]
  )


tp_color <-
  c("0" = ggsci::pal_d3()(n=5)[1],
    "30" = ggsci::pal_d3()(n=5)[2],
    "60" = ggsci::pal_d3()(n=5)[3],
    "120" = ggsci::pal_d3()(n=5)[4],
    "240" = ggsci::pal_d3()(n=5)[5])

subject_col <-
  c("#543005", "#5F3606", "#6A3D07", "#764408", "#814B09", "#8D520A", "#975C12", "#A26519",
    "#AC6F20", "#B77927", "#C08431", "#C79141", "#CD9F51", "#D4AC62", "#DAB972", "#DDC584",
    "#D8CD9A", "#D3D5AF", "#CEDDC4", "#C9E5DA", "#BFE7E1", "#B1E1D9", "#A2DBD2", "#94D5CB",
    "#85CFC3", "#76C6BA", "#67BBB0", "#57AFA6", "#48A49B", "#389991", "#2D8F87", "#22857D",
    "#177B73", "#0D7169", "#02675F", "#005E55", "#00554C", "#004D42", "#004439", "#003C30")

names(subject_col) <- 
  c("S1",  "S2",  "S3",  "S4",  "S5",  "S6",  "S7",  "S8",  "S9",  "S10", "S11", "S12", "S13", "S14",
    "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28",
    "S29", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40")


sex_color <- 
  c("M" = ggsci::pal_aaas()(n = 10)[10],
    "F" = ggsci::pal_aaas()(n = 10)[2])


ethnicity_color <- 
  c("A" = ggsci::pal_nejm()(n=8)[1],
    "B" = ggsci::pal_nejm()(n=8)[2],
    "C" = ggsci::pal_nejm()(n=8)[3],
    "H" = ggsci::pal_nejm()(n=8)[5])

match_data = function(sample_info1,
                      expression_data1,
                      variable_info1,
                      sample_info2,
                      expression_data2,
                      variable_info2,
                      tol = 0.25) {
  idx =
    sample_info1$accurate_time %>%
    purrr::map(function(x) {
      temp_idx = which(abs(difftime(
        x, sample_info2$accurate_time, units = "hours"
      )) < tol)
      
      temp_idx[which.min(abs(x - sample_info2$accurate_time[temp_idx]))]
      
    })
  
  
  ##match plot
  #  temp_data =
  #    1:length(idx) %>%
  #    # 1:10 %>%
  #    purrr::map(function(x) {
  #      if(length(idx[[x]]) > 0){
  #        data.frame("internal_omics" = sample_info1$accurate_time[x],
  #                   "wearable" = sample_info2$accurate_time[idx[[x]]],
  #                   group = x)
  #      }else{
  #        data.frame(internal_omics = sample_info1$accurate_time[x],
  #                   wearable = NA,
  #                   group = x)
  #      }
  #
  #    }) %>%
  #    do.call(rbind, .) %>%
  #    as.data.frame() %>%
  #    dplyr::mutate(group = as.character(group))
  #
  #  segment_data = temp_data
  #
  #  temp_data =
  #    temp_data %>%
  #  tidyr::pivot_longer(cols = -group, names_to = "class", values_to = "time")
  #
  #
  #  plot =
  #    ggplot() +
  #    geom_rect(
  #      mapping = aes(
  #        xmin = start,
  #        xmax = end,
  #        ymin = -Inf,
  #        ymax = Inf
  #      ),
  #      fill = "lightyellow",
  #      data = day_night_df,
  #      # data = day_night_df[1:1,],
  #      show.legend = FALSE
  #    ) +
  #      geom_segment(data = segment_data,
  #                   aes(x = internal_omics, y = "internal_omics",
  #                       xend = wearable, yend = "wearable",
  #                       color = group),
  #                   show.legend = FALSE) +
  #    geom_point(data = temp_data,
  #      aes(x = time, y = class, fill = group),
  #      show.legend = FALSE,
  #      shape = 21,
  #      size = 3
  #    ) +
  #      scale_x_datetime(
  #        breaks = scales::date_breaks("4 hour"),
  #        date_labels = "%a %H:%M",
  #        timezone = "America/Los_Angeles"
  #      ) +
  #      labs(x = "", y = "") +
  #      base_theme +
  #    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
  #          axis.line.x = element_blank(),
  #          # axis.ticks.x = element_blank(),
  #          panel.grid = element_blank(),
  #          panel.background = element_rect(fill = alpha("grey", 0.2)),
  #          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  # plot
  #  ggsave(plot,
  #         filename = "match_plot.pdf",
  #         width = 14,
  #         height = 7)
  
  remain_idx =
    lapply(idx, length) %>%
    unlist() %>%
    `>`(0) %>%
    which()
  
  expression_data1 =
    expression_data1[, remain_idx]
  
  sample_info1 =
    sample_info1[remain_idx, ]
  
  idx = idx[remain_idx]
  
  # temp_data =
  #   1:length(idx) %>%
  #   purrr::map(function(x) {
  #     temp = expression_data2[, idx[[x]], drop = FALSE]
  #     data.frame(time = sample_info1$accurate_time[x],
  #                value = as.numeric(temp[1, ]))
  #   }) %>%
  #   do.call(rbind, .) %>%
  #   as.data.frame() %>%
  #   dplyr::mutate(time1 = as.character(time))
  #
  # rsd =
  # temp_data %>%
  #   dplyr::group_by(time) %>%
  #   dplyr::summarise(rsd = sd(value)*100/mean(value)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(time1 = as.character(time))
  #
  # plot1 =
  #   ggplot() +
  #   geom_rect(
  #     mapping = aes(
  #       xmin = start,
  #       xmax = end,
  #       ymin = -Inf,
  #       ymax = Inf
  #     ),
  #     fill = "lightyellow",
  #     data = day_night_df,
  #     show.legend = FALSE
  #   ) +
  #   geom_segment(data = rsd,
  #                aes(x = time, xend = time,
  #                    y = 0, yend = rsd, color = time1),
  #                show.legend = FALSE) +
  #   geom_point(data = rsd,
  #     aes(x = time, y = rsd, fill = time1),
  #     show.legend = FALSE,
  #     shape = 21,
  #     size = 3
  #   ) +
  #     scale_x_datetime(
  #       breaks = scales::date_breaks("4 hour"),
  #       date_labels = "%a %H:%M",
  #       timezone = "America/Los_Angeles"
  #     ) +
  #     labs(x = "", y = "RSD") +
  #     scale_y_continuous(expand = expansion(mult = c(0.02,0.1))) +
  #     base_theme +
  #     theme(axis.text.x = element_blank(),
  #           axis.ticks.x = element_blank(),
  #           axis.line.x = element_blank(),
  #           panel.grid = element_blank(),
  #           panel.background = element_rect(fill = alpha("grey", 0.2)),
  #           plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  #
  # plot2 =
  #   ggplot() +
  #   geom_rect(
  #     mapping = aes(
  #       xmin = start,
  #       xmax = end,
  #       ymin = -Inf,
  #       ymax = Inf
  #     ),
  #     fill = "lightyellow",
  #     data = day_night_df,
  #     show.legend = FALSE
  #   ) +
  #   geom_boxplot(aes(time, value, color = time1),
  #                data = temp_data,
  #                show.legend = FALSE, outlier.shape = NA) +
  #   geom_jitter(aes(time, value, color = time1),
  #               data = temp_data,
  #               show.legend = FALSE,
  #               size = 3,
  #               alpha = 0.5) +
  #   labs(x = "", y = "CGM") +
  #   scale_x_datetime(
  #     breaks = scales::date_breaks("12 hour"),
  #     date_labels = "%a %H:%M",
  #     timezone = "America/Los_Angeles"
  #   ) +
  #   scale_y_continuous(expand = expansion(mult = c(0.02,0.1))) +
  #   base_theme +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
  #         axis.line.x = element_blank(),
  #         # axis.ticks.x = element_blank(),
  #         panel.grid = element_blank(),
  #         panel.background = element_rect(fill = alpha("grey", 0.2)),
  #         plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  # plot2
  #
  # library(patchwork)
  #
  # plot =
  # plot1 + plot2 + patchwork::plot_layout(ncol = 1, heights = c(1, 3))
  # plot
  # ggsave(plot,
  #        filename = "proteomics_vs_cgm.pdf",
  #        width = 14,
  #        height = 7)
  
  # ##rsd plot
  # plot =
  # seq(0, 30, by = 0.1) %>%
  #   purrr::map(function(x){
  #     sum(rsd$rsd < x)*100/nrow(rsd)
  #   }) %>%
  #   unlist() %>%
  #   data.frame(rsd = seq(0, 30, by = 0.1),
  #              percent = .) %>%
  #   ggplot(aes(x = rsd, y = percent)) +
  #   geom_vline(xintercept = c(10,15), color = "red") +
  #   geom_line(aes(group = 1)) +
  #   labs(x = "RSD", y = 'Cumulative percentage (%)') +
  #   base_theme
  #
  # ggsave(plot,
  #        filename = "rsd_plot.pdf",
  #        width = 7,
  #        height = 7)
  
  
  expression_data2 =
    idx %>%
    purrr::map(function(x) {
      temp = expression_data2[, x, drop = FALSE]
      if (ncol(temp) == 1) {
        return(temp[, 1])
      } else{
        apply(temp, 1, mean)
      }
    }) %>%
    do.call(cbind, .) %>%
    as.data.frame()
  
  colnames(expression_data2) = colnames(expression_data1)
  
  sample_info2 = sample_info1
  return(
    list(
      sample_info1 = sample_info1,
      variable_info1 = variable_info1,
      expression_data1 = expression_data1,
      sample_info2 = sample_info2,
      variable_info2 = variable_info2,
      expression_data2 = expression_data2
    )
  )
}




cor_local = function(x,
                     y,
                     time,
                     min_time_point = 5,
                     method = c("pearson", "spearman")) {
  method = match.arg(method)
  cor_value =
    purrr::map(
      1:length(x),
      .f = function(idx) {
        temp_idx = seq(idx - round(min_time_point / 2),
                       idx + round(min_time_point / 2),
                       1)
        if (temp_idx[1] < 0) {
          temp_idx = temp_idx + abs(temp_idx[1]) + 1
        }
        
        if (temp_idx[length(temp_idx)] > length(x)) {
          temp_idx = temp_idx - (temp_idx[length(temp_idx)] - length(x))
        }
        
        temp_cor = cor.test(x[temp_idx], y[temp_idx], method = method)
        c(temp_cor$estimate, temp_cor$p.value)
      }
    ) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  cor_value =
    cor_value %>%
    dplyr::mutate(time = time) %>%
    dplyr::rename(p = V2) %>%
    dplyr::select(time, dplyr::everything())
  cor_value$p.adjust = p.adjust(cor_value$p, method = "BH")
  cor_value
}



lag_cor =
  function(expression_data1,
           sample_info1,
           variable_info1,
           expression_data2,
           sample_info2,
           variable_info2,
           lag = -5,
           method = "pearson") {
    sample_info1_old = sample_info1
    sample_info1$accurate_time =
      sample_info1$accurate_time + lag * 60 * 60
    temp_data =
      match_data(
        sample_info1 = sample_info1,
        expression_data1 = expression_data1,
        variable_info1 = variable_info1,
        sample_info2 = sample_info2,
        variable_info2 = variable_info2,
        expression_data2 = expression_data2
      )
    
    cor_value =
      purrr::map(
        1:nrow(variable_info1),
        .f = function(idx) {
          temp =
            cor.test(
              as.numeric(temp_data$expression_data2[1,]),
              as.numeric(temp_data$expression_data1[idx,]),
              method = "pearson"
            )
          c(temp$estimate, temp$p.value)
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    
    colnames(cor_value) = c("cor", 'p')
    
    cor_value$p_adjust =
      p.adjust(cor_value$p, method = "BH")
    
    cor_value$variable_id = temp_data$variable_info1$variable_id
    cor_value$lag = lag
    cor_value
  }



#####
time_plot = function(x,
                     y = NULL,
                     x_color = "blue",
                     y_color = "red",
                     x_name = "protein",
                     y_name = "CGM",
                     y_axis_name = "CGM",
                     time,
                     day_night_df,
                     add_point = FALSE) {
  x = data.frame(time, value = as.numeric(x),
                 stringsAsFactors = FALSE)
  if (!is.null(y)) {
    y = data.frame(time, value = as.numeric(y), 
                   stringsAsFactors = FALSE)
  }
  
  plot =
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
      show.legend = FALSE
    ) +
    geom_line(aes(x = time,
                  y = value,
                  color = x_name),
              # color = x_color,
              data = x) +
    labs(x = "", y = y_axis_name) +
    scale_x_datetime(
      breaks = scales::date_breaks("12 hour"),
      date_labels = "%a %H:%M",
      timezone = "America/Los_Angeles"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    base_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 10
      ),
      axis.line.x = element_blank(),
      # axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = alpha("grey", 0.2)),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0,
        unit = "pt"
      )
    )
  
  if (add_point) {
    plot =
      plot +
      geom_point(aes(x = time,
                     y = value,
                     color = x_name),
                 shape = 16,
                 data = x)
  }
  
  
  if (!is.null(y)) {
    plot =
      plot +
      geom_line(aes(x = time,
                    y = value,
                    color = y_name),
                data = y)
    
    
    if (add_point) {
      plot =
        plot +
        geom_point(aes(x = time,
                       y = value,
                       color = y_name),
                   shape = 16,
                   data = y)
    }
    
    color = c(unname(x_color),
              unname(y_color))
    names(color) = c(x_name, y_name)
  } else{
    color = c(unname(x_color))
    names(color) = c(x_name)
  }
  
  plot +
    scale_color_manual(values = color) +
    guides(color = guide_legend(title = "")) +
    theme(legend.position = "top")
  
}


summary_annotation_result = function(object) {
  object =
    object %>%
    dplyr::filter(!is.na(Compound.name)) %>%
    dplyr::filter(Level != 3)
  
  temp_data =
    rbind(
      object %>%
        dplyr::select(fill = Level) %>%
        dplyr::mutate(class = "level") %>%
        dplyr::group_by(fill) %>%
        dplyr::mutate(n = dplyr::n()) %>%
        dplyr::distinct() %>%
        dplyr::mutate(fill = as.character(fill)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(rate = round(n * 100 / sum(n), 2)),
      object %>%
        dplyr::select(fill = Database) %>%
        dplyr::mutate(class = "database") %>%
        dplyr::group_by(fill) %>%
        dplyr::mutate(n = dplyr::n()) %>%
        dplyr::distinct() %>%
        dplyr::ungroup() %>%
        dplyr::mutate(rate = round(n * 100 / sum(n), 2))
    ) %>%
    dplyr::mutate(class = factor(class, levels = c("level", "database")))
  
  class_color =
    object %>%
    dplyr::select(Level, Database) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Level)
  
  temp_data$fill = factor(temp_data$fill,
                          levels = unique(c(
                            class_color$Level, class_color$Database
                          )))
  
  level_color = ggsci::pal_aaas()(n = 10)[1:length(unique(class_color$Level))]
  names(level_color) = unique(class_color$Level)
  
  database_color =
    purrr::map(unique(class_color$Level), function(x) {
      temp =
        class_color %>%
        dplyr::filter(Level == x) %>%
        dplyr::pull(Database)
      alpha = seq(1, 0.1, length.out = length(temp))
      database_color = purrr::map(alpha, function(y) {
        ggplot2::alpha(level_color[x], y)
      }) %>%
        unlist()
      names(database_color) = temp
      database_color
    }) %>%
    unlist()
  
  temp_data$new_name = paste(temp_data$fill,
                             " (",
                             temp_data$n,
                             ",",
                             temp_data$rate,
                             "%)",
                             sep = "")
  names(level_color) =
    temp_data$new_name[match(names(level_color), temp_data$fill)]
  
  names(database_color) =
    temp_data$new_name[match(names(database_color), temp_data$fill)]
  
  new_level =
    temp_data$new_name[match(levels(temp_data$fill), temp_data$fill)]
  
  temp_data$fill = temp_data$new_name
  
  temp_data$fill = factor(temp_data$fill, levels = new_level)
  
  ggplot(temp_data,
         aes(x = class, y = n, fill = fill)) +
    geom_col(color = "white", position = 'stack') +
    # scale_x_discrete(limits = c(" ", "level","database")) +
    scale_fill_manual(values = c(level_color, database_color)) +
    theme_void() +
    guides(fill = guide_legend(title = "", ncol = 1)) +
    coord_polar("y")
}


dtw_twoway_plot = function(object,
                           shift = 5,
                           accurate_time,
                           day_night_df,
                           add_point = TRUE,
                           color = c("blue", "red")) {
  shift = max(object$query) - min(object$reference)
  temp_data =
    data.frame(accurate_time,
               query = object$query,
               reference = object$reference) %>%
    dplyr::mutate(reference = reference + shift) %>%
    dplyr::mutate(index = 1:length(object$query)) %>%
    tidyr::pivot_longer(
      cols = -c(index, accurate_time),
      names_to = "class",
      values_to = "value"
    )
  
  connect_df =
    data.frame(index1 = object$index1,
               index2 = object$index2)
  
  connect_df$time1 = accurate_time[connect_df$index1]
  connect_df$time2 = accurate_time[connect_df$index2]
  connect_df$y1 = object$query[connect_df$index1]
  connect_df$y2 = object$reference[connect_df$index2] + shift
  
  names(color) = c("reference", "query")
  
  plot =
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
      show.legend = FALSE
    ) +
    geom_line(
      data = temp_data,
      aes(
        x = accurate_time,
        y = value,
        group = class,
        color = class
      ),
      show.legend = FALSE
    ) +
    scale_color_manual(values = color) +
    geom_segment(aes(
      x = time1,
      xend = time2,
      y = y1,
      yend = y2
    ),
    data = connect_df,
    color = "grey") +
    scale_x_datetime(
      breaks = scales::date_breaks("12 hour"),
      date_labels = "%a %H:%M",
      timezone = "America/Los_Angeles"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    base_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 10
      ),
      axis.line.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = alpha("grey", 0.2)),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0,
        unit = "pt"
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    labs(x = "", y = "Value")
  
  if (add_point) {
    plot =
      plot +
      geom_point(
        data = temp_data,
        aes(
          x = accurate_time,
          y = value,
          group = class,
          color = class
        ),
        shape = 16,
        show.legend = FALSE
      )
  }
  
  plot
  
}













library(ggforce)

parlDiag <- function(Parties,
                     shares,
                     cols = NULL,
                     repr = c("absolute", "proportion")) {
  repr = match.arg(repr)
  stopifnot(length(Parties) == length(shares))
  if (repr == "proportion") {
    stopifnot(sum(shares) == 1)
  }
  if (!is.null(cols)) {
    names(cols) <- Parties
  }
  
  # arc start/end in rads, last one reset bc rounding errors
  cc <-
    cumsum(c(-pi / 2, switch(
      repr,
      "absolute" = (shares / sum(shares)) * pi,
      "proportion" = shares * pi
    )))
  cc[length(cc)] <- pi / 2
  
  # get angle of arc midpoints
  meanAngles <-
    colMeans(rbind(cc[2:length(cc)], cc[1:length(cc) - 1]))
  
  # unit circle
  labelX <- sin(meanAngles)
  labelY <- cos(meanAngles)
  
  # prevent bounding box < y=0
  labelY <- ifelse(labelY < 0.015, 0.015, labelY)
  
  p <- ggplot() + theme_no_axes() + coord_fixed() +
    expand_limits(x = c(-1.3, 1.3), y = c(0, 1.3)) +
    theme(panel.border = element_blank()) +
    theme(legend.position = "none") +
    geom_arc_bar(aes(
      x0 = 0,
      y0 = 0,
      r0 = 0.5,
      r = 1,
      start = cc[1:length(shares)],
      end = c(cc[2:length(shares)], pi / 2),
      fill = Parties
    ), color = "white") +
    switch(is.null(cols) + 1, scale_fill_manual(values = cols), NULL) +
    switch(is.null(cols) + 1, scale_color_manual(values = cols), NULL) +
    # for label and line positions, just scale sin & cos to get in and out of arc
    geom_path(aes(
      x = c(0.9 * labelX, 1.15 * labelX),
      y = c(0.9 * labelY, 1.15 * labelY),
      group = rep(1:length(shares), 2)
    ),
    colour = "white",
    size = 2) +
    geom_path(aes(
      x = c(0.9 * labelX, 1.15 * labelX),
      y = c(0.9 * labelY, 1.15 * labelY),
      group = rep(1:length(shares), 2)
    ), size = 1) +
    ggrepel::geom_label_repel(
      aes(
        x = 1.15 * labelX,
        y = 1.15 * labelY,
        color = Parties,
        label = switch(
          repr,
          "absolute" = sprintf("%s\n%i", Parties, shares),
          "proportion" = sprintf("%s\n%i%%", Parties, round(shares *
                                                              100))
        )
      ),
      fontface = "bold",
      label.padding = unit(1, "points")
    ) +
    
    geom_point(aes(x = 0.9 * labelX, y = 0.9 * labelY),
               colour = "white",
               size = 2) +
    geom_point(aes(x = 0.9 * labelX, y = 0.9 * labelY)) +
    geom_text(aes(x = 0, y = 0, label = switch(
      repr,
      "absolute" = (sprintf("Total: %i internal moleculars", sum(shares))),
      "proportion" = ""
    )),
    fontface = "bold",
    size = 7)
  
  return(p)
}






######get the correlation matrix
get_cor_matrix = function(data, 
                          c, 
                          scale = TURE,
                          method = c("spearman", "pearson"),
                          which = c("median", "mean")){
  method = match.arg(method)
  which = match.arg(which)
  data = 
    data %>% 
    apply(1, function(x){
      (x - mean(x))/sd(x)}) %>% 
    t() %>% 
    as.data.frame()
  
  ###get the distance for each two cluster
  cluster = 
    sort(unique(c$cluster))
  
  cluster[-length(cluster)] %>%
    purrr::map(function(i) {
      cluster1 = which(c$cluster == i)
      purrr::map((i + 1):length(cluster),
                 .f = function(j) {
                   cluster2 = which(c$cluster == j)
                   data1 = data[cluster1, ]
                   data2 = data[cluster2, ]
                   all_cor =
                     get_cor_for_each_variable(data1 = data1,
                                               data2 = data2,
                                               method = method)
                   # if (which == "median") {
                   #   value = median(all_cor)
                   # } else{
                   #   value = mean(all_cor)
                   # }
                   data.frame(from = i, to = j, cor = unname(all_cor))
                 } 
      ) %>% 
        do.call(rbind, .) %>% 
        as.data.frame()
    }) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
} 

get_cor_for_each_variable =
  function(data1, data2, method = c("spearman", "pearson")) {
    method = match.arg(method)
    purrr::map(as.data.frame(t(data1)), .f = function(x){
      purrr::map(as.data.frame(t(data2)), .f = function(y){
        cor(x, y, method = method)
      }) %>% 
        unlist()
    }) %>% 
      unlist()
  }



modularity_plot = function(subnetworks){
  plot <- 
    ggplot(
      data.frame(index = 1:length(subnetworks$modularity),
                 modu = subnetworks$modularity, stringsAsFactors = FALSE),
      aes(index, modu) 
    ) +
    geom_vline(xintercept = which.max(subnetworks$modularity), 
               linetype = 2, colour = "#800000B2") + 
    labs(x = "Community analysis iteration", y = "Modularity") +
    geom_line(colour = "black") +
    # geom_point() +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 12))
  
  plot <-
    plot + 
    ggplot2::annotate(geom = "point", 
                      x = which.max(subnetworks$modularity),
                      y = max(subnetworks$modularity), 
                      size = 3, 
                      colour = "red") +
    annotate(geom = "text", 
             x = which.max(subnetworks$modularity),
             y = max(subnetworks$modularity), 
             label = paste("(",  which.max(subnetworks$modularity),
                           ",", 
                           max(subnetworks$modularity) %>% round(3),
                           ")"),
             size = 5,
             colour = "red"
    )
  
  plot
}


optimize_loess_span =
  function(x, y, span_range = seq(0.2, 0.6, 0.1)) {
    span_rmse =
      purrr::map(span_range, function(span) {
        temp_data =
          data.frame(x, y)
        
        prediction =
          purrr::map(
            2:(nrow(temp_data) - 1),
            .f = function(idx) {
              temp_result =
                loess(formula = y ~ x,
                      data = temp_data[-idx, ],
                      span = span)
              prediction =
                try(predict(object = temp_result,
                            newdata = temp_data[idx, -2, drop = FALSE]))
              
              if (class(prediction) == "try-error") {
                data.frame(real = temp_data$y[idx],
                           prediction = NA)
              } else{
                data.frame(real = temp_data$y[idx],
                           prediction = as.numeric(prediction))
              }
            }
          ) %>%
          dplyr::bind_rows()
        
        if (all(is.na(prediction$prediction))) {
          temp_rmse = NA
        } else{
          temp_rmse = sqrt(sum((
            prediction$real - prediction$prediction
          ) ^ 2) / nrow(prediction))
        }
        
        data.frame(span = span, rmse = temp_rmse)
      }) %>%
      dplyr::bind_rows()
    
    # span_rmse
    
    plot = 
      data.frame(x, y) %>% 
      ggplot(aes(x, y)) +
      geom_point(size = 5) +
      # geom_line() +
      base_theme
    
    span_rmse = 
      span_rmse %>% 
      dplyr::filter(!is.na(rmse))
    idx = which.min(span_rmse$rmse)
    # for(i in 1:nrow(span_rmse)){
    plot =
      plot +
      geom_smooth(
        se = FALSE,
        span = span_rmse$span[idx],
        color = ggsci::pal_lancet()(n = 9)[idx]
      )
    # }
    
    plot = 
      plot + 
      ggplot2::ggtitle(label = paste("Span: ", span_rmse$span[idx])) +
      theme(title = element_text(colour = ggsci::pal_lancet()(n = 9)[idx]))
    
    list(span_rmse, plot)
  }



match_time = function(time1, time2, tol = 1) {
  time1
  time2
  purrr::map(time1, function(x) {
    idx = which(abs(difftime(
      time1 = x,
      time2 = time2,
      units = "hour"
    )) < tol)
    if (length(idx) == 0) {
      return(NULL)
    } else{
      data.frame(time1 = x,
                 time2 = time2[idx])
      
    }
  }) %>% 
    dplyr::bind_rows()
} 


