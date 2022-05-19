
lagged_correlation = function(x,
                              y,
                              time1,
                              time2,
                              time_tol = 1,###unit is hour
                              step = 0.5/60,##unit is hour,
                              min_matched_sample = 10,
                              progressbar = TRUE,
                              all_idx = NULL
) {
  ##time_tol unit is hour
  ##step unit is hour 
  x = as.numeric(x) %>%
    scale() %>%
    as.numeric()
  
  y = as.numeric(y) %>%
    scale() %>%
    as.numeric()
  
  time_window1 = 
    seq(from = step/2, to = time_tol, by = step)
  
  time_window2 = 
    -rev(seq(from = step/2, to = time_tol, by = step))
  
  time_window = sort(c(time_window2, time_window1))
  
  temp_fun = 
    function(temp_idx,
             time_window,
             x,
             y,
             time1,
             time2
    ){
      idx =
        time1 %>%
        purrr::map(function(x) {
          diff_time =
            difftime(x, time2, units = "hours")
          which(diff_time > time_window[temp_idx] &
                  diff_time <= time_window[temp_idx + 1])
        })
    }
  
  bpparam = BiocParallel::MulticoreParam(workers = 10, 
                                         progressbar = progressbar)
  
  
  if(is.null(all_idx)){
    all_idx = 
      BiocParallel::bplapply(X = 1:(length(time_window) - 1), 
                             FUN = temp_fun, time_window = time_window,
                             x = x,
                             y = y,
                             time1 = time1,
                             time2 = time2, 
                             BPPARAM = bpparam)      
  }
  # lapply(all_idx, function(x){
  #   length(unlist(x))
  # }) %>%
  #   unlist() %>%
  #   plot()
  
  all_cor_result = 
    purrr::map(all_idx, function(idx){
      temp_y = 
        lapply(idx, function(x){
          mean(y[x])
        }) %>% 
        unlist()
      temp_x = x[which(!is.na(temp_y))]
      temp_y = temp_y[which(!is.na(temp_y))]
      if(length(temp_x) < min_matched_sample){
        return(NA)
      }else{
        tryCatch(expr = cor.test(temp_x, temp_y, method = "spearman"), 
                 error = function(na){return(NA)})    
      }
    }) 
  
  
  all_cor_p = 
    all_cor_result %>% 
    purrr::map(function(x){
      if(is.na(x)){
        return(NA)
      }else{
        x$p.value  
      }
    }) %>% 
    unlist()
  
  all_cor = 
    all_cor_result %>% 
    purrr::map(function(x){
      if(is.na(x)){
        return(NA)
      }else{
        x$estimate  
      }
    }) %>% 
    unlist() %>% 
    unname()
  
  which_max_idx = 
    which.max(abs(all_cor))
  
  max_idx = all_idx[[which_max_idx]]
  
  shift_time = 
    paste("(",
          paste(round(time_window[-length(time_window)] * 60, 2),
                round(time_window[-1] * 60, 2), sep = ','),
          "]", sep = ""
    )
  
  which_global_idx = 
    purrr::map(shift_time, function(x){
      x = 
        stringr::str_replace(x, "\\(", "") %>% 
        stringr::str_replace("\\]", "") %>% 
        stringr::str_split(",") %>% 
        `[[`(1) %>% 
        as.numeric()
      x[1] < 0 & x[2] > 0
    }) %>% 
    unlist() %>% 
    which()
  
  global_idx = all_idx[[which_global_idx]]
  
  global_cor = all_cor[which_global_idx]
  global_cor_p = all_cor_p[which_global_idx]
  
  return_result =
    list(
      x = x,
      time1 = time1,
      y = y,
      time2 = time2,
      idx = all_idx,
      all_cor = all_cor,
      all_cor_p = all_cor_p,
      all_cor_p = all_cor_p,
      shift_time = shift_time,
      which_max_idx = which_max_idx,
      which_global_idx = which_global_idx,
      max_idx = max_idx,
      max_cor = all_cor[which_max_idx],
      global_idx = global_idx,
      global_cor = global_cor
    )
}




lagged_alignment_plot = function(object,
                                 day_night_df,
                                 internal_omics_color = "#631879FF",
                                 wearable_color = "#E377C2FF",
                                 internal_omics_name = "cytokine_1",
                                 warable_name = "HR",
                                 which = c("global", "max"),
                                 x_limit = c(1,20),
                                 non_matched_point_size = 0.1,
                                 wearable_point_size = 1,
                                 internal_omics_point_size = 3,
                                 integrated = FALSE,
                                 add_connect_line = TRUE,
                                 add_point = TRUE) {
  if(!integrated){
    which = match.arg(which)
    
    if(which == "global"){
      idx = object$global_idx  
      shift = object$shift_time[object$which_global_idx]
      correlation = round(object$global_cor, 3)
    }else{
      idx = object$max_idx
      shift = object$shift_time[object$which_max_idx]
      correlation = round(object$max_cor, 3)
    }
    
    time1 = object$time1
    time2 = object$time2
    x = object$x
    y = object$y
    
    if(x_limit[2] > length(x)){
      x_limit[2] = length(x)
    }
    
    if(max(x) > min(y)){
      x = x - (max(x) - min(y))
    }
    
    ###get the segment data
    segment_data = 
      purrr::map(1:length(idx), function(i) {
        if (length(idx[[i]]) > 0) {
          data.frame(
            time1 = time1[i],
            x = x[i],
            time2 = time2[idx[[i]]],
            y = y[idx[[i]]]
          )
        }
      }) %>% 
      do.call(rbind, .) %>% 
      as.data.frame()
    
    x2 = data.frame(time = time1, value = x, class = "internal_omics")
    y2 = data.frame(time = time2, value = y, class = "wearable")
    
    x2$matched = 'NO'
    x2$matched[which(unlist(lapply(idx, length)) > 0)] =  "YES"
    
    y2$matched = 'NO'
    y2$matched[unique(unlist(idx))] = "YES"
    
    value = rbind(x2, y2)
    
    value = 
      value %>% 
      dplyr::mutate(matched = 
                      case_when(
                        matched == "YES" ~ class,
                        matched == "NO" ~ "NO",
                      ))
    
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
      geom_line(data = value,
                aes(x = time, 
                    y = value, 
                    group = class),
                color = "grey", 
                show.legend = FALSE) +
      # geom_point(data = value,
      #            aes(
      #              x = time,
      #              y = value,
      #              group = class,
      #              color = matched,
      #              size = matched
      #            ),
      #            shape = 16,
      #            show.legend = TRUE) +
      scale_size_manual(values = c("internal_omics" = internal_omics_point_size,
                                   "wearable" = wearable_point_size,
                                   "NO" = non_matched_point_size)) +
      scale_x_datetime(
        breaks = scales::date_breaks("4 hour"),
        date_labels = "%a %H:%M",
        timezone = "America/Los_Angeles", 
        limits = c(min(time1[x_limit[1]],
                       time2[x_limit[1]]
        ), 
        max(time1[x_limit[2]],
            time2[x_limit[2]]
        ))
        # expand = expansion(mult = c(0.1, 0.1))
      ) +
      labs(x = "", y = "", 
           title = paste("Shift window: ", 
                         shift, 
                         "; Correlation: ",
                         correlation,
                         sep = "")) +
      guides(color = guide_legend(title = "", 
                                  override.aes = list(size = 3)), 
             size = "none") +
      scale_color_manual(values = c("internal_omics" = unname(internal_omics_color),
                                    "wearable" = unname(wearable_color),
                                    'NO' = "grey"), 
                         labels = c("internal_omics" = internal_omics_name,
                                    "wearable" = warable_name,
                                    "NO" = "Non matched")
      ) +
      base_theme +
      theme(
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
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
        )
      )
    
    if(add_point){
      plot = 
        plot +
        geom_point(data = value,
                   aes(
                     x = time,
                     y = value,
                     group = class,
                     color = matched,
                     size = matched
                   ),
                   shape = 16,
                   show.legend = TRUE)
    }
    
    
    
    if(add_connect_line){
      plot = 
        plot +
        geom_segment(
          data = segment_data,
          aes(
            x = time1,
            y = x,
            xend = time2,
            yend = y
          ),
          color = internal_omics_color,
          show.legend = FALSE
        )
    }
    
    plot  
  }else{
    which = match.arg(which)
    
    if(which == "global"){
      idx = object$global_idx  
      shift = object$shift_time[object$which_global_idx]
      correlation = round(object$global_cor, 3)
    }else{
      idx = object$max_idx
      shift = object$shift_time[object$which_max_idx]
      correlation = round(object$max_cor, 3)
    }
    
    time1 = object$time1
    time2 = object$time2
    x = object$x
    y = object$y
    
    if(x_limit[2] > length(x)){
      x_limit[2] = length(x)
    }
    
    time2 = 
      purrr::map(idx, function(x){
        mean(time2[c(head(x,1), tail(x, 1))])
      }) %>% 
      unlist() %>% 
      lubridate::as_datetime(tz = "America/Los_Angeles")
    
    y = 
      purrr::map(idx, function(x){
        mean(y[x])
      }) %>% 
      unlist() 
    
    if(max(x, na.rm = TRUE) > min(y, na.rm = TRUE)){
      x = x - (max(x, na.rm = TRUE) - min(y, na.rm = TRUE))
    }
    
    ###get the segment data
    segment_data = 
      data.frame(time1, x, time2, y)
    
    x2 = data.frame(time = time1, value = x, class = "internal_omics")
    y2 = data.frame(time = time2, value = y, class = "wearable")
    
    x2$matched = 'NO'
    x2$matched[which(unlist(lapply(idx, length)) > 0)] =  "YES"
    
    y2$matched = 'YES'
    
    y2 = 
      y2 %>% 
      dplyr::filter(!is.na(time))
    
    value = rbind(x2, y2)
    
    value = 
      value %>% 
      dplyr::mutate(matched = 
                      case_when(
                        matched == "YES" ~ class,
                        matched == "NO" ~ "NO",
                      ))
    
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
      # geom_segment(
      #   data = segment_data,
      #   aes(
      #     x = time1,
      #     y = x,
      #     xend = time2,
      #     yend = y
      #   ),
      #   color = internal_omics_color,
      #   show.legend = FALSE
      # ) +
    geom_line(data = value,
              aes(x = time, 
                  y = value, 
                  group = class,
                  color = class),
              show.legend = FALSE) +
      geom_point(data = value,
                 aes(
                   x = time,
                   y = value,
                   group = class,
                   color = matched,
                   size = matched
                 ),
                 shape = 16,
                 show.legend = TRUE) +
      scale_size_manual(values = c("internal_omics" = internal_omics_point_size,
                                   "wearable" = wearable_point_size,
                                   "NO" = non_matched_point_size)) +
      scale_x_datetime(
        breaks = scales::date_breaks("4 hour"),
        date_labels = "%a %H:%M",
        timezone = "America/Los_Angeles", 
        limits = c(min(time1[x_limit[1]],
                       time2[x_limit[1]]
        ), 
        max(time1[x_limit[2]],
            time2[x_limit[2]]
        ))
        # expand = expansion(mult = c(0.1, 0.1))
      ) +
      labs(x = "", y = "", 
           title = paste("Shift window: ", 
                         shift, 
                         "; Correlation: ",
                         correlation,
                         sep = "")) +
      guides(color = guide_legend(title = "", 
                                  override.aes = list(size = 3)), 
             size = "none") +
      scale_color_manual(values = c("internal_omics" = unname(internal_omics_color),
                                    "wearable" = unname(wearable_color),
                                    'NO' = "grey"), 
                         labels = c("internal_omics" = internal_omics_name,
                                    "wearable" = warable_name,
                                    "NO" = "Non matched")
      ) +
      base_theme +
      theme(
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
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
        )
      )
    
    if(add_connect_line){
      plot = 
        plot +
        geom_segment(
          data = segment_data,
          aes(
            x = time1,
            y = x,
            xend = time2,
            yend = y
          ),
          color = internal_omics_color,
          show.legend = FALSE
        )
    }
    plot  
  }
  
}










lagged_sactter_plot = function(object,
                               internal_omics_name = "cytokine_1",
                               warable_name = "HR",
                               which = c("global", "max")) {
  if(!integrated){
    which = match.arg(which)
    
    if(which == "global"){
      idx = object$global_idx  
      shift = object$shift_time[object$which_global_idx]
      correlation = round(object$global_cor, 3)
    }else{
      idx = object$max_idx
      shift = object$shift_time[object$which_max_idx]
      correlation = round(object$max_cor, 3)
    }
    
    time1 = object$time1
    time2 = object$time2
    x = object$x
    y = object$y
    
    idx1 =
      lapply(idx, function(x) {
        length(x)
      }) %>%
      unlist() %>%
      `>`(0) %>%
      which()
    
    x2 = x[idx1]
    
    y2 = 
      lapply(idx, function(z){
        mean(y[z])
      }) %>% 
      unlist()
    
    y2 = y2[!is.na(y2)]
    
    x2 = data.frame(time = time1, value = x, class = "internal_omics")
    y2 = data.frame(time = time2, value = y, class = "wearable")
    
    x2$matched = 'NO'
    x2$matched[which(unlist(lapply(idx, length)) > 0)] =  "YES"
    
    y2$matched = 'NO'
    y2$matched[unique(unlist(idx))] = "YES"
    
    value = rbind(x2, y2)
    
    value = 
      value %>% 
      dplyr::mutate(matched = 
                      case_when(
                        matched == "YES" ~ class,
                        matched == "NO" ~ "NO",
                      ))
    
    plot =
      ggplot() +
      geom_line(data = value,
                aes(x = time, 
                    y = value, 
                    group = class),
                color = "grey", 
                show.legend = FALSE) +
      # geom_point(data = value,
      #            aes(
      #              x = time,
      #              y = value,
      #              group = class,
      #              color = matched,
      #              size = matched
      #            ),
      #            shape = 16,
      #            show.legend = TRUE) +
      scale_size_manual(values = c("internal_omics" = internal_omics_point_size,
                                   "wearable" = wearable_point_size,
                                   "NO" = non_matched_point_size)) +
      scale_x_datetime(
        breaks = scales::date_breaks("4 hour"),
        date_labels = "%a %H:%M",
        timezone = "America/Los_Angeles", 
        limits = c(min(time1[x_limit[1]],
                       time2[x_limit[1]]
        ), 
        max(time1[x_limit[2]],
            time2[x_limit[2]]
        ))
        # expand = expansion(mult = c(0.1, 0.1))
      ) +
      labs(x = "", y = "", 
           title = paste("Shift window: ", 
                         shift, 
                         "; Correlation: ",
                         correlation,
                         sep = "")) +
      guides(color = guide_legend(title = "", 
                                  override.aes = list(size = 3)), 
             size = "none") +
      scale_color_manual(values = c("internal_omics" = unname(internal_omics_color),
                                    "wearable" = unname(wearable_color),
                                    'NO' = "grey"), 
                         labels = c("internal_omics" = internal_omics_name,
                                    "wearable" = warable_name,
                                    "NO" = "Non matched")
      ) +
      base_theme +
      theme(
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
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
        )
      )
    
    if(add_point){
      plot = 
        plot +
        geom_point(data = value,
                   aes(
                     x = time,
                     y = value,
                     group = class,
                     color = matched,
                     size = matched
                   ),
                   shape = 16,
                   show.legend = TRUE)
    }
    
    
    
    if(add_connect_line){
      plot = 
        plot +
        geom_segment(
          data = segment_data,
          aes(
            x = time1,
            y = x,
            xend = time2,
            yend = y
          ),
          color = internal_omics_color,
          show.legend = FALSE
        )
    }
    
    plot  
  }else{
    which = match.arg(which)
    
    if(which == "global"){
      idx = object$global_idx  
      shift = object$shift_time[object$which_global_idx]
      correlation = round(object$global_cor, 3)
    }else{
      idx = object$max_idx
      shift = object$shift_time[object$which_max_idx]
      correlation = round(object$max_cor, 3)
    }
    
    time1 = object$time1
    time2 = object$time2
    x = object$x
    y = object$y
    
    if(x_limit[2] > length(x)){
      x_limit[2] = length(x)
    }
    
    time2 = 
      purrr::map(idx, function(x){
        mean(time2[c(head(x,1), tail(x, 1))])
      }) %>% 
      unlist() %>% 
      lubridate::as_datetime(tz = "America/Los_Angeles")
    
    y = 
      purrr::map(idx, function(x){
        mean(y[x])
      }) %>% 
      unlist() 
    
    if(max(x, na.rm = TRUE) > min(y, na.rm = TRUE)){
      x = x - (max(x, na.rm = TRUE) - min(y, na.rm = TRUE))
    }
    
    ###get the segment data
    segment_data = 
      data.frame(time1, x, time2, y)
    
    x2 = data.frame(time = time1, value = x, class = "internal_omics")
    y2 = data.frame(time = time2, value = y, class = "wearable")
    
    x2$matched = 'NO'
    x2$matched[which(unlist(lapply(idx, length)) > 0)] =  "YES"
    
    y2$matched = 'YES'
    
    y2 = 
      y2 %>% 
      dplyr::filter(!is.na(time))
    
    value = rbind(x2, y2)
    
    value = 
      value %>% 
      dplyr::mutate(matched = 
                      case_when(
                        matched == "YES" ~ class,
                        matched == "NO" ~ "NO",
                      ))
    
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
      # geom_segment(
      #   data = segment_data,
      #   aes(
      #     x = time1,
      #     y = x,
      #     xend = time2,
      #     yend = y
      #   ),
      #   color = internal_omics_color,
      #   show.legend = FALSE
      # ) +
    geom_line(data = value,
              aes(x = time, 
                  y = value, 
                  group = class,
                  color = class),
              show.legend = FALSE) +
      geom_point(data = value,
                 aes(
                   x = time,
                   y = value,
                   group = class,
                   color = matched,
                   size = matched
                 ),
                 shape = 16,
                 show.legend = TRUE) +
      scale_size_manual(values = c("internal_omics" = internal_omics_point_size,
                                   "wearable" = wearable_point_size,
                                   "NO" = non_matched_point_size)) +
      scale_x_datetime(
        breaks = scales::date_breaks("4 hour"),
        date_labels = "%a %H:%M",
        timezone = "America/Los_Angeles", 
        limits = c(min(time1[x_limit[1]],
                       time2[x_limit[1]]
        ), 
        max(time1[x_limit[2]],
            time2[x_limit[2]]
        ))
        # expand = expansion(mult = c(0.1, 0.1))
      ) +
      labs(x = "", y = "", 
           title = paste("Shift window: ", 
                         shift, 
                         "; Correlation: ",
                         correlation,
                         sep = "")) +
      guides(color = guide_legend(title = "", 
                                  override.aes = list(size = 3)), 
             size = "none") +
      scale_color_manual(values = c("internal_omics" = unname(internal_omics_color),
                                    "wearable" = unname(wearable_color),
                                    'NO' = "grey"), 
                         labels = c("internal_omics" = internal_omics_name,
                                    "wearable" = warable_name,
                                    "NO" = "Non matched")
      ) +
      base_theme +
      theme(
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
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
        )
      )
    
    if(add_connect_line){
      plot = 
        plot +
        geom_segment(
          data = segment_data,
          aes(
            x = time1,
            y = x,
            xend = time2,
            yend = y
          ),
          color = internal_omics_color,
          show.legend = FALSE
        )
    }
    plot  
  }
  
}




show_sample_matching = function(object, 
                                index = 1,
                                only_remain_matched = TRUE,
                                day_night_df = NULL,
                                add_text = FALSE){
  match_idx = object$idx[[index]]
  time1 = object$time1
  time2 = object$time2
  
  idx1 = 
    lapply(match_idx, length) %>% 
    unlist() %>% 
    `!=`(0) %>% 
    which()
  
  if (length(idx1) == 0) {
    return(NULL)
  }
  
  idx2 = match_idx[idx1]
  
  segment_data = 
    purrr::map2(.x = idx1, .y = idx2, function(x,y){
      data.frame(x = time1[x], 
                 xend = time2[y],
                 y = "time1",
                 yend = "time2")
    }) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
  
  if(only_remain_matched){
    temp_data =
      data.frame(x = c(time1[unique(idx1)], time2[unique(unlist(idx2))]), 
                 y = c(rep("time1", length(time1[unique(idx1)])),
                       rep("time2", length(time2[unique(unlist(idx2))]))))
  }else{
    temp_data =
      data.frame(x = c(time1, time2), 
                 y = c(rep("time1", length(time1)),rep("time2", length(time2))))    
  }
  
  if(!is.null(day_night_df)){
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
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend), 
                   data = segment_data) +
      geom_point(aes(x = x, y = y,
                     color = y), 
                 data = temp_data,
                 show.legend = FALSE) +
      scale_x_datetime(
        breaks = scales::date_breaks("12 hour"),
        date_labels = "%a %H:%M",
        timezone = "America/Los_Angeles"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid = element_blank()) +
      labs(x = "", y = "") +
      ggsci::scale_color_jama()
  }else{
    plot = 
      ggplot() +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend), 
                   data = segment_data) +
      geom_point(aes(x = x, y = y,
                     color = y), 
                 data = temp_data,
                 show.legend = FALSE) +
      scale_x_datetime(
        breaks = scales::date_breaks("12 hour"),
        date_labels = "%a %H:%M",
        timezone = "America/Los_Angeles"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid = element_blank()) +
      labs(x = "", y = "") +
      ggsci::scale_color_jama()
  }
  
  shift_time =
    object$shift_time[index]
  
  plot =   
    plot +
    ggtitle(label = paste("Shift time(min):",shift_time, sep = ""), 
            subtitle = paste(length(idx1), "matched samples"))
  
  if(add_text){
    plot = 
      plot +
      ggrepel::geom_text_repel(
        mapping = aes(x = x, y = y, 
                      label = paste(lubridate::hour(x), lubridate::minute(x), sep = ":")),
        data = temp_data,
        size = 2
      )
  }
  
  plot
  
}






fitpeaks <- function(y, pos) {
  names(y) <- NULL
  tabnames <- c("rt", "sd", "FWHM", "height", "area")
  noPeaksMat <- matrix(rep(NA, 5),
                       nrow = 1,
                       dimnames = list(NULL, tabnames)) %>% 
    as.data.frame()
  
  # on.edge <- sapply(pos,
  #                   function(x)
  #                     y[x + 1] == 0 | y[x - 1] == 0)
  # # browser()
  # if(lapply(on.edge, function(x){
  #   length(x) == 0
  # }) %>% 
  # unlist() %>% 
  # all){
  #   return(noPeaksMat)
  # }
  # 
  #   pos <- pos[!on.edge]
  
  
  if (length(pos) == 0){
    return(noPeaksMat)
  }
  
  fitpk <- function(xloc) {
    ## find all areas higher than half the current max
    peak.loc <- which(y > 0.2 * y[xloc])
    peak.loc.diff <- diff(peak.loc)
    boundaries <-
      c(0, 
        which(diff(peak.loc) != 1), 
        length(peak.loc) + 1)
    
    peaknrs <- rep(1:length(boundaries),
                   c(boundaries[1], diff(c(boundaries))))
    peaknrs[boundaries[-1]] <- NA
    current.peak <- peaknrs[peak.loc == xloc]
    current.peak <- current.peak[!is.na(current.peak)]
    if (length(current.peak) == 0)
      return(rep(NA, 5))
    
    ## only retain those points adjacent to the current max
    FWHM <- diff(range(peak.loc[peaknrs == current.peak],
                       na.rm = TRUE))
    pksd <- FWHM / (2 * sqrt(2 * log(2)))
    
    c(
      rt = xloc,
      sd = pksd,
      FWHM = FWHM,
      height = y[xloc],
      area = y[xloc] / dnorm(x = xloc, mean = xloc, sd = pksd)
    )
  }
  
  huhn <- t(sapply(pos, fitpk))
  colnames(huhn) <- tabnames
  
  huhn
}


evaluate_peak_quality = function(object,
                                 plot = TRUE){
  if(is.null(object)){
    return(list(score = 0, plot = NULL))
  }
  
  cor = 
    object$all_cor
  
  shift_time = 
    object$shift_time %>% 
    stringr::str_replace("\\(|\\]", "") %>% 
    stringr::str_replace("\\]", "") %>% 
    stringr::str_split("\\,") %>% 
    purrr::map(function(x){
      mean(as.numeric(x))
    }) %>% 
    unlist()
  
  
  ###loess fit
  cor[is.na(cor)] = 0
  span = 3/length(cor)
  if(span < 0.1){
    span = 0.1
  }
  loess_lm = loess(cor ~ shift_time, data = data.frame(cor, shift_time), 
                   span = span)
  new_shift_time = seq(min(shift_time), max(shift_time), 1)
  new_cor = 
    predict(object = loess_lm, 
            newdata = data.frame(shift_time = new_shift_time))
  
  # plot(shift_time, cor, pch = 16, cex = 3)
  # lines(new_shift_time, new_cor, col = "red", type = "b", pch = 2)
  # 
  ####positive fit
  if(any(new_cor < 0)){
    cor_positive = new_cor - min(new_cor)  
  }else{
    cor_positive = new_cor
  }
  
  pk.pos = which.max(cor_positive)
  pks_positive = fitpeaks(y = cor_positive, pos = pk.pos)
  
  if(!is.na(pks_positive[1,1])){
    fitted_y_positive = 
      dnorm(x = 1:length(cor_positive), 
            mean = as.numeric(pks_positive[1,"rt"]), 
            sd = as.numeric(pks_positive[1,"sd"])) * as.numeric(pks_positive[1,"area"])
    fitted_y_positive[is.na(fitted_y_positive)] = 0
    score_positive = cor(cor_positive,
                         fitted_y_positive, method = "spearman")  
  }else{
    score_positive = 0
  }
  
  if(is.na(score_positive)){
    score_positive = 0
  }
  
  ####negative fit
  if(any(new_cor > 0)){
    cor_negative = new_cor - max(cor)  
  }else{
    cor_negative = new_cor
  }
  
  pk.pos = which.max(-cor_negative)
  pks_negative = fitpeaks(y = -cor_negative, pos = pk.pos)
  
  if (!is.na(pks_negative[1, 1])) {
    fitted_y_negative =
      dnorm(
        x = 1:length(cor_negative),
        mean = as.numeric(pks_negative[1, "rt"]),
        sd = as.numeric(pks_negative[1, "sd"])
      ) * as.numeric(pks_negative[1, "area"])
    fitted_y_negative[is.na(fitted_y_negative)] = 0
    score_negative = cor(-cor_negative,
                         fitted_y_negative, method = "spearman")
  } else{
    score_negative = 0
  }
  
  if(is.na(score_negative)){
    score_negative = 0
  }
  
  positive_max_cor = max(object$all_cor)
  negative_max_cor = min(object$all_cor)
  
  if(positive_max_cor <= 0){
    score_positive = 0
  }
  
  if(negative_max_cor >= 0){
    score_negative = 0
  }
  
  score = max(c(score_positive, score_negative))
  
  if (score_positive >= score_negative) {
    true = new_cor
    fitted = fitted_y_positive
    if (any(new_cor < 0)) {
      fitted = fitted_y_positive + min(new_cor)
    } else{
      fitted = fitted_y_positive
    }
  }else{
    true = new_cor
    fitted = -fitted_y_negative
    if (any(new_cor > 0)) {
      fitted = -fitted_y_negative + max(new_cor) 
    } else{
      fitted = -fitted_y_negative
    }
  } 
  
  score
  
  if(plot){
    temp_data1 =   
      data.frame(new_shift_time, true, fitted)
    
    temp_data2 = 
      temp_data1 %>% 
      dplyr::filter(new_shift_time %in% shift_time)
    
    temp_data2$all_cor_p = object$all_cor_p
    
    temp_data2 = 
      temp_data2 %>% 
      dplyr::mutate(yesornot = 
                      case_when(
                        all_cor_p < 0.05 ~ "yes",
                        all_cor_p >= 0.05 ~ "no"
                      )) %>% 
      dplyr::mutate(log_p = -log(all_cor_p, 10))
    
    p = 
      ggplot() +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      geom_point(aes(x = new_shift_time, 
                     y = true, 
                     color = yesornot,
                     size = log_p), 
                 show.legend = FALSE, data = temp_data2) +
      scale_color_manual(values = c("yes" = "red", "no" = "grey")) +
      scale_size_continuous(range = c(2,5)) +
      geom_line(aes(new_shift_time, fitted, group = 1), data = temp_data1, color = "red") +
      theme_bw() +
      theme(panel.grid.minor = element_blank()) +
      labs(x = "Shift time (min)",
           y = "Pearsom correlation") 
  }else{
    p = NULL}
  
  
  list(score = score, 
       plot = p)
}









# cor(cor, fitted_y, method = "spearman")
