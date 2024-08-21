cellProximityBarplot_KC = function(CPscore, average=F){
  if (average){
    CPscore$original = CPscore %>% select(matches('original')) %>%
      rowMeans(na.rm=T)
  }
  first = CPscore %>% filter(grepl('KCs--',unified_int)) %>%
    separate(unified_int, c('receiver','sender'), sep='--')
  second = CPscore %>% filter(grepl('--KCs',unified_int)) %>%
    separate(unified_int, c('sender','receiver'), sep='--')
  table_mean_results_dc = rbind(first, second[second$sender!='KCs',])
  table_mean_results_dc$sender[table_mean_results_dc$sender=='Central Vein Endothelial cells'] = 'CVECs'
  table_mean_results_dc$sender[table_mean_results_dc$sender=='Portain Vein Endothelial cells'] = 'PVECs'
  table_mean_results_dc$sender[table_mean_results_dc$sender=='Lymphatic Endothelial cells'] = 'LECs'
  table_mean_results_dc = table_mean_results_dc %>% mutate(type_int = case_when(sender=='KCs'~'homo',sender!='KCs'~'hetero')) %>% dplyr::arrange(enrichm) %>% mutate(sender=factor(sender,levels=sender))
  ## filter to remove low number of cell-cell proximity interactions ##
  #table_mean_results_dc$sender = table_mean_results_dc$sender %>% as.character()
  #
  
  # data.table variables
  original <- simulations <- p_higher_orig <- p_lower_orig <- enrichm <- 
    type_int <- unified_int <- NULL
  
  # table_mean_results_dc_filter <- table_mean_results_dc[
  #   original >= min_orig_ints & simulations >= min_sim_ints, ]
  # table_mean_results_dc_filter <- table_mean_results_dc_filter[
  #   p_higher_orig <= p_val | p_lower_orig <= p_val, ]
  table_mean_results_dc_filter = table_mean_results_dc
  
  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::geom_bar(
    data = table_mean_results_dc_filter, 
    ggplot2::aes(x = sender, y = enrichm, fill = type_int), 
    stat = "identity", show.legend = FALSE)
  pl <- pl + ggplot2::coord_flip()
  pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::labs(y = "enrichment/depletion")
  pl
  
  bpl <- ggplot2::ggplot()
  bpl <- bpl + ggplot2::geom_bar(
    data = table_mean_results_dc_filter, 
    ggplot2::aes(x = sender, y = original, fill = type_int), 
    stat = "identity", show.legend = TRUE)
  bpl <- bpl + ggplot2::coord_flip()
  bpl <- bpl + ggplot2::theme_bw() + ggplot2::theme(
    axis.text.y = element_blank())
  bpl <- bpl + ggplot2::labs(y = "# of interactions")
  bpl
  
  combo_plot <- cowplot::plot_grid(
    pl, bpl, ncol = 2, rel_heights = c(1), 
    rel_widths = c(3, 1.5), align = "h")
  
  # output plot
  # return(GiottoVisuals::plot_output_handler(
  #   gobject = gobject,
  #   plot_object = combo_plot,
  #   save_plot = save_plot,
  #   return_plot = return_plot,
  #   show_plot = show_plot,
  #   default_save_name = default_save_name,
  #   save_param = save_param,
  #   else_return = NULL
  # ))
  return(combo_plot)
}

squidpy_cellProximityBarplot_KC = function(CPscore){
  CPscore$X = NULL
  CPscore = CPscore %>% mutate(original=count, enrichm=enrichment.depletion, sender=Sender)
  CPscore$sender[CPscore$sender=='Central Vein Endothelial cells'] = 'CVECs'
  CPscore$sender[CPscore$sender=='Portain Vein Endothelial cells'] = 'PVECs'
  CPscore$sender[CPscore$sender=='Lymphatic Endothelial cells'] = 'LECs'
  table_mean_results_dc = CPscore %>% mutate(type_int = case_when(sender=='KCs'~'homo',sender!='KCs'~'hetero')) %>% dplyr::arrange(enrichm) %>% mutate(sender=factor(sender,levels=sender))
  ## filter to remove low number of cell-cell proximity interactions ##
  
  # data.table variables
  original <- simulations <- p_higher_orig <- p_lower_orig <- enrichm <- 
    type_int <- unified_int <- NULL
  
  # table_mean_results_dc_filter <- table_mean_results_dc[
  #   original >= min_orig_ints & simulations >= min_sim_ints, ]
  # table_mean_results_dc_filter <- table_mean_results_dc_filter[
  #   p_higher_orig <= p_val | p_lower_orig <= p_val, ]
  table_mean_results_dc_filter = table_mean_results_dc
  
  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::geom_bar(
    data = table_mean_results_dc_filter, 
    ggplot2::aes(x = sender, y = enrichm, fill = type_int), 
    stat = "identity", show.legend = FALSE)
  pl <- pl + ggplot2::coord_flip()
  pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::labs(y = "enrichment/depletion")+
    guides(fill = "none")
  pl
  
  bpl <- ggplot2::ggplot()
  bpl <- bpl + ggplot2::geom_bar(
    data = table_mean_results_dc_filter, 
    ggplot2::aes(x = sender, y = original, fill = type_int), 
    stat = "identity", show.legend = TRUE)
  bpl <- bpl + ggplot2::coord_flip()
  bpl <- bpl + ggplot2::theme_bw() + ggplot2::theme(
    axis.text.y = element_blank())
  bpl <- bpl + ggplot2::labs(y = "# of interactions")+
    guides(fill = "none")
  bpl
  
  combo_plot <- cowplot::plot_grid(
    pl, bpl, ncol = 2, rel_heights = c(1), 
    rel_widths = c(3, 1.5), align = "h")
  
  # output plot
  # return(GiottoVisuals::plot_output_handler(
  #   gobject = gobject,
  #   plot_object = combo_plot,
  #   save_plot = save_plot,
  #   return_plot = return_plot,
  #   show_plot = show_plot,
  #   default_save_name = default_save_name,
  #   save_param = save_param,
  #   else_return = NULL
  # ))
  return(combo_plot)
}

interactionPie = function(CPscore){
  CPscore$original = CPscore %>% select(matches('original')) %>%
    rowMeans(na.rm=T)
  first = CPscore %>% filter(grepl('KCs--',unified_int)) %>%
    separate(unified_int, c('receiver','sender'), sep='--')
  second = CPscore %>% filter(grepl('--KCs',unified_int)) %>%
    separate(unified_int, c('sender','receiver'), sep='--')
  spatial_matrix = rbind(first, second[second$sender!='KCs',])
  spatial_matrix <- spatial_matrix %>% mutate(rank=rank(-original), sender = ifelse(rank<=5, sender, 'Others'))
  spatial_matrix$fraction = spatial_matrix$original / sum(spatial_matrix$original)
  spatial_matrix$ymax = cumsum(spatial_matrix$fraction)
  spatial_matrix$ymin = c(0, head(spatial_matrix$ymax, n=-1))
  spatial_matrix$labelPosition <- (spatial_matrix$ymax + spatial_matrix$ymin) / 2
  
  # Compute a good label
  #spatial_matrix$label <- paste0(spatial_matrix$sender, spatial_matrix$original)
  spatial_matrix$label <- spatial_matrix$original
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  # Make the plot
  fig = ggplot(spatial_matrix, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=sender)) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
    scale_color_manual(values=c25)+
    coord_polar(theta="y") +
    xlim(c(2, 4))
    #theme_void() +
    #theme(legend.position = "none")
  return(fig)
}
