scaled_prioritization_table = function(
    sender_receiver_info, sender_receiver_de, ligand_activities, 
    lr_condition_de = NULL, spatial_matrix = NULL,
    prioritizing_weights = c("de_ligand" = 1, "de_receptor" = 1,
                             "activity_scaled" = 2, "exprs_ligand" = 1,
                             "exprs_receptor" = 1,
                             "ligand_condition_specificity" = 0, 
                             "receptor_condition_specificity" = 0, 
                             "spatial_prioritization" = 0)){
  
  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
  
  # Ligand DE prioritization
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                  p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
    dplyr::arrange(-lfc_pval_ligand)
  
  # Receptor DE prioritization
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                  p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) )
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  # Ligand activity prioritization
  ligand_activity_prioritization = ligand_activities %>% select(test_ligand, aupr_corrected, rank) %>% rename(activity=aupr_corrected, ligand=test_ligand) %>%
    dplyr::mutate(activity_zscore = nichenetr::scaling_zscore(activity),
                  scaled_activity = scale_quantile_adapted(activity, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_zscore)
  
  
  # Cell-type specificity of expression of ligand:  per ligand scale across cell types
  ligand_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(sender, ligand, avg_ligand) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>%
    dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
  
  # Cell-type specificity of expression of receptor:  per receptor scale across cell types
  receptor_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(receiver, receptor, avg_receptor) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>%
    dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
  
  if (!is.null(spatial_matrix)){
    # Proximity (Spatial)
    #spatial_matrix = spatial_matrix$enrichm_res %>% separate(unified_int, c('sender','receiver'), sep='--')
    #proximity_prioritization = spatial_matrix %>% select(sender, receiver, enrichm, int_ranking) %>% dplyr::distinct() %>%
    #  dplyr::mutate(scaled_enrich = scale_quantile_adapted(enrichm)) %>% dplyr::arrange(-scaled_enrich)
    first = spatial_matrix %>% 
      filter(grepl('KCs--',unified_int)) %>%
      separate(unified_int, c('receiver','sender'), sep='--') %>%
      select(sender, enrichm)
    second = spatial_matrix %>% 
      filter(grepl('--KCs',unified_int)) %>%
      separate(unified_int, c('sender','receiver'), sep='--') %>%
      select(sender, enrichm)
    spatial_matrix = rbind(first, second[second$sender!='KCs',])
    proximity_prioritization = spatial_matrix %>% select(sender, enrichm) %>% dplyr::distinct() %>%
      dplyr::mutate(scaled_enrich = scale_quantile_adapted(enrichm)) %>% dplyr::arrange(-scaled_enrich)
    #dplyr::mutate(scaled_enrich = scale(enrichm, center=TRUE, scale=TRUE))
  } else {
    if (prioritizing_weights['spatial_proximity'] > 0) {
      stop("No spatial_matrix table given, yet the relevant weights are nonzero.\nEither set weights of 'spatial_proximity' to zero or provide spatial_matrix.")
    }
  }
  
  if (!is.null(lr_condition_de)){
    # Condition specificity of ligand (upregulation)
    ligand_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                    p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
    ligand_condition_prioritization = ligand_condition_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_ligand) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -ligand)
    
    # Condition specificity of receptor (upregulation)
    receptor_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                    p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor))
    receptor_condition_prioritization = receptor_condition_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                            scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_receptor) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -receptor)
    
  } else {
    if (any(prioritizing_weights[grep("specificity", names(prioritizing_weights))] > 0)) {
      stop("No lr_condition_de table given, yet the relevant weights are nonzero.\nEither set weights of 'ligand_condition_specificity' and 'receptor_condition_specificity' to zero or provide lr_condition_de.")
    }
  }
  
  weights <- prioritizing_weights
  # final group-based prioritization
  group_prioritization_tbl = sender_receiver_de %>%
    dplyr::inner_join(sender_receiver_info) %>%
    {if (weights["de_ligand"] > 0) dplyr::inner_join(., sender_ligand_prioritization) else (.)} %>%
    {if (weights["activity_scaled"] > 0) dplyr::inner_join(., ligand_activity_prioritization) else (.)} %>%
    {if (weights["de_receptor"] > 0) dplyr::inner_join(., receiver_receptor_prioritization) else (.)} %>%
    {if (weights["exprs_ligand"] > 0) dplyr::inner_join(., ligand_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["exprs_receptor"] > 0) dplyr::inner_join(., receptor_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["ligand_condition_specificity"] > 0) dplyr::inner_join(., ligand_condition_prioritization) else (.)} %>%
    {if (weights["receptor_condition_specificity"] > 0) dplyr::inner_join(., receptor_condition_prioritization) else (.)} %>%
    {if (weights["spatial_proximity"] > 0) dplyr::inner_join(., proximity_prioritization) else (.)}
  
  
  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = 2*weights["de_ligand"] + 2*weights["de_receptor"] + weights["activity_scaled"] + weights["exprs_ligand"] + weights["exprs_receptor"] + weights["ligand_condition_specificity"] + weights["receptor_condition_specificity"] + weights['spatial_proximity']
  group_prioritization_tbl = group_prioritization_tbl %>% rowwise() %>%
    dplyr::mutate(prioritization_score =
                    (
                      (prioritizing_weights["de_ligand"] * ifelse("scaled_lfc_ligand" %in% names(group_prioritization_tbl), scaled_lfc_ligand, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_lfc_receptor" %in% names(group_prioritization_tbl), scaled_lfc_receptor, 0)) +
                        (prioritizing_weights["de_ligand"] * ifelse("scaled_p_val_ligand_adapted" %in% names(group_prioritization_tbl), scaled_p_val_ligand_adapted, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_p_val_receptor_adapted" %in% names(group_prioritization_tbl), scaled_p_val_receptor_adapted, 0)) +
                        (prioritizing_weights["activity_scaled"] * ifelse("scaled_activity" %in% names(group_prioritization_tbl), scaled_activity, 0)) +
                        (prioritizing_weights["exprs_ligand"] * ifelse("scaled_avg_exprs_ligand" %in% names(group_prioritization_tbl), scaled_avg_exprs_ligand, 0)) +
                        (prioritizing_weights["exprs_receptor"] * ifelse("scaled_avg_exprs_receptor" %in% names(group_prioritization_tbl), scaled_avg_exprs_receptor, 0)) +
                        (prioritizing_weights["ligand_condition_specificity"] * ifelse("scaled_lfc_ligand_group" %in% names(group_prioritization_tbl), scaled_lfc_ligand_group, 0)) +
                        (prioritizing_weights["receptor_condition_specificity"] * ifelse("scaled_lfc_receptor_group" %in% names(group_prioritization_tbl), scaled_lfc_receptor_group, 0)) +
                        (prioritizing_weights["spatial_proximity"] * ifelse("scaled_enrich" %in% names(group_prioritization_tbl), scaled_enrich, 0))
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score) %>%
    ungroup()
  return(group_prioritization_tbl)
}

unscaled_prioritization_table = function(
    sender_receiver_info, sender_receiver_de, ligand_activities, 
    lr_condition_de = NULL, spatial_matrix = NULL,
    prioritizing_weights = c("de_ligand" = 1, "de_receptor" = 1,
                             "activity_scaled" = 2, "exprs_ligand" = 1,
                             "exprs_receptor" = 1,
                             "ligand_condition_specificity" = 0, 
                             "receptor_condition_specificity" = 0, 
                             "spatial_prioritization" = 0)){
  
  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
  
  # Ligand DE prioritization
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                  p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
    dplyr::arrange(-lfc_pval_ligand)
  
  # Receptor DE prioritization
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                  p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) )
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  # Ligand activity prioritization
  ligand_activity_prioritization = ligand_activities %>% select(test_ligand, aupr_corrected, rank) %>% rename(activity=aupr_corrected, ligand=test_ligand) %>%
    dplyr::mutate(activity_zscore = nichenetr::scaling_zscore(activity),
                  scaled_activity = scale_quantile_adapted(activity, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_zscore)
  
  
  # Cell-type specificity of expression of ligand:  per ligand scale across cell types
  ligand_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(sender, ligand, avg_ligand) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>%
    dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
  
  # Cell-type specificity of expression of receptor:  per receptor scale across cell types
  receptor_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(receiver, receptor, avg_receptor) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>%
    dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
  
  if (!is.null(spatial_matrix)){
  # Proximity (Spatial)
    #spatial_matrix = spatial_matrix$enrichm_res %>% separate(unified_int, c('sender','receiver'), sep='--')
    #proximity_prioritization = spatial_matrix %>% select(sender, receiver, enrichm, int_ranking) %>% dplyr::distinct() %>%
    #  dplyr::mutate(scaled_enrich = scale_quantile_adapted(enrichm)) %>% dplyr::arrange(-scaled_enrich)
    first = spatial_matrix %>% 
      filter(grepl('KCs--',unified_int)) %>%
      separate(unified_int, c('receiver','sender'), sep='--') %>%
      select(sender, enrichm)
    second = spatial_matrix %>% 
      filter(grepl('--KCs',unified_int)) %>%
      separate(unified_int, c('sender','receiver'), sep='--') %>%
      select(sender, enrichm)
    spatial_matrix = rbind(first, second[second$sender!='KCs',])
    proximity_prioritization = spatial_matrix %>% select(sender, enrichm) %>% dplyr::distinct() %>%
      dplyr::mutate(scaled_enrich = enrichm) %>% dplyr::arrange(-scaled_enrich)
    #dplyr::mutate(scaled_enrich = scale(enrichm, center=TRUE, scale=TRUE))
  } else {
    if (prioritizing_weights['spatial_proximity'] > 0) {
      stop("No spatial_matrix table given, yet the relevant weights are nonzero.\nEither set weights of 'spatial_proximity' to zero or provide spatial_matrix.")
    }
  }
  
  if (!is.null(lr_condition_de)){
    # Condition specificity of ligand (upregulation)
    ligand_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                    p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
    ligand_condition_prioritization = ligand_condition_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_ligand) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -ligand)
    
    # Condition specificity of receptor (upregulation)
    receptor_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                    p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor))
    receptor_condition_prioritization = receptor_condition_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                            scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_receptor) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -receptor)
    
  } else {
    if (any(prioritizing_weights[grep("specificity", names(prioritizing_weights))] > 0)) {
      stop("No lr_condition_de table given, yet the relevant weights are nonzero.\nEither set weights of 'ligand_condition_specificity' and 'receptor_condition_specificity' to zero or provide lr_condition_de.")
    }
  }
  
  weights <- prioritizing_weights
  # final group-based prioritization
  group_prioritization_tbl = sender_receiver_de %>%
    dplyr::inner_join(sender_receiver_info) %>%
    {if (weights["de_ligand"] > 0) dplyr::inner_join(., sender_ligand_prioritization) else (.)} %>%
    {if (weights["activity_scaled"] > 0) dplyr::inner_join(., ligand_activity_prioritization) else (.)} %>%
    {if (weights["de_receptor"] > 0) dplyr::inner_join(., receiver_receptor_prioritization) else (.)} %>%
    {if (weights["exprs_ligand"] > 0) dplyr::inner_join(., ligand_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["exprs_receptor"] > 0) dplyr::inner_join(., receptor_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["ligand_condition_specificity"] > 0) dplyr::inner_join(., ligand_condition_prioritization) else (.)} %>%
    {if (weights["receptor_condition_specificity"] > 0) dplyr::inner_join(., receptor_condition_prioritization) else (.)} %>%
    {if (weights["spatial_proximity"] > 0) dplyr::inner_join(., proximity_prioritization) else (.)}
  
  
  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = 2*weights["de_ligand"] + 2*weights["de_receptor"] + weights["activity_scaled"] + weights["exprs_ligand"] + weights["exprs_receptor"] + weights["ligand_condition_specificity"] + weights["receptor_condition_specificity"] + weights['spatial_proximity']
  group_prioritization_tbl = group_prioritization_tbl %>% rowwise() %>%
    dplyr::mutate(prioritization_score =
                    (
                      (prioritizing_weights["de_ligand"] * ifelse("scaled_lfc_ligand" %in% names(group_prioritization_tbl), scaled_lfc_ligand, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_lfc_receptor" %in% names(group_prioritization_tbl), scaled_lfc_receptor, 0)) +
                        (prioritizing_weights["de_ligand"] * ifelse("scaled_p_val_ligand_adapted" %in% names(group_prioritization_tbl), scaled_p_val_ligand_adapted, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_p_val_receptor_adapted" %in% names(group_prioritization_tbl), scaled_p_val_receptor_adapted, 0)) +
                        (prioritizing_weights["activity_scaled"] * ifelse("scaled_activity" %in% names(group_prioritization_tbl), scaled_activity, 0)) +
                        (prioritizing_weights["exprs_ligand"] * ifelse("scaled_avg_exprs_ligand" %in% names(group_prioritization_tbl), scaled_avg_exprs_ligand, 0)) +
                        (prioritizing_weights["exprs_receptor"] * ifelse("scaled_avg_exprs_receptor" %in% names(group_prioritization_tbl), scaled_avg_exprs_receptor, 0)) +
                        (prioritizing_weights["ligand_condition_specificity"] * ifelse("scaled_lfc_ligand_group" %in% names(group_prioritization_tbl), scaled_lfc_ligand_group, 0)) +
                        (prioritizing_weights["receptor_condition_specificity"] * ifelse("scaled_lfc_receptor_group" %in% names(group_prioritization_tbl), scaled_lfc_receptor_group, 0)) +
                        (prioritizing_weights["spatial_proximity"] * ifelse("scaled_enrich" %in% names(group_prioritization_tbl), scaled_enrich, 0))
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score) %>%
    ungroup()
  return(group_prioritization_tbl)
}

lr_diff_prioritization_table = function(sender_receiver_info, sender_receiver_de, ligand_activities, lr_condition_de = NULL, spatial_matrix_connected = NULL, spatial_matrix_paracrine=NULL,
                                        prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2, "exprs_ligand" = 1,"exprs_receptor" = 1,
                                                                 "ligand_condition_specificity" = 0, "receptor_condition_specificity"=0)){
  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
  
  # Ligand DE prioritization
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                  p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
    dplyr::arrange(-lfc_pval_ligand)
  
  # Receptor DE prioritization
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                  p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) )
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  # Ligand activity prioritization
  ligand_activity_prioritization = ligand_activities %>% select(test_ligand, aupr_corrected, rank) %>% rename(activity=aupr_corrected, ligand=test_ligand) %>%
    dplyr::mutate(activity_zscore = nichenetr::scaling_zscore(activity),
                  scaled_activity = scale_quantile_adapted(activity, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_zscore)
  
  
  # Cell-type specificity of expression of ligand:  per ligand scale across cell types
  ligand_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(sender, ligand, avg_ligand) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>%
    dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
  
  # Cell-type specificity of expression of receptor:  per receptor scale across cell types
  receptor_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(receiver, receptor, avg_receptor) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>%
    dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
  
  if (!is.null(spatial_matrix_connected)){
    library(CellChat)
    cellchat_table = CellChatDB.mouse
    cellchat_table = subsetDB(cellchat_table, search=c("Cell-Cell Contact"), key='annotation')
    cellchat_table = cellchat_table$interaction %>% select(ligand.symbol, receptor.symbol)
    cellchat_table = cellchat_table %>% separate_rows(receptor.symbol, sep=', ')
    cellchat_table = cellchat_table %>% rename(ligand=ligand.symbol, receptor=receptor.symbol)
    
    cellphone_table = read.table('data/interaction_input.csv', header=T, sep=',')
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'HLA-A', 'HLA_A'))
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'HLA-B', 'HLA_B'))
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'HLA-C', 'HLA_C'))
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'HLA-E', 'HLA_E'))
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'HLA-F', 'HLA_F'))
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'HLA-G', 'HLA_G'))
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'ERVH48-1', 'ERVH48_1'))
    cellphone_table = cellphone_table %>% mutate(interactors = str_replace(interactors, 'MT-RNR2', 'MT_RNR2'))
    cellphone_table = cellphone_table %>% filter(grepl('adhesion', classification, ignore.case=T)|grepl('adhesion', directionality, ignore.case=T))
    cellphone_table = cellphone_table %>% separate(interactors, c('ligand','receptor'),sep='-')
    cellphone_table = cellphone_table %>% separate_rows(receptor, sep='\\+')
    cellphone_table = cellphone_table %>% separate_rows(ligand, sep='\\+')
    cellphone_table = cellphone_table %>% select(ligand, receptor)
    cellphone_table = cellphone_table %>% mutate(ligand=str_to_title(ligand), receptor=str_to_title(receptor))
    
    lr_db_tbl = rbind(cellchat_table, cellphone_table)
    lr_db_tbl = unique(lr_db_tbl)
    
    first = spatial_matrix_connected %>% 
      filter(grepl('KCs--',unified_int)) %>%
      separate(unified_int, c('receiver','sender'), sep='--') %>%
      select(sender, enrichm)
    second = spatial_matrix_connected %>% 
      filter(grepl('--KCs',unified_int)) %>%
      separate(unified_int, c('sender','receiver'), sep='--') %>%
      select(sender, enrichm)
    spatial_matrix_connected = rbind(first, second[second$sender!='KCs',])
    
    first = spatial_matrix_paracrine %>% 
      filter(grepl('KCs--',unified_int)) %>%
      separate(unified_int, c('receiver','sender'), sep='--') %>%
      select(sender, enrichm)
    second = spatial_matrix_paracrine %>% 
      filter(grepl('--KCs',unified_int)) %>%
      separate(unified_int, c('sender','receiver'), sep='--') %>%
      select(sender, enrichm)
    spatial_matrix_paracrine = rbind(first, second[second$sender!='KCs',])
    
    connected_proximity_prioritization = spatial_matrix_connected %>% select(sender, enrichm) %>% dplyr::distinct() %>%
      dplyr::mutate(scaled_enrich = scale_quantile_adapted(enrichm)) %>% dplyr::arrange(-scaled_enrich)
    paracrine_proximity_prioritization = spatial_matrix_paracrine %>% select(sender, enrichm) %>% dplyr::distinct() %>%
      dplyr::mutate(scaled_enrich = scale_quantile_adapted(enrichm)) %>% dplyr::arrange(-scaled_enrich)
    #dplyr::mutate(scaled_enrich = scale(enrichm, center=TRUE, scale=TRUE))
  } else {
    if (prioritizing_weights['spatial_proximity'] > 0) {
      stop("No spatial_matrix table given, yet the relevant weights are nonzero.\nEither set weights of 'spatial_proximity' to zero or provide spatial_matrix.")
    }
  }
  
  if (!is.null(lr_condition_de)){
    # Condition specificity of ligand (upregulation)
    ligand_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                    p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
    ligand_condition_prioritization = ligand_condition_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_ligand) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -ligand)
    
    # Condition specificity of receptor (upregulation)
    receptor_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                    p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor))
    receptor_condition_prioritization = receptor_condition_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                            scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_receptor) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -receptor)
    
  } else {
    if (any(prioritizing_weights[grep("specificity", names(prioritizing_weights))] > 0)) {
      stop("No lr_condition_de table given, yet the relevant weights are nonzero.\nEither set weights of 'ligand_condition_specificity' and 'receptor_condition_specificity' to zero or provide lr_condition_de.")
    }
  }
  
  weights <- prioritizing_weights
  # final group-based prioritization
  group_prioritization_tbl = sender_receiver_de %>%
    dplyr::inner_join(sender_receiver_info) %>%
    {if (weights["de_ligand"] > 0) dplyr::inner_join(., sender_ligand_prioritization) else (.)} %>%
    {if (weights["activity_scaled"] > 0) dplyr::inner_join(., ligand_activity_prioritization) else (.)} %>%
    {if (weights["de_receptor"] > 0) dplyr::inner_join(., receiver_receptor_prioritization) else (.)} %>%
    {if (weights["exprs_ligand"] > 0) dplyr::inner_join(., ligand_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["exprs_receptor"] > 0) dplyr::inner_join(., receptor_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["ligand_condition_specificity"] > 0) dplyr::inner_join(., ligand_condition_prioritization) else (.)} %>%
    {if (weights["receptor_condition_specificity"] > 0) dplyr::inner_join(., receptor_condition_prioritization) else (.)}
  
  
  for (i in 1:nrow(group_prioritization_tbl)){
    ligand_i = group_prioritization_tbl$ligand[i]
    receptor_i = group_prioritization_tbl$receptor[i]
    sender_i = group_prioritization_tbl$sender[i]
    lr_subset = lr_db_tbl %>% filter(grepl(ligand_i, ligand, fixed=T)) %>% filter(grepl(receptor_i, receptor, fixed=T))
    if (nrow(lr_subset)>0 & sender_i %in% spatial_matrix_connected$sender){
      group_prioritization_tbl$scaled_enrich[i] = spatial_matrix_connected[spatial_matrix_connected$sender==sender_i,]$scaled_enrich
    }else{
      group_prioritization_tbl$scaled_enrich[i] = spatial_matrix_paracrine[spatial_matrix_paracrine$sender==sender_i,]$scaled_enrich
    }
  }
  #group_prioritization_tbl = group_prioritization_tbl[!is.na(group_prioritization_tbl$scaled_enrich),]
  
  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = 2*weights["de_ligand"] + 2*weights["de_receptor"] + weights["activity_scaled"] + weights["exprs_ligand"] + weights["exprs_receptor"] + weights["ligand_condition_specificity"] + weights["receptor_condition_specificity"] + weights['spatial_proximity']
  group_prioritization_tbl = group_prioritization_tbl %>% rowwise() %>%
    dplyr::mutate(prioritization_score =
                    (
                      (prioritizing_weights["de_ligand"] * ifelse("scaled_lfc_ligand" %in% names(group_prioritization_tbl), scaled_lfc_ligand, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_lfc_receptor" %in% names(group_prioritization_tbl), scaled_lfc_receptor, 0)) +
                        (prioritizing_weights["de_ligand"] * ifelse("scaled_p_val_ligand_adapted" %in% names(group_prioritization_tbl), scaled_p_val_ligand_adapted, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_p_val_receptor_adapted" %in% names(group_prioritization_tbl), scaled_p_val_receptor_adapted, 0)) +
                        (prioritizing_weights["activity_scaled"] * ifelse("scaled_activity" %in% names(group_prioritization_tbl), scaled_activity, 0)) +
                        (prioritizing_weights["exprs_ligand"] * ifelse("scaled_avg_exprs_ligand" %in% names(group_prioritization_tbl), scaled_avg_exprs_ligand, 0)) +
                        (prioritizing_weights["exprs_receptor"] * ifelse("scaled_avg_exprs_receptor" %in% names(group_prioritization_tbl), scaled_avg_exprs_receptor, 0)) +
                        (prioritizing_weights["ligand_condition_specificity"] * ifelse("scaled_lfc_ligand_group" %in% names(group_prioritization_tbl), scaled_lfc_ligand_group, 0)) +
                        (prioritizing_weights["receptor_condition_specificity"] * ifelse("scaled_lfc_receptor_group" %in% names(group_prioritization_tbl), scaled_lfc_receptor_group, 0)) +
                        (prioritizing_weights["spatial_proximity"] * ifelse("scaled_enrich" %in% names(group_prioritization_tbl), scaled_enrich, 0))
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score) %>%
    ungroup()
  return(group_prioritization_tbl)
}

squidpy_prioritization_table = function(
    sender_receiver_info, sender_receiver_de, ligand_activities, 
    lr_condition_de = NULL, spatial_matrix = NULL,
    prioritizing_weights = c("de_ligand" = 1, "de_receptor" = 1,
                             "activity_scaled" = 2, "exprs_ligand" = 1,
                             "exprs_receptor" = 1,
                             "ligand_condition_specificity" = 0, 
                             "receptor_condition_specificity" = 0, 
                             "spatial_prioritization" = 0)){
  
  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
  
  # Ligand DE prioritization
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                  p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
    dplyr::arrange(-lfc_pval_ligand)
  
  # Receptor DE prioritization
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
    dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                  p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) )
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)
  
  # Ligand activity prioritization
  ligand_activity_prioritization = ligand_activities %>% select(test_ligand, aupr_corrected, rank) %>% rename(activity=aupr_corrected, ligand=test_ligand) %>%
    dplyr::mutate(activity_zscore = nichenetr::scaling_zscore(activity),
                  scaled_activity = scale_quantile_adapted(activity, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_zscore)
  
  
  # Cell-type specificity of expression of ligand:  per ligand scale across cell types
  ligand_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(sender, ligand, avg_ligand) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>%
    dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)
  
  # Cell-type specificity of expression of receptor:  per receptor scale across cell types
  receptor_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(receiver, receptor, avg_receptor) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>%
    dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)
  
  if (!is.null(spatial_matrix)){
    # Proximity (Spatial)
    #spatial_matrix = spatial_matrix$enrichm_res %>% separate(unified_int, c('sender','receiver'), sep='--')
    #proximity_prioritization = spatial_matrix %>% select(sender, receiver, enrichm, int_ranking) %>% dplyr::distinct() %>%
    #  dplyr::mutate(scaled_enrich = scale_quantile_adapted(enrichm)) %>% dplyr::arrange(-scaled_enrich)
    # first = spatial_matrix %>% 
    #   filter(grepl('KCs--',unified_int)) %>%
    #   separate(unified_int, c('receiver','sender'), sep='--') %>%
    #   select(sender, enrichm)
    # second = spatial_matrix %>% 
    #   filter(grepl('--KCs',unified_int)) %>%
    #   separate(unified_int, c('sender','receiver'), sep='--') %>%
    #   select(sender, enrichm)
    # spatial_matrix = rbind(first, second[second$sender!='KCs',])
    proximity_prioritization = spatial_matrix %>% select(sender, enrichm) %>% dplyr::distinct() %>%
      dplyr::mutate(scaled_enrich = scale_quantile_adapted(enrichm)) %>% dplyr::arrange(-scaled_enrich)
    #dplyr::mutate(scaled_enrich = scale(enrichm, center=TRUE, scale=TRUE))
  } else {
    if (prioritizing_weights['spatial_proximity'] > 0) {
      stop("No spatial_matrix table given, yet the relevant weights are nonzero.\nEither set weights of 'spatial_proximity' to zero or provide spatial_matrix.")
    }
  }
  
  if (!is.null(lr_condition_de)){
    # Condition specificity of ligand (upregulation)
    ligand_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                    p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
    ligand_condition_prioritization = ligand_condition_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_ligand) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -ligand)
    
    # Condition specificity of receptor (upregulation)
    receptor_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                    p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor))
    receptor_condition_prioritization = receptor_condition_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                            scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_receptor) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -receptor)
    
  } else {
    if (any(prioritizing_weights[grep("specificity", names(prioritizing_weights))] > 0)) {
      stop("No lr_condition_de table given, yet the relevant weights are nonzero.\nEither set weights of 'ligand_condition_specificity' and 'receptor_condition_specificity' to zero or provide lr_condition_de.")
    }
  }
  
  weights <- prioritizing_weights
  # final group-based prioritization
  group_prioritization_tbl = sender_receiver_de %>%
    dplyr::inner_join(sender_receiver_info) %>%
    {if (weights["de_ligand"] > 0) dplyr::inner_join(., sender_ligand_prioritization) else (.)} %>%
    {if (weights["activity_scaled"] > 0) dplyr::inner_join(., ligand_activity_prioritization) else (.)} %>%
    {if (weights["de_receptor"] > 0) dplyr::inner_join(., receiver_receptor_prioritization) else (.)} %>%
    {if (weights["exprs_ligand"] > 0) dplyr::inner_join(., ligand_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["exprs_receptor"] > 0) dplyr::inner_join(., receptor_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["ligand_condition_specificity"] > 0) dplyr::inner_join(., ligand_condition_prioritization) else (.)} %>%
    {if (weights["receptor_condition_specificity"] > 0) dplyr::inner_join(., receptor_condition_prioritization) else (.)} %>%
    {if (weights["spatial_proximity"] > 0) dplyr::inner_join(., proximity_prioritization) else (.)}
  
  
  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = 2*weights["de_ligand"] + 2*weights["de_receptor"] + weights["activity_scaled"] + weights["exprs_ligand"] + weights["exprs_receptor"] + weights["ligand_condition_specificity"] + weights["receptor_condition_specificity"] + weights['spatial_proximity']
  group_prioritization_tbl = group_prioritization_tbl %>% rowwise() %>%
    dplyr::mutate(prioritization_score =
                    (
                      (prioritizing_weights["de_ligand"] * ifelse("scaled_lfc_ligand" %in% names(group_prioritization_tbl), scaled_lfc_ligand, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_lfc_receptor" %in% names(group_prioritization_tbl), scaled_lfc_receptor, 0)) +
                        (prioritizing_weights["de_ligand"] * ifelse("scaled_p_val_ligand_adapted" %in% names(group_prioritization_tbl), scaled_p_val_ligand_adapted, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_p_val_receptor_adapted" %in% names(group_prioritization_tbl), scaled_p_val_receptor_adapted, 0)) +
                        (prioritizing_weights["activity_scaled"] * ifelse("scaled_activity" %in% names(group_prioritization_tbl), scaled_activity, 0)) +
                        (prioritizing_weights["exprs_ligand"] * ifelse("scaled_avg_exprs_ligand" %in% names(group_prioritization_tbl), scaled_avg_exprs_ligand, 0)) +
                        (prioritizing_weights["exprs_receptor"] * ifelse("scaled_avg_exprs_receptor" %in% names(group_prioritization_tbl), scaled_avg_exprs_receptor, 0)) +
                        (prioritizing_weights["ligand_condition_specificity"] * ifelse("scaled_lfc_ligand_group" %in% names(group_prioritization_tbl), scaled_lfc_ligand_group, 0)) +
                        (prioritizing_weights["receptor_condition_specificity"] * ifelse("scaled_lfc_receptor_group" %in% names(group_prioritization_tbl), scaled_lfc_receptor_group, 0)) +
                        (prioritizing_weights["spatial_proximity"] * ifelse("scaled_enrich" %in% names(group_prioritization_tbl), scaled_enrich, 0))
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score) %>%
    ungroup()
  return(group_prioritization_tbl)
}
