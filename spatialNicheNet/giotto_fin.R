giotto_spatial_proximity = function(a1,a2,c2,d2,average=T){
  
  library(Giotto)
  library(Seurat)
  library(reshape2)
  library(tidyverse)
  
  # Python path
  instr <- createGiottoInstructions(python_path = 
    "C:/Users/yunseolp/Documents/.virtualenvs/r-reticulate/Scripts/python.exe")
  
  # A1-1
  Idents(a1) <- a1$annotationSave
  a1_giotto_obj <- createGiottoObject(
    GetAssayData(a1, layer="counts"),
    norm_expr = GetAssayData(a1),
    spatial_locs = a1@reductions[["spatial"]]@cell.embeddings,
    cell_metadata = a1@meta.data, instructions = instr
    )
  a1_giotto_obj <- createSpatialNetwork(gobject = a1_giotto_obj, 
                                        maximum_distance_delaunay = 200)
  a1_cell_proximities_deluanay <- cellProximityEnrichment(
    gobject = a1_giotto_obj,
    cluster_column = 'annotationSave',
    spatial_network_name = 'Delaunay_network',
    adjust_method = 'fdr',
    number_of_simulations = 1000,
    set_seed = TRUE,
    seed_number = 0
    )
  #A1-2
  Idents(a2) <- a2$annotationSave
  a2_giotto_obj <- createGiottoObject(
    GetAssayData(a2, layer="counts"),
    norm_expr = GetAssayData(a2),
    spatial_locs = a2@reductions[["spatial"]]@cell.embeddings,
    cell_metadata = a2@meta.data, instructions = instr
    )
  a2_giotto_obj <- createSpatialNetwork(gobject = a2_giotto_obj, 
                                        maximum_distance_delaunay = 200)
  a2_cell_proximities_deluanay <- cellProximityEnrichment(
    gobject = a2_giotto_obj,
    cluster_column = 'annotationSave',
    spatial_network_name = 'Delaunay_network',
    adjust_method = 'fdr',
    number_of_simulations = 1000,
    set_seed = TRUE,
    seed_number = 0
    )
  # C2
  Idents(c2) <- c2$annotationSave
  c2_giotto_obj <- createGiottoObject(
    GetAssayData(c2, layer="counts"),
    norm_expr = GetAssayData(c2),
    spatial_locs = c2@reductions[["spatial"]]@cell.embeddings,
    cell_metadata = c2@meta.data, instructions = instr
    )
  c2_giotto_obj <- createSpatialNetwork(gobject = c2_giotto_obj, 
                                        maximum_distance_delaunay = 200)
  c2_cell_proximities_deluanay <- cellProximityEnrichment(
    gobject = c2_giotto_obj,
    cluster_column = 'annotationSave',
    spatial_network_name = 'Delaunay_network',
    adjust_method = 'fdr',
    number_of_simulations = 1000,
    set_seed = TRUE,
    seed_number = 0
    )
  # D2-1
  Idents(d2) <- d2$annotationSave
  d2_giotto_obj <- createGiottoObject(
    GetAssayData(d2, layer="counts"),
    norm_expr = GetAssayData(d2),
    spatial_locs = d2@reductions[["spatial"]]@cell.embeddings,
    cell_metadata = d2@meta.data, instructions = instr
    )
  d2_giotto_obj <- createSpatialNetwork(gobject = d2_giotto_obj, 
                                        maximum_distance_delaunay = 200)
  d2_cell_proximities_deluanay <- cellProximityEnrichment(
    gobject = d2_giotto_obj,
    cluster_column = 'annotationSave',
    spatial_network_name = 'Delaunay_network',
    adjust_method = 'fdr',
    number_of_simulations = 1000,
    set_seed = TRUE,
    seed_number = 0
    )
  
  avg_delaunay = list(a1_cell_proximities_deluanay$enrichm_res, 
                      a2_cell_proximities_deluanay$enrichm_res,
                      c2_cell_proximities_deluanay$enrichm_res, 
                      d2_cell_proximities_deluanay$enrichm_res) %>%
    reduce(full_join, by='unified_int')
  avg_delaunay$enrichm = avg_delaunay %>% select(matches('enrichm')) %>%
    rowMeans(na.rm=T)
  
  return(c(a1_cell_proximities_deluanay, a2_cell_proximities_deluanay, 
         c2_cell_proximities_deluanay, d2_cell_proximities_deluanay,
         avg_delaunay))
}

giotto_diff = function(a1){
  library(Giotto)
  library(Seurat)
  library(reshape2)
  library(tidyverse)
  
  # Python path
  instr <- createGiottoInstructions(python_path = 
                                      "C:/Users/yunseolp/Documents/.virtualenvs/r-reticulate/Scripts/python.exe")
  
  # A1-1
  Idents(a1) <- a1$annotationSave
  a1_giotto_obj <- createGiottoObject(
    GetAssayData(a1, layer="counts"),
    norm_expr = GetAssayData(a1),
    spatial_locs = a1@reductions[["spatial"]]@cell.embeddings,
    cell_metadata = a1@meta.data, instructions = instr
  )
  a1_giotto_obj <- createSpatialNetwork(gobject = a1_giotto_obj, maximum_distance_delaunay = (10+15)/0.138, name='juxtacrine_network')
  a1_cell_proximities_juxtacrine <- cellProximityEnrichment(gobject = a1_giotto_obj,
                                                           cluster_column = 'annotationSave',
                                                           spatial_network_name = 'juxtacrine_network',
                                                           adjust_method = 'fdr',
                                                           number_of_simulations = 1000,
                                                           set_seed = TRUE,
                                                           seed_number = 0)
  a1_giotto_obj <- createSpatialNetwork(gobject = a1_giotto_obj, maximum_distance_delaunay = (20*25+15)/0.138, name='paracrine_network')
  a1_cell_proximities_paracrine <- cellProximityEnrichment(gobject = a1_giotto_obj,
                                                           cluster_column = 'annotationSave',
                                                           spatial_network_name = 'paracrine_network',
                                                           adjust_method = 'fdr',
                                                           number_of_simulations = 1000,
                                                           set_seed = TRUE,
                                                           seed_number = 0)
  a1_giotto_obj <- createSpatialNetwork(gobject = a1_giotto_obj, maximum_distance_delaunay = 15000, name='endocrine_network')
  a1_cell_proximities_endocrine <- cellProximityEnrichment(gobject = a1_giotto_obj,
                                                           cluster_column = 'annotationSave',
                                                           spatial_network_name = 'endocrine_network',
                                                           adjust_method = 'fdr',
                                                           number_of_simulations = 1000,
                                                           set_seed = TRUE,
                                                           seed_number = 0)
  
  return(list(a1_cell_proximities_juxtacrine, a1_cell_proximities_paracrine, 
              a1_cell_proximities_endocrine))
}