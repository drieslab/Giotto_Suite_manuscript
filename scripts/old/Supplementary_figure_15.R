######################## Download Visium Brain data #######################

# Download data from https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain


################################ Load Libraries ################################

library(terra)
library(patchwork)
library(Giotto)

######################## Instantiate Data/Save Locations #######################

my_seed_num <- 315
data_directory <- "data"
save_directory <- "results"
subcell_save_directory <- "pseudo_subcellular"

############################# Create Instructions ##############################

instrs <- createGiottoInstructions(show_plot = TRUE,
                                   save_plot = TRUE,
                                   save_dir = save_directory)

####################### New Functions for below analysis #######################

merge_expr_with_spat_locs <- function(gobject,
                                      gene_name){

  expr = getExpression(gobject = gobject)
  target_expr = expr[gene_name][]
  ids_expressing_gene = names(target_expr)

  target_expr = data.table::as.data.table(target_expr)
  data.table::setnames(target_expr, old = "target_expr", new = "count")
  target_expr$cell_ID = ids_expressing_gene

  all_spat_locs = getSpatialLocations(gobject = gobject)
  gene_spat_locs = all_spat_locs[][cell_ID == ids_expressing_gene]

  gene_spatial_expression = data.table::merge.data.table(x = target_expr,
                                                         y = gene_spat_locs,
                                                         by = "cell_ID")
  data.table::setnames(gene_spatial_expression, old = c("sdimx", "sdimy"), new = c("x","y"))

  return(gene_spatial_expression)
}


rasterize_expression_data <- function(gobject,
                                      gene_name,
                                      poly_info_to_use,
                                      num_row,
                                      num_col){

  gene_spatial_expression = merge_expr_with_spat_locs(gobject = gobject,
                                                      gene_name = gene_name)

  reference_polys = getPolygonInfo(gobject = gobject,
                                   polygon_name = poly_info_to_use)

  target_ext = terra::ext(reference_polys)

  fill_rast = terra::rast(target_ext,
                          nrow = num_row,
                          ncol = num_col)

  rast_locs = as.matrix(gene_spatial_expression[,.(x, y)])
  rast_vals = as.vector(gene_spatial_expression$count)

  rasterized_expr = terra::rasterize(x = rast_locs, # locations
                                     y = fill_rast, # raster to populate
                                     values = rast_vals) # values to use
  names(rasterized_expr) = "count"

  return(rasterized_expr)

}

############### Begin Analysis, or Reload From Previous Session ################

RERUN_FROM_ORIGINAL_DATA <- TRUE

if (RERUN_FROM_ORIGINAL_DATA){



  ############################# Create Giotto Object #############################

  v_brain = createGiottoVisiumObject(data_directory, gene_column_index = 2)

  ######### Filter, normalize, ID hvf, dimension reductions, clustering ##########


  # Subset to in tissue only
  cm = pDataDT(v_brain)
  in_tissue_barcodes = cm[in_tissue == 1]$cell_ID
  v_brain = subsetGiotto(v_brain, cell_ids = in_tissue_barcodes)

  # Filter
  v_brain = filterGiotto(gobject = v_brain,
                         expression_threshold = 1,
                         feat_det_in_min_cells = 50,
                         min_det_feats_per_cell = 1000,
                         expression_values = c('raw'))

  # Normalize
  v_brain = normalizeGiotto(gobject = v_brain,
                            scalefactor = 6000,
                            verbose = TRUE)

  # Add stats
  v_brain = addStatistics(gobject = v_brain)

  # ID HVF
  v_brain = calculateHVF(gobject = v_brain)
  fm = fDataDT(v_brain)
  hv_feats = fm[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID
  length(hv_feats) # Matt's 579 # Joss: 538

  # Dimension Reductions
  v_brain = runPCA(gobject = v_brain,
                   feats_to_use = hv_feats)

  v_brain = runUMAP(v_brain,
                    dimensions_to_use = 1:10,
                    n_neighbors = 15,
                    set_seed = TRUE,
                    seed_number = my_seed_num)

  # NN Network
  v_brain = createNearestNetwork(gobject = v_brain,
                                 dimensions_to_use = 1:10,
                                 k = 15)
  # Leiden Cluster
  v_brain = doLeidenCluster(gobject = v_brain,
                            resolution = 0.4,
                            n_iterations = 200,
                            seed_number = my_seed_num)

  # Spatial Network (kNN)
  v_brain <- createSpatialNetwork(gobject = v_brain,
                                  method = 'kNN',
                                  k = 5,
                                  maximum_distance_knn = 400,
                                  name = 'spatial_network')

  # Spatially Variable Features
  ranktest = binSpect(v_brain,
                      bin_method = 'rank',
                      calc_hub = TRUE,
                      hub_min_int = 5,
                      spatial_network_name = 'spatial_network',
                      seed = FALSE) #not able to provide a seed number, so do not set one

  # cluster the top 1500 spatial features into 20 clusters
  ext_spatial_features = ranktest[1:500,]$feats

  # Find pairwise distances between these top 1500 features
  spat_cor_netw_DT = detectSpatialCorFeats(v_brain,
                                           method = 'network',
                                           spatial_network_name = 'spatial_network',
                                           subset_feats = ext_spatial_features)


  spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT,
                                            name = 'spat_netw_clus',
                                            k = 20)

  heatmSpatialCorFeats(v_brain,
                       spatCorObject = spat_cor_netw_DT,
                       use_clus_name = 'spat_netw_clus',
                       heatmap_legend_param = list(title = NULL),
                       show_plot = TRUE,
                       save_plot = TRUE,
                       return_plot = FALSE)


  ########### Rasterize some subset of HVF or Spatially Variable Genes ###########

  fm_dt = getFeatureMetadata(v_brain, output = "data.table")
  fm_dt = fm_dt[hvf == "yes"]
  data.table::setorder(fm_dt, cols = -"mean_expr")


  fxn_out = rasterize_expression_data(gobject = v_brain,
                                      gene_name = "Plp1",
                                      poly_info_to_use = "cell",
                                      num_row = 50,
                                      num_col = 50)

  ############### Krige on the "cell" polys for some subset of hvf ###############

  finer_res_rast = terra::rast(ext(fxn_out), nrow = 100, ncol = 100)

  gene_spat_expr = merge_expr_with_spat_locs(gobject = v_brain,
                                             gene_name = "Plp1")

  count_model = gstat::gstat(id = "count",
                             formula = count~1,
                             locations = ~x+y,
                             data=gene_spat_expr,
                             nmax=7,
                             set=list(idp = .5))

  pred_gene = terra::interpolate(object = finer_res_rast,
                                 model = count_model,
                                 debug.level=0,
                                 index=1)


  ################ Add interpolated gene expression as gimage(s) #################

  pred_plp1_gimg = createGiottoLargeImage(pred_gene,
                                          name = "Plp1",
                                          extent = ext(pred_gene),
                                          use_rast_ext = TRUE)

  v_brain = addGiottoLargeImage(gobject = v_brain,
                                largeImages = list("Plp1" = pred_plp1_gimg))
} else {
  v_brain = loadGiotto("saveGiottoDir/")
}

######## Create/Load GiottoLargeImages with Predicted Gene Information #########

RECALCULATE_PREDICTED_GENES = TRUE

if (RECALCULATE_PREDICTED_GENES){
  tictoc::tic()
  g_limg_list = list()
  spat_var_genes = spat_cor_netw_DT$feat_order
  idx = 1

  for (gene in spat_var_genes){

    rast_gene_expr = rasterize_expression_data(gobject = v_brain,
                                               gene_name = gene,
                                               poly_info_to_use = "cell",
                                               num_row = 50,
                                               num_col = 50)

    finer_res_rast = terra::rast(ext(rast_gene_expr),
                                 nrow = 100,
                                 ncol = 100)

    gene_spat_expr = merge_expr_with_spat_locs(gobject = v_brain,
                                               gene_name = gene)

    count_model = gstat::gstat(id = "count",
                               formula = count~1,
                               locations = ~x+y,
                               data=gene_spat_expr,
                               nmax=7,
                               set=list(idp = .5))

    pred_gene = terra::interpolate(object = finer_res_rast,
                                   model = count_model,
                                   debug.level=0,
                                   index=1)

    pred_gene_gimg = createGiottoLargeImage(pred_gene,
                                            name = gene,
                                            extent = ext(pred_gene),
                                            use_rast_ext = TRUE,
                                            verbose = F)

    g_limg_list = append(g_limg_list, pred_gene_gimg)

    cat(paste0(as.character(idx),": ",gene,"\n"))

    idx = idx + 1

  }

  names(g_limg_list) = spat_var_genes
  tictoc::toc()

  ########################### SAVE PREDICTED EXPRESSION ##########################

  # Save rasterized predicted gene expression as tif files
  for (gene in spat_var_genes){

    f = file.path(paste0("predicted_expression_images/", gene, ".tif"))

    terra::writeRaster(g_limg_list[[gene]]@raster_object,
                       filename = f)

  }

  # save a .csv file with the name of the gene and the name of the file for streamlined loading
  link_gene_file = list()
  for( i in spat_var_genes){
    link_gene_file = append(link_gene_file, paste0(i,".tif"))
  }
  names(link_gene_file) = spat_var_genes

  data.table::fwrite(link_gene_file,
                     file = "predicted_expression_images/link_gene_file.csv")


} else{

  # Load link file:
  # format 2xN
  # gene_name    |    gene_file.tif

  genes_files = data.table::fread(file = "predicted_expression_images/link_gene_file.csv", header = FALSE)

  gene_names = as.vector(genes_files[1])
  file_names = as.list(paste0("predicted_expression_images/",genes_files[2]))
  names(file_names) = gene_names


  reloaded_g_limg_list = list()

  # Create Giotto Images for each .tif file
  for (file in file_names){
    x = createGiottoLargeImage(terra::rast(file[[1]]),
                               use_rast_ext = TRUE,
                               name = names(file),
                               verbose = FALSE)
    reloaded_g_limg_list = append(reloaded_g_limg_list, x)
    cat(paste0("Loaded ",file,"\n"))
  }

  names(reloaded_g_limg_list) = names(file_names)
  spat_var_genes = names(file_names)
  g_limg_list = reloaded_g_limg_list

}

####################### Setup Subcellular Giotto Object ########################

stardist_cell_poly_path = "segmentations/stardist_only_cell_bounds.geojson"
stardist_cell_gpoly = createGiottoPolygonsFromGeoJSON(GeoJSON = stardist_cell_poly_path,
                                                      name = "stardist_cell",
                                                      calc_centroids = TRUE)

bad_id = GiottoClass:::.identify_background_range_polygons(stardist_cell_gpoly[])
keep_these = spatIDs(stardist_cell_gpoly)[spatIDs(stardist_cell_gpoly) != bad_id]
stardist_cell_gpoly = stardist_cell_gpoly[keep_these]

stardist_cell_gpoly[] = subset(x = stardist_cell_gpoly[],
                               stardist_cell_gpoly[][["poly_ID"]] != bad_id)

stardist_cell_gpoly@spatVectorCentroids[] = subset(x = stardist_cell_gpoly@spatVectorCentroids[],
                                                   stardist_cell_gpoly@spatVectorCentroids[][["poly_ID"]] != bad_id)
stardist_cell_gpoly = flip(stardist_cell_gpoly)


sLocsObj = getSpatialLocations(gobject = v_brain, spat_unit = "cell")

vis_scale = Giotto:::.visium_read_scalefactors(json_path = paste0(data_directory,"/spatial/scalefactors_json.json"))

vis_gpolys = Giotto:::.visium_spot_poly(spatlocs = sLocsObj,
                                        json_scalefactors = vis_scale)

subcell_instrs = createGiottoInstructions(save_plot = TRUE,
                                          show_plot = TRUE,
                                          save_dir = subcell_save_directory)

####################### Create Subcellular Giotto Object #######################

subcell_v_brain = createGiottoObjectSubcellular(gpolygons = list("cell" = vis_gpolys,
                                                                 "stardist_cell" = stardist_cell_gpoly),
                                                largeImages = g_limg_list,
                                                instructions = subcell_instrs)

hires_gimg = getGiottoImage(gobject = v_brain,
                            image_type = "largeImage")

subcell_v_brain = addGiottoImage(gobject = subcell_v_brain,
                                 largeImages = list(hires_gimg))

############## Find overlap - make expression data for "cell" polys ############
tictoc::tic()
subcell_v_brain = calculateOverlapPolygonImages(gobject = subcell_v_brain,
                                                name_overlap = "rna",
                                                spatial_info = "cell",
                                                image_names = spat_var_genes)

subcell_v_brain = overlapImagesToMatrix(gobject = subcell_v_brain,
                                        poly_info = "cell",
                                        feat_info = "rna",
                                        aggr_function = "sum",
                                        image_names = spat_var_genes)
tictoc::toc() # 45.132 sec elapsed
######## Find overlap - make expression data for "stardist_cell" polys #########
tictoc::tic()
subcell_v_brain = calculateOverlapPolygonImages(gobject = subcell_v_brain,
                                                name_overlap = "rna",
                                                spatial_info = "stardist_cell",
                                                image_names = spat_var_genes)

subcell_v_brain = overlapImagesToMatrix(gobject = subcell_v_brain,
                                        poly_info = "stardist_cell",
                                        feat_info = "rna",
                                        aggr_function = "sum",
                                        image_names = spat_var_genes)
tictoc::toc() # 553.777 sec elapsed



###
initialize(subcell_v_brain)


############# Subcellular, Cell: Dimension reductions and Clustering ###########

subcell_v_brain = normalizeGiotto(gobject = subcell_v_brain,
                                  spat_unit = "cell",
                                  scalefactor = 6000,
                                  verbose = T)

# PCA
subcell_v_brain = runPCA(gobject = subcell_v_brain,
                         spat_unit = "cell",
                         set_seed = TRUE,
                         seed_number = my_seed_num)

# UMAP
subcell_v_brain = runUMAP(subcell_v_brain,
                          spat_unit = "cell",
                          dimensions_to_use = 1:10,
                          n_neighbors = 25,
                          min_dist = 0.5,
                          spread = 1,
                          seed_number = my_seed_num)

# NN Network
subcell_v_brain = createNearestNetwork(gobject = subcell_v_brain,
                                       spat_unit = "cell",
                                       dimensions_to_use = 1:10,
                                       k = 25)
RECALCULATE_CLUSTERS_CELL = TRUE
if(RECALCULATE_CLUSTERS_CELL){
  # Leiden Cluster
  subcell_v_brain = doLeidenCluster(gobject = subcell_v_brain,
                                    spat_unit = "cell",
                                    resolution = 0.15,
                                    n_iterations = 100,
                                    seed_number = my_seed_num)

  subcell_v_brain = doKmeans(gobject = subcell_v_brain,
                             spat_unit = "cell",
                             feat_type = "rna",
                             centers = 12,
                             seed_number = my_seed_num)
} else{

  # Load Saved Data, created with seed = 315
  og_cm = data.table::fread(paste0(subcell_save_directory, "og_cm_cluster_data.csv"))

  subcell_v_brain = addCellMetadata(gobject = subcell_v_brain,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    new_metadata = og_cm,
                                    by_column = TRUE,
                                    column_cell_ID = "cell_ID")

}



########### Subcellular, StarDist: Dimension reductions and Clustering #########

# Filter out any polygons with expression of 0
subcell_v_brain = filterGiotto(gobject = subcell_v_brain,
                               spat_unit = "stardist_cell",
                               expression_values = "raw",
                               expression_threshold = 1,
                               feat_det_in_min_cells = 0,
                               min_det_feats_per_cell = 1)


subcell_v_brain = normalizeGiotto(gobject = subcell_v_brain,
                                  spat_unit = "stardist_cell",
                                  scalefactor = 6000,
                                  verbose = T)

# PCA
subcell_v_brain = runPCA(gobject = subcell_v_brain,
                         spat_unit = "stardist_cell",
                         expression_values = "normalized",
                         feats_to_use = NULL,
                         set_seed = TRUE,
                         seed_number = my_seed_num)

# UMAP
subcell_v_brain = runUMAP(subcell_v_brain,
                          spat_unit = "stardist_cell",
                          dimensions_to_use = 1:10,
                          n_neighbors = 25,
                          min_dist = 0.01,
                          spread = 1,
                          seed_number = my_seed_num)

# NN Network
subcell_v_brain = createNearestNetwork(gobject = subcell_v_brain,
                                       spat_unit = "stardist_cell",
                                       dimensions_to_use = 1:10,
                                       k = 25)
RECALCULATE_CLUSTERS_STARDIST = TRUE
if(RECALCULATE_CLUSTERS_STARDIST){

  # Leiden Cluster
  subcell_v_brain = doLeidenCluster(gobject = subcell_v_brain,
                                    spat_unit = "stardist_cell",
                                    resolution = 0.2,
                                    n_iterations = 100,
                                    seed_number = my_seed_num)
  # K-Means
  tictoc::tic()
  subcell_v_brain = doKmeans(gobject = subcell_v_brain,
                             spat_unit = "stardist_cell",
                             feat_type = "rna",
                             centers = 12,
                             iter_max = 50,
                             nstart = 100,
                             seed_number = my_seed_num)
  tictoc::toc() # 481.145 sec elapsed

  sd_cm = subcell_v_brain@cell_metadata$stardist_cell$rna[]
  data.table::fwrite(sd_cm, file = paste0(subcell_save_directory,"sd_cm_cluster_data_v3.csv"))

} else{

  # Load Saved Data, created with seed = 315
  sd_cm = data.table::fread(paste0(subcell_save_directory, "sd_cm_cluster_data_v3.csv"))

  subcell_v_brain = addCellMetadata(gobject = subcell_v_brain,
                                    spat_unit = "stardist_cell",
                                    feat_type = "rna",
                                    new_metadata = sd_cm,
                                    by_column = TRUE,
                                    column_cell_ID = "cell_ID")
}



############################ Write out cluster data ############################

# og_cm = subcell_v_brain@cell_metadata$cell$rna[]
# sd_cm = subcell_v_brain@cell_metadata$stardist_cell$rna[]
#
# data.table::fwrite(og_cm, file = paste0(save_directory,"og_cm_cluster_data.csv"))
# data.table::fwrite(sd_cm, file = paste0(save_directory,"sd_cm_cluster_data.csv"))


############################ Visualize and Compare #############################

# COLORS FOR k=12
my_colors = getDistinctColors(12)
og_colors = my_colors

names(my_colors) = seq(12)
names(og_colors) = c(11, 9, 2, 5, 12, 1, 8, 10, 4, 6, 7, 3)

c_spat = spatPlot2D(subcell_v_brain,
                    spat_unit = "cell",
                    cell_color = "leiden_clus",
                    point_size = 1.25,
                    show_image = TRUE,
                    largeImage_name = "image",
                    show_legend = FALSE,
                    return_plot = TRUE) +
  ggplot2::labs(title = "Original Visium Spots, Leiden")

sd_spat = spatPlot2D(subcell_v_brain,
                     spat_unit = "stardist_cell",
                     cell_color = "leiden_clus",
                     point_size = 1,
                     show_image = TRUE,
                     largeImage_name = "image",
                     show_legend = FALSE,
                     return_plot = TRUE) +
  ggplot2::labs(title = "StarDist Cells, Leiden")


c_spat + sd_spat

c_spat_k = spatPlot2D(subcell_v_brain,
                      spat_unit = "cell",
                      cell_color = "kmeans",
                      point_size = 1.5,
                      show_image = TRUE,
                      largeImage_name = "image",
                      cell_color_code = og_colors,
                      show_legend = FALSE,
                      title = "Original Visium Spots, k=12",
                      return_plot = TRUE) #+
#  ggplot2::labs(title = "Original Visium Spots, k=12")

sd_spat_k = spatPlot2D(subcell_v_brain,
                       spat_unit = "stardist_cell",
                       cell_color = "kmeans",
                       point_size = 1,
                       show_image = TRUE,
                       largeImage_name = "image",
                       cell_color_code = my_colors,
                       show_legend = FALSE,
                       title = "StarDist Cells, k=12",
                       return_plot = TRUE) #+
#  ggplot2::labs(title = "StarDist Cells, k=12")

c_spat_k + sd_spat_k

c_isspat_k = spatInSituPlotPoints(gobject = subcell_v_brain,
                                  spat_unit = "cell",
                                  feat_type = "rna",
                                  feats = NULL,
                                  show_image = T,
                                  largeImage_name = "image",
                                  show_polygon = TRUE,
                                  polygon_feat_type = "cell",
                                  polygon_fill = "kmeans",
                                  polygon_fill_as_factor = TRUE,
                                  polygon_color = "white",
                                  polygon_line_size = 0.1,
                                  polygon_alpha = 0.8,
                                  polygon_fill_code = og_colors,
                                  return_plot = TRUE)+
  ggplot2::labs(title = "Original Visium Spots, k=12")

sd_isspat_k = spatInSituPlotPoints(gobject = subcell_v_brain,
                                   spat_unit = "stardist_cell",
                                   feat_type = "rna",
                                   feats = NULL,
                                   show_image = T,
                                   largeImage_name = "image",
                                   show_polygon = TRUE,
                                   polygon_feat_type = "stardist_cell",
                                   polygon_fill = "kmeans",
                                   polygon_fill_as_factor = TRUE,
                                   polygon_color = "white",
                                   polygon_line_size = 0.05,
                                   polygon_alpha = 0.9,
                                   polygon_fill_code = my_colors,
                                   return_plot = TRUE) +
  ggplot2::labs(title = "StarDist Cells, k=12")

c_isspat_k + sd_isspat_k

subset_scvb_cell = subsetGiottoLocs(gobject = subcell_v_brain,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    feat_type_ssub = "rna",
                                    x_min = 2500, x_max = 3500,
                                    y_min = -5000, y_max = -3000,
                                    poly_info = "cell",
                                    spat_loc_name = "raw")

subset_scvb_star = subsetGiottoLocs(gobject = subcell_v_brain,
                                    spat_unit = "stardist_cell",
                                    feat_type = "rna",
                                    feat_type_ssub = "rna",
                                    x_min = 2500, x_max = 3500,
                                    y_min = -5000, y_max = -3000,
                                    poly_info = "stardist_cell",
                                    spat_loc_name = "raw")

c_sub_isspat_k = spatInSituPlotPoints(gobject = subset_scvb_cell,
                                      spat_unit = "cell",
                                      feat_type = "rna",
                                      feats = NULL,
                                      show_image = T,
                                      largeImage_name = "image",
                                      show_polygon = TRUE,
                                      polygon_feat_type = "cell",
                                      polygon_fill = "kmeans",
                                      polygon_fill_as_factor = TRUE,
                                      polygon_color = "black",
                                      polygon_line_size = 0.15,
                                      polygon_alpha = 0.75,
                                      polygon_fill_code = og_colors,
                                      return_plot = TRUE) +
  ggplot2::labs(title = "Original Visium Spots, k=12")

sd_sub_isspat_k = spatInSituPlotPoints(gobject = subset_scvb_star,
                                       spat_unit = "stardist_cell",
                                       feat_type = "rna",
                                       feats = NULL,
                                       show_image = T,
                                       largeImage_name = "image",
                                       show_polygon = TRUE,
                                       polygon_feat_type = "stardist_cell",
                                       polygon_fill = "kmeans",
                                       polygon_fill_as_factor = TRUE,
                                       polygon_color = "black",
                                       polygon_line_size = 0.075,
                                       polygon_alpha = 0.75,
                                       polygon_fill_code = my_colors,
                                       return_plot = TRUE) +
  ggplot2::labs(title = "StarDist Cells, k=12")

c_sub_isspat_k + sd_sub_isspat_k

c_feat_pl = spatFeatPlot2D(gobject = subcell_v_brain,
                           spat_unit = "cell",
                           feats = "Pacs2",
                           point_size = 2,
                           show_legend = FALSE,
                           return_plot = TRUE)+
  ggplot2::labs(title = "Pacs2, Original")
# Note from Joss, gene Pantr1 doesn't exist, I had to replace it
sd_feat_pl = spatFeatPlot2D(gobject = subcell_v_brain,
                            spat_unit = "stardist_cell",
                            feats = "Pacs2",
                            point_size = 0.66,
                            show_legend = FALSE,
                            return_plot = TRUE) +
  ggplot2::labs(title = "Pacs2, Super Enhanced")

plotUMAP(subcell_v_brain,
         spat_unit = "cell",
         cell_color = "leiden_clus",
         save_plot = TRUE,
         save_param = list(save_name = 'umap_cell_leiden'))

plotUMAP(subcell_v_brain,
         spat_unit = "stardist_cell",
         cell_color = "leiden_clus",
         save_plot = TRUE,
         save_param = list(save_name = 'umap_stardist_cell_leiden'))

################################## Save Plots ##################################

ggplot2::ggsave(filename = "s15C_original_visium_k12.png",
                plot = c_isspat_k,
                device = "png",
                width = 5,
                height = 6)

ggplot2::ggsave(filename = "s15C_stardist_cells_k12.png",
                plot = sd_isspat_k,
                device = "png",
                width = 5,
                height = 6)

ggplot2::ggsave(filename = "s15C_subset_original_visium_k12.png",
                plot = c_sub_isspat_k,
                device = "png",
                width = 4,
                height = 5)

ggplot2::ggsave(filename = "s15C_subset_stardist_cells_k12.png",
                plot = sd_sub_isspat_k,
                device = "png",
                width = 4,
                height = 5)

ggplot2::ggsave(filename = "s15D_original_spatFeatPlot_Pantr1.png",
                plot = c_feat_pl,
                device = "png",
                width = 7,
                height = 5)

ggplot2::ggsave(filename = "s15D_stardist_spatFeatPlot_Pantr1.png",
                plot = sd_feat_pl,
                device = "png",
                width = 7,
                height = 5)



