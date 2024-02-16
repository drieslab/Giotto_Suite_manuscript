# 231108
library(data.table)
library(Giotto)

# Pastel palette from rcartocolor
pal10 = getColors('Pastel', n = 10, src = 'rcartocolor')
viv10 = getColors('Vivid', n = 10, src = 'rcartocolor')

# For convenience function gobject creation
saveDir = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/SCRIPTS/FIGURES/nanostring_supp/plots/'
instrs = createGiottoInstructions(save_plot = T,
                                  show_plot = F,
                                  return_plot = F,
                                  save_dir = saveDir)



# Test both listed methods of gobject creation
# ---------------- Convenience Function Gobject Creation ####

## provide path to nanostring folder
data_path = '/projectnb/rd-spat/DATA/Public_data/Nanostring/Lung12/Lung12-Flat_files_and_images/'

## create giotto cosmx object
fov_join = createGiottoCosMxObject(cosmx_dir = data_path,
                                   data_to_use = 'subcellular',
                                   instructions = instrs,
                                   cores = 8)


# NOTE: cell0 is those counts that are not within any cells






showGiottoImageNames(fov_join)

id_set = sprintf('%03d', 1:28)
new_names = paste0("fov", id_set)
image_names = paste0(new_names, '-image')


spatInSituPlotPoints(fov_join,
                     show_image = TRUE,
                     image_name = image_names,
                     feats = list('rna' = c('MMP2', 'VEGFA', 'IGF1R',
                                            'MKI67', 'EPCAM', 'KRT8')),
                     feats_color_code = viv10,
                     spat_unit = 'cell',
                     point_size = 0.01,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_alpha = 0.1,
                     polygon_line_size = 0.03,
                     save_param = list(base_height = 10,
                                       save_name = '1_inSituFeats'))
# 41 GB ram usage




spatPlot2D(gobject = fov_join,
           image_name = image_names,
           show_image = TRUE,
           point_shape = 'no_border',
           point_size = 0.01,
           point_alpha = 0.5,
           coord_fix_ratio = 1,
           save_param = list(base_height = 10,
                             save_name = '2_spatCentroids'))

# ---------------- Downstream Convenience Function ####

fov_join = calculateOverlapRaster(fov_join, feat_info = 'rna')
fov_join = calculateOverlapRaster(fov_join, feat_info = 'neg_probe')

fov_join = overlapToMatrix(fov_join, feat_info = 'rna')
fov_join = overlapToMatrix(fov_join, feat_info = 'neg_probe')

filterDistributions(fov_join,
                    plot_type = 'hist',
                    detection = 'cells',
                    method = 'sum',
                    feat_type = 'rna',
                    nr_bins = 100,
                    save_param = list(base_height = 3,
                                      save_name = '3.1_totalexpr'))

filterDistributions(fov_join,
                    plot_type = 'hist',
                    detection = 'cells',
                    method = 'sum',
                    feat_type = 'neg_probe',
                    nr_bins = 25,
                    save_param = list(base_height = 3,
                                      save_name = '3.2_totalnegprbe'))

# Extract data from giotto object
# combine cell data
morphometa = combineCellData(fov_join,
                             feat_type = 'rna')

# combine feature data
featmeta = combineFeatureData(fov_join,
                              feat_type = c('rna'))

# combine overlapping feature data
featoverlapmeta = combineFeatureOverlapData(fov_join, #  poly_ID NA???
                                            feat_type = c('rna'))



showGiottoExpression(fov_join)


# these density functions require overlap info
spatInSituPlotDensity(gobject = fov_join,
                      feats = c('MKI67', 'KRT8'),
                      cow_n_col = 2,
                      save_param = list(base_height = 10,
                                        save_name = '4_inSituDens'))




fov_join <- filterGiotto(gobject = fov_join,
                         feat_type = 'rna',
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 5)

# normalize
# standard method of normalization (log normalization based)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            norm_methods = 'standard',
                            verbose = TRUE)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'neg_probe',
                            norm_methods = 'standard',
                            library_size_norm = FALSE,
                            verbose = TRUE)

# new normalization method based on pearson correlations (Lause/Kobak et al. 2021)
# this normalized matrix is given the name 'pearson' using the update_slot param and will be used in the downstream steps
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            scalefactor = 5000,
                            verbose = TRUE,
                            norm_methods = 'pearson_resid',
                            update_slot = 'pearson')

showGiottoExpression(fov_join)

# add statistics based on 'raw' for features rna and negative probes
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'rna')
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'neg_probe')

# View cellular data (default is feat = 'rna')
showGiottoCellMetadata(fov_join)
# View feature data
showGiottoFeatMetadata(fov_join)



## VIEW NORM COUNT DATA
filterDistributions(fov_join,
                    detection = 'cells',
                    feat_type = 'rna',
                    expression_values = 'normalized',
                    method = 'sum',
                    nr_bins = 100,
                    save_param = list(base_height = 3,
                                      save_name = '5.1_rna_norm_total_hist'))

filterDistributions(fov_join,
                    detection = 'cell',
                    feat_type = 'neg_probe',
                    expression_values = 'normalized',
                    method = 'sum',
                    nr_bins = 20,
                    save_param = list(base_height = 3,
                                      save_name = '5.2_neg_norm_total_hist'))

spatPlot2D(gobject = fov_join,
           cell_color = 'total_expr',
           color_as_factor = FALSE,
           show_image = TRUE,
           image_name = image_names,
           point_size = 0.9,
           point_alpha = 0.75,
           save_param = list(base_height = 10,
                             save_name = '5.3_color_centroids'))

spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_color = 'gray',
                     polygon_line_size = 0.05,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE,
                     save_param = list(base_height = 10,
                                       save_name = '5.4_rna_color_polys'))

spatInSituPlotPoints(fov_join,
                     feat_type = 'neg_probe',
                     show_polygon = TRUE,
                     polygon_color = 'gray',
                     polygon_line_size = 0.05,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE,
                     save_param = list(base_height = 10,
                                       save_name = '5.5_neg_color_polys'))


# Dim red

fov_join = calculateHVF(fov_join,
                        method = 'var_p_resid',
                        expression_values = 'pearson',
                        save_param = list(base_height = 5,
                                          save_name = '6.1_pearson_HVF'))

# print HVFs
gene_meta = fDataDT(fov_join)
gene_meta[hvf == 'yes', feat_ID]

fov_join = runPCA(fov_join,
                  scale_unit = FALSE,
                  center = FALSE,
                  expression_values = 'pearson')

# screeplot uses the generated PCA. No need to specify expr values
screePlot(fov_join, ncp = 20, save_param = list(save_name = '6.2_screeplot'))

plotPCA(fov_join,
        cell_color = 'nr_feats', # (from log norm statistics)
        color_as_factor = FALSE,
        point_size = 0.1,
        point_shape = 'no_border',
        save_param = list(save_name = '6.3_PCA'))

fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:10,
                    n_threads = 4)

plotUMAP(gobject = fov_join, save_param = list(save_name = '6.4_UMAP'))


fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)

fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 0.05,
                            n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = fov_join,
         cell_color = 'leiden_clus',
         show_NN_network = T,
         point_size = 2.5)









































### ### ### ### ### ### ### ### ### ### ### ### ### ### 
fov_join = calculateOverlapRaster(fov_join, feat_info = 'neg_probe')
showGiottoSpatialInfo(fov_join)

fov_join = overlapToMatrix(fov_join, feat_info = 'neg_probe')

showGiottoExpression(fov_join)

fov_join <- addStatistics(gobject = fov_join, feat_type = 'neg_probe', expression_values = 'raw') # get cell stats for probes
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# View cellular data
pDataDT(fov_join)
# View rna data
fDataDT(fov_join)

cellmeta = pDataDT(fov_join, feat_type = 'rna')
hist(cellmeta$nr_feats, 100)
cellmeta_neg = pDataDT(fov_join, feat_type = 'neg_probe')
hist(cellmeta_neg$nr_feats, 20)

# filterDistributions(gobject = fov_join,feat_type = 'rna', method = 'mean', nr_bins = 100)






# typical way of calculating HVG
fov_join <- calculateHVF(gobject = fov_join,
                         HVFname = 'hvg_orig')

# new method based on variance of pearson residuals for each gene
fov_join <- calculateHVF(gobject = fov_join,
                         method = 'var_p_resid',
                         expression_values = 'pearson',
                         show_plot = T)

gene_meta = fDataDT(fov_join)
gene_meta[hvf == 'yes']

fov_join <- runPCA(gobject = fov_join,
                   expression_values = 'pearson',
                   scale_unit = F,
                   center = F)
screePlot(fov_join, ncp = 20)

plotPCA(fov_join,
        dim1_to_use = 1,
        dim2_to_use = 2)

fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:10,
                    n_threads = 4)
plotUMAP(gobject = fov_join,
         point_size = 0.1)


## ******
## 
fov_all = createGiottoCosMxObject(
  cosmx_dir = data_path,
  data_to_use = 'subcellular',
  instructions = instrs2,
  cores = 12
)


fov_all = calculateOverlapRaster(fov_all)

# Generate raw matrix
fov_all = overlapToMatrix(fov_all)


fov_all <- filterGiotto(
  gobject = fov_all,
  feat_type = 'rna', # NEW
  expression_threshold = 1,
  feat_det_in_min_cells = 5,
  min_det_feats_per_cell = 50
)

fov_all <- normalizeGiotto(
  gobject = fov_all,
  scalefactor = 5000,
  verbose = T,
  norm_methods = 'pearson_resid',
  update_slot = 'pearson'
)

fov_all = addStatistics(fov_all,
                        feat_type = 'rna',
                        expression_values = 'raw')


fov_all = adjustGiottoMatrix(gobject = fov_all,
                             expression_values = 'pearson',
                             covariate_columns = 'nr_feats', # from raw
                             update_slot = 'pearson_adj')


fov_all = calculateHVF(gobject = fov_all,
                       method = 'var_p_resid',
                       expression_values = 'pearson_adj',
                       show_plot = T)

fov_all = runPCA(gobject = fov_all,
                 expression_values = 'pearson_adj',
                 scale_unit = F,
                 center = F)

# new stats
fov_all = addStatistics(fov_all,
                        expression_values = 'pearson_adj',
                        feat_type = 'rna')

fm = get_feature_metadata(fov_all)

screePlot(fov_all,
          ncp = 20,
          expression_values = 'pearson_adj')

plotPCA(fov_all,
        cell_color = 'nr_feats', # pearson adj valuse
        color_as_factor = F,
        point_shape = 'no_border',
        dim1_to_use = 1,
        dim2_to_use = 2,
        point_size = 0.01)

fov_all = runUMAP(fov_all,
                  dimensions_to_use = 1:10,
                  expression_values = 'pearson_adj',
                  n_threads = 12)

plotUMAP(gobject = fov_all,
         point_size = 0.01)


# feat plots
dimFeatPlot2D(fov_all,
              point_shape = 'no_border',
              point_size = 0.01,
              expression_values = 'pearson_adj',
              feats = c('COL1A1','KRT8','CD8A','MS4A1'))

# no ncam
dimFeatPlot2D(fov_all,
              point_shape = 'no_border',
              point_size = 0.01,
              expression_values = 'pearson_adj',
              feats = c('CD8A','CD4','JCHAIN','HLA-DRA'))


dimPlot2D(fov_all,
          point_shape = 'no_border',
          point_size = 0.01,
          color_as_factor = T,
          cell_color = 'list_ID')



# stats

dimPlot2D(fov_all,
          point_shape = 'no_border',
          point_size = 0.01,
          color_as_factor = F,
          cell_color = 'total_expr')

dimPlot2D(fov_all,
          point_shape = 'no_border',
          point_size = 0.01,
          color_as_factor = F,
          cell_color = 'nr_feats')



# --------------


fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)
fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 0.05,
                            n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = fov_join,
         cell_color = 'leiden_clus',
         show_NN_network = T,
         point_size = 2.5)

fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)
fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 0.07,
                            n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = fov_join,
         cell_color = 'leiden_clus',
         show_NN_network = T,
         point_size = 2.5)

# visualize UMAP and spatial results
spatDimPlot2D(gobject = fov_join,
              show_image = T,
              image_name = image_names,
              cell_color = 'leiden_clus',
              cell_color_code = pal,
              spat_point_size = 2)

spatInSituPlotPoints(fov_join,
                     feats = list('rna' = c("MMP2", "VEGFA", "IGF1R",
                                            'CDH2', 'MKI67', 'EPCAM')),
                     point_size = 0.15,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.01,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = T,
                     polygon_fill_code = pal,
                     coord_fix_ratio = TRUE,
                     save_param = list(base_height = 2))

locs <-fov_join@spatial_locs$cell$raw

#subset a Giotto object based on spatial locations
smallfov <- subsetGiottoLocs(fov_join,
                             x_max = 800,
                             x_min = 507,
                             y_max = -158800,
                             y_min = -159600)

#extract all genes observed in new object
smallfeats <- smallfov@feat_metadata$cell$rna$feat_ID

#plot all genes
spatInSituPlotPoints(smallfov,
                     feats = list(smallfeats),
                     point_size = 0.15,
                     polygon_line_size = .1,
                     show_polygon = T,
                     polygon_color = 'white',
                     show_image = T,
                     image_name = image_names,
                     coord_fix_ratio = TRUE,
                     show_legend = FALSE)

# create spatial network based on physical distance of cell centroids
fov_join = createSpatialNetwork(gobject = fov_join,
                                minimum_k = 2,
                                maximum_distance_delaunay = 50)

# select features
feats = fov_join@feat_ID$rna
# perform Binary Spatial Extraction of genes - NOTE: Depending on your system this could take time
km_spatialgenes = binSpect(fov_join,
                           subset_feats = feats)

# visualize spatial expression of selected genes obtained from binSpect
spatFeatPlot2D(fov_join,
               expression_values = 'scaled',
               feats = km_spatialgenes$feats[1:10],
               cell_color_gradient = c('blue', 'white', 'red'),
               point_shape = 'border',
               point_border_stroke = 0.01,
               show_network = F,
               network_color = 'lightgrey',
               point_size = 1.2,
               cow_n_col = 2)

markers = findMarkers_one_vs_all(gobject = fov_join,
                                 method = 'gini',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_clus',
                                 min_feats = 1,
                                 rank_score = 2)
markers[, head(.SD, 5), by = 'cluster']



# violinplot
topgini_genes = unique(markers[, head(.SD, 2), by = 'cluster']$feats)
violinPlot(fov_join,
           feats = topgini_genes,
           cluster_column = 'leiden_clus',
           strip_position = 'right')

cluster_order = c(1, 2, 3, 4, 5, 6, 7, 8, 9 , 10)
plotMetaDataHeatmap(fov_join,
                    expression_values = 'scaled',
                    metadata_cols = c('leiden_clus'),
                    selected_feats = topgini_genes,
                    custom_cluster_order = cluster_order)

## add cell types ###
clusters_cell_types_lung = c('Normal Epithelial', 'Cancer', 'Stromal', 'Plasma Cells',
                             'Cytotoxic T Cells', 'Cancer Stem Cells',
                             'Macrophage', 'Memory B Cell', 'Memory B Cell')

names(clusters_cell_types_lung) = as.character(sort(cluster_order))
fov_join = annotateGiotto(gobject = fov_join,
                          annotation_vector = clusters_cell_types_lung,
                          cluster_column = 'leiden_clus')

plotUMAP(fov_join,
         cell_color = 'cell_types',
         point_size = 1.5)

spatDimPlot2D(gobject = fov_join,
              show_image = T,
              image_name = image_names,
              cell_color = 'cell_types',
              spat_point_size = 2)

spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE)

future::plan('multisession', workers = 4) # NOTE: Depending on your system this could take time

goi = findInteractionChangedFeats(gobject = fov_join,
                                  cluster_column = 'leiden_clus')

# Identify top ten interaction changed genes
goi$CPGscores[type_int == 'hetero']$feats[1:10]

# Visualize ICG expression
spatInSituPlotPoints(fov_join,
                     feats = list(goi$CPGscores[type_int == 'hetero']$feats[1:10]),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'black',
                     polygon_line_size = 0.1,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal,
                     coord_fix_ratio = TRUE)









rm(fov_join)
# Go to gobject list creation

# ---------------- Gobject List Creation ####

## provide path to nanostring folder
data_path = '/projectnb/rd-spat/DATA/Public_data/Nanostring/Lung12/Lung12-Flat_files_and_images/'

# load transcript coordinates
tx_coord_all = fread(paste0(data_path, 'Lung12_tx_file-002.csv')) #! This is probably the wrong name. Keeping it for now since it's on the drive as this name


# colnames(tx_coord_all)
# cat('\n')
# z planes
# tx_coord_all[, table(z)]
# cat('\n')
# Cell compartment
# tx_coord_all[, table(CellComp)]



all_IDs = tx_coord_all[, unique(target)]
# negative probe IDs
neg_IDs = all_IDs[grepl(pattern = 'NegPrb', all_IDs)]
cat('Negative Probe IDs\n')
neg_IDs
cat('\nFeature IDs\n')
feat_IDs = all_IDs[!all_IDs %in% neg_IDs]
length(feat_IDs)

# split detections
feat_coords_all = tx_coord_all[target %in% feat_IDs]
neg_coords_all = tx_coord_all[target %in% neg_IDs]

cat('\nFeatures: ', feat_coords_all[, .N], '\n',
    'NegProbes: ', neg_coords_all[, .N])



neg_points = createGiottoPoints(
  x = neg_coords_all[, .(target, x_global_px, y_global_px)]
)
plot(neg_points, point_size = 0.2, feats = neg_IDs)




#  load field of vision (fov) positions
fov_offset_file = fread(paste0(data_path, 'Lung12_fov_positions_file.csv'))

gobjects_list = list()

# select which FOV's you would like to work with
# the dataset includes 28, which is too much for most computers to handle at once. 
#For this example I am using 02, 03, and 04
id_set = c('02', '03', '04')

for(fov_i in 1:length(id_set)) {
  
  fov_id = id_set[fov_i]
  
  
  # 1. original composite image as png
  original_composite_image = paste0(data_path, 'CellComposite/CellComposite_F0', fov_id,'.jpg')
  
  # 2. input cell segmentation as mask file
  segmentation_mask = paste0(data_path, 'CellLabels/CellLabels_F0', fov_id, '.tif')
  
  # 3. input features coordinates + offset
  tx_coord = feat_coords_all[fov == as.numeric(fov_id)]
  tx_coord = tx_coord[,.(x_local_px, y_local_px, z, target)]
  colnames(tx_coord) = c('x', 'y', 'z', 'gene_id')
  tx_coord = tx_coord[,.(x, y, gene_id)]
  
  ng_coord = neg_coords_all[fov == as.numeric(fov_id)]
  ng_coord = ng_coord[,.(x_local_px, y_local_px, z, target)]
  colnames(ng_coord) = c('x', 'y', 'z', 'gene_id')
  ng_coord = ng_coord[,.(x, y, gene_id)]
  
  
  fovsubset = createGiottoObjectSubcellular(gpoints = list('rna' = tx_coord, 'neg_probe' = ng_coord),
                                            gpolygons = list('cell' = segmentation_mask),
                                            polygon_mask_list_params = list(mask_method = 'guess',
                                                                            flip_vertical = TRUE,
                                                                            flip_horizontal = FALSE,
                                                                            shift_horizontal_step = FALSE),
                                            instructions = instrs1)
  
  
  # centroids are now used to provide the spatial locations (centroid of each cell)
  fovsubset = addSpatialCentroidLocations(fovsubset,
                                          poly_info = 'cell')
  
  # create and add Giotto images
  composite = createGiottoLargeImage(raster_object = original_composite_image,
                                     negative_y = FALSE,
                                     name = 'composite')
  
  fovsubset = addGiottoImage(gobject = fovsubset,
                             largeImages = list(composite))
  
  
  fovsubset = convertGiottoLargeImageToMG(giottoLargeImage = composite,
                                          #mg_name = 'composite',
                                          gobject = fovsubset,
                                          return_gobject = TRUE)
  
  gobjects_list[[fov_i]] = fovsubset
  
  
}

new_names = paste0("fov0", id_set)

id_match = match(as.numeric(id_set), fov_offset_file$fov)
x_shifts = fov_offset_file[id_match]$x_global_px
y_shifts = fov_offset_file[id_match]$y_global_px

# Create Giotto object that includes all selected FOVs
fov_join = joinGiottoObjects(gobject_list = gobjects_list,
                             gobject_names = new_names,
                             join_method = 'shift',
                             x_shift = x_shifts,
                             y_shift = y_shifts)




# ---------------- Downstream List Creation ####

showGiottoImageNames(fov_join)

# Set up vector of image names
id_set = c('02', '03', '04')
new_names = paste0("fov0", id_set)
image_names = paste0(new_names, '-image')

spatInSituPlotPoints(fov_join,
                     show_image = TRUE,
                     image_name = image_names,
                     feats = list('rna' = c("MMP2", "VEGFA", "IGF1R",
                                            'CDH2', 'MKI67', 'EPCAM')),
                     spat_unit = 'cell',
                     point_size = 0.15,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.02,
                     coord_fix_ratio = TRUE,
                     background_color = NA)

spatPlot2D(gobject = fov_join,
           image_name = image_names,
           show_image = TRUE,
           point_size = 0.2,
           coord_fix_ratio = 1)

fov_join = calculateOverlapRaster(fov_join)

fov_join = overlapToMatrix(fov_join)

showGiottoExpression(fov_join)

# combine cell data
morphometa = combineCellData(fov_join,
                             feat_type = 'rna')

# combine feature data
featmeta = combineFeatureData(fov_join,
                              feat_type = c('rna'))

# combine overlapping feature data
featoverlapmeta = combineFeatureOverlapData(fov_join,
                                            feat_type = c('rna'))

# filter
fov_join <- filterGiotto(gobject = fov_join,
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 5)

# normalize
# standard method
fov_join <- normalizeGiotto(gobject = fov_join,
                            scalefactor = 5000,
                            verbose = T)

# new normalizaton method based on pearson correlations (Lause/Kobak et al. 2021)
# this normalized matrix is given the name 'pearson' and will be used in the downstream steps
fov_join <- normalizeGiotto(gobject = fov_join,
                            scalefactor = 5000,
                            verbose = T,
                            norm_methods = 'pearson_resid',
                            update_slot = 'pearson')
# add statistics
fov_join <- addStatistics(gobject = fov_join)

# View cellular data
pDataDT(fov_join)
# View rna data
fDataDT(fov_join)

cellmeta = pDataDT(fov_join, feat_type = 'rna')
hist(cellmeta$nr_feats, 100)

spatPlot2D(gobject = fov_join,
           cell_color = 'total_expr',
           color_as_factor = F,
           show_image = TRUE,
           image_name = image_names,
           point_size = 1.5,
           point_alpha = 0.75,
           coord_fix_ratio = T)

spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = F,
                     coord_fix_ratio = T)

# typical way of calculating HVG
fov_join <- calculateHVF(gobject = fov_join,
                         HVFname = 'hvg_orig')

# new method based on variance of pearson residuals for each gene
fov_join <- calculateHVF(gobject = fov_join,
                         method = 'var_p_resid',
                         expression_values = 'pearson',
                         show_plot = T)

gene_meta = fDataDT(fov_join)
gene_meta[hvf == 'yes']

fov_join <- runPCA(gobject = fov_join,
                   expression_values = 'pearson',
                   scale_unit = F,
                   center = F)
screePlot(fov_join, ncp = 20)

plotPCA(fov_join,
        dim1_to_use = 1,
        dim2_to_use = 2)

fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:10,
                    n_threads = 4)
plotUMAP(gobject = fov_join)

fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)
fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 0.05,
                            n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = fov_join,
         cell_color = 'leiden_clus',
         show_NN_network = T,
         point_size = 2.5)

fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)
fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 0.05,
                            n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = fov_join,
         cell_color = 'leiden_clus',
         show_NN_network = T,
         point_size = 2.5)

# visualize UMAP and spatial results
spatDimPlot2D(gobject = fov_join,
              show_image = T,
              image_name = image_names,
              cell_color = 'leiden_clus',
              spat_point_size = 2)

spatInSituPlotPoints(fov_join,
                     feats = list('rna' = c("MMP2", "VEGFA", "IGF1R",
                                            'CDH2', 'MKI67', 'EPCAM')),
                     point_size = 0.15,
                     show_polygon = TRUE,
                     polygon_color = 'white',
                     polygon_line_size = 0.01,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = TRUE)

locs <-fov_join@spatial_locs$cell$raw

#subset a Giotto object based on spatial locations
smallfov <- subsetGiottoLocs(fov_join,
                             x_max = 800,
                             x_min = 507,
                             y_max = -158800,
                             y_min = -159600)

#extract all genes observed in new object
smallfeats <- smallfov@feat_metadata$cell$rna$feat_ID

#plot all genes
spatInSituPlotPoints(smallfov,
                     feats = list(smallfeats),
                     point_size = 0.15,
                     polygon_line_size = .1,
                     show_polygon = T,
                     polygon_color = 'white',
                     show_image = T,
                     image_name = image_names,
                     coord_fix_ratio = TRUE,
                     show_legend = FALSE)

# create spatial network based on physical distance of cell centroids
fov_join = createSpatialNetwork(gobject = fov_join,
                                minimum_k = 2,
                                maximum_distance_delaunay = 50)

# select features
feats = fov_join@feat_ID$rna
# perform Binary Spatial Extraction of genes - NOTE: Depending on your system this could take time
km_spatialgenes = binSpect(fov_join,
                           subset_feats = feats)

# visualize spatial expression of selected genes obtained from binSpect
spatFeatPlot2D(fov_join,
               expression_values = 'scaled',
               feats = km_spatialgenes$feats[1:10],
               cell_color_gradient = c('blue', 'white', 'red'),
               point_shape = 'border',
               point_border_stroke = 0.01,
               show_network = F,
               network_color = 'lightgrey',
               point_size = 1.2,
               cow_n_col = 2)

markers = findMarkers_one_vs_all(gobject = fov_join,
                                 method = 'gini',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_clus',
                                 min_feats = 1,
                                 rank_score = 2)
markers[, head(.SD, 5), by = 'cluster']

# violinplot
topgini_genes = unique(markers[, head(.SD, 2), by = 'cluster']$feats)
violinPlot(fov_join,
           feats = topgini_genes,
           cluster_column = 'leiden_clus',
           strip_position = 'right')



### <-<-

dimPlot2D(fov_join, feats = c('KRT8', 'CD4', 'CD8A', 'JCHAIN', ))

### <-<-

cluster_order = c(1, 2, 3, 4, 5, 6, 7, 8, 9)
plotMetaDataHeatmap(fov_join,
                    expression_values = 'scaled',
                    metadata_cols = c('leiden_clus'),
                    selected_feats = topgini_genes,
                    custom_cluster_order = cluster_order)

## add cell types ###
clusters_cell_types_lung = c('Normal Epithelial', 'Cancer', 'Stromal', 'Plasma Cells',
                             'Cytotoxic T Cells', 'Cancer Stem Cells',
                             'Macrophage', 'Memory B Cell', 'Memory B Cell')

names(clusters_cell_types_lung) = as.character(sort(cluster_order))
fov_join = annotateGiotto(gobject = fov_join,
                          annotation_vector = clusters_cell_types_lung,
                          cluster_column = 'leiden_clus')

plotUMAP(fov_join,
         cell_color = 'cell_types',
         point_size = 1.5)

spatDimPlot2D(gobject = fov_join,
              show_image = T,
              image_name = image_names,
              cell_color = 'cell_types',
              spat_point_size = 2)

spatInSituPlotPoints(fov_join,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE)

future::plan('multisession', workers = 4) # NOTE: Depending on your system this could take time

goi = findInteractionChangedFeats(gobject = fov_join,
                                  cluster_column = 'leiden_clus')

# Identify top ten interaction changed genes
goi$CPGscores[type_int == 'hetero']$feats[1:10]

# Visualize ICG expression
spatInSituPlotPoints(fov_join,
                     feats = list(goi$CPGscores[type_int == 'hetero']$feats[1:10]),
                     point_size = 0.15,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'black',
                     polygon_line_size = 0.1,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     polygon_fill_code = pal,
                     coord_fix_ratio = TRUE)













