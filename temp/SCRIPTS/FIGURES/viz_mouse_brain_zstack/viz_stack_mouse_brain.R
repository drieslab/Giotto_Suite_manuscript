
library(Giotto)

# 1. create minivizgen2 dataset from mouse vizgen brain dataset ####
# ---------------------------------------------------------------- #

# see SCRIPTS/create_mini_viz_mousebrain

# 2. load giotto object ####
# ------------------------ #
minivizgen2 = loadGiotto(path_to_folder = 'DATA/mini7/')
instructions(minivizgen2, 'show_plot') = TRUE
instructions(minivizgen2, 'return_plot') = FALSE
instructions(minivizgen2, 'save_plot') = FALSE

spatInSituPlotPoints(minivizgen2,
                     feats = list('rna' = c("Htr1b", "Ackr1", "Epha7")),
                     feats_color_code = c("Htr1b" = 'green', 'Ackr1' = 'blue', 'Epha7' = 'red'),
                     show_image = TRUE,
                     largeImage_name = 'dapi_z0',
                     point_size = 0.7,
                     plot_method = 'ggplot',
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'z1',
                     polygon_color = 'yellow',
                     polygon_bg_color = 'yellow',
                     polygon_line_size = 0.2,
                     coord_fix_ratio = TRUE,
                     background_color = NA)




# 3. aggregate information to matrix: polygons and transcripts ####
# --------------------------------------------------------------- #

# we will use the z1 polygon information
# we can set a global option or specify this for each command
# options('giotto.spat_unit' = 'z1') # now you don't need to think about setting spat_unit each time

z_digits = 0:6
z_levels = paste0('z', z_digits)

for(z_digit in z_digits) {
  
  z_level = paste0('z', z_digit)
  
  cat('level ', z_level, '\n')
  
  minivizgen2 = calculateOverlapRaster(minivizgen2,
                                       spatial_info = z_level,
                                       feat_info = 'rna',
                                       feat_subset_column = 'z',
                                       feat_subset_ids = z_digit)
  
  minivizgen2 = overlapToMatrix(minivizgen2,
                                poly_info = z_level,
                                feat_info = 'rna',
                                name = 'raw')
}

minivizgen2 = aggregateStacks(gobject = minivizgen2,
                              spat_units = z_levels,
                              feat_type = 'rna',
                              values = 'raw',
                              summarize_expression = 'sum',
                              summarize_locations = 'mean',
                              new_spat_unit = 'aggregate')


showGiottoSpatialInfo(minivizgen2)
showGiottoFeatInfo(minivizgen2)


# 4. filter object on aggregated layer #####
# --------------------------------------- ##
minivizgen3 <- filterGiotto(gobject = minivizgen2,
                            spat_unit = 'aggregate',
                            expression_threshold = 1,
                            feat_det_in_min_cells = 3,
                            min_det_feats_per_cell = 25,
                            poly_info = c('aggregate'))

plot(minivizgen2@spatial_info$aggregate)
plot(minivizgen3@spatial_info$aggregate)

activeSpatUnit(gobject = minivizgen3) <- 'aggregate'



spatInSituPlotPoints(minivizgen3,
                     feats = list('rna' = c("Htr1b", "Ackr1", "Epha7")),
                     feats_color_code = c("Htr1b" = 'green', 'Ackr1' = 'blue', 'Epha7' = 'red'),
                     show_image = TRUE,
                     largeImage_name = 'dapi_z0',
                     point_size = 0.7,
                     plot_method = 'ggplot',
                     show_polygon = TRUE,
                     use_overlap = F,
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'yellow',
                     polygon_bg_color = 'yellow',
                     polygon_line_size = 0.2,
                     coord_fix_ratio = TRUE,
                     background_color = NA,
                     save_param = list(base_width = 7, base_height = 7))



# 5. normalize on aggregated layer #####
# ----------------------------------- ##

# rna data, default.
# other feature modalities can be processed and filtered in an anologous manner
minivizgen3 <- normalizeGiotto(gobject = minivizgen3,
                               spat_unit = 'aggregate',
                               scalefactor = 5000,
                               verbose = TRUE)
minivizgen3 <- addStatistics(gobject = minivizgen3, spat_unit = 'aggregate')
minivizgen3 <- normalizeGiotto(gobject = minivizgen3,
                               spat_unit = 'aggregate',
                               norm_methods = 'pearson_resid',
                               update_slot = 'pearson')

spatPlot2D(gobject = minivizgen3,
           spat_unit = 'aggregate',
           cell_color = 'total_expr',
           color_as_factor = F,
           largeImage_name = 'dapi_z0',
           show_image = TRUE,
           point_size = 3.5, 
           point_alpha = 0.5)

spatInSituPlotPoints(minivizgen3,
                     show_image = TRUE,
                     largeImage_name = 'dapi_z0',
                     show_polygon = TRUE,
                     spat_unit = 'aggregate',
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE)


# 6. highly variable genes ####
# ----------------------------- #

# typical way of calculating HVG
minivizgen3 <- calculateHVF(gobject = minivizgen3, spat_unit = 'aggregate', HVFname = 'hvg_orig')

# new method based on variance of pearson residuals for each gene
minivizgen3 <- calculateHVF(gobject = minivizgen3,
                            spat_unit = 'aggregate',
                            method = 'var_p_resid',
                            expression_values = 'pearson',
                            show_plot = T)



# 7. dimension reduction ####
# --------------------------- #

# ** 7.1 PCA ####

# we will run pca on the pre-scaled matrix from the pearson residual normalization
# if features are not specified it will automatically search for the hvf column in the feature metadata

?runPCA

showGiottoExpression(minivizgen3)

minivizgen3 <- runPCA(gobject = minivizgen3,
                      feats_to_use = NULL,
                      spat_unit = 'aggregate',
                      expression_values = 'normalized',
                      scale_unit = TRUE,
                      center = TRUE)

screePlot(minivizgen3,
          ncp = 20,
          spat_unit = 'aggregate')

showGiottoDimRed(minivizgen3)

plotPCA(minivizgen3,
        spat_unit = 'aggregate',
        dim_reduction_name = 'pca',
        dim1_to_use = 1,
        dim2_to_use = 2)


# ** 7.2 UMAP and TSN ####
minivizgen3 <- runUMAP(minivizgen3, dimensions_to_use = 1:12, n_threads = 4, spat_unit = 'aggregate')
plotUMAP(gobject = minivizgen3, spat_unit = 'aggregate')

minivizgen3 <- runtSNE(minivizgen3, dimensions_to_use = 1:12, spat_unit = 'aggregate')
plotTSNE(gobject = minivizgen3, spat_unit = 'aggregate')


# 8. graph-based clustering ####
# ---------------------------- #
minivizgen3 <- createNearestNetwork(gobject = minivizgen3, dimensions_to_use = 1:12, k = 10,
                                    spat_unit = 'aggregate')
minivizgen3 <- doLeidenClusterIgraph(gobject = minivizgen3,
                                     resolution = 0.05,
                                     n_iterations = 2000,
                                     spat_unit = 'aggregate')

pDataDT(minivizgen3)


# visualize UMAP / TSNE cluster results
plotUMAP(gobject = minivizgen3,
         spat_unit = 'aggregate',
         cell_color = 'leiden_clus',
         show_NN_network = T,
         point_size = 2.5)

spatInSituPlotPoints(minivizgen3,
                     show_polygon = TRUE,
                     spat_unit = 'aggregate',
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = T)

spatInSituPlotPoints(minivizgen3,
                     spat_unit = 'aggregate',
                     show_image = T,
                     largeImage_name = 'dapi_z0',
                     feats = list('rna' = c("Htr1b", "Ackr1", "Epha7")),
                     feats_color_code = c("Htr1b" = 'green', 'Ackr1' = 'blue', 'Epha7' = 'red'),
                     point_size = 0.75,
                     show_polygon = TRUE,
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = TRUE)


# 9. spatial network ####
# --------------------- #

# defaults to delaunay
minivizgen3 = createSpatialNetwork(minivizgen3, spat_unit = 'aggregate')
# kNN with k of 8 nearest
minivizgen3 = createSpatialNetwork(minivizgen3, spat_unit = 'aggregate', method = 'kNN', k = 8)
# create spatial weight matrix
minivizgen3 = createSpatialWeightMatrix(minivizgen3,
                                        spat_unit = 'aggregate',
                                        method = 'distance',
                                        spatial_network_to_use = 'kNN_network',
                                        return_gobject = TRUE)

pDataDT(minivizgen3, 'aggregate')

spatPlot(gobject = minivizgen3, spat_unit = 'aggregate', show_network = T,
         network_color = 'lightgray', spatial_network_name = 'Delaunay_network',
         point_size = 2.5, cell_color = 'leiden_clus')


## 9.1 spatial genes ####
km_spatialfeats = binSpect(minivizgen3, spat_unit = 'aggregate')

spatFeatPlot2D_single(minivizgen3,
                      spat_unit = 'aggregate',
                      expression_values = 'scaled',
                      feats = c('Slc17a6', 'Gpr182', 'Myh11', 'Gfap'),
                      point_shape = 'border',
                      point_border_stroke = 0.1,
                      show_network = F, network_color = 'lightgrey', point_size = 2.5,
                      cow_n_col = 2)


## 9.2 spatial co-expression ####
# here we use existing detectSpatialCorGenes function to calculate pairwise distances between genes
ext_spatial_genes = km_spatialfeats[1:200,]$feats
spat_cor_netw_DT = detectSpatialCorFeats(minivizgen3,
                                         spat_unit = 'aggregate',
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network',
                                         subset_feats = ext_spatial_genes)

# cluster and visualize spatial co-expression genes
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, name = 'spat_netw_clus', k = 7)

heatmSpatialCorFeats(minivizgen3,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus',
                     heatmap_legend_param = list(title = NULL),
                     save_param = list(base_height = 6, base_width = 8, units = 'cm'))

# create and visualize metafeatures
testweight = getBalancedSpatCoexpressionFeats(spat_cor_netw_DT, rank = 'weighted', maximum = 20)

minivizgen3 = createMetafeats(minivizgen3,
                              spat_unit = 'aggregate',
                              feat_clusters = testweight,
                              name = 'cluster_metagene')

spatCellPlot(minivizgen3,
             spat_unit = 'aggregate',
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = as.character(c(1:6)),
             point_size = 2, cow_n_col = 3, save_param = list(base_width = 15))



## 9.3 spatial structure ####
cell_proximities = cellProximityEnrichment(gobject = minivizgen3,
                                           spat_unit = 'aggregate',
                                           cluster_column = 'leiden_clus',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)
## barplot
cellProximityBarplot(gobject = minivizgen3,
                     CPscore = cell_proximities,
                     min_orig_ints = 5, min_sim_ints = 5)

## heatmap
cellProximityHeatmap(gobject = minivizgen3,
                     CPscore = cell_proximities,
                     order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5),
                     color_names = c('blue', 'white', 'red'))


# overlaps vs no overlaps
spatInSituPlotPoints(minivizgen3,
                     spat_unit = 'aggregate',
                     show_image = T,
                     largeImage_name = 'dapi_z0',
                     feats = list('rna' = ext_spatial_genes[1:50]),
                     show_legend = F,
                     point_size = 0.75,
                     show_polygon = TRUE,
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = 1)



spatInSituPlotPoints(minivizgen3,
                     spat_unit = 'aggregate',
                     show_image = T,
                     largeImage_name = 'dapi_z0',
                     feats = list('rna' = ext_spatial_genes[1:50]),
                     show_legend = F,
                     point_size = 0.3,
                     show_polygon = TRUE,
                     use_overlap = F,
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.2,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = 1,
                     plot_last = 'polygons')

showGiottoSpatialInfo(minivizgen3)


spatInSituPlotDensity(minivizgen3,
                      polygon_feat_type = 'aggregate',
                      feats = c("Htr1b", "Ackr1", "Epha7"),
                      alpha = 0.5, polygon_color = 'white')

spatInSituPlotHex(minivizgen3,
                  polygon_feat_type = 'aggregate',
                  feats = c("Htr1b", "Ackr1", "Epha7"))




# niche

pDataDT(minivizgen3)
showGiottoSpatNetworks(minivizgen3)


minivizgen3 = calculateSpatCellMetadataProportions(gobject = minivizgen3,
                                                   spat_unit = 'aggregate',
                                                   feat_type = 'rna',
                                                   metadata_column = 'leiden_clus',
                                                   spat_network = 'Delaunay_network')



prop_table = getSpatialEnrichment(minivizgen3, name = 'proportion', output = 'data.table')
prop_matrix = dt_to_matrix(prop_table)

prop_kmeans = kmeans(x = prop_matrix, centers = 8, iter.max = 100, nstart = 3)
prop_kmeansDT = data.table::data.table(cell_ID = names(prop_kmeans$cluster), niche = prop_kmeans$cluster)


minivizgen3 = addCellMetadata(minivizgen3, new_metadata = prop_kmeansDT, by_column = T, column_cell_ID = 'cell_ID')

spatPlot(gobject = minivizgen3, show_network = T,
         network_color = 'lightgray', spatial_network_name = 'Delaunay_network',
         point_size = 2.5, cell_color = 'niche')


cell_meta = pDataDT(minivizgen3)
source_IDs = cell_meta[niche == 3][['cell_ID']]
neighbor_cells = findNetworkNeighbors(minivizgen3,
                                      source_cell_ids = source_IDs, spatial_network_name = 'Delaunay_network')

minivizgen3 = addCellMetadata(minivizgen3, new_metadata = neighbor_cells, by_column = T, column_cell_ID = 'cell_ID')

spatPlot2D(gobject = minivizgen3, show_network = T,
           network_color = 'lightgray', spatial_network_name = 'Delaunay_network',
           point_size = 2.5, cell_color = 'nb_cells',
           cell_color_code = c(source = 'darkred', both = 'red', neighbor = 'lightblue', others = '#ececec'))



spatPlot(gobject = minivizgen3, show_network = T,
         network_color = 'lightgray', spatial_network_name = 'Delaunay_network',
         point_size = 2.5, cell_color = 'leiden_clus')


spatInSituPlotPoints(minivizgen3,
                     show_polygon = TRUE,
                     spat_unit = 'aggregate',
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'nb_cells',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = T)


# HMRF

# do HMRF with different betas on top 30 genes per spatial co-expression module
hmrf_folder = paste0('SCRIPTS/FIGURES/viz_stack_mouse_brain/','HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = TRUE)

showGiottoSpatNetworks(minivizgen3)


HMRF_spatial_genes = doHMRF(gobject = minivizgen3,
                            expression_values = 'scaled',
                            spatial_genes = names(testweight),
                            k = 3,
                            spatial_network_name = "kNN_network",
                            betas = c(0, 10, 5),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k20_scaled'))

minivizgen3 = addHMRF(gobject = minivizgen3,
                      HMRFoutput = HMRF_spatial_genes,
                      k = 3, betas_to_add = c(0, 10, 20, 30, 40),
                      hmrf_name = 'HMRF')

pDataDT(minivizgen3)

spatPlot2D(gobject = minivizgen3, cell_color = 'HMRF_k3_b.40')


spatInSituPlotPoints(minivizgen3,
                     show_polygon = TRUE,
                     spat_unit = 'aggregate',
                     polygon_feat_type = 'aggregate',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'HMRF_k3_b.40',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = T)



# 10. save Giotto object ####
# ------------------------- #













# image with all data
giottoLargeImage = minivizgen2@largeImages$dapi_z0
spatvector_poly = minivizgen2@spatial_info$z0@spatVector
spatvector_poly_centroids = minivizgen2@spatial_info$z0@spatVectorCentroids
spatvector_points = minivizgen2@feat_info$rna@spatVector

plot(giottoLargeImage)
terra::lines(spatvector_poly, col = 'white')
terra::points(spatvector_poly_centroids, col = 'yellow', cex = 0.6)
terra::points(spatvector_points, col = 'blue', cex = 0.2, alpha = 0.4)





# polygon differences #

# get polys
poly_list = minivizgen2@spatial_info[z_levels]

# set up raster
r = terra::rast(nrows = 1e3, ncols = 1e3, extent = ext(poly_list[[1]]))

poly_rasters = lapply(poly_list, function(p) {
  terra::rasterize(p@spatVector, r)
})

poly_rasterstack = Reduce('c', poly_rasters)
poly_rasterstack[is.na(poly_rasterstack)] = 0
poly_rs_sum = terra::app(poly_rasterstack, sum)

vcolors = viridisLite::turbo(8)
vcoltab = data.table::data.table(sum = 0:7, color = vcolors)
terra::coltab(poly_rs_sum) = vcoltab

plot(poly_rs_sum) # plot of multiple layers






## approach ##

# decide on the image resolution
featinfo = getFeatureInfo(minivizgen2, return_giottoPoints = F)
rasterimage = terra::rast(x = featinfo, nrows = 20, ncols = 20)
# this rasterizes the feature points into a 20x20 image
# The resolution is arbitrary, but requires that the resolution is high enough
# that it makes some visual sense, while also being low enough that the
# rasterization will aggregate counts

# create z specific feature information layers
z_layers = sort(unique(featinfo$z))
z_featinfo_list = lapply(z_layers, FUN = function(x) {
  layerinfo = featinfo[featinfo$z == x]
})

# rasterize each gene in all z-layers
genes = unique(featinfo$feat_ID)
rast_list = lapply(1:length(z_featinfo_list), FUN = function(i) {
  
  cat('for layer ', i, '\n')
  
  z_layer = z_featinfo_list[[i]]
  
  final_z_collection = list()
  
  for(gene_i in 1:length(genes)) {
    print(gene_i)
    gene = genes[gene_i]
    genevec = z_layer[z_layer$feat_ID == gene]
    generast = terra::rasterize(x = genevec, y = rasterimage, fun = 'sum', na.rm = TRUE)
    final_z_collection[[gene_i]] = generast
  }
  
  final_z_collection
  
})

# combine the genes together
new_savelist = lapply(1:length(genes), FUN = function(x) {
  
  print(x)
  
  # aggregate each gene across all rasterized layers
  comb_generast = NULL
  for(i_layer in 1:length(rast_list)) {
    generast = rast_list[[i_layer]][[x]]
    comb_generast = c(comb_generast, generast)
  }
  
  comb_generast = do.call('c', comb_generast)
  
})


# new_savelist
# each gene is converted into a spatraster with defined dimensions (rastermimage)
# each spatraster has 7 layers (each z-stack)

# example:
one_gene = new_savelist[[2]]
plot(one_gene)

?getColors
GiottoVisuals::showColorInstructions()
GiottoVisuals::getColors(n = 10)


plot_layers_gene = function(rasterized_feats,
                            feats_order = NULL,
                            feat = 'Gfap',
                            pal = 'Sunset',
                            n = 100,
                            rev = F) {
  
  
  feat_index = which(feats_order == feat)
  spatRasterGene = rasterized_feats[[feat_index]]
  
  min_value = min(unlist(lapply(1:terra::nlyr(spatRasterGene), FUN = function(x) {
    terra::minmax(spatRasterGene[[x]])[[1]]
  })))
  
  max_value = max(unlist(lapply(1:terra::nlyr(spatRasterGene), FUN = function(x) {
    terra::minmax(spatRasterGene[[x]])[[2]]
  })))
  
  mycolors = GiottoVisuals::getColors(pal = pal, n = n, rev = rev)
  plot(spatRasterGene, range = c(min_value, max_value), col = mycolors, main = feat)
  
}


plot_layers_gene(new_savelist, genes, 'Gfap')

plot_layers_gene(new_savelist, genes, 'Htr1a')

plot_layers_gene(new_savelist, genes, 'Cxcl12')



## example for all genes


# calculate the total sum
sum_spat_genes = lapply(X = new_savelist, function(x) {
  sum_gene = sum(x, na.rm = T)
})
sum_spat_genes_vector = do.call('c', sum_spat_genes)
sum_all = sum(sum_spat_genes_vector, na.rm = T)

mycolors = GiottoVisuals::getColors(pal = 'Sunset', n = 200, rev = T)

# set missing values to 0
sum_all[is.na(sum_all)] = 0
sum_range = terra::minmax(sum_all)



cell_poly = minivizgen2@spatial_info$z0
cell_centroids = minivizgen2@spatial_info$z0@spatVectorCentroids

# plot the total sum (figure 6C)
plot(sum_all, col = mycolors, range = sum_range, main = 'all layers')
plot(cell_poly, add = T)
plot(cell_centroids, add = T, cex = 0.5)



# plot the total sum for each layer (figure 6C continued)
vals = list()
for(i in 1:7) {
  sum_spat_genes_layer = lapply(X = new_savelist, function(x) {
    selected_layer = terra::subset(x = x, i)
  })
  sum_spat_genes_vector_layer = do.call('c', sum_spat_genes_layer)
  sum_spat_genes_vector_layer[is.na(sum_spat_genes_vector_layer)] = 0
  sum_all_layer = sum(sum_spat_genes_vector_layer)
  vals[[i]] = sum_all_layer |>
    terra::values() |>
    as.numeric()
  
  # plot
  plot(sum_all_layer, col = mycolors, main = paste0('layer = ', i))
  plot(cell_poly, add = T)
  plot(cell_centroids, add = T, cex = 0.5)
  
}

# plot as boxplot (figure 6E)
vals_DF = Reduce(cbind, vals)
colnames(vals_DF) = z_levels
boxplot(vals_DF)






# COV
cov_spat_genes = lapply(X = new_savelist, function(x) {
  std_gene = terra::stdev(x, na.rm = T)
  mean_gene = terra::mean(x, na.rm = T)
  cov_gene = std_gene/mean_gene
  cov_gene
})

cov_spat_genes_vector = do.call('c', cov_spat_genes)
cov_spat_genes_vector[is.na(cov_spat_genes_vector)] = 0
mean_cov_all = terra::mean(cov_spat_genes_vector, na.rm = T)


mycolors = GiottoVisuals::getColors(pal = 'Temps', n = 10, rev = F)

# plot mean cov
plot(mean_cov_all, col = mycolors, main = 'mean COV')
plot(minivizgen2@spatial_info$z0@spatVector , add = T)
plot(minivizgen2@spatial_info$z0@spatVectorCentroids, add = T, cex = 0.5)



sum_cov_all = sum(cov_spat_genes_vector, na.rm = T)

mycolors = GiottoVisuals::getColors(pal = 'Temps', n = 10, rev = F)

# plot total cov (figure 6F)
plot(sum_cov_all, col = mycolors, main = 'total COV')
plot(minivizgen2@spatial_info$z0@spatVector , add = T)
plot(minivizgen2@spatial_info$z0@spatVectorCentroids, add = T, cex = 0.5)






# spatial_diff = data.table(genes = genes,
#                           sum = sum(terra::values(sum_spat_genes[[1]], na.rm = T)),
#                           sd = sum(terra::values(sd_savelist[[1]], na.rm = T)))


spatial_diff = data.table::data.table(genes = genes, sum = unlist(sum_spat_genes), cov = unlist(cov_spat_genes))
spatial_diff[, sumt := sum(terra::values(sum[[1]], na.rm = T)), by = 1:nrow(spatial_diff)]
spatial_diff[, covt := sum(terra::values(cov[[1]], na.rm = T)), by = 1:nrow(spatial_diff)]

data.table::setorder(spatial_diff, -covt)
spatial_diff[1:50]

spatial_diff[sumt > 1500]

# (figure 6G)
library(ggplot2)
pl = ggplot()
pl = pl + geom_point(data = spatial_diff, aes(x = log(sumt+1), y = log(covt+1)))
pl = pl + theme_bw() + xlab('log(counts +1)') + ylab('log(COV+1)')
pl
# with Gad1 label
pl2 = pl + geom_point(data = spatial_diff[genes == 'Gad1'], aes(x = log(sumt+1), y = log(covt+1), color = 'red')) +
  geom_text(data = spatial_diff[genes == 'Gad1'], aes(x = log(sumt+1), y = log(covt+1), color = 'red'), label="Gad1", hjust=-0.2)



plot_layers_gene(new_savelist, genes, 'Cspg5')

# example of gene with differences in distribution across the z slices (figure 6H)
plot_layers_gene(new_savelist, genes, 'Gad1')





# gfap example
names(genes) = 1:length(genes)

selected_gene = 'Cspg5'
gene_index = as.numeric(names(genes[genes == selected_gene]))


high_cov_gene = 'Cspg5'
gene_index_high_cov = which(genes == high_cov_gene)


low_cov_gene = 'Gad1'
gene_index_low_cov = which(genes == low_cov_gene)

mycolors = GiottoVisuals::getColors(pal = 'Temps', n = 10, rev = F)
plot(c(cov_spat_genes[[gene_index_high_cov]], cov_spat_genes[[gene_index_low_cov]]), col = mycolors)
