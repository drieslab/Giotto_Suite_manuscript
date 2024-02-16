
library(Giotto)
library(GiottoDB)

# load DB materials
bInfo <- readRDS('DATA/stereoseq/backend_info.rds')
b <- reconnectBackend(bInfo)

s <- readRDS('DATA/stereoseq/ss_points.rds')
s <- initialize(s) # fix extent pointer


# get extent
e <- ext(s)

# radius is currently defined strangely (y top to y bottom of a single hex)
# single cells are expected to contain roughly 400 HDMIs
hexbin1 <- tessellate(e, shape = 'hexagon', radius = 200, name = 'hex200')
hexbin2 <- tessellate(e, shape = 'hexagon', radius = 50, name = 'hex50')


h1 <- dbvect(
  x = hexbin1[],
  db = b,
  remote_name = 'hex200',
  overwrite = TRUE
)
h2 <- dbvect(
  x = hexbin2[],
  db = b,
  remote_name = 'hex50',
  overwrite = TRUE
)


# 23 min, during the overlap calculations, it never went above 35GB.
# Then at the very end, it froze for a sec and then usage jumped to 53GB
# probably the finalization step where it checks for multiple overlaps of features
h1_overlaps = calculateOverlap(
  h1, s,
  n_per_chunk = 200, # 200 / 4290 geometries run per chunk)
  count_info_column = 'count',
  overwrite = TRUE
)

expected_cn <- unique(h1$poly_ID)
expected_rn <- s[] %>% 
  dplyr::distinct(feat_ID) %>% 
  dplyr::pull()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# REGENERATING


# fine resolution


# 1294.41 sec elapsed 21 min
tictoc::tic()
h2_overlaps = calculateOverlap(
  h2, s,
  n_per_chunk = 5000, # 5000 / 67595 geometries run per chunk
  count_info_column = 'count',
  remote_name = 'hex50_overlaps',
  overwrite = TRUE
)
tictoc::toc()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# 546.814 sec elapsed
tictoc::tic()
h1_M <- overlapToMatrix(
  h1_overlaps,
  col_names = expected_cn,
  row_names = expected_rn,
  count_info_column = 'count',
  write_step = c(1e5,1e3),
  filepath = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/hex200_rna_raw_matrix.h5',
  remote_name = "hex200_rna_raw",
  overwrite = TRUE)
tictoc::toc()

# fine resolution
# 1294.41 sec elapsed 21 min
tictoc::tic()
h2_overlaps = calculateOverlap(
  h2, s,
  n_per_chunk = 5000, # 5000 / 67595 geometries run per chunk
  count_info_column = 'count',
  remote_name = 'hex50_overlaps',
  overwrite = TRUE
)
tictoc::toc()

expected_cn <- unique(h2$poly_ID)
tictoc::tic()
h2_M <- overlapToMatrix(
  h2_overlaps,
  col_names = expected_cn,
  row_names = expected_rn,
  write_step = c(1e5,1e3),
  filepath = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/hex50_rna_raw_matrix.h5',
  remote_name = "hex50_rna_raw",
  overwrite = TRUE)
tictoc::toc()

x <- '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/hex50_rna_raw_matrix.h5'
h2_M <- HDF5Array::HDF5Array(x, name = 'hex50_rna_raw')

g <- addCellStatistics(g, spat_unit = 'hex200', expression_values = 'raw')
spatInSituPlotPoints(g, spat_unit = 'hex200', polygon_fill = 'total_expr',
                     polygon_line_size = 0, polygon_alpha = 1)
g <- addCellStatistics(g, expression_values = 'raw')
spatInSituPlotPoints(g, spat_unit = 'hex50', polygon_fill = 'total_expr',
                     polygon_alpha = 1, polygon_line_size = 0, polygon_fill_gradient_style = 's')








# AFTER RELOAD
### ######################################################################## ###

library(Giotto)
library(GiottoDB)


database_path <- '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/'

b = createBackend(
  drv = duckdb::duckdb(),
  dbdir = database_path
)

p <- getBackendPool(b)
rna_tbl <- tableBE(cPool = p, remote_name = 'rna')
s <- dbPointsProxy(data = rna_tbl)
e <- ext(s)



h1_M_filepath <- '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/hex200_rna_raw_matrix.h5'
h2_M_filepath <- '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/hex50_rna_raw_matrix.h5'

file.size(h1_M_filepath) %>%
  utils:::format.object_size("auto") # 32.5 Mb

file.size(h2_M_filepath) %>%
  utils:::format.object_size("auto") # 161 Mb on disk, 7.61 GB in mem as dense matrix, 1.53GB as dgCMatrix


h1_M <- HDF5Array::HDF5Array(h1_M_filepath, name = 'hex200_rna_raw')
h2_M <- HDF5Array::HDF5Array(h2_M_filepath, name = 'hex50_rna_raw')

exp1 <- createExprObj(
  h1_M,
  name = 'raw',
  spat_unit = 'hex200',
  feat_type = 'rna'
)

exp2 <- createExprObj(
  h2_M,
  name = 'raw',
  spat_unit = 'hex50',
  feat_type = 'rna'
)

hexbin1 = GiottoClass::tessellate(e, shape = 'hexagon', radius = 200, name = 'hex200')
hexbin2 = GiottoClass::tessellate(e, shape = 'hexagon', radius = 50, name = 'hex50')

g <- createGiottoObject(
  expression = list(exp1, exp2),
  spatial_info = list(hexbin1, hexbin2)
)

instructions(g, 'save_dir') <- '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/SCRIPTS/FIGURES/stereoseq_binning_deconv/plots/'
instructions(g, 'show_plot') <- FALSE
instructions(g, 'save_plot') <- TRUE
instructions(g, 'return_plot') <- FALSE

g <- addSpatialCentroidLocations(g, poly_info = list('hex50', 'hex200'))



# filtering

# roughly 3 - 4 GB (with just the HDF5matrices, no DB objects)
filterCombinations(
  gobject = g,
  spat_unit = 'hex50',
  expression_thresholds = c(1, 2),
  feat_det_in_min_cells = c(10, 100, 100, 100),
  min_det_feats_per_cell = c(100, 100, 300, 1000)
)

g <- filterGiotto(
  gobject = g,
  spat_unit = 'hex200',
  expression_threshold = 1,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 300,
  spat_unit_fsub = 'hex200' # limit the feats subsets to this spat unit
)

# Number of cells removed:  1039  out of  4290 
# Number of feats removed:  4158  out of  28133 

g <- filterGiotto(
  gobject = g,
  spat_unit = 'hex50',
  expression_threshold = 1,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 300,
  spat_unit_fsub = 'hex50' # limit the feats subsets to this spat unit
)

# Number of cells removed:  18644  out of  67595 
# Number of feats removed:  4138  out of  28133 

# 1.5 min
tictoc::tic()
g <- normalizeGiotto(
  gobject = g,
  spat_unit = 'hex50'
)
tictoc::toc()

# DelayedArray: 15.2 min elapsed
# dgeMatrix:    2.954 sec elapsed
tictoc::tic()
g <- addCellStatistics(g, spat_unit = 'hex50')
tictoc::toc()

spatInSituPlotPoints(g,
                     spat_unit = 'hex50',
                     polygon_fill = 'total_expr',
                     polygon_alpha = 1,
                     polygon_line_size = 0,
                     polygon_fill_gradient_style = 's',
                     save_param = list(
                       base_width = 10,
                       save_name = '01_filtered_total_expr'
                     ))

# calculate HVF on a subset so that it takes less time.
# Results should be comparable
set.seed(12345)
hex_ids_all <- spatIDs(g, 'hex50') # 48951
hex_ids_10p <- sample(hex_ids_all, length(hex_ids_all)/10) # 4895

# highly variable feats

# 28.268 sec elapsed
tictoc::tic()
g <- calculateHVF(
  gobject = g,
  spat_unit = 'hex50',
  cell_ids = hex_ids_10p
)
tictoc::toc()

# Use projection as a strategy for dealing with larger expression data.

# First define a number of expression values to use in the projection
n_25_percent = round(length(spatIDs(g, 'hex50')) * 0.25) # 25%

# The projection functions randomly sample cell IDs to use in the projection,
# so for the sake of this demonstration, we will set a seed
set.seed(12345)

tictoc::tic()
g <- runPCAprojection(
  gobject = g,
  spat_unit = 'hex50',
  expression_values = 'scaled',
  feats_to_use = 'hvf',
  name = 'pca.projection',
  random_subset = n_25_percent
)
tictoc::toc()
# DelayedArray:          did not finish after 5hrs or so.
# dgeMatrix:             5.4 min
# dgeMatrix 25% + proj:  49.615 sec elapsed





# Visualize Screeplot and PCA

screePlot(
  g,
  name = 'pca.projection',
  expression_values = 'normalized',
  spat_unit = 'hex50',
  ncp = 50,
  save_param = list(
    save_name = 'ss_hex50_screePlot')
)

plotPCA(
  g,
  spat_unit = 'hex50',
  dim_reduction_name = 'pca.projection',
  dim1_to_use = 1,
  dim2_to_use = 2,
  save_param = list(
    base_width = 10,
    save_name = 'ss_hex50_pca'
  )
)


# run the UMAP in a projected manner as well

# 25.003 sec elapsed
tictoc::tic()
# Worked on Matrix 1.6.3
g = runUMAPprojection(
  gobject = g,
  spat_unit = "hex50",
  dimensions_to_use = 1:25,
  name = "umap.projection",
  random_subset = n_25_percent
)
tictoc::toc()

# 33.633 sec elapsed
tictoc::tic()
g = createNearestNetwork(
  gobject = g,
  spat_unit = 'hex50',
  dim_reduction_to_use = 'pca',
  dim_reduction_name = 'pca.projection',
  dimensions_to_use = 1:25,
  k = 25,
  type = 'sNN',
  name = 'sNN.pca.projection'
)
tictoc::toc()



# We will again be using a projection approach for leiden clustering in order
# to speed things up.
# This is done by taking a subset of the total dataset, performing the
# clustering on the subset, then projecting those results to the full dataset.

set.seed(1234)
subset_IDs = sample(x = spatIDs(g, 'hex50'), size = n_25_percent)
subset_g = subsetGiotto(
  gobject = g,
  spat_unit = 'hex50',
  cell_ids = subset_IDs
)

# remake nn network since this was a random sampling, so connections are
# expected to be broken.
subset_g = createNearestNetwork(
  gobject = subset_g,
  spat_unit = 'hex50',
  dim_reduction_to_use = 'pca',
  dim_reduction_name = 'pca.projection',
  dimensions_to_use = 1:25,
  k = 25,
  type = 'sNN',
  name = 'sNN.pca.projection'
)

subset_g = doLeidenClusterIgraph(
  gobject = subset_g,
  spat_unit = "hex50",
  resolution_parameter = 0.55,
  name = "sub_leiden_clus",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca.projection",
) 

# preview clustering of subset
dimPlot2D(
  gobject = subset_g,
  spat_unit = 'hex50',
  dim_reduction_to_use = 'umap',
  dim_reduction_name = 'umap.projection',
  cell_color = 'sub_leiden_clus',
  save_param = list(
    save_name = 'ss_hex50_subset_leiden',
    base_width = 10
  )
)

# project clusterings back to full dataset
# 15.13 seconds
tictoc::tic()
g <- doClusterProjection(
  target_gobject = g,
  source_gobject = subset_g,
  spat_unit = "hex50",
  source_cluster_labels = "sub_leiden_clus",
  reduction_method = 'pca',
  reduction_name = 'pca.projection',
  prob = FALSE,
  knn_k = 25, 
  dimensions_to_use = 1:25
)
tictoc::toc()



# plot results
spatInSituPlotPoints(
  gobject = g,
  spat_unit = 'hex50',
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'knn_labels', # leiden clusterings projected by kNN
  polygon_alpha = 1,
  polygon_line_size = 0,
  save_param = list(
    base_width = 10,
    save_name = 'ss_hex50_leiden_clus'
  )
)








# pseudovisium

# create subset extent
brain_ext <- ext(1000, 3000, 9000, 12000)
visium_polygons = makePseudoVisium(
  extent = brain_ext,
  micron_size = 0.5, 
  name = "pseudo_vis"
)
plot(visium_polygons)

# vp <- dbvect(
#   x = visium_polygons[],
#   db = b,
#   remote_name = 'pseudo_vis',
#   overwrite = TRUE
# )

s_sub <- crop(s, brain_ext)


# These steps can be done in memory
s_sub_sv <- as.spatvector(s_sub)

# find overlapped features
vp_overlaps = calculateOverlap(
  visium_polygons[], s_sub_sv,
  count_info_column = 'count',
)

expected_cn <- unique(visium_polygons$poly_ID)
expected_rn <- s[] %>% 
  dplyr::distinct(feat_ID) %>% 
  dplyr::pull()

vp_M <- overlapToMatrix(
  vp_overlaps,
  col_names = expected_cn,
  row_names = expected_rn,
)

exp_vp <- createExprObj(
  expression_data = vp_M,
  spat_unit = 'pseudo_vis',
  feat_type = 'rna',
  provenance = 'pseudo_vis',
  name = 'raw'
)

g <- setExpression(g, exp_vp)
g <- setPolygonInfo(g, visium_polygons)




g <- addSpatialCentroidLocations(g, poly_info = 'pseudo_vis')

activeSpatUnit(g) <- 'pseudo_vis'

g <- filterGiotto(
  gobject = g,
  expression_threshold = 1,
  feat_det_in_min_cells = 1,
  min_det_feats_per_cell = 10,
  spat_unit_fsub = 'pseudo_vis'
)

g <- normalizeGiotto(
  gobject = g
)

# added NA values
g <- addStatistics(g, spat_unit = 'pseudo_vis')

# plot initial results
spatInSituPlotPoints(
  gobject = g,
  spat_unit = 'pseudo_vis',
  polygon_fill_as_factor = FALSE,
  polygon_fill = 'total_expr', # leiden clusterings projected by kNN
  polygon_alpha = 1,
  polygon_fill_gradient_style = 's',
  polygon_line_size = 0,
  save_param = list(
    base_width = 10,
    save_name = 'ss_pseudovis_total_expr'
  )
)


g <- calculateHVF(g)
g <- runPCA(g)
screePlot(g,
          save_param = list(
            save_name = 'pseudo_vis_scree'
          ))
g <- runUMAP(
  g,
  dimensions_to_use = 1:50,
  n_neighbors = 4
)

dimPlot2D(
  g,
  save_param = list(
    save_name = "pseudo_vis_umap"
  )
)

g <- createNearestNetwork(
  gobject = g,
  dimensions_to_use = 1:50,
  k = 4
)

g <- doLeidenCluster(
  gobject = g,
  resolution = 0.6
)

dimPlot2D(
  g,
  cell_color = 'leiden_clus',
  save_param = list(
    save_name = "pseudo_vis_leiden"
  ),
  point_size = 2
)

spatInSituPlotPoints(
  gobject = g,
  spat_unit = 'pseudo_vis',
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'leiden_clus', #
  polygon_alpha = 1,
  polygon_line_size = 0,
  save_param = list(
    base_width = 10,
    save_name = 'ss_pseudovis_leiden'
  )
)

saveGiotto(g, dir = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/',
           method = 'qs',
           overwrite = TRUE)





# load in single cell dataset to act as reference
giotto_SC_subset_filtered_noBlood = loadGiotto("/projectnb/rd-spat/HOME/ecruiz/duckdb/giotto_SC_subset_filtered_noBlood/")

pDataDT(giotto_SC_subset_filtered_noBlood)[, table(Class)]
#             Choroid plexus              Early lineage 
#                         85                       1900 
#            Ensheating cell                  Ependymal 
#                         63                         38 
#                 Fibroblast                     Immune 
#                        771                        397 
#    Neural Mesenchyme cells          Neural progenitor 
#                       1637                      15624 
# Olfactory ensheathing cell                 Oligo-like 
#                         29                       1526 
#               Pineal gland                Radial glia 
#                         11                       4063 
#       Subcommissural organ                  Undefined 
#                         12                         17 
#                   Vascular 
#                        719 

# find top 50 scran gene markers for each of the 15 cell types
markers_scran_sc = findMarkers_one_vs_all(gobject=giotto_SC_subset_filtered_noBlood,
                                          method="scran",
                                          expression_values="normalized",
                                          cluster_column = "Class",
                                          min_feats=3)
top_markers_sc <- markers_scran_sc[, head(.SD, 50), by="cluster"]

# create DWLS signmatrix
DWLS_matrix <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(
    giotto_SC_subset_filtered_noBlood,
    values = "normalized",
    output = "matrix"
  ),
  cell_type = pDataDT(giotto_SC_subset_filtered_noBlood)$Class,
  sign_gene = top_markers_sc$feats
)

g = runDWLSDeconv(
  gobject = g, 
  spat_unit = "pseudo_vis",
  sign_matrix = DWLS_matrix
)

spatDeconvPlot(
  g,
  spat_unit = "pseudo_vis",
  radius = 90,
  save_param = list(
    save_name = 'pseudo_vis_DWLS',
    base_width = 10
  )
)

saveGiotto(g, dir = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/',
           method = 'qs',
           overwrite = TRUE)




# issue with pseudovisium initialization of IDs and metadata
# strip those values then subset for hexbins
g@expression$pseudo_vis <- NULL
g@spatial_locs$pseudo_vis <- NULL
g@spatial_enrichment$pseudo_vis <- NULL
g@nn_network$pseudo_vis <- NULL
g@cell_metadata$pseudo_vis <- NULL
g@feat_metadata$pseudo_vis <- NULL
g@spatial_info$pseudo_vis <- NULL
g@dimension_reduction$cells$pseudo_vis <- NULL
activeSpatUnit(g) = 'hex50'

saveGiotto(g, dir = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq_no_pseudovis/',
           method = 'qs',
           overwrite = TRUE)




brain_ext <- ext(1000, 3000, 9000, 12000)

s_sub <- crop(s, brain_ext)
# pull into mem
s_sub_sv <- as.spatvector(s_sub)

g_brain_sub = giotto()

hexbin3 <- tessellate(brain_ext, shape = 'hexagon', shape_size = 20, name = 'hex10')

h3 <- dbvect(
  hexbin3[],
  db = b,
  remote_name = 'hex10',
  overwrite = TRUE
)


g_brain_sub <- setPolygonInfo(g_brain_sub, hexbin3)

# find overlapped features 
tictoc::tic()
hex3_overlaps = calculateOverlap(
  h3, s,
  remote_name = 'hex10_overlaps',
  count_info_column = 'count',
  verbose = 'debug'
)
tictoc::toc()

expected_cn <- unique(h1$poly_ID)
expected_rn <- s[] %>% 
  dplyr::distinct(feat_ID) %>% 
  dplyr::pull()


tictoc::tic()
h1_M <- overlapToMatrix(
  h1_overlaps,
  col_names = expected_cn,
  row_names = expected_rn,
  count_info_column = 'count',
  write_step = c(1e5,1e3),
  filepath = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq/hex10_rna_raw_matrix.h5',
  remote_name = "hex10_rna_raw",
  overwrite = TRUE)
tictoc::toc()






# 546.814 sec elapsed






hex3_M <- overlapToMatrix(
  col_names = hexbin3$poly_ID,
  row_names = expected_rn,
  hex3_overlaps
)

exp_hex3 <- createExprObj(
  expression_data = hex3_M,
  name = 'raw',
  spat_unit = 'hex5',
  feat_type = 'rna',
  provenance = 'hex5'
)

g_brain_sub <- setExpression(g_brain_sub, exp_hex3)


# now that it is a smaller region, lets re-filter and re-norm then calculate 
# HVFs again since the local spatial patterns are expected to be driven by 
# different genes

g_brain_sub_filter <- filterGiotto(
  gobject = g_brain_sub,
  spat_unit = 'hex5',
  expression_threshold = 1,
  feat_det_in_min_cells = 1,
  min_det_feats_per_cell = 10
)

g_brain_sub_filter <- normalizeGiotto(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex5'
)

g_brain_sub_filter <- addSpatialCentroidLocations(g_brain_sub_filter,
                                                  poly_info = 'hex5')

g_brain_sub_filter <- addStatistics(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex5'
)

spatInSituPlotPoints(
  g_brain_sub_filter,
  polygon_fill = 'total_expr',
  spat_unit = 'hex5',
  polygon_fill_gradient_style = 's',
  polygon_alpha = 1,
  save_param = list(
    save_name = 'brain_hex5_total_expr'
  )
)

g_brain_sub_filter <- createSpatialNetwork(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex5',
  method = 'kNN',
  k = 6
)

# saveRDS(out, file = 'DATA/stereoseq/bin20_binspect.rds')

# binspect genes
# out <- binSpect(
#   g_brain_sub_filter,
#   spat_unit = 'hex20',
#   spatial_network_name = 'kNN_network',
#   expression_values = 'normalized'
# )

# out[score > 168, svg := 'yes']
# data.table::setnames(out, old = 'feats', new = 'feat_ID')
# 
# g_brain_sub_filter <- addFeatMetadata(
#   g_brain_sub_filter,
#   spat_unit = 'hex20',
#   new_metadata = out[, c('feat_ID', 'svg')]
# )

# find 25% of cell_ids
hex5_ids <- spatIDs(g_brain_sub_filter, 'hex5')
set.seed(1234)
hex5_ids_25p <- sample(hex5_ids, length(hex5_ids)/25)

g_brain_sub_filter <- calculateHVF(
  gobject = g_brain_sub_filter,
  cell_ids = hex5_ids_25p,
  spat_unit = 'hex5'
)

n_25_percent_brain <- round(length(hex20_ids) * 0.25)

# 324.589 sec elapsed
tictoc::tic()
g_brain_sub_filter<- runPCAprojection(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex20',
  expression_values = 'scaled',
  feats_to_use = 'svg',
  name = 'pca.projection',
  random_subset = n_25_percent_brain
)
tictoc::toc()


screePlot(
  gobject = g_brain_sub_filter,
  name = 'pca.projection',
  expression_values = 'normalized',
  spat_unit = 'hex20',
  ncp = 100,
  save_param = list(
    save_name = 'brain_hex20_screePlot')
)

plotPCA(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex20',
  dim_reduction_name = 'pca.projection',
  dim1_to_use = 1,
  dim2_to_use = 2,
  save_param = list(
    base_width = 10,
    save_name = 'brain_hex20_pca'
  )
)


# run the UMAP in a projected manner as well

# 42.767 sec elapsed
tictoc::tic()
# Worked on Matrix 1.6.3
g_brain_sub_filter = runUMAPprojection(
  gobject = g_brain_sub_filter,
  dim_reduction_to_use = 'pca',
  dim_reduction_name = 'pca.projection',
  spat_unit = "hex20",
  dimensions_to_use = 1:10,
  n_neighbors = 10,
  name = "umap.projection",
  random_subset = n_25_percent_brain
)
tictoc::toc()

plotUMAP(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex20',
  dim_reduction_name = 'umap.projection',
  save_param = list(
    save_name = 'brain_hex20_umap'
  )
)


# 57.325 sec elapsed
tictoc::tic()
g_brain_sub_filter = createNearestNetwork(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex20',
  dim_reduction_to_use = 'pca',
  dim_reduction_name = 'pca.projection',
  dimensions_to_use = 1:10,
  k = 10,
  type = 'sNN',
  name = 'sNN.pca.projection'
)
tictoc::toc()




g_brain_sub_filter <- doLeidenCluster(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex20',
  nn_network_to_use = 'sNN',
  network_name = 'sNN.pca.projection',
  resolution = 0.03,
  name = 'brain_leiden'
)

plotUMAP(
  gobject = g_brain_sub_filter,
  dim_reduction_name = 'umap.projection',
  spat_unit = 'hex20',
  cell_color = 'brain_leiden',
  save_param = list(
    save_name = 'brain_hex20_umap_leiden'
  )
)

spatInSituPlotPoints(
  gobject = g_brain_sub_filter,
  spat_unit = 'hex20',
  polygon_fill = 'brain_leiden',
  polygon_fill_as_factor = TRUE,
  polygon_alpha = 1,
  polygon_line_size = 0,
  save_param = list(
    base_width = 10,
    save_name = 'brain_hex20_leiden'
  )
)


saveGiotto(g, dir = '/projectnb/rd-spat/MANUSCRIPTS/giotto_suite/DATA/stereoseq_no_pseudovis/',
           method = 'qs',
           overwrite = TRUE)




