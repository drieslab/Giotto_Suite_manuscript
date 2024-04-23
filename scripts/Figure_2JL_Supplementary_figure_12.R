
# This script relies on GiottoDB which is a Giotto module that improves
# scalability by performing calculations on-disk in the database. These
# `proxy` implementations are experimental and provide a way to store the data
# and perform simple queries in a way that only requires the base duckdb.
# More complex calculations are performed chunkwise in memory then moved back
# to an on-disk representation.
if (!requireNamespace("GiottoDB", quietly = TRUE)) {
  remotes::install_github("drieslab/GiottoDB")
}
if(!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
if (!requireNamespace("HDF5Array", quietly = TRUE)) {
  BiocManager::install(c("HDF5Array"))
}

library(Giotto)
library(GiottoDB)

# create the database that contains the gef information
# see [scripts/CREATE/create_gef_db.R]


outs_dir <- "scripts/CREATE/OUTS/stereoseq_E16.5_E2S6"


# load DB materials
bInfo <- readRDS(file.path(outs_dir, 'backend_info.rds'))
b <- dbBackend(bInfo)

# `s` is a dbPointsProxy, which is a spatial points representation that is 
# backed by a duckdb representation. `loadGDB()` is a specific load function
# that reads `dbData` inheriting objects back into memory and also repairs the
# database connection.
s <- loadGDB(file.path(outs_dir, 'ss_points.rds'))

# get extent
e <- ext(s)

# generate hexbin tessellations that span the entire extent of the points
# data in `s`.
# The diameter is calculated by the total x distance that each hexagon
# traverses.
# diameter 400 hexbins
hexbin400 <- tessellate(e, shape = 'hexagon', shape_size = 400, name = 'hex400') 
# diameter 100 hexbins (finer resolution)
hexbin100 <- tessellate(e, shape = 'hexagon', shape_size = 100, name = 'hex100')

# write the tesselations to the DB backend as dbPolygonProxy so that they can 
# interact with the points representations.
h400 <- dbvect(
  x = hexbin400[],
  db = b,
  remote_name = 'hex400',
  overwrite = TRUE
)
h100 <- dbvect(
  x = hexbin100[],
  db = b,
  remote_name = 'hex100',
  overwrite = TRUE
)



# calculate overlaps ####
# 
# These steps will take a while
# The calculateOverlap() method for dbPolygonProxy and dbPointsProxy objects
# will generate a new dbPointsProxy of overlaps information in a chunkwise
# manner. Spatial chunks of the data are retrieved from the database,
# loaded into terra and Giotto, processed, then the results are fed back into
# the database.
h400_overlaps <- calculateOverlap(
  h400, s,
  n_per_chunk = 400, # number of polys per chunk to process
  count_info_column = "count",
  overwrite = TRUE
)

# get the expected colnames and rownames as the unique poly_IDs and unique
# feat_IDs respectively. Without providing these, the matrix generated will
# only include bins (cols) and feats (rows) that were actually overlapped.
expected_cn <- spatIDs(hexbin400)
expected_rn <- s[] %>% # same for both resolutions
  dplyr::distinct(feat_ID) %>% 
  dplyr::pull()

# The overlaps are essentially a set of overlap relationships that have been
# calculated. From that, we can convert it into a set of counts information.
# These values are written stepwise to on-disk HDF5Matrix format.
h400_M <- overlapToMatrix(
  h400_overlaps,
  col_names = expected_cn,
  row_names = expected_rn,
  count_info_column = 'count',
  write_step = c(1e5,1e3),
  filepath = file.path(outs_dir, 'hex400_rna_raw_matrix.h5'),
  remote_name = "hex400_rna_raw",
  overwrite = TRUE
)


# repeat for a finer resolution
h100_overlaps = calculateOverlap(
  h100, s,
  n_per_chunk = 6000,
  count_info_column = 'count',
  overwrite = TRUE
)

h100_M <- overlapToMatrix(
  h100_overlaps,
  col_names = spatIDs(hexbin100),
  row_names = expected_rn,
  count_info_column = 'count',
  write_step = c(1e5,1e3),
  filepath = file.path(outs_dir, 'hex100_rna_raw_matrix.h5'),
  remote_name = "hex100_rna_raw",
  overwrite = TRUE
)

# * end of slow portion * #


# since this data is large in terms of number of "point" detections, but not
# that large when aggregated, we can simply convert the h5 back to Matrix
# formats in order to speed up the downstream steps.
# 
# We also load the values from disk here and remake other elements here in order 
# to make resuming of the script easier since the previous step can take a while.
e <- ext(2826, 29774, 2975, 19225) # extent of s from earlier

h400_M_path <- file.path(outs_dir, 'hex400_rna_raw_matrix.h5')
h100_M_path <- file.path(outs_dir, 'hex100_rna_raw_matrix.h5')
h400_M <- HDF5Array::HDF5Array(h400_M_path, name = 'hex400_rna_raw')
h100_M <- HDF5Array::HDF5Array(h100_M_path, name = 'hex100_rna_raw')
h400_Matrix <- Matrix::Matrix(as.matrix(h400_M), sparse = TRUE)
h100_Matrix <- Matrix::Matrix(as.matrix(h100_M), sparse = TRUE)
hexbin400 <- tessellate(e, shape = 'hexagon', shape_size = 400, name = 'hex400') 
hexbin100 <- tessellate(e, shape = 'hexagon', shape_size = 100, name = 'hex100')


exp1 <- createExprObj(
  h400_Matrix,
  name = 'raw',
  spat_unit = 'hex400',
  feat_type = 'rna',
  provenance = "hex400"
)

exp2 <- createExprObj(
  h100_Matrix,
  name = 'raw',
  spat_unit = 'hex100',
  feat_type = 'rna',
  provenance = "hex100"
)

# create Giotto object
g <- createGiottoObject(
  expression = list(exp1, exp2),
  spatial_info = list(hexbin400, hexbin100)
)

instructions(g, 'save_dir') <- "scripts/FIGURE/S12/"
instructions(g, 'show_plot') <- FALSE
instructions(g, 'save_plot') <- TRUE
instructions(g, 'return_plot') <- FALSE


# calculate and plot total raw expression
g <- addCellStatistics(g, spat_unit = "hex400", expression_values = 'raw')
g <- addCellStatistics(g, spat_unit = "hex100", expression_values = 'raw')

# FIGURE S12 A1 --------------------------------------------------------- ####
spatInSituPlotPoints(
  g, polygon_feat_type = 'hex400',
  polygon_fill = 'total_expr',
  polygon_alpha = 1, 
  polygon_line_size = 0, 
  polygon_fill_gradient_style = 's',
  save_param = list(
    save_name = "hex400_total_expr",
    save_format = "svg"
  )
)

# FIGURE S12 A2 --------------------------------------------------------- ####
spatInSituPlotPoints(
  g, polygon_feat_type = 'hex100', 
  polygon_fill = 'total_expr',
  polygon_alpha = 1,
  polygon_line_size = 0, 
  polygon_fill_gradient_style = "s",
  save_param = list(
    save_name = "hex100_total_expr",
    save_format = "svg"
  )
)




# filtering ####

filterCombinations(
  gobject = g,
  spat_unit = 'hex100',
  expression_thresholds = c(1, 2),
  feat_det_in_min_cells = c(10, 100, 100, 100),
  min_det_feats_per_cell = c(100, 100, 300, 1000),
  save_param = list(
    save_name = "filter_combinations"
  )
)
# Based on these results, keep the feat_det_in_min_cells param low so that
# the removal of features is not too aggressive. min_det_feats_per_cell is
# reasonable around 300. Raising the threshold to 2 removes features too
# aggressively.

g <- filterGiotto(
  gobject = g,
  spat_unit = 'hex400',
  expression_threshold = 1,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 300,
  spat_unit_fsub = 'hex400' # limit the feats subsets to this spat unit
)


g <- filterGiotto(
  gobject = g,
  spat_unit = 'hex100',
  expression_threshold = 1,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 300,
  spat_unit_fsub = 'hex100' # limit the feats subsets to this spat unit
)


# normalization ####
g <- normalizeGiotto(
  gobject = g,
  spat_unit = 'hex100'
)


# HVF detection ####
# calculate HVF on a 10% subset so that it takes less time.
# Results should be comparable
g <- calculateHVF(
  gobject = g,
  spat_unit = "hex100",
  random_subset = round(length(spatIDs(g, "hex100")) / 10),
  set_seed = TRUE,
  seed_number = 12345
)


# PCA projection ####
# Use PCA projection as a strategy for dealing with larger expression data.
n_25_percent <- round(length(spatIDs(g, 'hex100')) * 0.25)

g <- runPCAprojection(
  gobject = g,
  spat_unit = 'hex100',
  expression_values = 'scaled',
  feats_to_use = 'hvf',
  name = 'pca.projection',
  set_seed = TRUE,
  seed_number = 12345,
  random_subset = n_25_percent
)


# Visualize Screeplot and PCA

screePlot(
  g,
  name = 'pca.projection',
  expression_values = 'normalized',
  spat_unit = 'hex100',
  ncp = 50,
  save_param = list(
    save_name = 'ss_hex100_screePlot')
)

plotPCA(
  g,
  spat_unit = 'hex100',
  dim_reduction_name = 'pca.projection',
  dim1_to_use = 1,
  dim2_to_use = 2,
  save_param = list(
    base_width = 10,
    save_name = 'ss_hex100_pca'
  )
)

# UMAP projection ####
# run the UMAP in a projected manner as well
# Worked on Matrix 1.6.3
g = runUMAPprojection(
  gobject = g,
  spat_unit = "hex100",
  dim_reduction_to_use = 'pca',
  dim_reduction_name = "pca.projection",
  dimensions_to_use = 1:25,
  name = "umap.projection",
  random_subset = n_25_percent
)


g = createNearestNetwork(
  gobject = g,
  spat_unit = 'hex100',
  dim_reduction_to_use = 'pca',
  dim_reduction_name = 'pca.projection',
  dimensions_to_use = 1:25,
  k = 25,
  type = 'sNN',
  name = 'sNN.pca.projection'
)



# leiden clustering projection ####
# We will again be using a projection approach for leiden clustering in order
# to speed things up.
# This is done by taking a subset of the total dataset, performing the
# clustering on the subset, then projecting those results to the full dataset.

set.seed(1234)
subset_IDs = sample(x = spatIDs(g, 'hex100'), size = n_25_percent)
subset_g = subsetGiotto(
  gobject = g,
  spat_unit = 'hex100',
  cell_ids = subset_IDs
)

# remake nn network since this was a random sampling, so connections are
# expected to be broken.
subset_g = createNearestNetwork(
  gobject = subset_g,
  spat_unit = 'hex100',
  dim_reduction_to_use = 'pca',
  dim_reduction_name = 'pca.projection',
  dimensions_to_use = 1:25,
  k = 25,
  type = 'sNN',
  name = 'sNN.pca.projection'
)

subset_g = doLeidenClusterIgraph(
  gobject = subset_g,
  spat_unit = "hex100",
  resolution_parameter = 0.55,
  name = "sub_leiden_clus",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca.projection",
) 

# preview clustering of subset
dimPlot2D(
  gobject = subset_g,
  spat_unit = 'hex100',
  dim_reduction_to_use = 'umap',
  dim_reduction_name = 'umap.projection',
  cell_color = 'sub_leiden_clus',
  save_param = list(
    save_name = 'ss_hex100_subset_leiden',
    base_width = 10
  )
)

# project clusterings back to full dataset
g <- doClusterProjection(
  target_gobject = g,
  source_gobject = subset_g,
  spat_unit = "hex100",
  source_cluster_labels = "sub_leiden_clus",
  reduction_method = 'pca',
  reduction_name = 'pca.projection',
  prob = FALSE,
  knn_k = 25, 
  dimensions_to_use = 1:25
)



# FIGURE F2 J1 ---------------------------------------------------------- ####

# plot results
spatInSituPlotPoints(
  gobject = g,
  spat_unit = 'hex100',
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'knn_labels', # leiden clusterings projected by kNN
  polygon_alpha = 1,
  polygon_line_size = 0,
  save_param = list(
    base_width = 10,
    save_name = 'ss_hex100_leiden_clus',
    save_format = "svg"
  )
)






# brain plots ####
# Focus in on the brain area for some more focused examples

# create subset extent
brain_ext <- ext(3826, 5826, 11975, 14975)

# These steps can be done in memory because of the fewer feature detections
# in play

# polys
h100 <- tessellate(brain_ext, shape = 'hexagon', shape_size = 100, name = 'hex100') 
h50 <- tessellate(brain_ext, shape = 'hexagon', shape_size = 50, name = 'hex50') 
# points
s_sub <- crop(s, brain_ext) # spatially crop the feature detections
s_sub_sv <- as.spatvector(s_sub) # pull into memory

brain <- createGiottoObject(
  spatial_info = list(hex100 = h100, hex50 = h50),
  feat_info = createGiottoPoints(s_sub_sv)
)

instructions(brain, 'save_dir') <- "scripts/FIGURE/S12/"
instructions(brain, 'show_plot') <- TRUE
instructions(brain, 'save_plot') <- TRUE
instructions(brain, 'return_plot') <- FALSE

brain <- addSpatialCentroidLocations(brain, poly_info = "hex100")
brain <- addSpatialCentroidLocations(brain, poly_info = "hex50")

# FIGURE 2 J2 ---------------------------------------------------------- ####
# visualize raw transcripts
spatInSituPlotPoints(
  brain,
  show_polygon = FALSE,
  feats = list(rna = c(
    "Hmgb2",
    "Krt5",
    "Necab2",
    "Neurod6",
    "Pantr1",
    "Rnd2",
    "Trpm3"
  )),
  point_size = 0.75,
  use_overlap = FALSE,
  jitter = c(25,25),
  save_param = list(
    save_name = "brain_pnts",
    save_format = "svg"
  )
)


# Workflow for filtering, normalization, and highly variable gene detection
prep_workflow <- function(
    gobject,
    spat_unit,
    expression_threshold,
    feat_det_in_min_cells,
    min_det_feats_per_cell
) {
  activeSpatUnit(gobject) <- spat_unit
  
  gobject <- calculateOverlap(
    gobject,
    spatial_info = spat_unit,
    count_info_column = "count",
    verbose = FALSE
  )
  
  # expected rownames and colnames are automatically detected for
  # the gobject overlapToMatrix method.
  gobject <- overlapToMatrix(
    gobject,
    poly_info = spat_unit,
    count_info_column = "count",
    verbose = FALSE
  )
  
  gobject <- filterGiotto(
    gobject = gobject,
    expression_threshold = expression_threshold,
    feat_det_in_min_cells = feat_det_in_min_cells,
    min_det_feats_per_cell = min_det_feats_per_cell,
    spat_unit_fsub = spat_unit, # limit the feats subsets to this spat unit
    verbose = FALSE
  )
  
  gobject <- normalizeGiotto(
    gobject = gobject, verbose = FALSE
  )
  
  gobject <- calculateHVF(gobject, save_plot = FALSE, show_plot = FALSE)
  
  return(gobject)
}

# workflow for dimension reduction and leiden clustering
leiden_workflow <- function(
    gobject,
    spat_unit,
    run_prep = TRUE,
    ..., # passed to prep_workflow
    pca_dims_to_use,
    pca_feats = "hvf",
    k_neighbors = 4,
    leiden_res
) {
  activeSpatUnit(gobject) <- spat_unit
  
  if (run_prep) {
    gobject <- prep_workflow(
      gobject = gobject,
      spat_unit = spat_unit,
      ...
    )
  }
  
  gobject <- runPCA(
    gobject, feats_to_use = pca_feats, verbose = FALSE
  )
  
  gobject <- runUMAP(
    gobject,
    dimensions_to_use = pca_dims_to_use,
    n_neighbors = 4,
    verbose = FALSE
  )
  
  gobject <- createNearestNetwork(
    gobject = gobject,
    dimensions_to_use = pca_dims_to_use,
    k = k_neighbors,
    verbose = FALSE
  )
  
  gobject <- doLeidenCluster(
    gobject = gobject,
    resolution = leiden_res
  )
  
  return(gobject)
}

dim_spat_plots <- function(
    spat_unit,
    color,
    categorical = TRUE, 
    dim_point_size = 1,
    name
  ) {
  dimPlot2D(
    gobject = brain,
    spat_unit = spat_unit,
    point_size = dim_point_size,
    cell_color = color,
    save_param = list(
      save_name = sprintf("brain_%s_umap", name),
      save_format = "svg"
    )
  )
  
  spatInSituPlotPoints(
    gobject = brain,
    polygon_feat_type = spat_unit,
    polygon_fill = color,
    polygon_fill_as_factor = categorical,
    polygon_alpha = 1,
    save_param = list(
      save_name = sprintf("brain_%s_spat", name),
      save_format = "svg"
    ),
    verbose = FALSE
  )
}

# FIGURE S12 B1 ---------------------------------------------------------- ####
# hex 100
brain <- leiden_workflow(
    gobject = brain,
    spat_unit = "hex100",
    expression_threshold = 1,
    feat_det_in_min_cells = 10,
    min_det_feats_per_cell = 300,
    pca_dims_to_use = seq(25),
    k_neighbors = 25,
    leiden_res = 1.7
) 

dim_spat_plots(
  spat_unit = "hex100",
  color = "leiden_clus",
  name = "hex100"
)

# hex50
brain <- leiden_workflow(
  gobject = brain,
  spat_unit = "hex50",
  expression_threshold = 1,
  feat_det_in_min_cells = 1,
  min_det_feats_per_cell = 10,
  pca_dims_to_use = seq(10),
  k_neighbors = 50,
  leiden_res = 1
) 

dim_spat_plots(
  spat_unit = "hex50",
  color = "leiden_clus",
  name = "hex50"
)

# Clustering is not very clear at this binning resolution. This is due to the
# noisiness of the RNA expression being captured coming through as bins become
# smaller. This noisiness disturbs unsupervised clustering.
# 
# Since the clustering should be highly spatially organized, we can use only
# the genes that are spatially variable when generating the PCA space that
# UMAP and leiden are performed on in order to get around the noisiness of
# the data.


# Spatially Variable Genes ####

# Detection of spatially patterned genes via `binSpect()`.
# A spatial network of 6 nearest neighbors is convenient since hexbins have
# 6 equidistant neighbors.

brain <- createSpatialNetwork(
  gobject = brain,
  spat_unit = 'hex50',
  method = 'kNN',
  k = 6
)

# will take a little time
svgs <- binSpect(
  brain,
  spat_unit = 'hex50',
  spatial_network_name = 'kNN_network',
  expression_values = 'normalized'
)
data.table::setorder(svgs, -score) # ensure ordered by descending score
data.table::setnames(svgs, old = "feats", new = "feat_ID")

# svgs is now ordered by how spatially patterned the expression of these
# features is. However, while spatial patterns are defined by expression,
# the number genes that follow that similar spatial pattern or
# "spatial coexpression module" may not be the same within each module.
# We must first find the coexpression modules formed by these spatial genes
# then select equal numbers of genes from each module to formed a balanced
# set of genes that can be used to construct the PCA space without biasing
# for any particular spatial pattern.

# top 1000 spatial genes
top_svgs = svgs[seq(1000), feat_ID]

# here we use existing detectSpatialCorGenes function to calculate pairwise 
# distances between genes while smoothening the gene expression matrix using
# spatial network neighbors to deal with noisy expression.
spatcor_obj = detectSpatialCorFeats(
  brain,
  spat_unit = "hex50",
  method = 'network',
  spatial_network_name = 'kNN_network',
  subset_feats = top_svgs
)

# cluster spatial genes into 20 modules
spatcor_obj = clusterSpatialCorFeats(
  spatcor_obj, 
  name = 'spat_netw_clus', 
  k = 15
)


# FIGURE S12 C1 ---------------------------------------------------------- ####

# visualize clusters
heatmSpatialCorFeats(
  brain,
  spatCorObject = spatcor_obj,
  use_clus_name = 'spat_netw_clus',
  heatmap_legend_param = list(title = NULL),
  save_param = list(
    save_name = "brain_spatcor_heatmap",
    save_format = "svg"
  )
)

# create a metagene for each of the modules
module_genes_dt <- showSpatialCorFeats(
  spatCorObject = spatcor_obj, 
  use_clus_name = 'spat_netw_clus', 
  show_top_feats = 1
)

module_genes <- module_genes_dt$clus
names(module_genes) <- module_genes_dt$feat_ID

brain <- createMetafeats(
  brain,
  spat_unit = "hex50",
  feat_clusters = module_genes, 
  name = 'spatcor_modules'
)

# FIGURE S12 C2-5 --------------------------------------------------------- ####

module_plot <- function(modules) {
  modules <- as.character(modules)
  for (m_i in seq_along(modules)) {
    spatInSituPlotPoints(
      brain,
      polygon_feat_type = "hex50",
      spat_enr_names = 'spatcor_modules',
      polygon_fill = modules[m_i],
      background = "black",
      polygon_alpha = 1,
      polygon_fill_gradient_style = "sequential",
      save_param = list(
        save_name = paste0("hex50_spatcor_modules_", m_i)
      ),
      verbose = FALSE
    )
  }
}

module_plot(c(1, 8, 13, 10))

# pull the top 30 representative genes for each spatial module
# total of 200 balanced SVGs
balanced_spatmodule_feats <- getBalancedSpatCoexpressionFeats(
  spatCorObject = spatcor_obj,
  maximum = 20,
  rank = "weighted"
)

bsm_feats_dt <- data.table::data.table(
  feat_ID = names(balanced_spatmodule_feats),
  svg = "yes"
)

# append SVG information into metdata to be used during PCA construction
brain <- addFeatMetadata(
  brain,
  spat_unit = 'hex50',
  by_column = TRUE,
  new_metadata = bsm_feats_dt
)

# FIGURE S12 B2 ---------------------------------------------------------- ####
# repeat with svgs for PCA space generation
brain <- leiden_workflow(
  gobject = brain,
  spat_unit = "hex50",
  run_prep = FALSE,
  pca_feats = "svg",
  pca_dims_to_use = seq(25),
  k_neighbors = 20,
  leiden_res = 1.2
) 

dim_spat_plots(
  spat_unit = "hex50",
  color = "leiden_clus",
  name = "hex50_svg"
)


# There is a lot of spatial patterns being found, but the clustering
# is a bit noisy. We can further refine these results via niche
# clustering.



# niche clustering ####

# Calculate the enrichment of a leiden annotation for any particular
# hexbin based on the leiden cluster assignment of itself and its 
# neighbors in the previously generated kNN k = 6 spatial network, 
# producing a niche-based spatial enrichment. Then use kmeans to select
# the final niche to belong to.
brain <- calculateSpatCellMetadataProportions(
  brain,
  spat_unit = "hex50",
  spat_network = "kNN_network",
  metadata_column = "leiden_clus",
  name = "leiden_niche",
  return_gobject = TRUE
)

prop_table = getSpatialEnrichment(
  brain, spat_unit = 'hex50', name = 'leiden_niche', output = 'data.table'
)
# convert the data.table to a sparse Matrix with row and colnames
# here we use a utility function to perform the operation
prop_matrix = GiottoUtils::dt_to_matrix(prop_table)


# These enrichments are essentially a measure of how many cells of each
# leiden cluster exist in the local region
# 
# Using kmeans, we can classify each cell by its niche leiden cluster proportions
set.seed(12345) # set seed for kmeans
prop_kmeans = kmeans(x = prop_matrix, centers = 15, iter.max = 100, nstart = 3)
prop_kmeansDT = data.table::data.table(
  cell_ID = names(prop_kmeans$cluster), 
  niche = prop_kmeans$cluster
)

# add kmeans clustering of niches to cell metadata
brain <- addCellMetadata(
  brain,
  spat_unit = "hex50" , 
  new_metadata = prop_kmeansDT, 
  by_column = TRUE,
  column_cell_ID = "cell_ID"
)

# FIGURE F2 J3 ---------------------------------------------------------- ####

spatInSituPlotPoints(
  brain,
  polygon_feat_type = "hex50",
  polygon_fill = "niche",
  polygon_fill_as_factor = TRUE,
  polygon_alpha = 1,
  save_param = list(
    save_name = "brain_hex50_spat_svg_niche",
    save_format = "svg"
  )
)




# pseudovisium ####

# create pseudovisium polygons. Micron size taken from StereoSeq.
# The StereoSeq dataset is scaled based the spatial DNB grid array used to
# capture mRNA. The center to center distance between two DNBs is roughly 
# 500nm or 0.5 micron.
visium_polygons = makePseudoVisium(
  extent = brain_ext,
  micron_size = 0.5,
  name = "pseudo_vis"
)

brain <- setGiotto(brain, visium_polygons)
brain <- addSpatialCentroidLocations(brain, poly_info = "pseudo_vis")


# FIGURE S12 D1 ---------------------------------------------------------- ####
spatInSituPlotPoints(
  brain,
  feats = list(rna = c(
    "Hmgb2",
    "Krt5",
    "Necab2",
    "Neurod6",
    "Pantr1",
    "Rnd2",
    "Trpm3"
  )),
  point_size = 0.3,
  use_overlap = FALSE,
  polygon_feat_type = "pseudo_vis",
  polygon_color = "white",
  polygon_alpha = 0.4,
  jitter = c(25,25),
  save_param = list(
    save_name = "brain_pts_pvis",
    save_format = "svg"
  )
)

# FIGURE S12 D2 ---------------------------------------------------------- ####
brain <- leiden_workflow(
  gobject = brain,
  spat_unit = "pseudo_vis",
  expression_threshold = 1,
  feat_det_in_min_cells = 1,
  min_det_feats_per_cell = 10,
  pca_dims_to_use = seq(30),
  k_neighbors = 4,
  leiden_res = 1
) 

dim_spat_plots(
  spat_unit = "pseudo_vis",
  dim_point_size = 3,
  color = "leiden_clus",
  name = "pvis"
)





# DWLS ####

# load in single cell dataset to act as reference
scmb_link <- "https://zenodo.org/records/10960725/files/scmb_blood_filtered_gobject.RDS"
dest_dir <- "scripts/CREATE/EXT_DATA/sc_mousebrain"
dir.create(dest_dir, recursive = TRUE)
utils::download.file(
  url = scmb_link, 
  destfile = file.path(dest_dir, "gobject.RDS")
)

scmb = loadGiotto(dest_dir)

# print the types of cell classifications available
pDataDT(scmb)[, table(Class)]

# find top 50 scran gene markers for each of the 15 cell types
markers_scran_sc = findMarkers_one_vs_all(
  gobject = scmb,
  method = "scran",
  expression_values = "normalized",
  cluster_column = "Class",
  min_feats = 3
)

top_markers_sc <- markers_scran_sc[, head(.SD, 50), by = "cluster"]

# create DWLS signmatrix
DWLS_matrix <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(
    scmb,
    values = "normalized",
    output = "matrix"
  ),
  cell_type = pDataDT(scmb)$Class,
  sign_gene = top_markers_sc$feats
)

brain = runDWLSDeconv(
  gobject = brain, 
  spat_unit = "pseudo_vis",
  sign_matrix = DWLS_matrix
)

# FIGURE F2 L ----------------------------------------------------------- ####
spatDeconvPlot(
  brain,
  spat_unit = "pseudo_vis",
  radius = 90,
  save_param = list(
    save_name = 'pseudo_vis_DWLS',
    base_width = 10
  )
)

