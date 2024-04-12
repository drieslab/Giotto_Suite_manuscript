
## %######################################################%##
#                                                          #
####     Supplementary figure 3 Nanostring dataset      ####
#                                                          #
## %######################################################%##


############################## LOAD LIBRARIES ##################################

library(terra)
library(patchwork)
library(ggplot2)
library(Giotto)
# Custom color palettes from rcartocolor

# pal10 = rcartocolor::carto_pal(n = 10, name = 'Pastel')
pal10 <- c(
  "#66C5CC", "#F6CF71", "#F89C74", "#DCB0F2", "#87C55F",
  "#9EB9F3", "#FE88B1", "#C9DB74", "#8BE0A4", "#B3B3B3"
)

# viv10 = rcartocolor::carto_pal(n = 10, name = 'Vivid')
viv10 <- c(
  "#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0",
  "#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#A5AA99"
)

############################ Set seed, Initialize ###############################


# set working directory to results folder for plot saving
results_folder <- "results_nanostring/"
my_seed_num <- 315
set.seed(my_seed_num)
my_python_path <- NULL # alternatively, "/local/python/path/python" if desired.

instrs <- createGiottoInstructions(
  save_dir = results_folder,
  save_plot = FALSE,
  show_plot = TRUE,
  return_plot = TRUE,
  python_path = my_python_path
)

## provide path to nanostring folder
data_path <- "Lung12/Lung12-Flat_files_and_images/"

########################### Create Giotto Object ###############################

## create giotto cosmx object
fov_join <- createGiottoCosMxObject(
  cosmx_dir = data_path,
  data_to_use = "subcellular", # only subcellular
  FOVs = c(
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28
  ),
  instructions = instrs
)

showGiottoFeatInfo(fov_join)
showGiottoSpatialInfo(fov_join)


id_set <- c(
  "01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
  "14", "15", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28"
)

new_names <- paste0("fov0", id_set)

## Set up vector of image names
image_names <- paste0(new_names, "-image")

######################### Create Expression Matrix #############################

# Find the feature points overlapped by polygons. This overlap information is then
# returned to the relevant giottoPolygon object's overlaps slot.
fov_join <- calculateOverlapRaster(fov_join, feat_info = "rna")
fov_join <- calculateOverlapRaster(fov_join, feat_info = "neg_probe")

# Convert the overlap information into a cell by feature expression matrix which
# is then stored in the Giotto object's expression slot
fov_join <- overlapToMatrix(fov_join, feat_info = "rna")
fov_join <- overlapToMatrix(fov_join, feat_info = "neg_probe")

showGiottoExpression(fov_join)

############################# Combine Metadata #################################

# combine cell data
morphometa <- combineCellData(fov_join,
                              feat_type = "rna"
)

# combine feature data
featmeta <- combineFeatureData(fov_join,
                               feat_type = c("rna")
)

# combine overlapping feature data
featoverlapmeta <- combineFeatureOverlapData(fov_join,
                                             feat_type = c("rna")
)

################### Filtering, Normalization, and Stats ########################

# filter
fov_join <- filterGiotto(
  gobject = fov_join,
  feat_type = "rna",
  expression_threshold = 1,
  feat_det_in_min_cells = 5,
  min_det_feats_per_cell = 5
)

# normalize
# standard method of normalization (log normalization based)
fov_join <- normalizeGiotto(
  gobject = fov_join,
  feat_type = "rna",
  norm_methods = "standard",
  verbose = TRUE
)
fov_join <- normalizeGiotto(
  gobject = fov_join,
  feat_type = "neg_probe",
  norm_methods = "standard",
  library_size_norm = FALSE,
  verbose = TRUE
)



showGiottoExpression(fov_join)

# add statistics based on log normalized values for features rna and negative probes
fov_join <- addStatistics(
  gobject = fov_join,
  expression_values = "normalized",
  feat_type = "rna"
)
fov_join <- addStatistics(
  gobject = fov_join,
  expression_values = "normalized",
  feat_type = "neg_probe"
)

# View cellular data (default is feat = 'rna')
showGiottoCellMetadata(fov_join)
# View feature data
showGiottoFeatMetadata(fov_join)


############### Dimension Reductions, NN Networks, Clustering ##################

# PCA
fov_join <- runPCA(fov_join,
                   scale_unit = FALSE,
                   center = FALSE,
                   expression_values = "normalized"
)

# Generate UMAP from PCA
fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:10,
                    seed_number = my_seed_num,
                    set_seed = TRUE
)

# Shared Nearest Neighbor Network
fov_join <- createNearestNetwork(
  gobject = fov_join,
  dimensions_to_use = 1:10,
  k = 10
)
# Leiden Clustering
fov_join <- doLeidenCluster(
  gobject = fov_join,
  resolution = 0.15,
  n_iterations = 100,
  seed_number = my_seed_num,
  set_seed = TRUE
)


################## Spatially Variable Gene Identification ######################

# create spatial network based on physical distance of cell centroids
fov_join <- createSpatialNetwork(
  gobject = fov_join,
  minimum_k = 2,
  maximum_distance_delaunay = 50
)


# Perform Binary Spatial Extraction of genes
# NOTE: Depending on your system this could take time
rank_spatialgenes <- binSpect(fov_join,
                              bin_method = "rank",
                              spatial_network_name = "Delaunay_network",
                              calc_hub = TRUE,
                              hub_min_int = 5
)

# Extract the top 500 spatially variable genes
top_spat_var_genes <- rank_spatialgenes$feats[1:500]

# Find spatially correlated features
scn_DT <- detectSpatialCorFeats(fov_join,
                                method = "network",
                                spatial_network_name = "Delaunay_network",
                                subset_feats = top_spat_var_genes
)

scn_DT <- clusterSpatialCorFeats(scn_DT,
                                 name = "spat_netw_clus",
                                 k = 6
)

# extract an equal amount of features from each cluster
cluster_genes <- getBalancedSpatCoexpressionFeats(scn_DT,
                                                  maximum = 30
)

# Create a metagene representation of the extracted features from above
fov_join <- createMetafeats(fov_join,
                            feat_clusters = cluster_genes,
                            name = "cluster_metagene"
)

svg_names <- names(cluster_genes)


####################### Spatially Informed Clustering ##########################

# Provide spatially variable gene names to the feats_to_use argument
# and give this PCA a new name
fov_join <- runPCA(
  gobject = fov_join,
  feats_to_use = svg_names,
  name = "custom_pca",
  set_seed = TRUE,
  seed_number = my_seed_num
)


fov_join <- runUMAP(
  gobject = fov_join,
  dim_reduction_name = "custom_pca",
  dimensions_to_use = 1:10,
  n_neighbors = 5,
  name = "custom_UMAP",
  set_seed = TRUE,
  seed_number = my_seed_num
)

fov_join <- createNearestNetwork(
  gobject = fov_join,
  dim_reduction_name = "custom_pca",
  dimensions_to_use = 1:10,
  name = "custom_NN"
)

fov_join <- doLeidenClusterIgraph(
  gobject = fov_join,
  network_name = "custom_NN",
  resolution_parameter = 0.3,
  n_iterations = 100,
  name = "custom_leiden",
  set_seed = TRUE,
  seed_number = my_seed_num
)


################## Visualize Spatially Informed Clustering #####################

pdf("spatPlot_custom_clustering.pdf", width = 20, height = 10)
sp_pl <- spatPlot2D(fov_join,
                    cell_color = "custom_leiden",
                    point_size = 4,
                    show_plot = TRUE,
                    save_plot = FALSE,
                    return_plot = TRUE
)
dev.off()

leiden_colors <- getDistinctColors(length(unique(getCellMetadata(fov_join)[]$custom_leiden)))
names(leiden_colors) <- sort(unique(getCellMetadata(fov_join)[]$custom_leiden))

sisp_pl <- spatInSituPlotPoints(fov_join,
                                show_polygon = TRUE,
                                polygon_color = "white",
                                polygon_line_size = 0.05,
                                polygon_fill = "custom_leiden",
                                polygon_fill_code = leiden_colors,
                                polygon_fill_as_factor = TRUE,
                                show_plot = TRUE,
                                return_plot = TRUE,
                                save_plot = FALSE
)

ggplot2::ggsave(
  filename = "s3_nanostring_spatInSituPlot_custom_clustering.png",
  plot = sisp_pl,
  device = "png",
  dpi = "retina",
  width = 11,
  height = 11,
  units = "in"
)


# Subset and visualize polygons with transcripts

ROI <- subsetGiottoLocs(
  gobject = fov_join,
  spat_unit = "cell",
  feat_type = "rna",
  spat_loc_name = "raw",
  x_min = 9000,
  x_max = 12000,
  y_min = -138000,
  y_max = -120000,
  return_gobject = TRUE
)

roi_sisp_pl <- spatInSituPlotPoints(ROI,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    spat_loc_name = "raw",
                                    feats = list(
                                      "RPL21",
                                      "IGHG1",
                                      "MALAT1",
                                      "KRT19",
                                      "COL1A1",
                                      "CD74"
                                    ),
                                    point_size = 0.25,
                                    show_polygon = TRUE,
                                    use_overlap = FALSE,
                                    polygon_color = "white",
                                    polygon_line_size = 0.05,
                                    polygon_alpha = 0.33,
                                    polygon_fill = "custom_leiden",
                                    polygon_fill_code = leiden_colors,
                                    polygon_fill_as_factor = TRUE,
                                    show_legend = FALSE,
                                    show_plot = TRUE,
                                    return_plot = TRUE,
                                    save_plot = FALSE
)


ggplot2::ggsave(
  filename = "s3_nanostring_ROI_spatInSituPlot_custom_clustering.png",
  plot = roi_sisp_pl,
  device = "png",
  dpi = "retina",
  width = 6,
  height = 6,
  units = "in"
)
