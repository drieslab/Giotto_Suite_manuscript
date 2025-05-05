## %######################################################%##
#                                                          #
####      Supplementary figure 3 SeqScope dataset        ####
#                                                          #
## %######################################################%##

############################## Download dataset  ###############################
## Download
## - 1st-Seq coordinate information for liver dataset (MiSeq-DraI-100pM-mbcore-RD2-revHDMIs-pos-uniq.txt)
## - 1st-Seq raw MiSeq FASTQ output file for liver dataset (DraI-100pM-mbcore-RD2.fastq.gz)

## Download the Seq-Scope pipeline (https://github.com/leeju-umich/Cho_Xi_Seqscope.git)
## for creating the coordinates file

## Run the following command to create spatialcoordinates.txt
## bash extractCoord.sh DraI-100pM-mbcore-RD2.fastq MiSeq-DraI-100pM-mbcore-RD2-revHDMIs-pos-uniq.txt 20

## Download the expression matrix files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5212844
## You will need barcodes.tsv, features.tsv, and matrix.mtx. Place the three files in your data directory

library(Giotto)
library(data.table)
library(tidyverse)

library(terra)
library(Matrix)
library(scattermore)

############################## Set instructions  ###############################
python_path <- NULL

results_folder <- "results"

instructions <- createGiottoInstructions(
    save_dir = results_folder,
    save_plot = TRUE,
    show_plot = TRUE,
    python_path = python_path
)

dataDir <- "data/raw/"

############################## Create the object  ##############################

# read raw data

## expression matrix
expr_matrix <- Giotto::get10Xmatrix(path_to_data = dataDir,
                                    gene_column_index = 2)

## Spatial coordinates
spatial_coords <- data.table::fread("data/spatialcoordinates.txt")
colnames(spatial_coords) <- c("HDMI", "Lane", "Tile", "X", "Y")

# Subset expression and spatial info by tile
spatial_coords_tile <- spatial_coords[Tile == "2104"]

expr_matrix_tile <- expr_matrix[, as.character(colnames(expr_matrix)) %in% spatial_coords_tile$HDMI]

# convert expression matrix to minimal data.table object
matrix_tile_dt <- data.table::as.data.table(Matrix::summary(expr_matrix_tile))
genes <- expr_matrix_tile@Dimnames[[1]]
samples <- expr_matrix_tile@Dimnames[[2]]
matrix_tile_dt[, gene := genes[i]]
matrix_tile_dt[, hdmi := samples[j]]

# merge data.table matrix and spatial coordinates to create input for Giotto Polygons
gpoints <- data.table::merge.data.table(matrix_tile_dt, spatial_coords_tile,
                                        by.x = "hdmi", by.y = "HDMI")

gpoints <- gpoints[, .(hdmi, X, Y, gene, x)]
colnames(gpoints) <- c("hdmi", "y", "x", "gene", "counts")

# check total counts per hdmi
gpoints_aggr <- gpoints[, sum(counts), by = .(hdmi, x, y)]
colnames(gpoints_aggr) <- c("hdmi", "x", "y", "total_counts")
data.table::setorder(gpoints_aggr, -total_counts)

pl <- ggplot() +
    geom_point(data = gpoints_aggr[total_counts > 8 & total_counts < 1000,],
               aes(x = x, y = y, color = total_counts), size = 0.5) +
    scale_color_gradient2(midpoint = 10, low = "blue", mid = "yellow", high = "red")
pl

### Use Segmentation from Qupath
geojson_path <- "data/2104_seg.geojson"
test <- terra::vect(geojson_path)

cellpolygons <- test[test$objectType != "annotation"]
plot(cellpolygons)

final_polygons <- createGiottoPolygon(cellpolygons)

original_points <- createGiottoPoints(x = gpoints[, .(x, y, gene, hdmi, counts)])
original_feat_ext <- ext(original_points@spatVector)

# convert polygon to spatRaster to change extent to that of original points
final_spatraster <- polygon_to_raster(polygon = final_polygons@spatVector)
ext(final_spatraster$raster) <- original_feat_ext
final_polygons@spatVector <- as.polygons(final_spatraster$raster)
final_polygons@spatVector$poly_ID <- final_spatraster$ID_vector[final_polygons@spatVector$poly_i]

# flip and shift
final_polygons@spatVector <- flip(final_polygons@spatVector)
yshift <- terra::ymin(original_feat_ext) - terra::ymax(original_feat_ext)
final_polygons@spatVector <- terra::shift(final_polygons@spatVector, dy = -yshift)

plot(final_polygons)


# add giotto points class
gpoints_subset <- gpoints[hdmi %in% gpoints_aggr[total_counts > 1]$hdmi]

# multiply rows with multiple counts and add jitter
gpoints_extra <- gpoints_subset[counts > 1]
gpoints_extra <- gpoints_extra[, rep(counts, counts), by = .(hdmi, gene, x, y)]
gpoints_extra <- rbind(gpoints_extra[, .(hdmi, gene, x, y)], gpoints_subset[counts == 1, .(hdmi, gene, x, y)])
jitter_x <- runif(nrow(gpoints_extra), min = -3, max = 3)
jitter_y <- runif(nrow(gpoints_extra), min = -3, max = 3)
gpoints_extra[, x := x + jitter_x]
gpoints_extra[, y := y + jitter_y]

gpoints_aggr <- gpoints_aggr[total_counts > 8 & total_counts < 1000]

seqscope <- createGiottoObjectSubcellular(
    gpoints = list(rna = gpoints_extra[, .(x, y, gene, hdmi)]),
    gpolygons = list(cell = final_polygons),
    instructions = instructions
)

# Overlap to Polygon information
seqscope <- calculateOverlapRaster(seqscope)
seqscope <- overlapToMatrix(seqscope)
seqscope <- addSpatialCentroidLocations(seqscope,
                                        poly_info = "cell")

#################################### Filtering  ################################
seqscope <- filterGiotto(
    gobject = seqscope,
    expression_threshold = 1,
    feat_det_in_min_cells = 10,
    min_det_feats_per_cell = 50
)


################################ Normalization  ################################
seqscope <- normalizeGiotto(gobject = seqscope,
                            scalefactor = 100,
                            verbose = TRUE)

################################# Statistics  ##################################
seqscope <- addStatistics(gobject = seqscope)

############################## Dimension reduction  ############################
seqscope <- calculateHVF(gobject = seqscope,
                         HVFname = "hvg_orig")

seqscope <- runPCA(
    gobject = seqscope,
    expression_values = "normalized",
    feats_to_use = "hvg_orig",
    scale_unit = TRUE,
    center = TRUE
)

################################### Clustering  ################################

seqscope <- runUMAP(seqscope,
                    dimensions_to_use = 1:100)

seqscope <- createNearestNetwork(gobject = seqscope,
                                 dimensions_to_use = 1:100,
                                 k = 10)

seqscope <- doLeidenCluster(gobject = seqscope,
                            resolution = 0.45,
                            n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = seqscope,
         cell_color = "leiden_clus",
         show_NN_network = FALSE,
         point_size = 3.5
)

spatInSituPlotPoints(seqscope,
    show_polygon = TRUE,
    polygon_color = "white",
    polygon_line_size = 0.1,
    polygon_fill = "leiden_clus",
    polygon_fill_as_factor = TRUE,
    show_plot = TRUE,
    coord_fix_ratio = TRUE,
    save_param = list(base_width = 10)
)

################################### Subset  ####################################

seqscope_sub <- subsetGiottoLocs(seqscope,
                                 x_max = 15000,
                                 x_min = 10000,
                                 y_max = 15000,
                                 y_min = 10000
)

pdt <- pDataDT(seqscope_sub)

seqscope_sub <- subsetGiotto(seqscope,
                             spat_unit = "cell",
                             feat_type = "rna",
                             cell_ids = pdt$cell_ID)

seqscope_sub <- filterGiotto(gobject = seqscope_sub,
                             expression_threshold = 2,
                             feat_det_in_min_cells = 5,
                             min_det_feats_per_cell = 10
)

fmeta <- fDataDT(seqscope_sub)
setorder(fmeta, -total_expr)

spatInSituPlotPoints(seqscope_sub,
                     show_legend = FALSE,
                     show_image = FALSE,
                     show_plot = TRUE,
                     feats = list("rna" = fmeta$feat_ID[21:30]),
                     spat_unit = "cell",
                     point_size = 1,
                     show_polygon = TRUE,
                     use_overlap = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_alpha = 0.1,
                     polygon_line_size = 0.5,
                     coord_fix_ratio = TRUE,
                     background_color = "black",
                     save_param = list(base_width = 10)
)
