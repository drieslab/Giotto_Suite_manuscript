# Vizgen Mouse Brain Receptor Map mini dataset
# This mini dataset contains polygons and transcript detection data for all 7
# different layers.


# Setup ####

library(Giotto)

GiottoUtils::package_check(
  pkg_name = c("rhdf5", "qs", "future.apply"),
  repository = c("Bioc:rhdf5", "CRAN:qs", "CRAN:future.apply")
)


## filepaths ####

## ------------------------ EDIT THESE ------------------------- ##
setwd("[PATH TO Giotto_Suite_manuscript REPO]")
data_dir <- "[PATH TO MERSCOPE MOUSE BRAIN MAP SLICE1 REP1 DATASET]"
## ------------------------------------------------------------- ##

source("scripts/SOURCE/vizgen_hdf5.R") # utils for working with H5TileProxy

# output directories
H5TP_dir <- "scripts/CREATE/OUTS/viz_mb_H5TP"
out_dir <- "scripts/CREATE/OUTS/"

# vizgen MERSCOPE dataset paths
img_dir <- file.path(data_dir, "images/")
poly_dir <- file.path(data_dir, "cell_boundaries/")

# image paths
tfs <- file.path(img_dir, "datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_images_micron_to_mosaic_pixel_transform.csv")

dapi <- file.path(img_dir, "datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_images_mosaic_DAPI_z0.tif")


## H5TileProxy ####
# Generate an H5TileProxy that provides a spatial index of where the polygonal
# information is within the Vizgen HDF5 files folder. (time consuming)
# --------------------------------------------------------------------------- #
#                               **DO NOT RUN**
#                                   start
#
#              Pre-made H5TileProxy already exists in DATA folder

future::plan(future::multisession)
options("giotto.logdir" = H5TP_dir)

fovIndexVizgenHDF5(
  poly_dir = poly_dir,
  out_dir = H5TP_dir
)

#                                    end
#                               **DO NOT RUN**
# --------------------------------------------------------------------------- #



## prepare crop extent ####
# mini extent used to subset the data for a quick example
crop_ext <- terra::ext(6700, 7200, 3100, 3600)



# load data ####

## load image ####
dapi_img <- createMerscopeLargeImage(
  image_file = dapi,
  transforms_file = tfs,
  name = "dapi_z0"
)[[1]]

# Not necessary for making a dataset, but done so that the saved object is smaller
dapi_z0_crop_img <- crop(dapi_img, crop_ext)



## load polygons ####
# load the H5TileProxy for this dataset that contains a manifest and spatial
# index of the cell boundary polygons for this dataset (provided)
mb_vizproxy <- qs::qread(file.path(H5TP_dir, "vizH5Proxy.qs"))
mb_vizproxy <- H5TPmigratePaths(mb_vizproxy, new_dir = poly_dir)






# load and pre-process data ####

## query terra SpatVectors from H5TileProxy ####
# set up tileproxy filters as a named character list
# these are used to filter on col V2 of the H5TileProxy manifest
# See mb_vizproxy@manifest
filter_list <- paste0("zIndex_", 0:6)
names(filter_list) <- paste0("z", 0:6)

# perform ext query to get the polygons ####
spatvector_poly_list <- lapply(
  filter_list,
  function(token_filter) {
    mb_vizproxy@filter$V2 <- token_filter
    H5TPqueryPolys(mb_vizproxy, dapi_img@raster_object, crop_ext)
  }
)

polys_list <- lapply(seq_along(spatvector_poly_list), function(sv_i) {
  giottoPolygon(
    spatVector = spatvector_poly_list[[sv_i]],
    name = paste0("z", (sv_i - 1)),
    unique_ID_cache = unique(spatvector_poly_list[[sv_i]]$poly_ID)
  )
})
names(polys_list) <- paste0("z", 0:6)



## load the transcript detections points data ####
tx_dt <- data.table::fread(
  file.path(data_dir, "datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_detected_transcripts_S1R1.csv")
)
tx_dt[, c("V1", "x", "y") := NULL]
# needs the x and y as colnames for downstream
data.table::setnames(tx_dt, c("global_x", "global_y", "global_z"), c("x", "y", "z"))

tx_dt_sub <- ext_query_points(
  xy = tx_dt,
  SpatRaster = dapi_img,
  crop_ext = crop_ext
)

gpoints_list <- createGiottoPoints(
  tx_dt_sub,
  feat_type = c("rna", "blanks"),
  split_keyword = list("Blank")
)





# make sure the data is all in the right place ####
plot(dapi_z0_crop_img, col = grDevices::hcl.colors(100, "viridis"))
plot(polys_list$z0, add = TRUE, col = "cyan", alpha = 0.1)
plot(sample(gpoints_list$rna@spatVector, 1e4), add = TRUE, cex = 0.3, col = "yellow", alpha = 0.3)




# create giotto object ####
g <- createGiottoObject(
  feat_info = gpoints_list,
  spatial_info = polys_list,
  largeImages = list(dapi_z0_crop_img)
)

# calculate centroids for all layers ####
g <- addSpatialCentroidLocations(
  gobject = g,
  poly_info = paste0("z", 0:6),
  provenance = as.list(paste0("z", 0:6)), # provide as vector for multiple sets
  # of spatial info to be assigned as provenance. Provide as list of character
  # vectors for when multiple centroid calculations are being performed at the
  # same time
  return_gobject = TRUE
)

# save mini7 object ####
saveGiotto(
  gobject = g,
  foldername = "mini7",
  dir = out_dir,
  method = "qs",
  overwrite = TRUE,
  image_filetype = "COG"
)
