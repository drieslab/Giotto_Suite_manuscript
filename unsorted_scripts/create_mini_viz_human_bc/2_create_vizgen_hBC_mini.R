# Vizgen Human FFPE Breast Cancer mini dataset
# This mini dataset contains both Vizgen provided cell polygon annotations and
# hand annotated nuclear annotations.


# Setup ####

# set working directory to project
setwd("PATH TO Giotto_Suite_manuscript REPO")

# path to Vizgen dataset output directory
data_dir = "PATH TO DATASET"


# Generate an H5TileProxy that provides a spatial index of where the polygonal
# information is within the Vizgen HDF5 files folder
#
#                               **DO NOT RUN**
#              Pre-made H5TileProxy already exists in DATA folder

poly_dir = file.path(data_dir, 'cell_boundaries/')
out_dir = 'aux_data/viz_human_breast_cancer_h5tileproxy'

source("aux_scripts/source/vizgen_hdf5.R")

fovIndexVizgenHDF5(
  poly_dir = poly_dir,
  out_dir = out_dir
)







# prepare crop extent ####
# mini extent used to subset the data for a quick example
crop_ext = terra::ext(5500, 6e3, 3500, 4e3)


# filepaths ####
out_dir = 'DATA/'

# dataset directory
data_dir = '/projectnb/rd-spat/DATA/Public_data/vizgen_ffpe/vizgen_human_breast_cancer/Patient1/'
img_dir = paste0(data_dir, 'images/')
poly_dir = paste0(data_dir, 'cell_boundaries/')


# image paths
tfs = paste0(img_dir, 'micron_to_mosaic_pixel_transform.csv')

b1 = paste0(img_dir, 'mosaic_Cellbound1_z0.tif')
b2 = paste0(img_dir, 'mosaic_Cellbound2_z0.tif')
b3 = paste0(img_dir, 'mosaic_Cellbound3_z0.tif')

dapi = paste0(img_dir, 'mosaic_DAPI_z0.tif')
polyT = paste0(img_dir, 'mosaic_PolyT_z0.tif')


# nuclei segmentation .geojson filepath
# Hand annotated for this mini object using QuPath. Needed for comparisons
# between cellular compartments
nuclei = 'DATA/vizgen_ffpe_BC_mini_subset/viz_mini_DAPI.geojson'





# load data ####

# load images
gimg_list <- createMerscopeLargeImage(
  image_file = c(dapi, polyT, b1, b2, b3),
  transforms_file = tfs,
  name = c("dapi_z0", "polyT_z0", "cellbound1_z0",
           "cellbound2_z0", "cellbound3_z0")
)

gimg_crop_list <- lapply(gimg_list, function(x) {
  crop(x, crop_ext)
})



# load polygons
# load [nuclear] polygons .geojson (external)
nuc_poly = createGiottoPolygonsFromGeoJSON(GeoJSON = nuclei, name = 'nucleus')
ext(nuc_poly) = terra::ext(crop_ext)
nuc_poly = terra::flip(nuc_poly, y0 = mean(c(crop_ext[3], crop_ext[4])))
nuc_poly@spatVector$objectType = NULL

# load [cell boundary] polygons (provided)
# load needed functions for working with Vizgen HDF5 datasets
source('SCRIPTS/SOURCE/vizgen_hdf5.R')
# load the H5TileProxy for this dataset that contains a manifest and spatial
# index of the cell boundary polygons for this dataset
hBC_vizproxy = qs::qread('DATA/viz_human_breast_cancer_h5tileproxy/vizH5Proxy.qs')

cell_sv = ext_query_polys(hBC_vizproxy, gimg_list[[1]]@raster_object, crop_ext)
cell_poly = giottoPolygon(spatVector = cell_sv,
                          unique_ID_cache = unique(cell_sv$poly_ID),
                          name = 'cell')



# load transcripts
tx_dt = data.table::fread(paste0(data_dir, 'detected_transcripts.csv'),
                          drop = c('x', 'y')) # remove local coordinates
data.table::setnames(tx_dt, old = c('global_x', 'global_y'), new = c('x', 'y'))
tx_sub = ext_query_points(tx_dt, gimg_list[[1]]@raster_object, crop_ext)

# create giotto points
gpoints_list = createGiottoPoints(
  tx_sub,
  feat_type = c('rna', 'blanks'),
  split_keyword = list('Blank')
)





# check spatial overlap ####
plot(gimg_crop_list[[1]]) # dapi
plot(cell_poly, add = TRUE, col = '#1E88E5', alpha = 0.3)
plot(nuc_poly, add = TRUE, col = '#FFC107', alpha = 0.3)
plot(sample(gpoints_list$rna@spatVector, 1e5), col = 'magenta', add = TRUE, alpha = 0.3, cex = 0.1)
# good overlap




# create giotto object ####
g = createGiottoObject(
  spatial_info = list(
    'cell' = cell_poly,
    'nucleus' = nuc_poly
  ),
  feat_info = gpoints_list,
  largeImages = gimg_crop_list
)


# generate spatial locations (centroids) ####
g = addSpatialCentroidLocations(
  gobject = g,
  poly_info = c('cell', 'nucleus'),
  provenance = list('cell', 'nucleus'),
  return_gobject = TRUE
)





# save the giotto object ####
saveGiotto(gobject = g,
           dir = out_dir,
           foldername = 'vizgen_ffpe_BC_mini',
           image_filetype = 'COG',
           overwrite = TRUE,
           method = 'qs')





