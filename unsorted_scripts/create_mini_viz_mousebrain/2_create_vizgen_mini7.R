# Vizgen Mouse Brain Receptor Map mini dataset
# This mini dataset contains polygons and transcript detection data for all 7
# different layers.


# Setup ####

library(Giotto)

# project paths
# outputs directory
out_dir = 'DATA/'

# dataset directory
mb_dir = '/projectnb2/rd-spat/DATA/Public_data/vizgen_mouse_brain/Slice1_rep1/'
img_dir = paste0(mb_dir, 'images/')
poly_dir = paste0(mb_dir, 'cell_boundaries/')

# Load needed functions for working with Vizgen HDF5 datasets
source('SCRIPTS/SOURCE/vizgen_hdf5.R')

# load z0 dapi image from the dataset
dapi_full_mb = create_vizgen_large_image(
  image_file = paste0(img_dir, 'datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_images_mosaic_DAPI_z0.tif'),
  transforms_file = paste0(img_dir, 'datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_images_micron_to_mosaic_pixel_transform.csv')
)

# Load pre-created H5TileProxy for the mouse brain receptor map dataset
mb_vizproxy = qs::qread(file = paste0('DATA/viz_mouse_brain_h5tileproxy/vizH5Proxy.qs'))





# define subset crop extent ####
crop_ext = terra::ext(6700, 7200, 3100, 3600)




# load and pre-process data ####

## query terra SpatVectors from H5TileProxy ####
# set up tileproxy filters as a named character list
# these are used to filter on col V2 of the H5TileProxy manifest
# See mb_vizproxy@manifest
filter_list = paste0('zIndex_', 0:6)
names(filter_list) = paste0('z', 0:6)

# perform ext query to get the polygons ####
spatvector_poly_list = lapply(filter_list,
       function(token_filter) {
         mb_vizproxy@filter$V2 = token_filter
         ext_query_polys(mb_vizproxy, dapi_full_mb@raster_object, crop_ext)
       })

polys_list = lapply(seq_along(spatvector_poly_list), function(sv_i) {
    giottoPolygon(spatVector = spatvector_poly_list[[sv_i]],
                  name = paste0('z', (sv_i - 1)),
                  unique_ID_cache = unique(spatvector_poly_list[[sv_i]]$poly_ID))
  })
names(polys_list) = paste0('z', 0:6)



## crop image ####
# Not necessary for making a dataset, but done so that the saved object is smaller
dapi_z0_crop = terra::crop(dapi_full_mb@raster_object, crop_ext) 
dapi_z0_crop_gimg = createGiottoLargeImage(raster_object = dapi_z0_crop,
                                           name = 'dapi_z0',
                                           use_rast_ext = TRUE)



## load the transcript detections points data ####
tx_dt = data.table::fread(
  paste0(mb_dir, 'datasets_mouse_brain_map_BrainReceptorShowcase_Slice1_Replicate1_detected_transcripts_S1R1.csv')
)
tx_dt[, c('V1', 'x', 'y') := NULL]
# needs the x and y as colnames for downstream
data.table::setnames(tx_dt, c('global_x', 'global_y', 'global_z'), c('x', 'y', 'z'))

tx_dt_sub = ext_query_points(
  xy = tx_dt,
  SpatRaster = dapi_full_mb,
  crop_ext = crop_ext
)

gpoints_list = createGiottoPoints(
  tx_dt_sub,
  feat_type = c('rna', 'blanks'),
  split_keyword = list('Blank')
)





# make sure the data is all in the right place ####
plot(dapi_z0_crop_gimg, col = grDevices::hcl.colors(100, 'viridis'))
plot(polys_list$z0, add = TRUE, col = 'cyan', alpha = 0.1)
plot(sample(gpoints_list$rna@spatVector, 1e4), add = TRUE, cex = 0.3, col = 'yellow', alpha = 0.3)




# create giotto object ####
g = createGiottoObject(feat_info = gpoints_list,
                       spatial_info = polys_list,
                       largeImages = list(dapi_z0_crop_gimg))

# calculate centroids for all layers ####
g = addSpatialCentroidLocations(
  gobject = g,
  poly_info = paste0('z', 0:6),
  provenance = as.list(paste0('z', 0:6)), # provide as vector for multiple sets 
  # of spatial info to be assigned as provenance. Provide as list of character
  # vectors for when multiple centroid calculations are being performed at the 
  # same time
  return_gobject = TRUE
)

# save mini7 object ####
saveGiotto(gobject = g,
           foldername = 'mini7',
           dir = out_dir,
           method = 'qs',
           overwrite = TRUE,
           image_filetype = 'COG')


