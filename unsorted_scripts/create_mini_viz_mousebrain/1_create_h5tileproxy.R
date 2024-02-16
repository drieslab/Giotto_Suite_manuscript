# Generate an H5TileProxy that is provides a spatial index of where the polygonal
# information is within the Vizgen HDF5 files folder
#
#                               **DO NOT RUN**
#              Pre-made H5TileProxy already exists in DATA folder


# dataset directory
mb_dir = '/projectnb2/rd-spat/DATA/Public_data/vizgen_mouse_brain/Slice1_rep1/'

poly_dir = paste0(mb_dir, 'cell_boundaries/')
out_dir = 'DATA/viz_mouse_brain_h5tileproxy/'

# Load needed functions for working with Vizgen HDF5 datasets
source('SCRIPTS/SOURCE/vizgen_hdf5.R')

# polygons ####
# create H5TileProxy
mb_vizproxy = h5TileProxy(root = 'featuredata',
                          files = poly_dir,
                          parser = parse_array)

token_filter = list(V3 = 'p_0',
                    V2 = 'zIndex_0')
mb_vizproxy@filter = token_filter

qs::qsave(mb_vizproxy, file = paste0(out_dir, 'vizH5Proxy.qs'))
