##%######################################################%##
#                                                          #
####               Supplementary figure 9               ####
#                                                          #
##%######################################################%##


library(Giotto)

# set seed for this analysis
my_seed_num = 315
set.seed(my_seed_num)

# create instructions
results_folder = "/Users/rubendries/Google Drive/Shared drives/Dries_lab/Project/Lab/Giotto_dev/testprojects/Ruben/GiottoSuite_manuscript_IMC/"
my_python_path = NULL
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  return_plot = FALSE,
                                  python_path = my_python_path)


## 1. create subcellular Giotto object ####
# --------------------------------------- #

# create polygons from geojson file
prefix = "/Users/rubendries/Google Drive/Shared drives/Dries_lab/Datasets/Public_data/Spatial/Multiplexing_protein/IMC/20210614_LN_panorama_test_wd/TIFF_ROI_001 (ROIs, 20210614_LN_panorama_test.mcd)/"
maskfile = paste0(prefix, "193Ir.geojson")
polygons = createGiottoPolygonsFromGeoJSON(GeoJSON = maskfile, calc_centroids = T)

# image input files
MX1 = paste0(prefix,'MX1_2375((2706))Yb172.tif')
Ki67 = paste0(prefix,'Ki-67_142((3020))Pt198.tif')
CD20 = paste0(prefix,'CD20_36((3317))Nd148.tif')
CD69 = paste0(prefix,'CD69_2257((2445))Gd156.tif')
DNA1 = paste0(prefix,'191Ir.tif')
Actin = paste0(prefix,'Cleaved_93((346))Lu175.tif')
FoxP3 = paste0(prefix,'FOXP3_115((3001))Dy163.tif')
MMP9 = paste0(prefix,'MMP9_2241((2912))Gd158.tif')
CXCL13 = paste0(prefix,'CXCL13 _2316((2458))Dy164.tif')
CD45 = paste0(prefix,'CD45RO_2014((3019))Dy162.tif')
vimentin = paste0(prefix,'Vimenti_655((3468))Pt196.tif')



testg = createGiottoObjectSubcellular(gpolygons = list('cell' = polygons),
                                      polygon_mask_list_params = list(flip_vertical = TRUE,
                                                                      flip_horizontal = FALSE,
                                                                      shift_vertical_step = TRUE,
                                                                      shift_horizontal_step = FALSE),
                                      largeImages = list('MX1' = MX1,
                                                         'Ki67' = Ki67,
                                                         'CD20' = CD20,
                                                         'CD69' = CD69,
                                                         'DNA1' = DNA1,
                                                         'Actin' = Actin,
                                                         'FoxP3' = FoxP3,
                                                         'MMP9' = MMP9,
                                                         'CXCL13' = CXCL13,
                                                         'CD45' = CD45,
                                                         'vimentin' = vimentin),
                                      largeImages_list_params = list(negative_y = FALSE,
                                                                     extent = NULL,
                                                                     use_rast_ext = FALSE,
                                                                     image_transformations = NULL,
                                                                     xmax_bound = NULL,
                                                                     xmin_bound = NULL,
                                                                     ymax_bound = NULL,
                                                                     ymin_bound = NULL,
                                                                     scale_factor = 1,
                                                                     verbose = TRUE),
                                      instructions = instrs)



## 2. create rescaled polygons as another polygon information layer ####
# -------------------------------------------------------------------- #
testg = rescalePolygons(gobject = testg,
                        poly_info = 'cell',
                        name = 'smallcell',
                        fx = 0.75, fy = 0.75,
                        calculate_centroids = T)

testg = rescalePolygons(gobject = testg,
                        poly_info = 'cell',
                        name = 'largecell',
                        fx = 1.25, fy = 1.25,
                        calculate_centroids = T)

showGiottoSpatialInfo(testg)






## 3. calculate overlap for original (cell) and rescaled (smallcell) polygons ####
# ------------------------------------------------------------------------------ #

## cell
testg = calculateOverlapPolygonImages(testg,
                                      spatial_info = 'cell',
                                      image_names = c('MX1', 'Ki67', 'CD20', 'CD69',
                                                      'DNA1', 'Actin', 'FoxP3', 'MMP9',
                                                      'CXCL13',  'CD45', 'vimentin'))

testg = overlapImagesToMatrix(testg,
                              poly_info = 'cell',
                              aggr_function = 'sum', # will sum the intensity values per pixel, use 'mean' to get mean intensity
                              image_names = c('MX1', 'Ki67', 'CD20', 'CD69',
                                              'DNA1', 'Actin', 'FoxP3', 'MMP9',
                                              'CXCL13',  'CD45', 'vimentin'))


## small cell
testg = calculateOverlapPolygonImages(testg,
                                      spatial_info = 'smallcell',
                                      image_names = c('MX1', 'Ki67', 'CD20', 'CD69',
                                                      'DNA1', 'Actin', 'FoxP3', 'MMP9',
                                                      'CXCL13',  'CD45', 'vimentin'))

testg = overlapImagesToMatrix(testg,
                              poly_info = 'smallcell',
                              aggr_function = 'sum', # will sum the intensity values per pixel, use 'mean' to get mean intensity
                              image_names = c('MX1', 'Ki67', 'CD20', 'CD69',
                                              'DNA1', 'Actin', 'FoxP3', 'MMP9',
                                              'CXCL13',  'CD45', 'vimentin'))

## large cell
testg = calculateOverlapPolygonImages(testg,
                                      spatial_info = 'largecell',
                                      image_names = c('MX1', 'Ki67', 'CD20', 'CD69',
                                                      'DNA1', 'Actin', 'FoxP3', 'MMP9',
                                                      'CXCL13',  'CD45', 'vimentin'))

testg = overlapImagesToMatrix(testg,
                              poly_info = 'largecell',
                              aggr_function = 'sum', # will sum the intensity values per pixel, use 'mean' to get mean intensity
                              image_names = c('MX1', 'Ki67', 'CD20', 'CD69',
                                              'DNA1', 'Actin', 'FoxP3', 'MMP9',
                                              'CXCL13',  'CD45', 'vimentin'))



## plot raw data images

# for spatial unit = cell
spatInSituPlotPoints(testg,
                     show_image = T,
                     largeImage_name = 'DNA1',
                     spat_unit = 'cell',
                     polygon_feat_type = "cell",
                     show_polygon = T,
                     polygon_color = 'red',
                     polygon_alpha = 0,
                     polygon_line_size = 0.1,
                     polygon_fill = NULL,
                     feat_type = "protein",
                     feats = NULL)

spatInSituPlotPoints(testg,
                     show_image = T,
                     largeImage_name = 'vimentin',
                     spat_unit = 'cell',
                     polygon_feat_type = "cell",
                     show_polygon = T,
                     polygon_color = 'red',
                     polygon_alpha = 0,
                     polygon_line_size = 0.1,
                     polygon_fill = NULL,
                     feat_type = "protein",
                     feats = NULL)


# for spatial unit = small cell
spatInSituPlotPoints(testg,
                     show_image = T,
                     largeImage_name = 'vimentin',
                     spat_unit = 'smallcell',
                     polygon_feat_type = "smallcell",
                     show_polygon = T,
                     polygon_color = 'red',
                     polygon_alpha = 0,
                     polygon_line_size = 0.1,
                     polygon_fill = NULL,
                     feat_type = "protein",
                     feats = NULL)

spatInSituPlotPoints(testg,
                     show_image = T,
                     largeImage_name = 'vimentin',
                     spat_unit = 'largecell',
                     polygon_feat_type = "largecell",
                     show_polygon = T,
                     polygon_color = 'red',
                     polygon_alpha = 0,
                     polygon_line_size = 0.1,
                     polygon_fill = NULL,
                     feat_type = "protein",
                     feats = NULL)

## 4. filter giotto object and add statistics ####
# --------------------------------------------- #
# check for cells with zero overlaps
testg = filterGiotto(gobject = testg,
                     spat_unit = 'smallcell',
                     poly_info = 'smallcell',
                     feat_type = 'protein',
                     expression_threshold = 0,
                     feat_det_in_min_cells = 10,
                     min_det_feats_per_cell = 2)

testg = filterGiotto(gobject = testg,
                     spat_unit = 'cell',
                     poly_info = 'cell',
                     feat_type = 'protein',
                     expression_threshold = 0,
                     feat_det_in_min_cells = 10,
                     min_det_feats_per_cell = 2)

testg = filterGiotto(gobject = testg,
                     spat_unit = 'largecell',
                     poly_info = 'largecell',
                     feat_type = 'protein',
                     expression_threshold = 0,
                     feat_det_in_min_cells = 10,
                     min_det_feats_per_cell = 2)

testg = addStatistics(testg, spat_unit = 'smallcell', expression_values = 'raw')
testg = addStatistics(testg, spat_unit = 'cell', expression_values = 'raw')
testg = addStatistics(testg, spat_unit = 'largecell', expression_values = 'raw')


## 5. normalize data using simple pearson residuals ####
# ---------------------------------------------------- #

testg = normalizeGiotto(testg, spat_unit = 'smallcell', norm_methods = 'pearson_resid')
testg = normalizeGiotto(testg, spat_unit = 'cell', norm_methods = 'pearson_resid')
testg = normalizeGiotto(testg, spat_unit = 'largecell', norm_methods = 'pearson_resid')

showGiottoExpression(testg)


## 6. dimension reduction and clustering ####
# ----------------------------------------- #

## PCA
testg = runPCA(testg, spat_unit = 'smallcell', expression_values = 'scaled',
               scale_unit = F, center = F, ncp = 20,
               set_seed = TRUE, seed_number = my_seed_num)

testg = runPCA(testg, spat_unit = 'cell', expression_values = 'scaled',
               scale_unit = F, center = F, ncp = 20,
               set_seed = TRUE, seed_number = my_seed_num)

testg = runPCA(testg, spat_unit = 'largecell', expression_values = 'scaled',
               scale_unit = F, center = F, ncp = 20,
               set_seed = TRUE, seed_number = my_seed_num)



# UMAP
testg = runUMAP(gobject = testg, spat_unit = 'smallcell',
                dimensions_to_use = 1:5,
                set_seed = TRUE, seed_number = my_seed_num)

testg = runUMAP(gobject = testg, spat_unit = 'cell',
                dimensions_to_use = 1:5,
                set_seed = TRUE, seed_number = my_seed_num)

testg = runUMAP(gobject = testg, spat_unit = 'largecell',
                dimensions_to_use = 1:5,
                set_seed = TRUE, seed_number = my_seed_num)


# nearest network & clustering
testg = createNearestNetwork(testg, spat_unit = 'smallcell', dimensions_to_use = 1:5)
testg = createNearestNetwork(testg, spat_unit = 'cell', dimensions_to_use = 1:5)
testg = createNearestNetwork(testg, spat_unit = 'largecell', dimensions_to_use = 1:5)

# K-means
testg = doKmeans(testg,
                 spat_unit = "smallcell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = 5,
                 iter_max = 1000,
                 nstart = 1000,
                 set_seed = TRUE,
                 seed_number = my_seed_num)

testg = doKmeans(testg,
                 spat_unit = "cell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = 5,
                 iter_max = 1000,
                 nstart = 1000,
                 set_seed = TRUE,
                 seed_number = my_seed_num)

testg = doKmeans(testg,
                 spat_unit = "largecell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = 5,
                 iter_max = 1000,
                 nstart = 1000,
                 set_seed = TRUE,
                 seed_number = my_seed_num)


## 7. compare k-means clustering results for spatial units cell and smallcell ####
# ------------------------------------------------------------------------------ #
library(data.table)

# combine metadata from spatial units 'cell' and 'smallcell'
small_cell_meta = getCellMetadata(testg, spat_unit = 'smallcell', output = 'data.table')
setnames(small_cell_meta, old = 'kmeans', new = 'km_smallcell')

cell_meta = getCellMetadata(testg, spat_unit = 'cell', output = 'data.table')
setnames(cell_meta, old = 'kmeans', new = 'km_cell')

comb_meta = data.table::merge.data.table(small_cell_meta[,.(cell_ID, km_smallcell)],
                                         cell_meta[,.(cell_ID, km_cell)], by = 'cell_ID')

# create quantitative tables for visual inspection
counttable = table(comb_meta$km_smallcell, comb_meta$km_cell)
proptable = prop.table(counttable)



# identifies most stable cluster combinations
comb_meta[, km_comp := paste0(km_smallcell, '-', km_cell), by = .I]
summary_comb_meta = comb_meta[, .N, by = .(km_cell, km_comp)]
summary_comb_meta[, max_n := max(N), by = .(km_cell)]
stable_combos = summary_comb_meta[N == max_n][['km_comp']]

# annotate cell IDs that are stable or switch cluster
comb_meta[, km_stable := ifelse(km_comp %in% stable_combos, 'stable', 'switch')]
comb_meta_stable = comb_meta[, .(cell_ID, km_stable)]

# add metadata information back for spatial units 'cell' and 'smallcell'
testg = addCellMetadata(testg, new_metadata = comb_meta_stable,
                        by_column = T, column_cell_ID = 'cell_ID',
                        spat_unit = 'smallcell')
testg = addCellMetadata(testg, new_metadata = comb_meta_stable,
                        by_column = T, column_cell_ID = 'cell_ID',
                        spat_unit = 'cell')

pDataDT(testg)


## 7. Plotting and Visualizations ####
# ---------------------------------- #

## * 7.1. plot original cells ####

# definte cell colors
my_colors = getDistinctColors(5)
order_cell = as.numeric(unlist(lapply(stable_combos, FUN = function(x) {
  strsplit(x, split = '-')[[1]][2]
})))
names(my_colors) = order_cell

upl <- plotUMAP(testg,
                spat_unit = "cell",
                cell_color = "kmeans",
                cell_color_code = my_colors,
                point_size = 2,
                return_plot = T) +
  ggplot2::labs(title = "Original UMAP and Kmeans Clusters")
upl

ggplot2::ggsave(filename = paste0(results_folder, "/s9_umap_cell.png"),
                plot = upl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")

# spatplot kmeans
mypl <- spatInSituPlotPoints(gobject = testg,
                             spat_unit = "cell",
                             polygon_feat_type = "cell",
                             show_polygon = T,
                             feat_type = "protein",
                             feats = NULL,
                             polygon_fill = 'kmeans',
                             polygon_fill_code = my_colors,
                             polygon_fill_as_factor = TRUE,
                             polygon_line_size = 0.05,
                             polygon_color = 'white',
                             coord_fix_ratio = 1,
                             return_plot = T,
                             background_color = "black")
mypl <- mypl + ggplot2::labs(title = "Original KMeans Clusters")
mypl

ggplot2::ggsave(filename = paste0(results_folder, "/s9_spatinsitu_cell.png"),
                plot = mypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")

# spatplot switch
mypl <- spatInSituPlotPoints(gobject = testg,
                             spat_unit = "cell",
                             polygon_feat_type = "cell",
                             show_polygon = T,
                             feat_type = "protein",
                             feats = NULL,
                             polygon_fill = 'km_stable',
                             polygon_fill_code = c('gray', 'red'),
                             polygon_fill_as_factor = TRUE,
                             polygon_line_size = 0.05,
                             polygon_color = 'white',
                             coord_fix_ratio = 1,
                             return_plot = T,
                             background_color = "black")

mypl <- mypl + ggplot2::labs(title = "Switches cell")
mypl

ggplot2::ggsave(filename = paste0(results_folder, "/s9_spatinsitu_cell_switches.png"),
                plot = mypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")



## * 7.2. plot smaller cells ####

# definte smallcell colors
my_colors_small = getDistinctColors(5)
order_small_cell = as.numeric(unlist(lapply(stable_combos, FUN = function(x) {
  strsplit(x, split = '-')[[1]][1]
})))
names(my_colors_small) = order_small_cell

# umap
upl2 <- plotUMAP(testg,
                 spat_unit = "smallcell",
                 cell_color = "kmeans",
                 cell_color_code = my_colors_small,
                 point_size = 2,
                 return_plot = T) +
  ggplot2::labs(title = "Rescaled UMAP and KMeans Clusters")
upl2

ggplot2::ggsave(filename = paste0(results_folder, "/s9_umap_smallcell.png"),
                plot = upl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")

# kmeans spatplot
mypl2 <- spatInSituPlotPoints(gobject = testg,
                              spat_unit = "smallcell",
                              polygon_feat_type = "smallcell",
                              show_polygon = T,
                              feat_type = "protein",
                              feats = NULL,
                              polygon_fill = 'kmeans',
                              polygon_fill_code = my_colors_small,
                              polygon_fill_as_factor = TRUE,
                              polygon_line_size = 0.05,
                              polygon_color = 'white',
                              coord_fix_ratio = 1,
                              return_plot = T,
                              background_color = "black")

mypl2 <- mypl2 + ggplot2::labs(title = "Rescaled KMeans Clusters")
mypl2

ggplot2::ggsave(filename = paste0(results_folder, "/s9_spatinsitu_smallcell.png"),
                plot = mypl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")


# switch spatplot
mypl2 <- spatInSituPlotPoints(gobject = testg,
                              spat_unit = "smallcell",
                              polygon_feat_type = "smallcell",
                              show_polygon = T,
                              feat_type = "protein",
                              feats = NULL,
                              polygon_fill = 'km_stable',
                              polygon_fill_code = c('gray', 'red'),
                              polygon_fill_as_factor = TRUE,
                              polygon_line_size = 0.05,
                              polygon_color = 'white',
                              coord_fix_ratio = 1,
                              return_plot = T,
                              background_color = "black")
mypl2 <- mypl2 + ggplot2::labs(title = "Switches smallcell")
mypl2

ggplot2::ggsave(filename = paste0(results_folder, "/s9_spatinsitu_smallcell_switches.png"),
                plot = mypl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")




## 8. Identify spatial network neighbors ####
# ----------------------------------------- #

# create spatial network
testg = createSpatialNetwork(testg)
showGiottoSpatNetworks(testg)

# cell cell IDs from cells that switched cluster
cell_metadata = getCellMetadata(testg,
                                spat_unit = 'cell',
                                output = 'data.table')
switch_cell_IDs = cell_metadata[km_stable == 'switch'][['cell_ID']]

# find neighbor cell IDs
find_switch_neighbors = findNetworkNeighbors(gobject = testg, spat_unit = 'cell',
                                             spatial_network_name = 'Delaunay_network',
                                             source_cell_ids = switch_cell_IDs)

# add to metadata for cell
testg = addCellMetadata(testg, spat_unit = 'cell',
                        new_metadata = find_switch_neighbors,
                        by_column = T, column_cell_ID = 'cell_ID')

# add to metadata for small cell
testg = addCellMetadata(testg, spat_unit = 'smallcell',
                        new_metadata = find_switch_neighbors,
                        by_column = T, column_cell_ID = 'cell_ID')

# visualize
switch_nb_colors = c('gray', 'orange', 'blue', 'red');
names(switch_nb_colors) = c('others', 'neighbor', 'source', 'both')


# cell
mypl <- spatInSituPlotPoints(gobject = testg,
                             spat_unit = "cell",
                             polygon_feat_type = "cell",
                             show_polygon = T,
                             feat_type = "protein",
                             feats = NULL,
                             polygon_fill = 'nb_cells',
                             polygon_fill_code = switch_nb_colors,
                             polygon_fill_as_factor = TRUE,
                             polygon_line_size = 0.05,
                             polygon_color = 'white',
                             coord_fix_ratio = 1,
                             return_plot = T,
                             background_color = "black")

mypl <- mypl + ggplot2::labs(title = "neighbors cell")
mypl

ggplot2::ggsave(filename = paste0(results_folder, "/s9_spatinsitu_cell_nbs.png"),
                plot = mypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")


# small cell
mypl <- spatInSituPlotPoints(gobject = testg,
                             spat_unit = "smallcell",
                             polygon_feat_type = "smallcell",
                             show_polygon = T,
                             feat_type = "protein",
                             feats = NULL,
                             polygon_fill = 'nb_cells',
                             polygon_fill_code = switch_nb_colors,
                             polygon_fill_as_factor = TRUE,
                             polygon_line_size = 0.05,
                             polygon_color = 'white',
                             coord_fix_ratio = 1,
                             return_plot = T,
                             background_color = "black")

mypl <- mypl + ggplot2::labs(title = "neighbors cell")
mypl

ggplot2::ggsave(filename = paste0(results_folder, "/s9_spatinsitu_smallcell_nbs.png"),
                plot = mypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")



## barplot
cell_metadata = getCellMetadata(testg,
                                spat_unit = 'cell',
                                output = 'data.table')
cell_metadata[, nb_cells := factor(nb_cells, levels = c('source', 'both', 'neighbor', 'others'))]
cell_metadata[, (.N/1679)*100, by = .(nb_cells)]

# neighbors
pl = ggplot()
pl = pl + geom_bar(data = cell_metadata, aes(x = '', fill = nb_cells),
                   position = 'fill')
pl = pl + ggplot2::scale_fill_manual(values = mycolors)
pl = pl + theme_classic()
pl

ggplot2::ggsave(filename = paste0(results_folder, "/s9_barplot_nbs.png"),
                plot = pl,
                device = "png",
                width = 3,
                height = 6,
                dpi = "retina")

# switches
cell_metadata[, km_stable := factor(km_stable, levels = c('switch', 'stable'))]

pl = ggplot()
pl = pl + geom_bar(data = cell_metadata, aes(x = '', fill = km_stable),
                   position = 'fill')
pl = pl + ggplot2::scale_fill_manual(values = c(stable = 'gray', switch = 'red'))
pl = pl + theme_classic()
pl

ggplot2::ggsave(filename = paste0(results_folder, "/s9_barplot_switches.png"),
                plot = pl,
                device = "png",
                width = 3,
                height = 6,
                dpi = "retina")


## 9. zoomed-in regions ####
# ------------------------ #

## for spatial unit = smallcell
sub_sc_testg = subsetGiottoLocs(testg,
                                spat_unit = "smallcell",
                                spat_loc_name = "raw",
                                x_max = 500,
                                x_min = 400,
                                y_max = 100,
                                y_min = 0,
                                poly_info = 'smallcell')

# kmeans
smypl2 <- spatInSituPlotPoints(gobject = sub_sc_testg,
                               spat_unit = "smallcell",
                               polygon_feat_type = "smallcell",
                               show_polygon = T,
                               feat_type = "protein",
                               feats = NULL,
                               polygon_fill = 'kmeans',
                               polygon_fill_code = my_colors_small,
                               polygon_fill_as_factor = TRUE,
                               polygon_line_size = 0.2,
                               polygon_color = 'white',
                               coord_fix_ratio = 1,
                               return_plot = T,
                               background_color = "black")

smypl2 <- smypl2 + ggplot2::labs(title = "kmeans smallcell")
smypl2

ggplot2::ggsave(filename = paste0(results_folder, "./s9_spatinsitu_subset_smallcell_kmeans.png"),
                plot = smypl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")

# switches
smypl2 <- spatInSituPlotPoints(gobject = sub_sc_testg,
                               spat_unit = "smallcell",
                               polygon_feat_type = "smallcell",
                               show_polygon = T,
                               feat_type = "protein",
                               feats = NULL,
                               polygon_fill = 'km_stable',
                               polygon_fill_code = c('gray', 'red'),
                               polygon_fill_as_factor = TRUE,
                               polygon_line_size = 0.2,
                               polygon_color = 'white',
                               coord_fix_ratio = 1,
                               return_plot = T,
                               background_color = "black")

smypl2 <- smypl2 + ggplot2::labs(title = "switches smallcell")
smypl2

ggplot2::ggsave(filename = paste0(results_folder, "./s9_spatinsitu_subset_smallcell_switches.png"),
                plot = smypl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")


# neighbors
smypl2 <- spatInSituPlotPoints(gobject = sub_sc_testg,
                               spat_unit = "smallcell",
                               polygon_feat_type = "smallcell",
                               show_polygon = T,
                               feat_type = "protein",
                               feats = NULL,
                               polygon_fill = 'nb_cells',
                               polygon_fill_code = switch_nb_colors,
                               polygon_fill_as_factor = TRUE,
                               polygon_line_size = 0.2,
                               polygon_color = 'white',
                               coord_fix_ratio = 1,
                               return_plot = T,
                               background_color = "black")

smypl2 <- smypl2 + ggplot2::labs(title = "neighbors smallcell")
smypl2

ggplot2::ggsave(filename = paste0(results_folder, "./s9_spatinsitu_subset_smallcell_neighbors.png"),
                plot = smypl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")






## for spatial unit = cell
sub_c_testg = subsetGiottoLocs(testg,
                               spat_unit = "cell",
                               spat_loc_name = "raw",
                               x_max = 500,
                               x_min = 400,
                               y_max = 100,
                               y_min = 0,
                               poly_info = 'cell')

smypl <- spatInSituPlotPoints(gobject = sub_c_testg,
                              spat_unit = "cell",
                              polygon_feat_type = "cell",
                              show_polygon = T,
                              feat_type = "protein",
                              feats = NULL,
                              polygon_fill = 'kmeans',
                              polygon_fill_code = my_colors,
                              polygon_fill_as_factor = TRUE,
                              polygon_line_size = 0.2,
                              polygon_color = 'white',
                              coord_fix_ratio = 1,
                              return_plot = T,
                              background_color = "black")

smypl <- smypl2 + ggplot2::labs(title = "kmeans cell")
smypl

ggplot2::ggsave(filename = paste0(results_folder, "./s9_spatinsitu_subset_cell_kmeans.png"),
                plot = smypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")


# switches
smypl <- spatInSituPlotPoints(gobject = sub_sc_testg,
                              spat_unit = "cell",
                              polygon_feat_type = "cell",
                              show_polygon = T,
                              feat_type = "protein",
                              feats = NULL,
                              polygon_fill = 'km_stable',
                              polygon_fill_code = c('gray', 'red'),
                              polygon_fill_as_factor = TRUE,
                              polygon_line_size = 0.2,
                              polygon_color = 'white',
                              coord_fix_ratio = 1,
                              return_plot = T,
                              background_color = "black")

smypl <- smypl + ggplot2::labs(title = "switches cell")
smypl

ggplot2::ggsave(filename = paste0(results_folder, "./s9_spatinsitu_subset_cell_switches.png"),
                plot = smypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")

# neighbors
smypl <- spatInSituPlotPoints(gobject = sub_c_testg,
                              spat_unit = "cell",
                              polygon_feat_type = "cell",
                              show_polygon = T,
                              feat_type = "protein",
                              feats = NULL,
                              polygon_fill = 'nb_cells',
                              polygon_fill_code = switch_nb_colors,
                              polygon_fill_as_factor = TRUE,
                              polygon_line_size = 0.2,
                              polygon_color = 'white',
                              coord_fix_ratio = 1,
                              return_plot = T,
                              background_color = "black")


smypl <- smypl2 + ggplot2::labs(title = "neighbors cell")
smypl

ggplot2::ggsave(filename = paste0(results_folder, "./s9_spatinsitu_subset_cell_neighbors.png"),
                plot = smypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")





## 10. Sankey plot ####
# ------------------- #

alluvial_dt = comb_meta[, .N, by = .(km_cell, km_smallcell)]
alluvial_dt[, km_cell := as.factor(km_cell)]
alluvial_dt[, km_smallcell := as.factor(km_smallcell)]

library(ggalluvial)
# see https://corybrunson.github.io/ggalluvial/
# see https://cheatography.com/seleven/cheat-sheets/ggalluvial/

pl = ggplot(data = alluvial_dt, aes(axis1 = km_cell, axis2 = km_smallcell, y = N))
pl = pl + geom_alluvium(aes(fill = km_cell))
pl = pl + geom_stratum()
pl = pl + geom_text(stat = "stratum", aes(label = after_stat(stratum)))
pl = pl + scale_fill_manual(values = my_colors)
pl = pl + scale_x_discrete(limits = c("km_cell", "km_smallcell"), expand = c(.1, .05))
pl = pl + theme_classic()
pl


ggplot2::ggsave(filename = paste0(results_folder, "./s9_sankeyplot.png"),
                plot = pl,
                device = "png",
                width = 8,
                height = 4,
                dpi = "retina")




## 11. niche clustering ####
# ------------------------ #
testg = calculateSpatCellMetadataProportions(gobject = testg,
                                             spat_unit = 'cell',
                                             spat_network = 'Delaunay_network',
                                             metadata_column = 'kmeans', name = 'km_niche')
showGiottoSpatEnrichments(testg)

# get spatial enrichment data
prop_table = getSpatialEnrichment(testg, name = 'km_niche', output = 'data.table')
prop_matrix = dt_to_matrix(prop_table)

# perform kmeans on spatial enrichment data
prop_kmeans = kmeans(x = prop_matrix, centers = 5, iter.max = 1000, nstart = 100)
prop_kmeansDT = data.table(cell_ID = names(prop_kmeans$cluster), niche = prop_kmeans$cluster)

# add results back to giotto object
testg = addCellMetadata(testg, new_metadata = prop_kmeansDT,
                        by_column = T, column_cell_ID = 'cell_ID')

# visualize nich clustering
mypl <- spatInSituPlotPoints(gobject = testg,
                             spat_unit = "cell",
                             polygon_feat_type = "cell",
                             show_polygon = T,
                             feat_type = "protein",
                             feats = NULL,
                             polygon_fill = 'niche',
                             polygon_fill_as_factor = TRUE,
                             polygon_line_size = 0.05,
                             polygon_color = 'white',
                             coord_fix_ratio = 1,
                             return_plot = T,
                             background_color = "black")


mypl <- mypl + ggplot2::labs(title = "niche clustering on cell kmeans")
mypl

ggplot2::ggsave(filename = paste0(results_folder, "./s9_spatinsitu_cell_niches.png"),
                plot = mypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")