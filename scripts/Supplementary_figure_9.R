##%######################################################%##
#                                                          #
####               Supplementary figure 9               ####
#                                                          #
##%######################################################%##


library(Giotto)

my_seed_num = 315
set.seed(my_seed_num)
setwd("C:/Users/matth/Documents/images_presentations/Giotto_Suite_Manuscript/s9/")

## DATA FROM
#   https://data.mendeley.com/datasets/ncfgz5xxyb/1
# VIA:
#   https://www.nature.com/articles/s41592-022-01692-z#data-availability

## SCC PATH /projectnb/rd-spat/DATA/Public_data/Spatial/Multiplexing_protein/IMC/20210614_LN_panorama_test/README.txt


reload_IMC <- function(){
  results_folder = "C:/Users/matth/Documents/IMC/"
  my_python_path = NULL

  instrs = createGiottoInstructions(save_dir = results_folder,
                                    save_plot = FALSE,
                                    show_plot = TRUE,
                                    return_plot = FALSE,
                                    python_path = my_python_path)


  prefix = "C:/Users/matth/Documents/data_storage/IMC_ROI_001/"

  maskfile = paste0(prefix, "193Ir.geojson")

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
  vimenti = paste0(prefix,'Vimenti_655((3468))Pt196.tif')

  mask = terra::vect(maskfile)

  testg = createGiottoObjectSubcellular(gpolygons = list('cell' = spatVector_to_dt(mask)),
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
                                                           'vimenti' = vimenti),
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
  return(testg)
}

testg = reload_IMC()

## 1. create subset ####
# -------------------- #
testg = subsetGiottoLocsSubcellular(testg, x_max = 500, x_min = 0, y_max = 500, y_min = 0, poly_info = 'cell')
testg = addSpatialCentroidLocations(testg, poly_info = 'cell')

## 2. create rescaled polygons as another polygon information layer ####
# -------------------------------------------------------------------- #
testg = rescalePolygons(gobject = testg,
                        poly_info = 'cell',
                        name = 'smallcell',
                        fx = 0.75, fy = 0.75, calculate_centroids = T)
testg = rescalePolygons(gobject = testg,
                        poly_info = 'cell',
                        name = 'largecell',
                        fx = 1.25, fy = 1.25, calculate_centroids = T)

showGiottoSpatialInfo(testg)

## 3. calculate overlap for original (cell) and rescaled (smallcell) polygons ####
# ------------------------------------------------------------------------------ #

## cell
testg = calculateOverlapPolygonImages(testg,
                                      spatial_info = 'cell',
                                      image_names = c('MX1', 'Ki67', 'CD20', 'CD69','DNA1', 'Actin', 'FoxP3', 'MMP9', 'CXCL13',  'CD45', 'vimenti'))

testg = overlapImagesToMatrix(testg,
                              poly_info = 'cell',
                              aggr_function = 'sum', # will sum the intensity values per pixel, use 'mean' to get mean intensity
                              image_names = c('MX1', 'Ki67', 'CD20', 'CD69','DNA1', 'Actin', 'FoxP3', 'MMP9', 'CXCL13',  'CD45', 'vimenti'))


## small cell
testg = calculateOverlapPolygonImages(testg,
                                      spatial_info = 'smallcell',
                                      image_names = c('MX1', 'Ki67', 'CD20', 'CD69','DNA1', 'Actin', 'FoxP3', 'MMP9', 'CXCL13',  'CD45', 'vimenti'))

testg = overlapImagesToMatrix(testg,
                              poly_info = 'smallcell',
                              aggr_function = 'sum', # will sum the intensity values per pixel, use 'mean' to get mean intensity
                              image_names = c('MX1', 'Ki67', 'CD20', 'CD69','DNA1', 'Actin', 'FoxP3', 'MMP9', 'CXCL13',  'CD45', 'vimenti'))

## large cell
testg = calculateOverlapPolygonImages(testg,
                                      spatial_info = 'largecell',
                                      image_names = c('MX1', 'Ki67', 'CD20', 'CD69','DNA1', 'Actin', 'FoxP3', 'MMP9', 'CXCL13',  'CD45', 'vimenti'))

testg = overlapImagesToMatrix(testg,
                              poly_info = 'largecell',
                              aggr_function = 'sum', # will sum the intensity values per pixel, use 'mean' to get mean intensity
                              image_names = c('MX1', 'Ki67', 'CD20', 'CD69','DNA1', 'Actin', 'FoxP3', 'MMP9', 'CXCL13',  'CD45', 'vimenti'))


## 4. filter giotto object and add statistics ####
# --------------------------------------------- #
# check for cells with zero overlaps
testg = filterGiotto(gobject = testg,
                       spat_unit = 'cell',
                       poly_info = 'cell',
                       feat_type = 'protein',
                       expression_threshold = 0,
                       feat_det_in_min_cells = 10,
                       min_det_feats_per_cell = 2)

testg = filterGiotto(gobject = testg,
                       spat_unit = 'smallcell',
                       poly_info = 'smallcell',
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

testg = addStatistics(testg, spat_unit = 'cell', expression_values = 'raw')
testg = addStatistics(testg, spat_unit = 'smallcell', expression_values = 'raw')
testg = addStatistics(testg, spat_unit = 'largecell', expression_values = 'raw')


## 5. normalize data using simple pearson residuals ####
# ---------------------------------------------------- #
testg = normalizeGiotto(testg, spat_unit = 'cell', norm_methods = 'pearson_resid')
testg = normalizeGiotto(testg, spat_unit = 'smallcell', norm_methods = 'pearson_resid')
testg = normalizeGiotto(testg, spat_unit = 'largecell', norm_methods = 'pearson_resid')

showGiottoExpression(testg)


## 6. dimension reduction and clustering ####
# ----------------------------------------- #

## PCA
testg = runPCA(testg, spat_unit = 'cell', expression_values = 'scaled',
               scale_unit = F, center = F, ncp = 20,
               set_seed = TRUE, seed_number = my_seed_num)

testg = runPCA(testg, spat_unit = 'smallcell', expression_values = 'scaled',
               scale_unit = F, center = F, ncp = 20,
               set_seed = TRUE, seed_number = my_seed_num)

testg = runPCA(testg, spat_unit = 'largecell', expression_values = 'scaled',
               scale_unit = F, center = F, ncp = 20,
               set_seed = TRUE, seed_number = my_seed_num)



# UMAP
testg = runUMAP(gobject = testg, spat_unit = 'cell',
                dimensions_to_use = 1:5,
                set_seed = TRUE, seed_number = my_seed_num)

testg = runUMAP(gobject = testg, spat_unit = 'smallcell',
                dimensions_to_use = 1:5,
                set_seed = TRUE, seed_number = my_seed_num)

testg = runUMAP(gobject = testg, spat_unit = 'largecell',
                dimensions_to_use = 1:5,
                set_seed = TRUE, seed_number = my_seed_num)


# nearest network & clustering
testg = createNearestNetwork(testg, spat_unit = 'cell', dimensions_to_use = 1:5)
testg = createNearestNetwork(testg, spat_unit = 'smallcell', dimensions_to_use = 1:5)
testg = createNearestNetwork(testg, spat_unit = 'largecell', dimensions_to_use = 1:5)


testg = doLeidenCluster(testg,
                        spat_unit = 'cell',
                        feat_type = 'protein',
                        n_iterations = 100,
                        resolution = 0.5,
                        set_seed = TRUE,
                        seed_number = my_seed_num)

testg = doLeidenCluster(testg,
                        spat_unit = 'smallcell',
                        feat_type = 'protein',
                        n_iterations = 100,
                        resolution = 0.6,
                        set_seed = TRUE,
                        seed_number = my_seed_num)
testg = doLeidenCluster(testg,
                        spat_unit = 'largecell',
                        feat_type = 'protein',
                        n_iterations = 100,
                        resolution = 0.45,
                        set_seed = TRUE,
                        seed_number = my_seed_num)

# K-means
testg = doKmeans(testg,
                 spat_unit = "cell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = 7,
                 set_seed = TRUE,
                 seed_number = my_seed_num)

testg = doKmeans(testg,
                 spat_unit = "smallcell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = 7,
                 set_seed = TRUE,
                 seed_number = my_seed_num)

testg = doKmeans(testg,
                 spat_unit = "largecell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = 7,
                 set_seed = TRUE,
                 seed_number = my_seed_num)


## 7. Plotting and Visualizations ####
# ---------------------------------- #
my_colors = getDistinctColors(7)
names(my_colors) = unique(getCellMetadata(testg, output = "data.table")$kmeans)

mypl <- spatInSituPlotPoints(gobject = testg,
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

mypl <- mypl + ggplot2::labs(title = "Original KMeans Clusters")
mypl

ggplot2::ggsave(filename = "./s9B_original_in_situ.png",
                plot = mypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")


mypl2 <- spatInSituPlotPoints(gobject = testg,
                              spat_unit = "smallcell",
                              polygon_feat_type = "smallcell",
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

mypl2 <- mypl2 + ggplot2::labs(title = "Rescaled KMeans Clusters")
mypl2

ggplot2::ggsave(filename = "./s9B_rescaled_in_situ.png",
                plot = mypl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")

upl <- plotUMAP(testg, spat_unit = "cell",
                cell_color = "kmeans",
                point_size = 2,
                return_plot = T) +
       ggplot2::labs(title = "Original UMAP and Kmeans Clusters")
upl

ggplot2::ggsave(filename = "./s9B_original_UMAP.png",
                plot = upl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")

upl2 <- plotUMAP(testg,
                 spat_unit = "smallcell",
                 cell_color = "kmeans",
                 point_size = 2,
                 return_plot = T) +
        ggplot2::labs(title = "Rescaled UMAP and KMeans Clusters")
upl2

ggplot2::ggsave(filename = "./s9B_rescaled_UMAP.png",
                plot = upl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")




sub_c_testg = subsetGiottoLocs(testg,
                               spat_unit = "cell",
                               spat_loc_name = "raw",
                               x_max = 500,
                               x_min = 400,
                               y_max = 100,
                               y_min = 0,
                               poly_info = 'cell')
sub_sc_testg = subsetGiottoLocs(testg,
                                spat_unit = "smallcell",
                                spat_loc_name = "raw",
                                x_max = 500,
                                x_min = 400,
                                y_max = 100,
                                y_min = 0,
                                poly_info = 'smallcell')


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

smypl

ggplot2::ggsave(filename = "./s9B_original_in_situ_subset.png",
                plot = smypl,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")


smypl2 <- spatInSituPlotPoints(gobject = sub_sc_testg,
                               spat_unit = "smallcell",
                               polygon_feat_type = "smallcell",
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

smypl2

ggplot2::ggsave(filename = "./s9B_rescaled_in_situ_subset.png",
                plot = smypl2,
                device = "png",
                width = 6,
                height = 6,
                dpi = "retina")


## 8. Polygon Resizing and Visualizations ####
# ------------------------------------------ #


testg = showPolygonSizeInfluence(gobject = testg,
                                 spat_unit = "cell",
                                 alt_spat_unit = "smallcell",
                                 feat_type = "protein",
                                 clus_name = "kmeans",
                                 verbose = T)

switch_plot = showPolygonSizeInfluence(gobject = testg,
                                       spat_unit = "cell",
                                       alt_spat_unit = "smallcell",
                                       feat_type = "protein",
                                       clus_name = "kmeans",
                                       return_plot = TRUE,
                                       verbose = FALSE)

switch_plot <- switch_plot + ggplot2::labs(title = "Polygon Resizing Effects on KMeans Clusters")
switch_plot

ggplot2::ggsave(filename = "./s9C_in_situ_switches.png",
                plot = switch_plot,
                device = "png",
                width = 5,
                height = 5,
                dpi = "retina")

showCellProportionSwitchedPie(gobject = testg,
                              spat_unit = "cell",
                              feat_type = "protein")

ggplot2::ggsave(filename = "./s9D_switch_pie.png",
                plot = last_plot(),
                device = "png",
                width = 4,
                height = 4,
                dpi = "retina")

# Manual save required
showCellProportionSwitchedSanKey(gobject = testg,
                                 spat_unit = "cell",
                                 alt_spat_unit = "smallcell",
                                 feat_type = "protein")




# testg = showPolygonSizeInfluence(gobject = testg,
#                                  spat_unit = "cell",
#                                  alt_spat_unit = "largecell",
#                                  feat_type = "protein",
#                                  clus_name = "kmeans",
#                                  verbose = T)
#
# showCellProportionSwitchedPie(gobject = testg, spat_unit = "cell", feat_type = "protein")
#
# showCellProportionSwitchedSanKey(gobject = testg,
#                                  spat_unit = "cell",
#                                  alt_spat_unit = "largecell",
#                                  feat_type = "protein")

