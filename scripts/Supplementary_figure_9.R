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
results_folder = "scripts/FIGURE/S9/"
dir.create(results_folder, showWarnings = FALSE)
my_python_path = NULL
instrs = createGiottoInstructions(
  save_dir = results_folder,
  save_plot = TRUE,
  show_plot = TRUE,
  return_plot = FALSE,
  python_path = my_python_path
)


## 1. create subcellular Giotto object ####
# --------------------------------------- #

# create polygons from geojson file
img_out = "scripts/CREATE/OUTS/imc_tiff/"
maskfile = "scripts/CREATE/EXT_DATA/193Ir.geojson" # (pre-made)
polygons = createGiottoPolygonsFromGeoJSON(
  GeoJSON = maskfile, 
  calc_centroids = TRUE
)
# spatially orient polygons to fit spatial information
e <- ext(polygons)
polygons <- flip(polygons)
ext(polygons) <- e

DNA1 = paste0(img_out,'191Ir.tif')

# image input files
img_list <- list.files(img_out, full.names = TRUE)
# further panel info here:
# https://data.mendeley.com/datasets/ncfgz5xxyb/1/files/76410bf6-5b40-4664-9eca-0a8f9c594e81
names(img_list) <- c(
  "xenon1",      "xenon2",          "DNA1",       "DNA2",             "AICDA",
  "Caveolin-1",  "CCL19",           "CCL21",      "CCR7",             "CD163",
  "CD20",        "CD209",           "CD274",      "CD278",            "CD3",
  "CD303",       "CD31",            "CD40L",      "CD45RO",           "CD56",
  "CD69",        "CD8a",            "CD9",        "Cleaved_Caspase3", "CXCL13",
  "FOXP3",       "TNFRSF18",        "Granzyme B", "HIF-1a",           "H3K9Ac",
  "HLA-DR",      "IgM",             "Ki67",       "LEF1",             "MMP9",
  "MX1",         "Myeloperoxidase", "SMA",        "Vimentin",         "XBP1"
)


testg = createGiottoObjectSubcellular(
  gpolygons = list('cell' = polygons),
  largeImages = img_list[4:length(img_list)], # everything but xenon & DNA1 channels
  largeImages_list_params = list(
    negative_y = FALSE,
    extent = NULL,
    use_rast_ext = FALSE,
    image_transformations = NULL,
    xmax_bound = NULL,
    xmin_bound = NULL,
    ymax_bound = NULL,
    ymin_bound = NULL,
    scale_factor = 1,
    verbose = TRUE
  ),
  instructions = instrs
)



## 2. create rescaled polygons as another polygon information layer ####
# -------------------------------------------------------------------- #
testg = rescalePolygons(
  gobject = testg,
  poly_info = 'cell',
  name = 'smallcell',
  fx = 0.75, fy = 0.75,
  calculate_centroids = TRUE
)

# larger cells are also possible.
testg = rescalePolygons(
  gobject = testg,
  poly_info = 'cell',
  name = 'largecell',
  fx = 1.25, fy = 1.25,
  calculate_centroids = TRUE
)

showGiottoSpatialInfo(testg)






## 3. calculate overlap for original (cell) and rescaled (smallcell) polygons ####
# ------------------------------------------------------------------------------ #

## cell
testg = calculateOverlap(
  testg,
  spatial_info = 'cell',
  name_overlap = "protein",
  image_names = names(img_list)[5:length(img_list)] # ignore DNA channel
)

testg = overlapToMatrix(
  testg,
  poly_info = 'cell',
  feat_info = "protein",
  type = "intensity",
  aggr_function = 'sum', # will sum the intensity values per pixel, use 'mean' to get mean intensity
)


## small cell
testg = calculateOverlap(
  testg,
  spatial_info = 'smallcell',
  image_names = names(img_list)[5:length(img_list)] # ignore DNA channel
)

testg = overlapToMatrix(
  testg,
  poly_info = 'smallcell',
  feat_info = "protein",
  type = "intensity",
  aggr_function = 'sum' # will sum the intensity values per pixel, use 'mean' to get mean intensity
)

## large cell
testg = calculateOverlap(
  testg,
  spatial_info = 'largecell',
  image_names = names(img_list)[5:length(img_list)] # ignore DNA channel
)

testg = overlapToMatrix(
  testg,
  poly_info = 'largecell',
  feat_info = "protein",
  type = "intensity",
  aggr_function = 'sum' # will sum the intensity values per pixel, use 'mean' to get mean intensity
)

showGiottoExpression(testg)


## plot raw data images

# for spatial unit = cell, plot DNA2 and the segmentations that were created
# based on it in QuPath
svg(file.path(results_folder, "cell_overlap.svg"))
plot(testg@largeImages$DNA2, max_intensity = 30)
plot(testg@spatial_info$cell, add = TRUE, border = "red", lwd = 0.7)
dev.off()

spatInSituPlotPoints(
  testg,
  show_image = T,
  largeImage_name = 'Vimentin',
  spat_unit = 'cell',
  polygon_feat_type = "cell",
  show_polygon = T,
  polygon_color = 'red',
  polygon_alpha = 1,
  polygon_line_size = 0.1,
  polygon_fill = NULL,
  feat_type = "protein",
  feats = NULL,
  save_param = list(
    save_name = "cell_vimentin"
  )
)


# for spatial unit = small cell
spatInSituPlotPoints(
  testg,
  show_image = T,
  largeImage_name = 'Vimentin',
  spat_unit = 'smallcell',
  polygon_feat_type = "smallcell",
  show_polygon = T,
  polygon_color = 'red',
  polygon_alpha = 1,
  polygon_line_size = 0.1,
  polygon_fill = NULL,
  feat_type = "protein",
  feats = NULL,
  save_param = list(
    save_name = "smallcell_vimentin"
  )
)

spatInSituPlotPoints(
  testg,
  show_image = T,
  largeImage_name = 'Vimentin',
  spat_unit = 'largecell',
  polygon_feat_type = "largecell",
  show_polygon = T,
  polygon_color = 'red',
  polygon_alpha = 1,
  polygon_line_size = 0.1,
  polygon_fill = NULL,
  feat_type = "protein",
  feats = NULL,
  save_param = list(
    save_name = "largecell_vimentin"
  )
)

## 4. filter giotto object and add statistics ####
# --------------------------------------------- #
# check for cells with zero overlaps
testg = filterGiotto(
  gobject = testg,
  spat_unit = 'smallcell',
  poly_info = 'smallcell',
  feat_type = 'protein',
  expression_threshold = 0,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 2
)

testg = filterGiotto(
  gobject = testg,
  spat_unit = 'cell',
  poly_info = 'cell',
  feat_type = 'protein',
  expression_threshold = 0,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 2
)

testg = filterGiotto(
  gobject = testg,
  spat_unit = 'largecell',
  poly_info = 'largecell',
  feat_type = 'protein',
  expression_threshold = 0,
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 2
)

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

screePlot(testg, spat_unit = "cell", save_param = list(save_name = "cell_scree")) # 10PCs
pcs <- seq(10)

# UMAP
testg = runUMAP(gobject = testg, spat_unit = 'smallcell',
                dimensions_to_use = pcs,
                set_seed = TRUE, seed_number = my_seed_num)

testg = runUMAP(gobject = testg, spat_unit = 'cell',
                dimensions_to_use = pcs,
                set_seed = TRUE, seed_number = my_seed_num)

testg = runUMAP(gobject = testg, spat_unit = 'largecell',
                dimensions_to_use = pcs,
                set_seed = TRUE, seed_number = my_seed_num)


# nearest network & clustering
testg = createNearestNetwork(testg, spat_unit = 'smallcell', dimensions_to_use = pcs)
testg = createNearestNetwork(testg, spat_unit = 'cell', dimensions_to_use = pcs)
testg = createNearestNetwork(testg, spat_unit = 'largecell', dimensions_to_use = pcs)

# K-means
k <- 9
testg = doKmeans(testg,
                 spat_unit = "smallcell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = k,
                 iter_max = 1000,
                 nstart = 1000,
                 set_seed = TRUE,
                 seed_number = my_seed_num)

testg = doKmeans(testg,
                 spat_unit = "cell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = k,
                 iter_max = 1000,
                 nstart = 1000,
                 set_seed = TRUE,
                 seed_number = my_seed_num)

testg = doKmeans(testg,
                 spat_unit = "largecell",
                 feat_type = "protein",
                 expression_values = "scaled",
                 centers = k,
                 iter_max = 1000,
                 nstart = 1000,
                 set_seed = TRUE,
                 seed_number = my_seed_num)


## 7. compare k-means clustering results for spatial units cell and smallcell ####
# ------------------------------------------------------------------------------ #
library(data.table)

# combine metadata from spatial units 'cell' and 'smallcell'
small_cell_meta = pDataDT(testg, spat_unit = 'smallcell')
setnames(small_cell_meta, old = 'kmeans', new = 'km_smallcell')

cell_meta = pDataDT(testg, spat_unit = 'cell')
setnames(cell_meta, old = 'kmeans', new = 'km_cell')

comb_meta = data.table::merge.data.table(
  small_cell_meta[,.(cell_ID, km_smallcell)],
  cell_meta[,.(cell_ID, km_cell)], by = 'cell_ID'
)

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

# define cell colors
my_colors = getDistinctColors(k)
order_cell = as.numeric(unlist(lapply(stable_combos, FUN = function(x) {
  strsplit(x, split = '-')[[1]][2]
})))
names(my_colors) = order_cell


# spatplot kmeans
spatInSituPlotPoints(
  gobject = testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "S9B1_spat_cell",
    save_format = "svg"
  )
)

plotUMAP(
  testg,
  spat_unit = "cell",
  cell_color = "kmeans",
  cell_color_code = my_colors,
  point_size = 1,
  point_shape = "no_border",
  show_center_label = FALSE,
  save_param = list(
    save_name = "S9B2_umap_cell",
    save_format = "svg"
  )
)


# spatplot switch
spatInSituPlotPoints(
  gobject = testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "cell_switches",
    save_format = "svg"
  )
)




## * 7.2. plot smaller cells ####

# definte smallcell colors
my_colors_small = getDistinctColors(k)
order_small_cell = as.numeric(unlist(lapply(stable_combos, FUN = function(x) {
  strsplit(x, split = '-')[[1]][1]
})))
names(my_colors_small) = order_small_cell

# kmeans spatplot
spatInSituPlotPoints(
  gobject = testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "S9C1_spat_smallcell",
    save_format = "svg"
  )
)

# umap
plotUMAP(
  testg,
  spat_unit = "smallcell",
  cell_color = "kmeans",
  cell_color_code = my_colors_small,
  point_size = 2,
  point_shape = "no_border",
  show_center_label = FALSE,
  save_param = list(
    save_name = "S9C2_umap_smallcell",
    save_format = "svg"
  )
)


# switch spatplot
spatInSituPlotPoints(
  gobject = testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "smallcell_switches",
    save_format = "svg"
  )
)




## 8. Identify spatial network neighbors ####
# ----------------------------------------- #

# Show whether cells have switched kmeans cluster between cell and smallcell
# spat_units and the spatial neighbors/niches that may be affected by the
# switch.

# create spatial network
testg = createSpatialNetwork(testg)
showGiottoSpatNetworks(testg)

# cell cell IDs from cells that switched cluster
cell_metadata = pDataDT(testg, spat_unit = 'cell')
switch_cell_IDs = cell_metadata[km_stable == 'switch', cell_ID]

# find neighbor cell IDs
# Find relationships between all cells and the cells that switched.
# - neighbor: spatial neighbor of a switched cell
# - source: switched cell
# - both: a source and a neighbor of another source
# - other: not a neighbor or source
find_switch_neighbors = findNetworkNeighbors(
  gobject = testg, 
  spat_unit = 'cell',
  spatial_network_name = 'Delaunay_network',
  source_cell_ids = switch_cell_IDs
)

# add to metadata for cell
testg = addCellMetadata(
  testg, 
  spat_unit = 'cell',
  new_metadata = find_switch_neighbors,
  by_column = T, 
  column_cell_ID = 'cell_ID'
)

# add to metadata for small cell
testg = addCellMetadata(
  testg, 
  spat_unit = 'smallcell',
  new_metadata = find_switch_neighbors,
  by_column = T, 
  column_cell_ID = 'cell_ID'
)

# visualize
switch_nb_colors = c('gray', 'orange', 'blue', 'red')
names(switch_nb_colors) = c('others', 'neighbor', 'source', 'both')


# cell
spatInSituPlotPoints(
  gobject = testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "S9D_cell_switch_nb",
    save_format = "svg"
  )
)




# small cell
spatInSituPlotPoints(
  gobject = testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "smallcell_switch_nb",
    save_format = "svg"
  )
)



## barplot
cell_metadata = pDataDT(testg, spat_unit = 'cell')
cell_metadata[, nb_cells := factor(nb_cells, levels = c('source', 'both', 'neighbor', 'others'))]
cell_metadata[, (.N/1679)*100, by = .(nb_cells)]

# neighbors
library(ggplot2)

# switches
cell_metadata[, km_stable := factor(km_stable, levels = c('switch', 'stable'))]

pl = ggplot()
pl = pl + geom_bar(data = cell_metadata, aes(x = '', fill = km_stable),
                   position = 'fill')
pl = pl + ggplot2::scale_fill_manual(values = c(stable = 'gray', switch = 'red'))
pl = pl + theme_classic()
pl

ggplot2::ggsave(filename = paste0(results_folder, "/S9E1_barplot_switches.svg"),
                plot = pl,
                device = "svg",
                width = 3,
                height = 6)

pl = ggplot()
pl = pl + geom_bar(data = cell_metadata, aes(x = '', fill = nb_cells),
                   position = 'fill')
pl = pl + ggplot2::scale_fill_manual(values = switch_nb_colors)
pl = pl + theme_classic()
pl

ggplot2::ggsave(filename = paste0(results_folder, "/S9E2_barplot_nbs.svg"),
                plot = pl,
                device = "svg",
                width = 3,
                height = 6)



## 9. zoomed-in regions ####
# ------------------------ #

## for spatial unit = smallcell
sub_sc_testg = subsetGiottoLocs(
  testg,
  spat_unit = "smallcell",
  spat_loc_name = "raw",
  x_max = 500,
  x_min = 400,
  y_max = 100,
  y_min = 0,
  poly_info = 'smallcell'
)

# kmeans
spatInSituPlotPoints(
  gobject = sub_sc_testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "S9C3_mini_spat_smallcell",
    save_format = "svg"
  )
)

# switches
spatInSituPlotPoints(
  gobject = sub_sc_testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "mini_smallcell_switches",
    save_format = "svg"
  )
)


# neighbors
spatInSituPlotPoints(
  gobject = sub_sc_testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "mini_smallcell_nb",
    save_format = "svg"
  )
)




## for spatial unit = cell
sub_c_testg <- subsetGiottoLocs(
  testg,
  spat_unit = "cell",
  spat_loc_name = "raw",
  x_max = 500,
  x_min = 400,
  y_max = 100,
  y_min = 0,
  poly_info = 'cell'
)

spatInSituPlotPoints(
  gobject = sub_c_testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "S9B3_mini_spat_cell",
    save_format = "svg"
  )
)


# switches
spatInSituPlotPoints(
  gobject = sub_c_testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "mini_cell_switches",
    save_format = "svg"
  )
)


# neighbors
spatInSituPlotPoints(
  gobject = sub_c_testg,
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
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "mini_cell_nb",
    save_format = "svg"
  )
)





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


ggplot2::ggsave(filename = paste0(results_folder, "./S9F_sankeyplot.svg"),
                plot = pl,
                device = "svg",
                width = 8,
                height = 4)




## 11. niche clustering ####
# ------------------------ #
testg <- calculateSpatCellMetadataProportions(
  gobject = testg,
  spat_unit = 'cell',
  spat_network = 'Delaunay_network',
  metadata_column = 'kmeans', 
  name = 'km_niche'
)
showGiottoSpatEnrichments(testg)

# get spatial enrichment data
prop_table = getSpatialEnrichment(testg, name = 'km_niche', output = 'data.table')
prop_matrix = GiottoUtils::dt_to_matrix(prop_table)

# perform kmeans on spatial enrichment data
prop_kmeans = kmeans(x = prop_matrix, centers = k, iter.max = 1000, nstart = 100)
prop_kmeansDT = data.table(cell_ID = names(prop_kmeans$cluster), niche = prop_kmeans$cluster)

# add results back to giotto object
testg = addCellMetadata(testg, new_metadata = prop_kmeansDT,
                        by_column = T, column_cell_ID = 'cell_ID')

# visualize nich clustering
spatInSituPlotPoints(
  gobject = testg,
  spat_unit = "cell",
  polygon_feat_type = "cell",
  show_polygon = T,
  feat_type = "protein",
  feats = NULL,
  polygon_fill = 'niche',
  polygon_fill_as_factor = TRUE,
  polygon_line_size = 0.05,
  polygon_color = 'white',
  polygon_alpha = 1,
  background_color = "black",
  save_param = list(
    save_name = "cell_niches.png",
    save_format = "svg"
  )
)
