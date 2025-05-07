##%##########################################################%##
#                                                              #
####          Figure 2 B - G                                   #
#                                                              #
##%##########################################################%##



library(Giotto)

# set working directory to project
setwd("PATH TO Giotto_Suite_manuscript REPO")

# create dataset from vizgen FFPE breast cancer dataset ####

# see scripts/CREATE/create_vizgen_hBC_mini.R




# load giotto object ####
python_path <- NULL
g = loadGiotto(path_to_folder = 'scripts/CREATE/OUTS/vizgen_ffpe_BC_mini/',
               python_path = python_path)
instructions(g, 'save_dir') = 'scripts/FIGURE/F2/'
instructions(g, 'show_plot') = TRUE
instructions(g, 'return_plot') = FALSE
instructions(g, 'save_plot') = TRUE




##### Fig 2 C -----------------------------------------------------start ####

# explore contents
force(g)

# plot transcripts and image information
tx_points <- getFeatureInfo(g, return_giottoPoints = TRUE)
dapi_img <- getGiottoImage(g, image_type = "largeImage", name = "dapi_z0")

nrow(tx_points) # number of transcripts within this area

# plot a rasterized representation of the points information
# without density-based colorization, it displays more clearly where the
# points spatially cover.
plot(tx_points) #                                     [F2 C1]
# With the density, the visualization better shows how the transcripts
# recapitulate the morphology.
plot(tx_points, dens = TRUE)


# plot a resampling of the image
plot(dapi_img)
distGiottoImage(giottoLargeImage = dapi_img)
plot(dapi_img, max_intensity = 15000) #               [F2 C2]
##### Fig 2 C --------------------------------------------------------end ###


# process data ####
# only cell and nucleus spat units

giotto_data_process = function(spat_unit) {
  
  g |>
    # aggregate feature counts per polygon
    calculateOverlapRaster(spatial_info = spat_unit) |>
    # transfer aggregate values to expression count matrix
    overlapToMatrix(poly_info = spat_unit) |>
    # filter expression matrix.
    filterGiotto(
      spat_unit = spat_unit,
      expression_threshold = 1,
      feat_det_in_min_cells = 3,
      min_det_feats_per_cell = 25
    ) |>
    # run standard library and log normalization of expression matrix
    normalizeGiotto(spat_unit = spat_unit,
                    norm_methods = 'standard') |>
    # new normalization method based on pearson correlations (Lause/Kobak et al. 2021)
    # this normalized matrix is given the name 'pearson' using the update_slot param
    # results only used in hvf detection and pca generation
    normalizeGiotto(spat_unit = spat_unit,
                    scalefactor = 5000,
                    norm_methods = 'pearson_resid',
                    update_slot = 'pearson')
}

activeFeatType(g) = 'rna'
g = giotto_data_process(spat_unit = 'cell')
g = giotto_data_process(spat_unit = 'nucleus')



# Since there are only 500 genes to consider, we will skip highly variable 
# feature detection and directly use all 500 for dimension reduction.

# dimension reduction ####

# generate PCA
giotto_pca = function(spat_unit) {
  g = runPCA(gobject = g,
             feats_to_use = NULL,
             spat_unit = spat_unit,
             # using pearson normalized values
             expression_values = 'pearson')
}

g = giotto_pca('cell')
g = giotto_pca('nucleus')

screePlot(g, spat_unit = 'cell', ncp = 50)
screePlot(g, spat_unit = 'nucleus', ncp = 50)

# Use similar dims of PCA
dims_to_use = list(
  'cell' = 1:15,
  'nucleus' = 1:15
)

# generate UMAP
g = runUMAP(gobject = g,
            spat_unit = 'cell',
            n_neighbors = 10,
            dimensions_to_use = dims_to_use[['cell']])
g = runUMAP(gobject = g,
            spat_unit = 'nucleus',
            dimensions_to_use = dims_to_use[['nucleus']])

# plot UMAPs
plotUMAP(g, spat_unit = 'cell')
plotUMAP(g, spat_unit = 'nucleus')



# leiden clustering
giotto_leiden_clus = function(spat_unit, res) {
  g |>
    createNearestNetwork(
      spat_unit = spat_unit,
      dimensions_to_use = dims_to_use[[spat_unit]]
    ) |>
    doLeidenCluster(
      spat_unit = spat_unit,
      resolution = res,
      n_iterations = 100
    )
}

# same resolution used for leiden
g = giotto_leiden_clus('cell', res = 0.7)
g = giotto_leiden_clus('nucleus', res = 0.7)

# Know from having run this already:
# 7 clusters for cell and nucleus

# get arbitrary colors to use

cell_color_codes <- getColors('distinct', n = 7, src = "rcartocolor")
nuc_color_codes <- getColors('vivid', n = 7, src = 'rcartocolor')




##### Fig 2 E - F -------------------------------------------------start ####
# Plot leiden clusters on UMAP
plotUMAP(g, spat_unit = 'cell',                       # [F2 E1]
         cell_color = 'leiden_clus',
         show_center_label = FALSE,
         point_size = 1,
         point_shape = "no_border",
         cell_color_code = cell_color_codes,
         save_param = list(
           save_name = "E1",
           save_format = "svg"
         ))
plotUMAP(g, spat_unit = 'nucleus',                    # [F2 F1]
         cell_color = 'leiden_clus',
         cell_color_code = nuc_color_codes,
         show_center_label = FALSE,
         point_size = 1,
         point_shape = "no_border",
         save_param = list(
           save_name = "F1",
           save_format = "svg"
         ))


# Spatially plot the leiden clusters
spatInSituPlotPoints(
  g,
  polygon_feat_type = 'cell',                         # [F2 E2]
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'leiden_clus',
  polygon_line_size = 0.1,
  polygon_fill_code = cell_color_codes,
  polygon_alpha = 1,
  save_param = list(
    save_name = "E2",
    save_format = "svg"
  )
)
spatInSituPlotPoints(
  g,
  polygon_feat_type = 'nucleus',                      # [F2 F2]
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'leiden_clus',
  polygon_line_size = 0.1,
  polygon_fill_code = nuc_color_codes,
  polygon_alpha = 1,
  save_param = list(
    save_name = "F2",
    save_format = "svg"
  )
)
##### Fig 2 E - F ----------------------------------------------------end ###






# niche neighborhoods ####

# niche detection
# Find spatially local proportions for each of the leiden clustering annotations
# based on the cell and its neighbors in the delaunay network.

# In order to do this, we will need a spatial network
# create spatial network
g = createSpatialNetwork(
  gobject = g,
  spat_unit = 'cell',
  method = 'Delaunay'
)
g = createSpatialNetwork(
  gobject = g,
  spat_unit = 'nucleus',
  method = 'Delaunay'
)




# Calculate niche composition
# This is performed by looking at the neighbors of each cell and finding what
# proportion of them fall within each of categories of a metadata column. The
# metadata column that we are selecting here for the niche calculation is the
# leiden clustering.
# These results are then stored in a new metadata column called 'proportion'
g = calculateSpatCellMetadataProportions(
  gobject = g,
  spat_unit = 'cell',
  feat_type = 'rna',
  spat_network = 'Delaunay_network',
  metadata_column = 'leiden_clus',
  name = 'proportion'
)


# visualize niche-level enrichment for leiden cluster 3 in spat_unit "cell"
spatPlot2D(
  gobject = g,
  spat_unit = 'cell',
  spat_enr_names = 'proportion',
  color_as_factor = FALSE,
  gradient_style = 'sequential',
  cell_color = '3'
)


# get proportions matrix
# Enrichments information such as niche enrichment are categorized as and stored
# within the spatial enrichments slot of the Giotto object.
# Spatial enrichments are represented as data.table under the hood.
prop_table = getSpatialEnrichment(
  g, spat_unit = 'cell', name = 'proportion', output = 'data.table'
)
# convert the data.table to a sparse Matrix with row and colnames
# here we use a utility function to perform the operation
prop_matrix = GiottoUtils::dt_to_matrix(prop_table)



# These enrichments are essentially a measure of how many cells of each
# leiden cluster exist in the local region
# 
# Using kmeans, we can classify each cell by its niche leiden cluster proportions
set.seed(12345) # set seed for kmeans
prop_kmeans = kmeans(x = prop_matrix, centers = 6, iter.max = 100, nstart = 3)
prop_kmeansDT = data.table::data.table(
  cell_ID = names(prop_kmeans$cluster), 
  niche = prop_kmeans$cluster
)

# add kmeans clustering of niches to cell metadata
g = addCellMetadata(
  g, spat_unit = 'cell' , 
  new_metadata = prop_kmeansDT, 
  by_column = TRUE,
  column_cell_ID = 'cell_ID'
)

# Spatially visualize the niches
spatPlot(gobject = g, 
         show_network = TRUE,
         network_color = 'lightgray', 
         spatial_network_name = 'Delaunay_network',
         point_size = 2.5, 
         cell_color = 'niche')




# Splitting niche type info to neighborhoods
# The produced niche classifications are classified by niche composition. To
# turn these annotations into unique spatial neighborhoods, we split the niche
# classifications based on whether they are spatially touching (as shown by
# the previous plot which also has a Delaunay spatial network displayed
# alongside the niche plot.)
# 
# spatialSplitCluster will take any cluster column, then split those cluster
# annotations into several new ones depending on if they are spatially touching
# here the results are stored in a new metadata column called neighborhood.
g = spatialSplitCluster(g, cluster_col = 'niche', split_clus_name = 'neighborhood')

pDataDT(g)[, max(neighborhood)]
# 148 niches generated

# Visualize as points and spatial network
spatPlot(g,
    show_network = TRUE,
    network_color = '#493944',
    spatial_network_name = 'Delaunay_network',
    point_size = 2.5,
    cell_color = 'neighborhood',
    show_legend = FALSE,
    cell_color_code = getColors('distinct', 200),
    background = '#1F0A00'
)

# Visualize as polygons
spatInSituPlotPoints(g,
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'neighborhood',
  show_legend = FALSE,
  polygon_alpha = 1
)








# HMRF domains ####

# Use HMRF results as new spatial unit
# Start by getting a set of balanced spatial coexpression genes to use with
# HMRF.
km_spatialfeats = binSpect(g, spat_unit = 'cell')
ext_spatial_genes = km_spatialfeats[1:200,]$feats

spat_cor_netw_DT = detectSpatialCorFeats(
  g,
  spat_unit = 'cell',
  method = 'network',
  spatial_network_name = 'Delaunay_network',
  subset_feats = ext_spatial_genes
)

# cluster and visualize spatial co-expression genes
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, 
                                          name = 'spat_netw_clus', 
                                          k = 7)


testweight = getBalancedSpatCoexpressionFeats(spat_cor_netw_DT,
                                              rank = 'weighted',
                                              maximum = 20)


# run HMRF
hmrf_folder = file.path('scripts/FIGURE/S4','HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = TRUE)

# kNN with k of 8 nearest
g = createSpatialNetwork(g, spat_unit = 'cell', method = 'kNN', k = 8)

HMRF_spatial_genes = doHMRF(
  gobject = g,
  expression_values = 'scaled',
  spatial_genes = names(testweight),
  k = 3,
  spatial_network_name = "kNN_network",
  betas = c(0, 10, 5),
  output_folder = file.path(hmrf_folder, "Spatial_genes/SG_topgenes_k20_scaled")
)

g = addHMRF(gobject = g,
            HMRFoutput = HMRF_spatial_genes,
            k = 3, 
            betas_to_add = c(0, 10, 20, 30, 40),
            hmrf_name = 'HMRF')

pDataDT(g)

spatPlot2D(gobject = g, cell_color = 'HMRF_k3_b.40')

spatInSituPlotPoints(g,
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = 'HMRF_k3_b.40',
                     polygon_alpha = 1)

# spatially split clustering
g = spatialSplitCluster(g, cluster_col = 'HMRF_k3_b.40', 
                        split_clus_name = 'hmrf_split')

spatInSituPlotPoints(g,
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = 'hmrf_split',
                     polygon_alpha = 1)






# sankey plot
# for different levels

# get polygons
## cell level
cell <- getPolygonInfo(g, polygon_name = 'cell', return_giottoPolygon = TRUE)
cell[] <- terra::makeValid(cell[])
## nucleus level
nuc <- getPolygonInfo(g, polygon_name = 'nucleus', return_giottoPolygon = TRUE)
nuc[] <- terra::makeValid(nuc[])

# group of spatvector poly as multipolygon
meta <- pDataDT(g)

# workflow function to apply the same operations to two spatial units
multipoly_from_groupings <- function(
    group_col = NULL,
    cell_meta = meta,
    poly_to_use = cell,
    poly_name = 'multipoly',
    prefix = NULL
) {
  
  groups = cell_meta[, c('cell_ID', group_col), with = FALSE]
  if (!is.null(prefix)) {
    groups[, (group_col) := paste0(paste0(prefix, '_'), get(group_col))]
  }
  data.table::setnames(groups, c('poly_ID', 'group_ID'))
  
  combineToMultiPolygon(
    x = poly_to_use,
    groups = groups,
    name = poly_name
  )
}


# multipolygons from neighborhood and domain clusters ####

## neighborhood level
neighborhood = multipoly_from_groupings(
  group_col = 'neighborhood',
  poly_name = 'cell_neighborhoods',
  prefix = 'neighborhood'
)

## domain level
domain = multipoly_from_groupings(
  group_col = 'hmrf_split',
  poly_name = 'cell_domains',
  prefix = 'domain'
)

# find intersects using multipolygons ####

# workflow function for finding intersects and returning data.table reports
calc_intersect_dt = function(a, b) {
  # take the poly_ID col only
  terra::intersect(a[,1][], b[,1][])[,1:2] |>
    terra::values() |>
    data.table::setDT() |>
    data.table::setnames(new = c('a_ID', 'b_ID')) |>
    data.table::setkeyv(c('a_ID', 'b_ID'))
}


# preview cross-spat_unit relationships
cell_nuc = calc_intersect_dt(cell, nuc)
sankeyPlot(cell_nuc, nodePadding = 0)

neighbor_cell = calc_intersect_dt(neighborhood, cell)
sankeyPlot(neighbor_cell, nodePadding = 0)

domain_neighbor = calc_intersect_dt(domain, neighborhood)
sankeyPlot(domain_neighbor, nodePadding = 0)


# each of these sets of relationships have too many links to be comfortably
# plotted with sankey.

# To make a clearer plot:
# for neighborhoods, pick only those that have 5 or more cells in the
# neighborhood.
# This can be done using the group_n attribute which is generated during
# multipolygon generation
neighbor_sub = neighborhood[neighborhood$group_n >= 10]
# That cuts down the number of neighborhoods to show from 148 to 36

domain_neighbor_sub = calc_intersect_dt(domain, neighbor_sub)

sankeyPlot(domain_neighbor_sub)

# select a large cell neighborhood to plot
neighbor_sub_2 = neighborhood[neighborhood$group_n == max(neighborhood$group_n)]
neighbor_cell_sub = calc_intersect_dt(neighbor_sub_2, cell)

sankeyPlot(neighbor_cell_sub)

# find subset of cells and the nucleus they map to
cell_ID_sub = neighbor_cell_sub$b_ID
cell_sub = cell[cell_ID_sub]

cell_nuc_sub = calc_intersect_dt(cell_sub, nuc)

sankeyPlot(cell_nuc_sub)

# combine relations to plot
sankey_dt = data.table::rbindlist(list(domain_neighbor,
                                       neighbor_cell_sub,
                                       cell_nuc_sub))
sankeyPlot(sankey_dt)



# figure 2B -------------------------------------------------------------####
# there are too many labels for easy viewing, so we will hide less relevant ones
sankeyPlot(
  sankey_dt,
  focus_names = c(paste0('domain_', c(3)),
                  unique(neighbor_cell_sub$a_ID)),
  unfocused_color = FALSE,
  fontSize = 23,
  nodeWidth = 50
)



# FIGURE EXPORTS (C,D,G) ####

fig_dir <- "scripts/FIGURE/S4/"

# Fig 2 C1 ####
# transcript level
tx = getFeatureInfo(g, return_giottoPoints = TRUE)
svg(file.path(fig_dir, "C1.svg"))
plot(tx)
dev.off()

# Fig 2 C2 ####
# plot dapi (intensity info)
distGiottoImage(g, image_name = 'dapi_z0', image_type = 'largeImage')
plotGiottoImage(g, image_name = 'dapi_z0',
                image_type = 'largeImage',
                largeImage_max_intensity = 15000)

# Fig 2 D ####
# nucleus level # Fig 2 D1
svg(file.path(fig_dir, "D1.svg"))
plot(nuc, col = getDistinctColors(nrow(nuc)), background = 'black', alpha = 0.4)
nuc_sub = nuc[cell_nuc_sub$b_ID]
plot(nuc_sub, col = getDistinctColors(nrow(nuc)), background = 'black', alpha = 0.7, border = 'white', add = TRUE, lwd = 0.5)
dev.off()

# cell level # Fig 2 D2
svg(file.path(fig_dir, "D2.svg"))
plot(cell, col = getDistinctColors(nrow(cell)), background = 'black', alpha = 0.4)
plot(cell_sub, col = getDistinctColors(nrow(cell)), background = 'black', alpha = 0.7, border = 'white', add = TRUE, lwd = 0.5)
dev.off()

# neighborhood level # Fig 2 D3
svg(file.path(fig_dir, "D3.svg"))
plot(neighborhood, col = getDistinctColors(nrow(neighborhood)), background = 'black', alpha = 0.4)
plot(neighbor_sub_2, col = getDistinctColors(nrow(neighborhood)), background = 'black', alpha = 0.7, border = 'white', add = TRUE, lwd = 0.5)
dev.off()

# domain level # Fig 2 D4
svg(file.path(fig_dir, "D4.svg"))
plot(domain, col = getDistinctColors(nrow(domain)), background = 'black', alpha = 0.4)
domain_sub = domain[paste0('domain_', c(3))]
plot(domain_sub, col = getDistinctColors(nrow(domain)), background = 'black', alpha = 0.7, border = 'white', add = TRUE, lwd = 0.5)
dev.off()

# compare cross spat_unit UMAP

# nuc_leiden = sankeySet(spat_unit = 'nucleus', feat_type = 'rna', col = 'leiden_clus')
# cell_leiden = sankeySet(spat_unit = 'cell', feat_type = 'rna', col = 'leiden_clus')
#
# leiden_compare_plan = nuc_leiden + cell_leiden
# sankeyRelate(leiden_compare_plan) = c(0,1)
#
# sankeyPlot(g, leiden_compare_plan)


nuc_leiden = pDataDT(g, spat_unit = 'nucleus')
cell_leiden = pDataDT(g, spat_unit = 'cell')

# previously the only intersects given were a subset of the cells
cell_nuc = calc_intersect_dt(cell, nuc)

data.table::setnames(cell_nuc, new = c('cell_ID', 'nuc_ID'))
data.table::setnames(nuc_leiden, old = c('cell_ID', 'leiden_clus'), new = c('nuc_ID', 'leiden_nuc'))
data.table::setnames(cell_leiden, old = c('leiden_clus'), new = c('leiden_cell'))

cell_nuc_leiden = merge(cell_nuc, nuc_leiden, by = 'nuc_ID', all.y = TRUE)
cell_nuc_leiden = merge(cell_nuc_leiden, cell_leiden, by = 'cell_ID', all.y = TRUE)
cell_nuc_leiden = cell_nuc_leiden[, c('leiden_cell', 'leiden_nuc')]
# change labels so that node names are not overlapped
cell_nuc_leiden[, leiden_cell := paste0('cell_leiden_', leiden_cell)]
cell_nuc_leiden[, leiden_nuc := paste0('nuc_leiden_', leiden_nuc)]

# Fig 2 G ####
sankeyPlot(cell_nuc_leiden,
           nodePadding = 30,
           fontSize = 15,
           nodeWidth = 10)

# no labels
sankeyPlot(cell_nuc_leiden,
           nodePadding = 30,
           fontSize = 0,
           nodeWidth = 10)

