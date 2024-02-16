
library(Giotto)

# 1. create g dataset from vizgen FFPE breast cancer dataset ####
# ---------------------------------------------------------------- #

# see SCRIPTS/create_mini_viz_human_bc

# 2. load giotto object ####
# ------------------------ #
g = loadGiotto(path_to_folder = 'DATA/vizgen_ffpe_BC_mini/')
instructions(g, 'save_dir') = 'SCRIPTS/FIGURES/viz_ffpe_bc_hierarchical/plots/'
instructions(g, 'show_plot') = TRUE
instructions(g, 'return_plot') = FALSE
instructions(g, 'save_plot') = FALSE


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

cell_color_codes <- getColors('distinct', n = 7)
nuc_color_codes <- getColors('vivid', n = 7, src = 'rcartocolor')

# supplemental 4
plotUMAP(g, spat_unit = 'cell',
         cell_color = 'leiden_clus',
         label_size = 0,
         point_size = 4,
         cell_color_code = cell_color_codes
         )
plotUMAP(g, spat_unit = 'nucleus',
         cell_color = 'leiden_clus',
         cell_color_code = nuc_color_codes,
         label_size = 0,
         point_size = 4)


# Spatially plot the leiden clusters
spatInSituPlotPoints(
  g,
  polygon_feat_type = 'cell',
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'leiden_clus',
  polygon_line_size = 0.1,
  polygon_fill_code = cell_color_codes,
  polygon_alpha = 1
)
spatInSituPlotPoints(
  g,
  polygon_feat_type = 'nucleus',
  polygon_fill_as_factor = TRUE,
  polygon_fill = 'leiden_clus',
  polygon_line_size = 0.1,
  polygon_fill_code = nuc_color_codes,
  polygon_alpha = 1
)







# niche detection ####
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



# Calculate niche composition --------------------------------------------- #
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
  color_as_factor = F,
  gradient_style = 's',
  cell_color = '3'
)


# get the enrichments as a proportions matrix
prop_table = getSpatialEnrichment(g, spat_unit = 'cell',
                                  name = 'proportion',
                                  output = 'data.table')
prop_matrix = dt_to_matrix(prop_table)

# These enrichments are essentially a measure of how many cells of each
# leiden cluster exist in the local region
# 
# Using kmeans, we can classify each cell by its niche leiden cluster proportions
set.seed(12345) # set seed for kmeans
prop_kmeans = kmeans(x = prop_matrix, centers = 6, iter.max = 100, nstart = 3)
prop_kmeansDT = data.table::data.table(cell_ID = names(prop_kmeans$cluster), niche = prop_kmeans$cluster)
g = addCellMetadata(g, spat_unit = 'cell' , new_metadata = prop_kmeansDT, by_column = T, column_cell_ID = 'cell_ID')


# Spatially visualize the niches
spatPlot(gobject = g, show_network = TRUE,
         network_color = 'lightgray', spatial_network_name = 'Delaunay_network',
         point_size = 2.5, cell_color = 'niche')




# Splitting niche type info to neighborhoods --------------------------------- #
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
spatPlot(gobject = g,
         show_network = T,
         network_color = '#493944',
         spatial_network_name = 'Delaunay_network',
         point_size = 2.5,
         cell_color = 'neighborhood',
         show_legend = F,
         cell_color_code = getColors('distinct', 200),
         background = '#1F0A00')

# Visualize as polygons
spatInSituPlotPoints(g,
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = 'neighborhood',
                     polygon_alpha = 1, 
                     show_legend = FALSE)



# Use neighborhood as new spatial unit ####
# making a spat unit like this requires a set of cell_IDs and also a spatial
# network to base the spatial interactions on.







# Use hmrf as new spatial unit

km_spatialfeats = binSpect(g, spat_unit = 'cell')
ext_spatial_genes = km_spatialfeats[1:200,]$feats

spat_cor_netw_DT = detectSpatialCorFeats(g,
                                         spat_unit = 'cell',
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network',
                                         subset_feats = ext_spatial_genes)

# cluster and visualize spatial co-expression genes
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, name = 'spat_netw_clus', k = 7)


testweight = getBalancedSpatCoexpressionFeats(spat_cor_netw_DT,
                                              rank = 'weighted',
                                              maximum = 20)



hmrf_folder = paste0('SCRIPTS/FIGURES/viz_ffpe_bc_hierarchical//','HMRF/')
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
  output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k20_scaled')
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
g = spatialSplitCluster(g, cluster_col = 'HMRF_k3_b.40', split_clus_name = 'hmrf_split')

spatInSituPlotPoints(g,
                     polygon_fill_as_factor = TRUE,
                     polygon_fill = 'hmrf_split',
                     polygon_alpha = 1)






# sankey plot
# for different levels

# get polygons
## cell level
cell = getPolygonInfo(g, polygon_name = 'cell', return_giottoPolygon = TRUE)
## nucleus level
nuc = getPolygonInfo(g, polygon_name = 'nucleus', return_giottoPolygon = TRUE)

# group of spatvector poly as multipolygon
meta = pDataDT(g)

multipoly_from_groupings = function(group_col = NULL,
                                    cell_meta = meta,
                                    poly_to_use = cell,
                                    poly_name = 'multipoly',
                                    prefix = NULL) {
  
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

# find intersects
calc_intersect_dt = function(a, b) {
  # take the poly_ID col only
  terra::intersect(a[,1][], b[,1][])[,1:2] |>
    terra::values() |>
    data.table::setDT() |>
    data.table::setnames(new = c('a_ID', 'b_ID')) |>
    data.table::setkeyv(c('a_ID', 'b_ID'))
}

# cell_nuc = calc_intersect_dt(cell, nuc)
# neighbor_cell = calc_intersect_dt(neighborhood, cell)
domain_neighbor = calc_intersect_dt(domain, neighborhood)

sankeyPlot(domain_neighbor, nodePadding = 0)

# each of these sets of relationships have too many links to be comfortably
# plotted with sankey. This is demonstrated with the domain vs neighbor plot
# above

# for neighborhoods, pick only those that have 5 or more cells in the
# neighborhood.
# This can be done using the group_n attribute which is generated during 
# multipolygon generation
neighbor_sub = neighborhood[neighborhood$group_n >= 10]
# That cuts down the number of neighborhoods to show from 184 to 44

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

# figure 2B ####
# there are too many labels for easy viewing, so we will hide less relevant ones
sankeyPlot(
  sankey_dt,
  focus_names = c(paste0('domain_', c(3)),
                  unique(neighbor_cell_sub$a_ID)),
  unfocused_color = FALSE,
  fontSize = 23,
  nodeWidth = 50
)



# (supplemental 4) ####

# transcript level
tx = getFeatureInfo(g, return_giottoPoints = TRUE)
plot(tx)

# plot dapi (intensity info)
distGiottoImage(g, image_name = 'dapi_z0', image_type = 'largeImage')
plotGiottoImage(g, image_name = 'dapi_z0',
                image_type = 'largeImage',
                largeImage_max_intensity = 15000)

# nucleus level
plot(nuc, col = getDistinctColors(nrow(nuc)), background = 'black', alpha = 0.4)
nuc_sub = nuc[cell_nuc_sub$b_ID]
plot(nuc_sub, col = getDistinctColors(nrow(nuc)), background = 'black', alpha = 0.7, border = 'white', lwd = 3, add = TRUE)

# cell level
plot(cell, col = getDistinctColors(nrow(cell)), background = 'black', alpha = 0.4)
plot(cell_sub, col = getDistinctColors(nrow(cell)), background = 'black', alpha = 0.7, border = 'white', lwd = 3, add = TRUE)


# neighborhood level
plot(neighborhood, col = getDistinctColors(nrow(neighborhood)), background = 'black', alpha = 0.4)

plot(neighbor_sub_2, col = getDistinctColors(nrow(neighborhood)), background = 'black', alpha = 0.7, border = 'white', add = TRUE, lwd = 3)

# domain level
plot(domain, col = getDistinctColors(nrow(domain)), background = 'black', alpha = 0.4)
domain_sub = domain[paste0('domain_', c(3))]
plot(domain_sub, col = getDistinctColors(nrow(domain)), background = 'black', alpha = 0.7, border = 'white', lwd = 3, add = TRUE)


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

sankeyPlot(cell_nuc_leiden,
           nodePadding = 30,
           fontSize = 15,
           nodeWidth = 10)

  # no labels
sankeyPlot(cell_nuc_leiden,
           nodePadding = 60,
           fontSize = 0,
           nodeWidth = 10)
