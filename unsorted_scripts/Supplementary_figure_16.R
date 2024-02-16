## --------------------------------------------------------------------------------------
# Load the object after running the code from Supplementary figure 3 DBiT-seq 

library(Giotto)

giottoObject = loadGiotto("giottoObject/")

## --------------------------------------------------------------------------------------
# Spatial grid
giottoObject <- createSpatialGrid(gobject = giottoObject,
                                  sdimx_stepsize = 1,
                                  sdimy_stepsize = 2,
                                  minimum_padding = 0)

showGiottoSpatGrids(giottoObject)

spatPlot2D(giottoObject, 
           cell_color = 'leiden_clus', 
           show_grid = T,
           grid_color = 'red', 
           spatial_grid_name = 'spatial_grid')

## --------------------------------------------------------------------------------------
# Spatial network
giottoObject <- createSpatialNetwork(gobject = giottoObject,
                                     method = 'kNN', 
                                     k = 6,
                                     maximum_distance_knn = 5,
                                     name = 'spatial_network')

showGiottoSpatNetworks(giottoObject)

spatPlot2D(gobject = giottoObject,  
           show_network= T,
           network_color = 'blue', 
           spatial_network_name = 'spatial_network')

## --------------------------------------------------------------------------------------
# Spatial genes
## rank binarization
ranktest = binSpect(giottoObject, 
                    bin_method = 'rank',
                    calc_hub = T, 
                    hub_min_int = 5,
                    spatial_network_name = 'spatial_network')

spatFeatPlot2D(giottoObject, 
               expression_values = 'scaled',
               feats = ranktest$feats[1:6], 
               cow_n_col = 2, 
               point_size = 1.5)

## --------------------------------------------------------------------------------------
# Spatial co-expression modules

# cluster the top 500 spatial genes into 20 clusters
ext_spatial_genes = ranktest[1:1500,]$feats

# here we use existing detectSpatialCorGenes function to calculate pairwise distances between genes (but set network_smoothing=0 to use default clustering)
spat_cor_netw_DT = detectSpatialCorFeats(giottoObject,
                                         method = 'network',
                                         spatial_network_name = 'spatial_network',
                                         subset_feats = ext_spatial_genes)

# cluster spatial genes
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, 
                                          name = 'spat_netw_clus', 
                                          k = 20)

# visualize clusters
heatmSpatialCorFeats(giottoObject,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus',
                     heatmap_legend_param = list(title = NULL),
                     save_param = list(base_height = 6, 
                                       base_width = 8, 
                                       units = 'cm'))

## --------------------------------------------------------------------------------------
# 4. rank spatial correlated clusters and show genes for selected clusters
netw_ranks = rankSpatialCorGroups(giottoObject,
                                  spatCorObject = spat_cor_netw_DT, 
                                  use_clus_name = 'spat_netw_clus',
                                  save_plot = FALSE)

top_netw_spat_cluster = showSpatialCorFeats(spat_cor_netw_DT, 
                                            use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, 
                                            show_top_feats = 1)

## --------------------------------------------------------------------------------------
# 5. create metagene enrichment score for clusters
cluster_genes_DT = showSpatialCorFeats(spat_cor_netw_DT, 
                                       use_clus_name = 'spat_netw_clus', 
                                       show_top_feats = 1)
cluster_genes = cluster_genes_DT$clus 
names(cluster_genes) = cluster_genes_DT$feat_ID

giottoObject = createMetafeats(giottoObject, 
                               feat_clusters = cluster_genes, 
                               name = 'cluster_metagene')

spatCellPlot(giottoObject,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1, 
             cow_n_col = 5, 
             save_param = list(base_width = 15))

## --------------------------------------------------------------------------------------
# Spatially informed clusters
# top 30 genes per spatial co-expression cluster
table(spat_cor_netw_DT$cor_clusters$spat_netw_clus)
coexpr_dt = data.table::data.table(genes = names(spat_cor_netw_DT$cor_clusters$spat_netw_clus),
                                   cluster = spat_cor_netw_DT$cor_clusters$spat_netw_clus)
data.table::setorder(coexpr_dt, cluster)
top30_coexpr_dt = coexpr_dt[, head(.SD, 30) , by = cluster]
my_spatial_genes <- top30_coexpr_dt$genes

giottoObject <- runPCA(gobject = giottoObject,
                       feats_to_use = my_spatial_genes,
                       name = 'custom_pca')

giottoObject <- runUMAP(giottoObject, 
                        dim_reduction_name = 'custom_pca', 
                        dimensions_to_use = 1:20,
                        name = 'custom_umap')

giottoObject <- createNearestNetwork(gobject = giottoObject,
                                     dim_reduction_name = 'custom_pca',
                                     dimensions_to_use = 1:20, 
                                     k = 5,
                                     name = 'custom_NN')

giottoObject <- doLeidenCluster(gobject = giottoObject, 
                                network_name = 'custom_NN',
                                resolution = 0.15, 
                                n_iterations = 1000,
                                name = 'custom_leiden')

cell_meta = pDataDT(giottoObject)

cell_clusters = unique(cell_meta$custom_leiden)

selected_colors = getDistinctColors(length(cell_clusters))
names(selected_colors) = cell_clusters

spatPlot2D(giottoObject, 
           cell_color = 'custom_leiden', 
           cell_color_code = selected_colors, 
           point_size = 3.5,
           axis_text = 14,
           axis_title = 18,
           legend_text = 14)


