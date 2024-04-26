

############################# Preprocess and Load ##############################

library(terra)
library(patchwork)
library(ggplot2)
library(tictoc)
library(Giotto)

results_folder = './co_cluster_results/final_PDFs/'
setwd(results_folder)

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE)

# ** SET PATH TO FOLDER CONTAINING XENIUM DATA **
xenium_folder = './Xenium/'
aligned_folder = './Xenium/Aligned_Xe_rep1/'
images_folder = './Xenium/Images/'
# files (SUBCELLULAR): (tutorial focuses on working with these files)
cell_bound_path = paste0(aligned_folder, 'cell_boundaries_aligned.csv.gz')
nuc_bound_path = paste0(aligned_folder, 'nucleus_boundaries_aligned.csv.gz')
tx_path = paste0(aligned_folder, 'tx_aligned.csv.gz')
IF_CD20_path = paste0(aligned_folder, 'CD20_aligned.csv.gz')
IF_DAPI_path = paste0(aligned_folder, 'DAPI_aligned.csv.gz')
IF_HER2_path = paste0(aligned_folder, 'HER2_aligned.csv.gz')
feat_meta_path = paste0(xenium_folder, 'cell_feature_matrix/features.tsv.gz') # (also used in aggregate)
HE_img_path = paste0(images_folder,'Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_compressed.png')

#--- load HE image
Xenium_HE_img = imager::load.image(HE_img_path)

#--- load features metadata
# (make sure cell_feature_matrix folder is unpacked)
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('feat_ID','feat_name','feat_type')

# find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])

feat_types_IDs = lapply(
  feat_types, function(type) feature_dt[feat_type == type, unique(feat_name)]
)
names(feat_types_IDs) = feat_types


#---- Use Aligned Coordinate system
# load polygons
cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(cellPoly_dt,
                     old = c('cell_id', 'aligned_x', 'aligned_y'),
                     new = c('poly_ID', 'x', 'y'))
data.table::setnames(nucPoly_dt,
                     old = c('cell_id', 'aligned_x', 'aligned_y'),
                     new = c('poly_ID', 'x', 'y'))

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)

#----------transcript
tx_dt = data.table::fread(tx_path)

data.table::setnames(x = tx_dt,
                     old = c('feature_name', 'aligned_x', 'aligned_y'),
                     new = c('feat_ID', 'x', 'y'))
cat('Transcripts info available:\n ', paste0('"', colnames(tx_dt), '"'), '\n',
    'with', tx_dt[,.N], 'unfiltered detections\n')

# filter by qv (Phred score)
tx_dt_filtered = tx_dt[qv >= 30]
cat('and', tx_dt_filtered[,.N], 'filtered detections\n\n')

# separate detections by feature type
tx_dt_types = lapply(
  feat_types_IDs, function(types) tx_dt_filtered[feat_ID %in% types]
)

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], 'detections: ', tx_dt_types[[x]][,.N], '\n')
}))

# IF
IF_dt_CD20 = data.table::fread(IF_CD20_path)
IF_dt_DAPI = data.table::fread(IF_DAPI_path)
IF_dt_HER2 = data.table::fread(IF_HER2_path)

IF_dt_HER2$IF_molecule = rep("HER2",nrow(IF_dt_HER2))
IF_dt_CD20$IF_molecule = rep("CD20",nrow(IF_dt_CD20))
IF_dt_DAPI$IF_molecule = rep("DAPI",nrow(IF_dt_DAPI))

IF_dt = rbind(IF_dt_CD20,IF_dt_HER2,IF_dt_DAPI)

data.table::setnames(x = IF_dt,
                     old = c('IF_molecule', 'aligned_x', 'aligned_y','x'),
                     new = c('feat_ID', 'x', 'y','intensity'))

tx_dt_types$IF = IF_dt[intensity> 0.1, ]

gpoints_list = lapply(
  tx_dt_types, function(x) createGiottoPoints(x = x)
)


gene_groups_path = paste0(xenium_folder,"Xenium_FFPE_Human_Breast_Cancer_Rep1_gene_groups.csv")
feat_groups = data.table::fread(gene_groups_path, header = FALSE)
colnames(feat_groups) = c("feature", "cell_type")

################### Aligned polygon loading: Baysor ############################

baysor_dt = data.table::fread("./xenium/segmentations_aligned/Baysor_polygons_aligned.csv",
                              drop = c(1,2,3,4)) #drop index, old x, and old y columns

data.table::setnames(baysor_dt,
                     old = c("aligned_x", "aligned_y"),
                     new = c("x", "y"))

baysor_poly = GiottoClass:::dt_to_spatVector_polygon(baysor_dt)

baysor_poly$poly_ID = paste0("polygon_", 1:length(baysor_poly))

baysor_gpoly = giottoPolygon(spatVector = baysor_poly,
                             unique_ID_cache = names(baysor_poly),
                             name = "Baysor")


################### Co-Clustering Analysis: Object Creation ####################

xen_cell = createGiottoObjectSubcellular(
  gpoints = list(rna = gpoints_list$`Gene Expression`),
  gpolygons = list(cell = gpoly_cells),
  instructions = instrs
)

# Cell RNA
xen_cell = calculateOverlapRaster(xen_cell,
                                  spatial_info = 'cell',
                                  feat_info = 'rna')

xen_cell = overlapToMatrix(xen_cell,
                           poly_info = 'cell',
                           feat_info = 'rna',
                           name = 'raw')

cID_and_type_dt = data.table::fread("./xenium/FFPE_Human_Breast_Cancer/Rep1/Xenium_cell_IDs_and_types.csv")
cID_and_type_dt$Barcode = as.character(cID_and_type_dt$Barcode)


#prep for merge into joined gobject
xcm = getCellMetadata(xen_cell, spat_unit = "cell", output = "data.table")
xcm = data.table::merge.data.table(xcm, cID_and_type_dt, by.x = "cell_ID", by.y = "Barcode")

# Convert cell_IDs to integers for sorting so that
# the IDs align with the 10X cell types
xcm$cell_ID = as.integer(xcm$cell_ID)
setorder(xcm, cell_ID)
xcm$cell_ID = paste0("cell-", seq_along(xcm$cell_ID))

xen_baysor = createGiottoObjectSubcellular(
  gpoints = list(rna = gpoints_list$`Gene Expression`),
  gpolygons = list(cell = baysor_gpoly),
  instructions = instrs
)

# Baysor RNA
xen_baysor = calculateOverlapRaster(xen_baysor,
                                    spatial_info = "cell",
                                    feat_info = "rna")

xen_baysor = overlapToMatrix(xen_baysor,
                             name = "raw",
                             poly_info = "cell",
                             feat_info = "rna")

xen_baysor = addSpatialCentroidLocations(xen_baysor,
                                         poly_info = "cell",
                                         spat_loc_name = "cell")

#prep for merge into joined gobject
xbm = getCellMetadata(xen_baysor, spat_unit = "cell", output = "data.table")
xbm$cell_ID = paste0("Baysor-polygon_", seq_along(xbm$cell_ID))
xbm$Cluster = "BAYSOR UNDEFINED"

xm = rbind(xcm, xbm)

REMAKE_JOINT_OBJECT = FALSE

if(REMAKE_JOINT_OBJECT){


  # JOIN

  join_xen = joinGiottoObjects(gobject_list = list(xen_cell, xen_baysor),
                               gobject_names = c("cell", "Baysor"),
                               join_method = "no_change")

  HE_gobj = createGiottoLargeImage(terra::rast(HE_img_path),
                                   negative_y = FALSE,
                                   name = "HE_image")

  join_xen = addGiottoLargeImage(join_xen, list(HE_gobj))

  join_xen = filterGiotto(gobject = join_xen,
                          spat_unit = 'cell',
                          poly_info = 'cell',
                          expression_threshold = 1,
                          feat_det_in_min_cells = 3,
                          min_det_feats_per_cell = 5)

  join_xen = normalizeGiotto(gobject = join_xen,
                             spat_unit = 'cell',
                             scalefactor = 5000,
                             verbose = T)

  join_xen = addStatistics(join_xen, spat_unit = "cell")

  # Add cell typing information
  join_xen = addCellMetadata(join_xen,
                             spat_unit = "cell",
                             feat_type = "rna",
                             new_metadata = xm,
                             by_column = TRUE,
                             column_cell_ID = "cell_ID")


  ######### Co-Clustering Analysis: Dimension Reduction and Clustering ###########

  tictoc::tic()
  join_xen = runPCA(gobject = join_xen,
                    spat_unit = 'cell',
                    expression_values = 'scaled',
                    feats_to_use = NULL,
                    scale_unit = F,
                    center = F)
  tictoc::toc()

  #save.image(file = "joined_xenium_cell_baysor.RData")
  #load("joined_xenium_cell_baysor.RData")

  # Subset to 25% of all cells


  subset_25_percent = round(length(get_cell_id(join_xen)) * 0.25) # 25% 96,329
  # 10% --> 38,582

  my_seed_number = 315
  set.seed(my_seed_number)
  subset_IDs = sample(x = join_xen@cell_ID$cell, size = subset_25_percent)

  sub_join_xen = subsetGiotto(join_xen,
                              cell_ids = subset_IDs)


  tictoc::tic()
  join_xen = runUMAPprojection(join_xen,
                               spat_unit = "cell",
                               feat_type = "rna",
                               dimensions_to_use = 1:25,
                               n_threads = determine_cores(),
                               n_neighbors = 30,
                               min_dist = 0.1,
                               name = "umap.projection",
                               n_epochs = 200,
                               spread = 1,
                               random_subset = subset_25_percent,
                               set_seed = TRUE,
                               seed_number = my_seed_number)
  tictoc::toc()



  tictoc::tic()
  sub_join_xen = createNearestNetwork(sub_join_xen,
                                      dimensions_to_use = 1:25,
                                      k = 30)
  tictoc::toc()

  tictoc::tic()
  sub_join_xen = doLeidenClusterIgraph(gobject = sub_join_xen,
                                       spat_unit = "cell",
                                       feat_type = "rna",
                                       resolution_parameter = 0.55,
                                       n_iterations = 100,
                                       name = "sub_leiden_clus",
                                       set_seed = TRUE,
                                       seed_number = my_seed_number)
  tictoc::toc()

  tictoc::tic()
  # " classification is decided by majority vote, with ties broken at random. "
  join_xen = doClusterProjection(target_gobject = join_xen,
                                 source_gobject = sub_join_xen,
                                 spat_unit = "cell",
                                 feat_type = "rna",
                                 source_cluster_labels = "sub_leiden_clus",
                                 prob = FALSE,
                                 knn_k = 30,
                                 dimensions_to_use = 1:25)
  tictoc::toc()

  saveGiotto(gobject = join_xen, foldername = "joined_gobject", dir = results_folder)
} else{

  join_xen = loadGiotto("/projectnb/rd-spat/HOME/mobrien2/xenium/co_cluster_results/joined_gobject/")

}

plotUMAP(join_xen,
         dim_reduction_name = "umap.projection",
         cell_color = "knn_labels",
         save_param = list(base_width = 15,
                           base_height = 15),
         save_plot = F)

##################### Co-Clustering Analysis: Cell Typing ######################

### MANUAL ANNOTATION OF LEIDEN CLUSTERS

leiden_clusters = sort(as.integer(unique(join_xen@cell_metadata$cell$rna[]$knn_labels)))
leiden_clusters = as.character(leiden_clusters)

ct_c1 = "Stromal_and_Tumor"
ct_c2 = "Macrophages"
ct_c3 = "Myoepithelium_ACTA2+"
ct_c4 = "Endothelial"
ct_c5 = "DCIS"
ct_c6 = "Invasive_Tumor"
ct_c7 = "CD8+_T_Cells"
ct_c8 = "Stromal"
ct_c9 = "Mast_Cells"
ct_c10 = "CD4+_T_Cells"
ct_c11 = "B_Cells"
ct_c12 = "Dendritic_Cells"
ct_c13 = "Myoepithelium_KRT15+"

names(leiden_clusters) = c(ct_c1, ct_c2, ct_c3, ct_c4, ct_c5, ct_c6, ct_c7,
                           ct_c8, ct_c9, ct_c10, ct_c11, ct_c12, ct_c13)


ct_c_idx_map = data.table::data.table(leiden_clusters, names(leiden_clusters))
names(ct_c_idx_map)[[2]] = "mapped_type"

jcm = getCellMetadata(join_xen, output = "data.table")
jcm = merge(jcm, ct_c_idx_map, by.x = "knn_labels", by.y = "leiden_clusters")
new_jcm = jcm[, .(cell_ID, mapped_type)]
join_xen = addCellMetadata(join_xen,
                           spat_unit = "cell",
                           feat_type = "rna",
                           new_metadata = new_jcm,
                           by_column = TRUE,
                           column_cell_ID = "cell_ID")

plotUMAP(join_xen,
         dim_reduction_name = "umap.projection",
         cell_color = "mapped_type",
         save_param = list(base_width = 15,
                           base_height = 15),
         save_plot = F)

########################### Color setup for plotting ###########################

map_cell_meta = getCellMetadata(join_xen, output = "data.table")

my_colors = getDistinctColors(length(unique(map_cell_meta$mapped_type)))

names(my_colors) = unique(map_cell_meta$mapped_type)

################# Co-Clustering Analysis: Cell Type Bar Plots ##################

### ORIGINAL
map_og_cell_meta = map_cell_meta[map_cell_meta[, list_ID == "cell"]]


map_og_cell_types = unique(map_og_cell_meta$mapped_type)
map_og_ct_freq_dt = data.table::data.table(table(map_og_cell_meta$mapped_type))
colnames(map_og_ct_freq_dt) = c("cell_type", "num_cells")

og_total_cells = sum(map_og_ct_freq_dt$num_cells)

for ( i in map_og_cell_types){
  nullvar = map_og_ct_freq_dt[cell_type == i, perc := num_cells/og_total_cells * 100]
}

# # Pie Chart
# pl_og = ggplot2::ggplot(as.data.frame(map_og_ct_freq_dt), aes(x="", y=perc, fill = cell_type)) +
#   geom_bar(stat="identity", width = 1) +
#   coord_polar("y", start = 0) +
#   scale_fill_manual(values = my_colors) +
#   theme_void() +
#   labs(title = paste("Original Cell Types (", as.character(og_total_cells), " Cells)"))

### BAYSOR
map_br_cell_meta = map_cell_meta[map_cell_meta[, list_ID == "Baysor"]]


map_br_cell_types = unique(map_br_cell_meta$mapped_type)
map_br_ct_freq_dt = data.table::data.table(table(map_br_cell_meta$mapped_type))
colnames(map_br_ct_freq_dt) = c("cell_type", "num_cells")

br_total_cells = sum(map_br_ct_freq_dt$num_cells)

for ( i in map_br_cell_types){
  nullvar = map_br_ct_freq_dt[cell_type == i, perc := num_cells/br_total_cells * 100]
}

# # Pie Chart
# pl_br = ggplot2::ggplot(as.data.frame(map_br_ct_freq_dt), aes(x="", y=perc, fill = cell_type)) +
#   geom_bar(stat="identity", width = 1) +
#   coord_polar("y", start = 0) +
#   scale_fill_manual(values = my_colors) +
#   theme_void() +
#   labs(title = paste("Original Cell Types (", as.character(br_total_cells), " Cells)"))


map_br_ct_freq_dt[,segmentation := "Baysor"]
map_og_ct_freq_dt[,segmentation := "Original"]

br_for_chi = map_br_ct_freq_dt
og_for_chi = map_og_ct_freq_dt

setorder(map_br_ct_freq_dt, cols = -perc)
setorder(map_og_ct_freq_dt, cols = -perc)

dt_ct <- rbind(map_br_ct_freq_dt,
               map_og_ct_freq_dt)


test_dt = dt_ct
dec_order = unique(test_dt$cell_type)
rm(test_dt)

setorder(dt_ct, cols = -segmentation)
dt_ct$segmentation = factor(dt_ct$segmentation, c("Original", "Baysor"))

### Percentage Bar Plots
pdf("./s8E_percentage_cell_type_by_segmentation.pdf", width = 14, height = 6)
ggplot2::ggplot(dt_ct, aes( x = factor(cell_type, levels = dec_order), y = perc, fill = segmentation)) +
  ggplot2::geom_bar(stat = "identity",
                    position = position_dodge()) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = element_text(angle = 45),
                 axis.ticks.length.x =unit(1.5, "cm")) +
  ggplot2::labs(title = "Cell Types by Segmentation",
                x = "Cell Type",
                y = "Percentage of Segmented Cells")
dev.off()

ggplot2::ggsave(filename = "s8E_percentage_cell_type_by_segmentation.png",
                plot = last_plot(),
                device = "png",
                width = 14,
                height = 6,
                dpi = "retina")

### Bar Plot (stacked)
ggplot2::ggplot(dt_ct, aes( x = factor(cell_type, levels = dec_order),
                            y = num_cells, fill = segmentation_type)) +
  ggplot2::geom_bar(stat = "identity",
                    position = position_stack()) +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = element_text(angle = 45),
                 axis.ticks.length.x =unit(1.5, "cm")) +
  ggplot2::labs(title = "Cell Types by Segmentation",
                x = "Cell Type",
                y = "Number of Segmented Cells")

### Bar Plot (stacked) seg type x axis
ggplot2::ggplot(dt_ct, aes( x = factor(segmentation_type),
                            y = num_cells, fill = cell_type)) +
  ggplot2::geom_bar(stat = "identity",
                    position = position_stack()) +
  ggplot2::theme_classic() +
  ggplot2::scale_fill_manual(values = my_colors) +
  ggplot2::theme(axis.text.x = element_text(angle = 45),
                 axis.ticks.length.x =unit(1.5, "cm")) +
  ggplot2::labs(title = "Cell Types by Segmentation",
                x = "Cell Type",
                y = "Number of Segmented Cells")



### Bar Plot (stacked) seg type x axis, percentage y axis
pdf("./s8F_stacked_percentage_cell_type_by_segmentation.pdf", width = 8, height = 6)
ggplot2::ggplot(dt_ct, aes( x = factor(segmentation_type),
                            y = perc, fill = cell_type)) +
  ggplot2::geom_bar(stat = "identity",
                    position = position_stack()) +
  ggplot2::theme_classic() +
  ggplot2::scale_fill_manual(values = my_colors) +
  ggplot2::theme(axis.text.x = element_text(size = 10,
                                            angle = 45),
                 axis.ticks.length.x =unit(1.5, "cm")) +
  ggplot2::labs(title = "Cell Types by Segmentation",
                x = "Segmentation",
                y = "Percentage of Segmented Cells") +
  ggplot2::guides(fill=guide_legend(title="Cell Types"))
dev.off()


################# Co-Clustering Analysis: Visualize Cell Types #################

both_umaps = plotUMAP(join_xen,
                      dim_reduction_name = "umap.projection",
                      group_by = "list_ID",
                      cell_color = "mapped_type",
                      point_size = 1.25,
                      cell_color_code = my_colors,
                      save_plot = FALSE,
                      return_plot = TRUE)

ggplot2::ggsave(filename = "s8C_UMAPs_each_segmentation.pdf",
                plot = both_umaps,
                device = "pdf",
                width = 40,
                height = 20)

ggplot2::ggsave(filename = "s8C_UMAPs_each_segmentation.png",
                plot = both_umaps,
                device = "png",
                width = 30,
                height = 15,
                dpi = "retina")



both_spatPlots = spatPlot2D(join_xen,
                            group_by = "list_ID",
                            cell_color = "mapped_type",
                            cell_color_code = my_colors,
                            point_size = 0.7,
                            legend_text = 14,
                            axis_text = 10,
                            axis_title = 10,
                            save_param = list(base_width = 40,
                                              base_height = 15),
                            save_plot = F,
                            return_plot = T)

ggplot2::ggsave(filename = "s8B_spatPlot_cell_types_each_segmentation.pdf",
                plot = both_spatPlots,
                device = "pdf",
                width = 40,
                height = 15)

ggplot2::ggsave(filename = "s8B_spatPlot_cell_types_each_segmentation.png",
                plot = both_spatPlots,
                device = "png",
                width = 30,
                height = 15)

########################## Zoom To Stromal: ORIGINAL ##########################

# Pull out metadata from joint object pertinent to the original segmentation

cell_merge = jcm[grepl("cell-", jcm$cell_ID)]
cell_merge$list_ID = NULL
cell_merge$Cluster = NULL
xcm = data.table::merge.data.table(xcm, cell_merge, by = "cell_ID")

merge_in_cell = getCellMetadata(xen_cell, output = "data.table")
merge_in_cell[, join_ID := paste0("cell-", cell_ID)]
xcm = data.table::merge.data.table(merge_in_cell, xcm, by.x = 'join_ID', by.y = "cell_ID")
xcm$join_ID = NULL
setorder(xcm,cell_ID)

# Add this metadata to the original segmentation Giotto Object

xen_cell = addCellMetadata(xen_cell,
                           spat_unit = "cell",
                           feat_type = "rna",
                           new_metadata = xcm,
                           by_column = TRUE,
                           column_cell_ID = "cell_ID")

# Subset the original segmentation object

cell_ROI = subsetGiottoLocs(xen_cell,
                            spat_unit = "cell",
                            feat_type = "rna",
                            x_min = 1450,
                            x_max = 1650,
                            y_min = 1400,
                            y_max = 1600,
                            return_gobject = T)

roi_feats = list("ERBB2", "LUM", "CEACAM6")
roi_colors = getDistinctColors(length(roi_feats))
names(roi_colors) = roi_feats

roi_polys = unique(cell_ROI@cell_metadata$cell$rna[]$mapped_type)
roi_poly_colors = getDistinctColors(length(roi_polys))
names(roi_poly_colors) = roi_polys

cell_roi_pl = spatInSituPlotPoints(cell_ROI,
                                   spat_unit = "cell",
                                   feat_type = "rna",
                                   feats = roi_feats,
                                   feats_color_code = roi_colors,
                                   point_size = 0.1,
                                   background_color = "black",
                                   use_overlap = FALSE,
                                   show_polygon = TRUE,
                                   polygon_color = "white",
                                   polygon_line_size = 0.001,
                                   polygon_fill_as_factor = T,
                                   polygon_fill = "mapped_type",
                                   polygon_alpha = 0.5,
                                   polygon_fill_code = my_colors,
                                   plot_last = "points",
                                   show_legend = FALSE,
                                   save_plot = F,
                                   return_plot = T) +
  ggplot2::labs(title = "Original Segmentation")

pdf("./s8B_cell_seg_ROI.pdf", width = 9, height = 9)
cell_roi_pl
dev.off()

# ggplot2::ggsave(paste0(results_folder,"/cell_roi_stromal.png"),
#                 plot = cell_roi_pl,
#                 units = "in",
#                 width = 9,
#                 height = 9)

########################### Zoom To Stromal: BAYSOR ############################

# Pull out metadata from joint object pertinent to the Baysor segmentation

baysor_merge = jcm[grepl("Baysor-", jcm$cell_ID)]
baysor_merge$list_ID = NULL
baysor_merge$Cluster = NULL
xbm = data.table::merge.data.table(xbm, baysor_merge, by = "cell_ID")

merge_in_baysor = getCellMetadata(xen_baysor, output = "data.table")
merge_in_baysor[, join_ID := paste0("Baysor-", cell_ID)]
xbm = data.table::merge.data.table(merge_in_baysor, xbm, by.x = 'join_ID', by.y = "cell_ID")
xbm$join_ID = NULL
setorder(xbm,cell_ID)

# Add this metadata to the Baysor segmentation Giotto Object

xen_baysor = addCellMetadata(xen_baysor,
                             spat_unit = "cell",
                             feat_type = "rna",
                             new_metadata = xbm,
                             by_column = TRUE,
                             column_cell_ID = "cell_ID")

# Subset the Baysor segmentation object

baysor_ROI = subsetGiottoLocs(xen_baysor,
                              spat_unit = "cell",
                              feat_type = "rna",
                              x_min = 1450,
                              x_max = 1650,
                              y_min = 1400,
                              y_max = 1600,
                              return_gobject = T)

roi_feats = list("ERBB2", "LUM", "CEACAM6")
roi_colors = getDistinctColors(length(roi_feats))
names(roi_colors) = roi_feats

roi_polys = unique(cell_ROI@cell_metadata$cell$rna[]$mapped_type)
roi_poly_colors = getDistinctColors(length(roi_polys))
names(roi_poly_colors) = roi_polys

baysor_roi_pl = spatInSituPlotPoints(baysor_ROI,
                                     spat_unit = "cell",
                                     feat_type = "rna",
                                     feats = roi_feats,
                                     feats_color_code = roi_colors,
                                     point_size = 0.25,
                                     background_color = "black",
                                     use_overlap = FALSE,
                                     show_polygon = TRUE,
                                     polygon_color = "white",
                                     polygon_line_size = 0.001,
                                     polygon_fill_as_factor = T,
                                     polygon_fill = "mapped_type",
                                     polygon_alpha = 0.5,
                                     polygon_fill_code = my_colors,
                                     plot_last = "points",
                                     show_legend = FALSE,
                                     save_plot = F,
                                     return_plot = T) +
  ggplot2::labs(title = "Baysor Segmentation")

pdf("./s8B_baysor_seg_ROI.pdf", width = 9, height = 9)
baysor_roi_pl
dev.off()

# ggplot2::ggsave(paste0(results_folder,"/baysor_roi_stromal.png"),
#                 plot = baysor_roi_pl,
#                 units = "in",
#                 width = 9,
#                 height = 9)

merged_pl = cell_roi_pl + baysor_roi_pl

pdf("./both_segmentations_ROI.pdf", width = 9, height = 9)
merged_pl
dev.off()


# ggplot2::ggsave(paste0(results_folder,"/both_seg_roi_stromal.png"),
#                 plot = merged_pl,
#                 units = "in",
#                 width = 18,
#                 height = 9)


#################### Polygon Size Difference Quantification ####################

cell_areas = terra::expanse(gpoly_cells@spatVector)
baysor_areas = terra::expanse(baysor_gpoly@spatVector)

cell_area_dt = data.table::as.data.table(cell_areas)
baysor_area_dt = data.table::as.data.table(baysor_areas)

cell_area_dt$segmentation = "Original"
baysor_area_dt$segmentation = "Baysor"

names(cell_area_dt)[[1]] = names(baysor_area_dt)[[1]] = "polygonal_area"

polygonal_areas = rbind(baysor_area_dt, cell_area_dt)

# flip order for parallel structure
polygonal_areas$segmentation = factor(polygonal_areas$segmentation, c("Original", "Baysor"))


vi_pl = ggpubr::ggviolin(polygonal_areas, x = "segmentation", y = "polygonal_area", fill = "segmentation",
                         palette = c("#00AFBB", "#FC4E07"),
                         xlab = "Segmentation Method",
                         ylab = "Polygon Area",
                         add = "boxplot", add.params = list(size = 0.25),
                         #draw_quantiles = c(0.25,0.5,0.75),
                         ) +
  ggpubr::stat_compare_means(comparisons = list(c("Original", "Baysor")), label = "p.signif") + # Add significance levels
  ggpubr::stat_compare_means(label.x = 1.5,label.y = 125)

pdf("./s8D_violin_plot_poly_area.pdf", width = 5, height = 5)
vi_pl
dev.off()


# > mean(cell_areas)
# [1] 14.19084
# > mean(baysor_areas)
# [1] 4.876179


