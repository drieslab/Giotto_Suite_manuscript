## %######################################################%##
#                                                          #
####      Supplementary figure 8 Multi-Segmentation        ####
#                                                          #
## %######################################################%##

## Download data
## Original Xenium data can be downloaded from 10X: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## This example script download orginal Xenium data to './Xenium/'
## Use the Baysor_run.sh to run Baysor and Place the Baysor output "segmentation_polygons.json" as "baysor.json" in './Segmentation_dir/'

############################# Preprocess and Load ##############################
## Assign the data directory
Xenium_dir <- paste0(data_dir,'/Xenium/')
Segmentation_dir <- paste0(data_dir,'/Segmentation_dir/')
library(Giotto)


# Load transcripts and default segmentations
x <- importXenium(paste0(Xenium_dir,'/outs/'))
x$qv <- 20 # default
tx <- x$load_transcripts(split_keyword =  list(
    c("BLANK"),
    c("NegControlCodeword"),
    c("NegControlProbe", "antisense")
))
cell <- x$load_polys()

#Flip vertically
aff <- affine() |> flip(direction = 'vertical')
cell <- affine(cell,aff)
tx_pts <- affine(tx[[1]]$rna,aff)
tx_pts
######################################################################## Multi_segmentation #########################################################################
# Segmentations Wrappers
IF_xen <- read10xAffineImage(file = paste0(Xenium_dir, "/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.ome.tif"),
                             imagealignment_path = paste0(Xenium_dir,"/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_imagealignment.csv"),
                             micron = 0.2125)

IF_xen <- IF_xen |> flip(direction = "vertical")
IF_xen <- IF_xen@funs$realize_magick(size = prod(dim(IF_xen)))
zoom <- terra::ext(6000,7000,3800,5100)
img_cropped <- terra::crop(IF_xen,y = zoom)
temp_rast <- terra::rast(ext = zoom, nrow = 1300, ncol = 1000)
IF_rast <- terra::resample(img_cropped@raster_object, temp_rast, method = "bilinear")
terra::writeRaster(IF_rast,filename = 'IF_cropped.tif',overwrite = T)

doMesmerSegmentation(
    Image_dir = 'IF_cropped.tif',
    mask_output = 'segmentation/mesmer.tif', 
    python_env = 'giotto_segmentation',
    Nucleus_channel = 3,
    Memberane_channel = 2,
    pixel_per_micron = 1)

doCellposeSegmentation(python_env = 'giotto_segmentation',
                       image_dir = 'IF_cropped.tif',
                       model_name = 'cyto3',
                       channel_1 = 2,
                       channel_2 = 3,
                       mask_output = 'segmentation/cellpose.tiff')


doStardistSegmentation(Image_dir = 'IF_cropped.tif',
                       python_env = 'giotto_segmentation',
                       mask_output = 'segmentation/stardist.tiff',
                       model_name = '2D_demo',
                       nuclei_channel = 3)

.rs.restartR()
stardist_poly <- createGiottoPolygonsFromMask('segmentation/stardist.tiff',
                                              shift_vertical_step = F,
                                              shift_horizontal_step = F,
                                              calc_centroids = T)
stardist_poly <- buffer(stardist_poly,5)


mesmer_poly <- createGiottoPolygonsFromMask('segmentation/mesmer.tif',
                                            shift_vertical_step = F,
                                            shift_horizontal_step = F,
                                            calc_centroids = T)

cellpose_poly <- createGiottoPolygonsFromMask('segmentation/cellpose.tiff',
                                              shift_vertical_step = F,
                                              shift_horizontal_step = F,
                                              calc_centroids = T)


affine_matrix <- matrix(c(1, 0, 6000, 
                          0, 1, 3800, 
                          0, 0, 1), 
                        nrow = 3, byrow = TRUE)
gpoly_mesmer_aligned <- affine(mesmer_poly,affine_matrix)
gpoly_cellpose_aligned <- affine(cellpose_poly,affine_matrix)
gpoly_stardist_aligned <- affine(stardist_poly,affine_matrix)

# BAYSOR
# No need to register as it was generated using transcript coordinates
Baysor_JSON = paste0(Segmentation_dir,'/baysor.json')
json_data = jsonlite::fromJSON(Baysor_JSON)
baysor_poly_df = json_data$geometries
baysor_poly_dt_list = list()
for (i in 1:(nrow(baysor_poly_df))){ 
    tmp_dt <- data.table::data.table()
    tmp_dt$x = json_data$geometries$coordinates[[i]][,,1]
    tmp_dt$y = json_data$geometries$coordinates[[i]][,,2]
    tmp_dt$geom = rep(i,length(json_data$geometries$coordinates[[i]][,,1]))
    tmp_dt$part = rep(0,length(json_data$geometries$coordinates[[i]][,,1]))
    tmp_dt$hole = rep(0,length(json_data$geometries$coordinates[[i]][,,1]))
    baysor_poly_dt_list[[i]] <- tmp_dt
}
baysor_poly_dt = data.table::rbindlist(baysor_poly_dt_list)

data.table::setnames(baysor_poly_dt,
                     old = 'geom',
                     new = 'poly_ID')
Baysor_gpoly = createGiottoPolygonsFromDfr(baysor_poly_dt,name = 'Baysor',calc_centroids = T)

####################################################################### Figure 8A #####################################################################
tx_ROI <- crop(tx_pts,zoom)
cell_ROI <- crop(cell,zoom)
Baysor_ROI <- crop(Baysor_gpoly,zoom)
mesmer_ROI <- crop(gpoly_mesmer_aligned,zoom)
cellpose_ROI <- crop(gpoly_cellpose_aligned,zoom)
stardist_ROI <- crop(gpoly_stardist_aligned,zoom)
## Create Giotto Xenium Object
xen_ROI <- createGiottoObjectSubcellular(gpoints = list('rna' = tx_ROI),
                                         gpolygons = list('cell' = cell_ROI,
                                                          'Baysor' = Baysor_ROI,
                                                          'CellPose' = cellpose_ROI,
                                                          'StarDist' = stardist_ROI,
                                                          'mesmer' = mesmer_ROI))

spatInSituPlotPoints(xen_ROI,
                     show_image = FALSE,
                     feats = list('rna' = c("TACSTD2", "CXCR4", "ITGAX")),
                     feats_color_code = c("TACSTD2" = 'green',
                                          'CXCR4' = 'blue',
                                          'ITGAX' = 'red'),
                     point_size = 0.02,
                     polygon_alpha = 0.5,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'CellPose',
                     polygon_color = 'pink',
                     polygon_line_size = 0.25,
                     coord_fix_ratio = TRUE)

spatInSituPlotPoints(xen_ROI,
                     show_image = FALSE,
                     feats = list('rna' = c("TACSTD2", "CXCR4", "ITGAX")),
                     feats_color_code = c("TACSTD2" = 'green',
                                          'CXCR4' = 'blue',
                                          'ITGAX' = 'red'),
                     point_size = 0.05,
                     polygon_alpha = 0.3,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'pink',
                     polygon_line_size = 0.25,
                     coord_fix_ratio = TRUE)

spatInSituPlotPoints(xen_ROI,
                     show_image = FALSE,
                     feats = list('rna' = c("TACSTD2", "CXCR4", "ITGAX")),
                     feats_color_code = c("TACSTD2" = 'green',
                                          'CXCR4' = 'blue',
                                          'ITGAX' = 'red'),
                     point_size = 0.05,
                     polygon_alpha = 0.3,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'StarDist',
                     polygon_color = 'pink',
                     polygon_line_size = 0.25,
                     coord_fix_ratio = TRUE)

spatInSituPlotPoints(xen_ROI,
                     show_image = FALSE,
                     feats = list('rna' = c("TACSTD2", "CXCR4", "ITGAX")),
                     feats_color_code = c("TACSTD2" = 'green',
                                          'CXCR4' = 'blue',
                                          'ITGAX' = 'red'),
                     point_size = 0.05,                     
                     polygon_alpha = 0.3,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'Baysor',
                     polygon_color = 'pink',
                     polygon_line_size = 0.25,
                     coord_fix_ratio = TRUE)



spatInSituPlotPoints(xen_ROI,
                     show_image = FALSE,
                     feats = list('rna' = c("TACSTD2", "CXCR4", "ITGAX")),
                     feats_color_code = c("TACSTD2" = 'green',
                                          'CXCR4' = 'blue',
                                          'ITGAX' = 'red'),
                     point_size = 0.05,                     
                     polygon_alpha = 0.3,
                     show_polygon = TRUE,
                     use_overlap = FALSE,
                     polygon_feat_type = 'mesmer',
                     polygon_color = 'pink',
                     polygon_line_size = 0.25,
                     coord_fix_ratio = TRUE)

####################################################################### Figure 8B #####################################################################
# Preprocess
xen_cell <- createGiottoObjectSubcellular(gpoints = list('rna' = tx_pts),
                                          gpolygons = list('cell' = cell))

xen_cell = calculateOverlapRaster(xen_cell,
                                  spatial_info = 'cell',
                                  feat_info = 'rna')

xen_cell = overlapToMatrix(xen_cell,
                           poly_info = 'cell',
                           feat_info = 'rna',
                           name = 'raw')

xen_baysor <- createGiottoObjectSubcellular(gpoints = list('rna' = tx_pts),
                                            gpolygons = list('cell' = Baysor_gpoly))
xen_baysor <- subsetGiottoLocs(xen_baysor,x_min = 0,y_min = 0,x_max = 9000, y_max = 7000) # remove invalid geoms
xen_baysor = calculateOverlap(xen_baysor,
                                    spatial_info = "cell",
                                    feat_info = "rna")

xen_baysor = overlapToMatrix(xen_baysor,
                             name = "raw",
                             poly_info = "cell",
                             feat_info = "rna")



# Joined process
join_xen = joinGiottoObjects(gobject_list = list(xen_cell, xen_baysor),
                             gobject_names = c("cell", "Baysor"),
                             join_method = "no_change")


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



my_seed_number = 2024
set.seed(my_seed_number)
subset_25_percent = round(nrow(getCellMetadata(join_xen)) * 0.25) 
subset_IDs = sample(x = join_xen@cell_ID$cell, size = subset_25_percent)

join_xen = runPCAprojection(gobject = join_xen,
                            spat_unit = 'cell',
                            expression_values = 'scaled',
                            feats_to_use = NULL,
                            scale_unit = F,
                            center = F)

join_xen = runUMAPprojection(join_xen,
                             spat_unit = "cell",
                             feat_type = "rna",
                             dim_reduction_name = 'pca.projection',
                             dimensions_to_use = 1:25,
                             n_neighbors = 30,
                             min_dist = 0.15,
                             name = "umap.projection",
                             n_epochs = 200,
                             spread = 1,
                             random_subset = subset_25_percent,
                             set_seed = TRUE,
                             seed_number = my_seed_number)


join_xen = createNearestNetwork(join_xen,
                                dim_reduction_name = 'pca.projection',
                                dimensions_to_use = 1:25,
                                k = 30)

join_xen = doLeidenClusterIgraph(gobject = join_xen,
                                 spat_unit = "cell",
                                 feat_type = "rna",
                                 resolution = 0.5, 
                                 n_iterations = 100,
                                 name = "leiden_clus",
                                 set_seed = TRUE,
                                 seed_number = my_seed_number) 
dimPlot2D(join_xen,
          dim_reduction_name = "umap.projection",
          cell_color = "leiden_clus")




fin_clusters = sort(as.integer(unique(join_xen@cell_metadata$cell$rna[]$leiden_clus)))
fin_clusters = as.character(fin_clusters)

ct_c1 = "Stromal_and_Tumor"
ct_c2 = "Invasive_Tumor"
ct_c3 = "Stromal_and_Tumor"
ct_c4 = "Stromal"
ct_c5 = "DCIS"
ct_c6 = "Macrophages"
ct_c7 = "B_Cells"
ct_c8 = "Mast_Cells"
ct_c9 = "CD8+_T_Cells"
ct_c10 = "Myoepithelium_ACTA2+"
ct_c11 = "Endothelial"
ct_c12 = "Myoepithelium_KRT15+"
ct_c13 = "CD4+_T_Cells"
ct_c14 = "CD4+_T_Cells"
ct_c15 = "Stromal"
ct_c16 = "Dendritic_Cells"
ct_c17 = "Invasive_Tumor"

annotation = c(ct_c1, ct_c2, ct_c3, ct_c4, ct_c5, ct_c6, ct_c7,
                           ct_c8, ct_c9, ct_c10, ct_c11, ct_c12, ct_c13,ct_c14,ct_c15,ct_c16,ct_c17)
names(annotation) = 1:17

join_xen <- annotateGiotto(join_xen,
                           annotation_vector = annotation,
                           cluster_column = 'leiden_clus',
                           name = 'mapped_type')


map_cell_meta = getCellMetadata(join_xen, output = "data.table")
my_colors = getDistinctColors(length(unique(map_cell_meta$mapped_type)))
names(my_colors) = unique(map_cell_meta$mapped_type)


####################################################################### Figure 8B #####################################################################
original_cell_meta = map_cell_meta[map_cell_meta$list_ID == 'cell']
original_cell_meta$cell_ID = sub(".*-", "", original_cell_meta$cell_ID)
xen_cell <- subsetGiotto(xen_cell,cell_ids = original_cell_meta$cell_ID)
xen_cell <- addCellMetadata(xen_cell,
                            spat_unit = 'cell',
                            feat_type = 'rna',
                            original_cell_meta,
                            by_column = TRUE,
                            column_cell_ID = "cell_ID")
spatPlot2D(xen_cell,
           cell_color = "mapped_type",
           cell_color_code = my_colors,point_border_stroke = 0,
           point_size = 1e-5,
           background_color = 'black',
           title = 'cell')

baysor_cell_meta = map_cell_meta[map_cell_meta$list_ID == 'Baysor']
baysor_cell_meta$cell_ID = sub(".*-", "", baysor_cell_meta$cell_ID)
xen_baysor <- addCellMetadata(xen_baysor,
                            spat_unit = 'cell',feat_type = 'rna',
                            baysor_cell_meta,
                            by_column = TRUE,
                            column_cell_ID = "cell_ID")
spatPlot2D(xen_baysor,
           cell_color = "mapped_type",
           cell_color_code = my_colors,point_border_stroke = 0,
           point_size = 1e-5,
           background_color = 'black',
           title = 'cell')

###########################------------- Zoomed In
roi_feats = list("ERBB2", "LUM", "CEACAM6")
roi_colors = getDistinctColors(length(roi_feats))
names(roi_colors) = roi_feats

baysor_zoom <- subsetGiottoLocs(xen_baysor,
                                x_min = 3900,
                                x_max = 4500,
                                y_min = 4600,
                                y_max = 5200)

spatInSituPlotPoints(baysor_zoom,
                     spat_unit = "cell",
                     feat_type = "rna",
                     feats = list('rna' = c("ERBB2", "LUM", "CEACAM6")),
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
                     save_plot = F)

cell_zoom <- subsetGiottoLocs(xen_cell,
                              x_min = 3900,
                              x_max = 4500,
                              y_min = 4600,
                              y_max = 5200)

spatInSituPlotPoints(cell_zoom,
                     spat_unit = "cell",
                     feat_type = "rna",
                     feats = list('rna' = c("ERBB2", "LUM", "CEACAM6")),
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
                     save_plot = F)
####################################################################### Figure 8C #####################################################################
dimPlot2D(join_xen,
          dim_reduction_name = "umap.projection",
          group_by = 'list_ID',
          show_center_label  = F,
          cell_color = "mapped_type",point_shape = 'no_border',
          cell_color_code = my_colors)

####################################################################### Figure 8D #####################################################################
cell_areas = terra::expanse(cell@spatVector)
baysor_areas = terra::expanse(Baysor_gpoly@spatVector)

cell_area_dt = data.table::as.data.table(cell_areas)
baysor_area_dt = data.table::as.data.table(baysor_areas)

cell_area_dt$segmentation = "Original"
baysor_area_dt$segmentation = "Baysor"

names(cell_area_dt)[[1]] = names(baysor_area_dt)[[1]] = "polygonal_area"

polygonal_areas = rbind(baysor_area_dt, cell_area_dt)

# flip order for parallel structure
polygonal_areas$segmentation = factor(polygonal_areas$segmentation, c("Original", "Baysor"))


ggpubr::ggviolin(polygonal_areas, x = "segmentation", y = "polygonal_area", fill = "segmentation",
                         palette = c("#00AFBB", "#FC4E07"),
                         xlab = "Segmentation Method",
                         ylab = "Polygon Area",
                         add = "boxplot", add.params = list(size = 0.25)) +
    ggpubr::stat_compare_means(comparisons = list(c("Original", "Baysor")), label = "p.signif") +
    ggpubr::stat_compare_means(label.x = 1.5,label.y = 125)

####################################################################### Figure 8E #####################################################################
# Calculate counts and percentages
library(tidyverse)
df <- as.data.frame(map_cell_meta) 
df$list_ID[df$list_ID == 'cell'] <- 'original'

counts <- df %>%
    group_by(list_ID, mapped_type) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    group_by(list_ID) %>%
    mutate(total_count = sum(count),
           percentage = (count / total_count) * 100) %>%
    ungroup() %>%
    group_by(list_ID) %>%
    mutate(mapped_type = reorder(mapped_type, -percentage))
ggplot(counts, aes(x = mapped_type, y = percentage, fill = list_ID)) +
    geom_bar(stat = "identity",
                      position = position_dodge()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45),
                   axis.ticks.length.x =unit(1.5, "cm")) +
    labs(title = "Cell Types by Segmentation",
                  x = "Cell Type",
                  y = "Percentage of Segmented Cells")

####################################################################### Figure 8F #####################################################################
ggplot2::ggplot(counts, aes( x = factor(list_ID),
                            y = percentage, fill = mapped_type)) +
    ggplot2::geom_bar(stat = "identity",
                      position = position_stack()) +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(values = my_colors) +
    ggplot2::theme(axis.text.x = element_text(size = 10),
                   axis.ticks.length.x =unit(1.5, "cm")) +
    ggplot2::labs(title = "Cell Types by Segmentation",
                  x = "Segmentation",
                  y = "Percentage of Segmented Cells") +
    ggplot2::guides(fill=guide_legend(title="Cell Types"))





