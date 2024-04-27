## %######################################################%##
#                                                          #
####      Supplementary figure 7 Xeinium Co-Register        ####
#                                                          #
## %######################################################%##

## Download data
## Original Xenium data can be downloaded from 10X: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## Registered intermediate files
## Registered files can be downloaded from Zenodo in order to reproduce the figure: 10.5281/zenodo.11075079
## These files include Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_compressed.png(target coordinate system), CD20_aligned.csv.gz, DAPI_aligned.csv.gz, HER2_aligned.csv.gz, cell_boundaries_aligned.csv.gz, nucleus_boundaries_aligned.csv.gz, Xenium_Rep1_STalign_to_VisiumHE.csv

############################# Preprocess and Load ##############################

library(terra)
library(Giotto)

results_folder = '/projectnb/rd-spat/HOME/mobrien2/xenium/T_vs_IF_results'
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
# --- load IF Image
Xenium_IF_img = imager::load.image(IF_img_path)

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


#################### Raster Manipulation Function  #############################

replaceRasterNAWithZero <- function(layered_raster = NULL){

  num_layers = dim(layered_raster)[[3]]

  for (layer in seq(num_layers)){
    lay_rast = layered_raster[[layer]]
    idx_na = is.na(terra::values(lay_rast))
    terra::values(lay_rast)[idx_na] = 0
    layered_raster[[layer]] = lay_rast
  }
  return(layered_raster)
}

################### Aligned polygon loading: Baysor ############################

baysor_dt = data.table::fread("/projectnb/rd-spat/HOME/mobrien2/xenium/segmentations_aligned/Baysor_polygons_aligned.csv",
                              drop = c(1,2,3,4)) #drop index, old x, and old y columns

data.table::setnames(baysor_dt,
                     old = c("aligned_x", "aligned_y"),
                     new = c("x", "y"))

baysor_poly = dt_to_spatVector_polygon(baysor_dt)

baysor_poly$poly_ID = paste0("polygon_", 1:length(baysor_poly))

baysor_gpoly = giottoPolygon(spatVector = baysor_poly,
                             unique_ID_cache = names(baysor_poly),
                             name = "Baysor")

################ Aligned polygon loading: StarDist Nuclei ######################

sd_nuc_dt = data.table::fread("/projectnb/rd-spat/HOME/mobrien2/xenium/segmentations_aligned/aligned_stardist_nuclei.csv",
                              drop = c(1,4,5,9))

data.table::setnames(sd_nuc_dt, old = c("aligned_x", "aligned_y"), new = c("x", "y"))

sd_nuc_poly = dt_to_spatVector_polygon(sd_nuc_dt)

names(sd_nuc_poly)[[1]] = "poly_ID"

sd_nuc_gpoly = giottoPolygon(spatVector = sd_nuc_poly,
                             unique_ID_cache = names(sd_nuc_poly)[[1]],
                             name = "StarDist_Nuclei")

RELOAD_XENIUM_MULTIOMIC_GOBJECT = FALSE
if (RELOAD_XENIUM_MULTIOMIC_GOBJECT){
  ############################# Create Giotto Object #############################

  # Create Xenium_obj
  xen = createGiottoObjectSubcellular(
    gpoints = list(rna = gpoints_list$`Gene Expression`,
                   IF = gpoints_list$IF),
    gpolygons = list(cell = gpoly_cells,
                     Baysor = baysor_gpoly,
                     nucleus = gpoly_nucs,
                     StarDist_nuclei = sd_nuc_gpoly),
    instructions = instrs
  )

  HE_gobj = createGiottoLargeImage(terra::rast(HE_img_path),
                                   negative_y = FALSE)

  xen = addGiottoLargeImage(xen,list(HE_gobj))

  ################ Create, Filter, Normalize Expression Matrices #################

  # Cell RNA
  xen = calculateOverlapRaster(xen,
                               spatial_info = 'cell',
                               feat_info = 'rna')

  xen = overlapToMatrix(xen,
                        poly_info = 'cell',
                        feat_info = 'rna',
                        name = 'raw')

  # Nucleus RNA
  xen = calculateOverlapRaster(xen,
                               spatial_info = 'nucleus',
                               feat_info = 'rna')

  xen = overlapToMatrix(xen,
                        poly_info = 'nucleus',
                        feat_info = 'rna',
                        name = 'raw')

  # Cell IF
  xen = calculateOverlapRaster(xen,
                               spatial_info = 'cell',
                               feat_info = 'IF')

  xen = overlapToMatrix(xen,
                        poly_info = 'cell',
                        feat_info = 'IF',
                        name = 'raw')

  # StarDist Nuclei RNA
  xen = calculateOverlapRaster(xen,
                               spatial_info = "StarDist_nuclei",
                               feat_info = "rna")

  xen = overlapToMatrix(xen,
                        name = "StarDist_nuclei",
                        poly_info = "StarDist_nuclei",
                        feat_info = "rna")

  # Baysor RNA
  xen = calculateOverlapRaster(xen,
                               spatial_info = "Baysor",
                               feat_info = "rna")

  xen = overlapToMatrix(xen,
                        name = "Baysor",
                        poly_info = "Baysor",
                        feat_info = "rna")


  # Add spatlocs, filter, normalize

  xen = filterGiotto(gobject = xen,
                     spat_unit = 'cell',
                     poly_info = 'cell',
                     expression_threshold = 1,
                     feat_det_in_min_cells = 3,
                     min_det_feats_per_cell = 5)

  xen = normalizeGiotto(gobject = xen,
                        spat_unit = 'cell',
                        scalefactor = 5000,
                        verbose = T)

  xen = filterGiotto(gobject = xen,
                     spat_unit = 'nucleus',
                     poly_info = 'nucleus',
                     expression_threshold = 1,
                     feat_det_in_min_cells = 3,
                     min_det_feats_per_cell = 5)

  xen = normalizeGiotto(gobject = xen,
                        spat_unit = 'nucleus',
                        scalefactor = 5000,
                        verbose = T)

  xen = addSpatialCentroidLocations(xen,
                                    poly_info = "StarDist_nuclei",
                                    spat_loc_name = "StarDist_nuclei")

  xen = filterGiotto(xen,
                     spat_unit = "StarDist_nuclei",
                     feat_type = "rna",
                     expression_values = "StarDist_nuclei",
                     poly_info = "StarDist_nuclei",
                     expression_threshold = 1,
                     feat_det_in_min_cells = 3,
                     min_det_feats_per_cell = 5)

  xen = normalizeGiotto(xen,
                        expression_values = "StarDist_nuclei",
                        spat_unit = "StarDist_nuclei",
                        feat_type = "rna",
                        norm_methods = "standard",
                        scalefactor = 5000)

  xen = addSpatialCentroidLocations(xen,
                                    poly_info = "Baysor",
                                    spat_loc_name = "Baysor")

  xen = filterGiotto(xen,
                     spat_unit = "Baysor",
                     feat_type = "rna",
                     expression_values = "Baysor",
                     poly_info = "Baysor",
                     expression_threshold = 1,
                     feat_det_in_min_cells = 3,
                     min_det_feats_per_cell = 5)

  xen = normalizeGiotto(xen,
                        expression_values = "Baysor",
                        spat_unit = "Baysor",
                        feat_type = "rna",
                        norm_methods = "standard",
                        scalefactor = 5000)

  # Add Statistics to Metadata
  xen = addStatistics(xen, expression_values = 'raw')
  xen = addStatistics(xen, expression_values = 'raw', spat_unit = "nucleus")
  xen = addStatistics(xen, expression_values = 'raw', spat_unit = "cell", feat_type = "IF")
  xen = addStatistics(xen, expression_values = 'StarDist_nuclei', spat_unit = "StarDist_nuclei", feat_type = "rna")
  xen = addStatistics(xen, expression_values = 'Baysor', spat_unit = "Baysor", feat_type = "rna")

  #saveGiotto(gobject = xen, foldername = "gobjects")
} else{
  xen = loadGiotto("/projectnb/rd-spat/HOME/mobrien2/xenium/T_vs_IF_results/gobjects")
}

####################### Transcript vs IF: ERBB2 and HER2 #######################

# ERBB2 is translated into HER2

################# ERBB2 setup

TX_dt_ERBB2 = tx_dt_types$`Gene Expression`[tx_dt_types$`Gene Expression`$feat_ID == "ERBB2"]
TX_dt_ERBB2 = TX_dt_ERBB2[, .(feat_ID, x, y)]
TX_dt_ERBB2[, expr_source := "RNA"]

ERBB2_vect = terra::vect(TX_dt_ERBB2, geom = c('x','y'))
ERBB2_ext = ext(ERBB2_vect)

ERBB2_rast = terra::rast(ext = ERBB2_ext, nrow = 500, ncol = 700)
ERBB2_rast = terra::rasterize(ERBB2_vect, ERBB2_rast, fun = sum)


################# HER2 setup

save_original_HER2 = IF_dt_HER2

data.table::setnames(x = IF_dt_HER2,
                     old = c('IF_molecule', 'aligned_x', 'aligned_y','x'),
                     new = c('feat_ID', 'x', 'y','intensity'))

IF_dt_HER2 = IF_dt_HER2[, .(feat_ID, x , y, intensity)]
IF_dt_HER2[, expr_source := "IF"]

HER2_vect = terra::vect(IF_dt_HER2, geom = c('x','y'))
HER2_ext = ext(HER2_vect)

fill_rast = terra::rast(ext = ERBB2_ext, nrow = 500, ncol = 700)
HER2_rast = terra::rasterize(as.matrix(IF_dt_HER2[, .(x, y)]), y = fill_rast, values = as.vector(IF_dt_HER2$intensity))

################# Combine and Plot

combined_rast = c(HER2_rast, ERBB2_rast)
names(combined_rast) = c("HER2", "ERBB2")
plot(combined_rast)

final_combined_rast = replaceRasterNAWithZero(combined_rast)

her2_erbb2_gimg = createGiottoLargeImage(combined_rast, name = "ERBB2_HER2", ext = ext(combined_rast))
xen = addGiottoLargeImage(xen, list(her2_erbb2_gimg))

pearson_R = terra::layerCor(final_combined_rast, fun = "pearson")
# 0.5995345

dt = as.data.table(final_combined_rast)
dt_ERBB2 = dt$ERBB2
dt_HER2 = dt$HER2

lin_model = lm(formula = dt_ERBB2 ~ dt_HER2,
               data = as.data.frame(dt[, .(ERBB2, HER2)]))

plot(x = values(final_combined_rast$HER2),
     y = values(final_combined_rast$ERBB2),
     xlab = "HER2 Intensity",
     ylab = "ERBB2 Counts",
     main = "Correlation of HER2 and ERBB2 (R = 0.6239516)",
     pch = 20,
     font.main = 50,
     font.lab = 18,
     cex = 0.50,
     col = "darkgoldenrod")
abline(lin_model, col = "coral1", lwd = 2)


# View all available colors: GiottoVisuals::pal_names

png("s7B_her2_erbb2_r.png",width = 1500, height = 750)

par(mfrow = c(1,2))
plot(combined_rast$HER2,
     col = GiottoVisuals::getColors(pal = "viridis"),
     background = "black",
     plg = list(loc = "bottom",
                title = "Intensity",
                title.cex = 2.5,
                ext = c(750,2250,0,30),
                cex = 2,
                at = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)))

plot(combined_rast$ERBB2,
     col = GiottoVisuals::getColors(pal = "viridis"),
     background = "black",
     plg = list(loc = "bottom",
                title = "Counts",
                ext = c(750,2250,0,30),
                title.cex = 2.5,
                cex = 2))

mtext(sprintf("HER2 and ERBB2 (r = %s)",
              round(pearson_R$pearson[[2]],2)),
      side = 3,
      line = -4,
      cex = 4,
      outer = TRUE)

dev.off()

her2_glimg = createGiottoLargeImage(HER2_rast, name = "HER2_img", extent = HER2_ext)
xen = addGiottoLargeImage(gobject = xen,
                          largeImages = list(her2_glimg))

# ## Alternative, Giotto-native visualizations
# spatInSituPlotDensity(gobject = xen,
#                       feat_type = "rna",
#                       feats = c("ERBB2"),
#                       save_plot = F)
#
# spatInSituPlotHex(gobject = xen,
#                   feat_type = "rna",
#                   feats = c("ERBB2"),
#                   save_plot = F)

####################### Transcript vs IF: MS4A1 and CD20 #######################

# MS4A1 is translated into CD20

################# MS4A1 setup

TX_dt_MS4A1 = tx_dt_types$`Gene Expression`[tx_dt_types$`Gene Expression`$feat_ID == "MS4A1"]
TX_dt_MS4A1 = TX_dt_MS4A1[, .(feat_ID, x, y)]
TX_dt_MS4A1[, expr_source := "RNA"]

MS4A1_vect = terra::vect(TX_dt_MS4A1, geom = c('x','y'))
MS4A1_ext = ext(MS4A1_vect)

MS4A1_rast = terra::rast(ext = MS4A1_ext, nrow = 500, ncol = 700)
MS4A1_rast = terra::rasterize(MS4A1_vect, MS4A1_rast, fun = sum)

################# CD20 setup

save_original_CD20 = IF_dt_CD20

data.table::setnames(x = IF_dt_CD20,
                     old = c('IF_molecule', 'aligned_x', 'aligned_y','x'),
                     new = c('feat_ID', 'x', 'y','intensity'))

IF_dt_CD20 = IF_dt_CD20[, .(feat_ID, x , y, intensity)]
IF_dt_CD20[, expr_source := "IF"]

CD20_vect = terra::vect(IF_dt_CD20, geom = c('x','y'))
cropped_CD20_vect = terra::crop(CD20_vect, MS4A1_ext)
CD20_ext = ext(CD20_vect)


fill_rast_2 = terra::rast(ext = MS4A1_ext, nrow = 500, ncol = 700)
CD20_rast = terra::rasterize(as.matrix(IF_dt_CD20[, .(x, y)]), y = fill_rast_2, values = as.vector(IF_dt_CD20$intensity))

################# Combine and Plot

combined_rast_2 = c(CD20_rast, MS4A1_rast)
names(combined_rast_2) = c("CD20", "MS4A1")
plot(combined_rast_2)

final_combined_rast_2 = replaceRasterNAWithZero(layered_raster = combined_rast_2)

pearson_R_2 = terra::layerCor(final_combined_rast_2, fun = "pearson")
# 0.2297314

dt2 = as.data.table(final_combined_rast_2)
dt2_MS4A1 = dt2$MS4A1
dt2_CD20 = dt2$CD20

lin_model_2 = lm(formula = dt2_MS4A1 ~ dt2_CD20,
               data = as.data.frame(dt2[, .(MS4A1, CD20)]))

plot(x = terra::values(final_combined_rast_2$CD20),
     y = terra::values(final_combined_rast_2$MS4A1))
abline(lin_model_2)

################# Save Plots

png("s7B_cd20_ms4a1_r_both_rasters.png",width = 1400, height = 700)

par(mfrow = c(1,2))
plot(combined_rast_2$CD20,
     col = GiottoVisuals::getColors(pal = "viridis"),
     background = "black",
     plg = list(loc = "bottom",
                title = "Intensity",
                title.cex = 2.5,
                ext = c(750,2250,0,30),
                cex = 2))

plot(combined_rast_2$MS4A1,
     col = GiottoVisuals::getColors(pal = "viridis"),
     background = "black",
     plg = list(loc = "bottom",
                title = "Counts",
                title.cex = 2.5,
                ext = c(750,2250,0,30),
                cex = 2))

mtext(sprintf("CD20 and MS4A1 (r = %s)",
              round(pearson_R_2$pearson[[2]],2)),
      side = 3,
      line = -4,
      cex = 4,
      outer = TRUE)

dev.off()

png("s7B_cd20_raster_ms4a1_counts_r.png",width = 1400, height = 700)

par(mfrow = c(1,2))
plot(combined_rast_2$CD20,
     col = GiottoVisuals::getColors(pal = "viridis"),
     background = "black",
     plg = list(loc = "bottom",
                title = "Intensity",
                title.cex = 2.5,
                ext = c(750,2250,0,30),
                cex = 2))

plot(MS4A1_vect,
     ext = MS4A1_ext,
     cex = 0.33,
     col = "cyan",
     background = "black",
     alpha = 0.75)

mtext(sprintf("CD20 and MS4A1 (r = %s)",
              round(pearson_R_2$pearson[[2]],2)),
      side = 3,
      line = -4,
      cex = 4,
      outer = TRUE)

dev.off()


subset_plot_ext = ext(850, 1050, 820, 1160)


png("s7B_cd20_raster_ms4a1_counts_subset_r.png", width = 1000, height = 1000)
par(mfrow = c(1,2))
plot(combined_rast_2$CD20,
     ext = subset_plot_ext,
     col = GiottoVisuals::getColors(pal = "viridis"),
     background = "black",
     plg = list(loc = "bottom",
                title = "Intensity",
                title.cex = 2.5,
                cex = 2))

plot(MS4A1_vect,
     ext = subset_plot_ext,
     cex = 1.5,
     col = "cyan",
     background = "black",
     alpha = 0.75)

dev.off()



####################### Visualize Transcripts and IF  ##########################
par(mfrow = c(1,1))



cropped_HER2_rast = terra::crop(HER2_rast, ext(2100, 2230, 1200, 1300))
cropped_ERBB2_vect = terra::crop(ERBB2_vect, ext(2100, 2230, 1200, 1300))

her2_glimg_cropped = giottoLargeImage(name = "HER2",
                                      raster_object = cropped_HER2_rast,
                                      extent = ext(cropped_HER2_rast),
                                      overall_extent = ext(cropped_HER2_rast),
                                      min_intensity = cropped_HER2_rast@ptr$range_min,
                                      max_intensity = cropped_HER2_rast@ptr$range_max,
                                      is_int = F)

xen = addGiottoLargeImage(gobject = xen, largeImages = list(her2_glimg_cropped), spat_loc_name = "raw")


roi3 = subsetGiottoLocs(xen,
                        spat_unit = "cell",
                        feat_type = "rna",
                        x_min = 2100, x_max = 2230 , y_min = 1200, y_max = 1300)

spatInSituPlotPoints(roi3,
                     show_image = T,
                     gimage = c(HE_gobj, her2_glimg_cropped),
                     largeImage_name = c("image","HER2"),
                     feats = list('rna' = c("ERBB2")),
                     feat_type = c('rna'),
                     point_size = 0.25,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     use_overlap = F,
                     polygon_alpha = 0.15,
                     polygon_color = "white",
                     polygon_bg_color = "white",
                     polygon_fill_code = "white",
                     polygon_fill_as_factor = T,
                     polygon_line_size = 0.1,
                     background_color = "black",
                     plot_last = c("points"),
                     feats_color_code = c("red"),
                     feat_shape_code = c(16),
                     coord_fix_ratio = T)


####### FIGURE  2F
terra::plot(HER2_rast,
            col = "blue")

terra::plot(ERBB2_vect,
            col = "red", cex = 0.5, alpha = 0.6, add = T)


## Subset
terra::plot(her2_glimg_cropped@raster_object,
            col = GiottoVisuals::getColors(pal = "blues", rev = T)[50:100])
terra::plot(roi3@spatial_info$cell@spatVector,
            lwd = 0.5,
            add = T)
terra::plot(cropped_ERBB2_vect, col = "red", cex = 0.5, alpha = 0.6, add = T)

### Similar to Figure2F but with Baysor Polygons
# terra::plot(her2_glimg_cropped@raster_object,
#             col = GiottoVisuals::getColors(pal = "blues", rev = T)[50:100])
# terra::plot(roi3@spatial_info$Baysor@spatVector,
#             lwd = 0.5,
#             add = T)
# terra::plot(cropped_ERBB2_vect, col = "red", cex = 0.5, alpha = 0.4, add = T)





