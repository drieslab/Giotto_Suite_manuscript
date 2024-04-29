## %######################################################%##
#                                                          #
####      Supplementary figure 8 Multi-Segmentation        ####
#                                                          #
## %######################################################%##

## Download data
## Original Xenium data can be downloaded from 10X: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## This example script download orginal Xenium data to './Xenium/'
## Registered intermediate files
## Registered files can be downloaded from Zenodo in order to reproduce the figure: 10.5281/zenodo.11075079
## These files include Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_compressed.png(target coordinate system), CD20_aligned.csv.gz, DAPI_aligned.csv.gz, HER2_aligned.csv.gz, cell_boundaries_aligned.csv.gz, nucleus_boundaries_aligned.csv.gz, Xenium_Rep1_STalign_to_VisiumHE.csv"
## Place the downloaded aligned files into './Xenium/Aligned_Xe_rep1/' to use the script

############################# Preprocess and Load ##############################

library(terra)
library(patchwork)
library(ggplot2)
library(tictoc)
library(Giotto)

results_folder = './xenium/all_segmentations/'
setwd(results_folder)

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = FALSE,
                                  show_plot = TRUE)

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
# Xenium_HE_img = imager::load.image(HE_img_path) # this was never used

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

baysor_dt = data.table::fread("./Baysor_polygons_aligned.csv",
                              drop = c(1,2,3,4)) #drop index, old x, and old y columns

data.table::setnames(baysor_dt,
                     old = c("aligned_x", "aligned_y"),
                     new = c("x", "y"))

baysor_dt$poly_ID = paste0("polygon_", baysor_dt$geom)

baysor_gpoly = createGiottoPolygonsFromDfr(baysor_dt,
                                          name = "Baysor",
                                          calc_centroids = TRUE)

#################### Aligned polygon loading: CellPose #########################

cellpose_dt = data.table::fread("./aligned_cellpose_v3.csv",
                                drop = c(1,4,5,7))

data.table::setnames(cellpose_dt,
                     old = c("aligned_x", "aligned_y"),
                     new = c("x", "y"))

cellpose_dt$poly_ID = paste0("polygon_", cellpose_dt$geom)

cellpose_gpoly = createGiottoPolygonsFromDfr(cellpose_dt,
                                             name = "CellPose",
                                             calc_centroids = TRUE)

#################### Aligned polygon loading: StarDist #########################

stardist_dt = data.table::fread("./aligned_stardist_v3.csv",
                                drop = c(1,4,5,7,8))

data.table::setnames(stardist_dt,
                     old = c("aligned_x", "aligned_y"),
                     new = c("x", "y"))

stardist_dt$poly_ID = paste0("polygon_", stardist_dt$geom)

stardist_gpoly = createGiottoPolygonsFromDfr(stardist_dt,
                                             name = "CellPose",
                                             calc_centroids = TRUE)

################### Co-Clustering Analysis: Object Creation ####################

xen_multi = createGiottoObjectSubcellular(
  gpoints = list(rna = gpoints_list$`Gene Expression`),
  gpolygons = list(cell = gpoly_cells,
                   Baysor = baysor_gpoly,
                   CellPose = cellpose_gpoly,
                   StarDist = stardist_gpoly),
  instructions = instrs
)


# Zoom to ROI

xen_ROI = subsetGiottoLocs(gobject = xen_multi,
                           spat_unit = c("cell", "Baysor", "CellPose","StarDist"),
                           spat_loc_name = c("raw"),
                           feat_type = "rna",
                           poly_info = list(cell = "cell",
                                            Baysor = "Baysor",
                                            CellPose = "CellPose",
                                            StarDist = "StarDist"),
                           x_min = 2050,
                           x_max = 2350,
                           y_min = 1200,
                           y_max = 1500,
                           verbose = T)

og_ROI = spatInSituPlotPoints(xen_ROI,
                              show_image = FALSE,
                              feats = list('rna' = c("TACSTD2",
                                                     "CXCR4",
                                                     "ITGAX")),
                              feats_color_code = c("TACSTD2" = 'green',
                                                   'CXCR4' = 'blue',
                                                   'ITGAX' = 'red'),
                              point_size = 0.05,
                              show_polygon = TRUE,
                              use_overlap = FALSE,
                              background_color = "black",
                              polygon_feat_type = 'cell',
                              polygon_color = 'white',
                              polygon_line_size = 0.075,
                              show_legend = FALSE,
                              coord_fix_ratio = TRUE,
                              return_plot = TRUE) +
  ggplot2::labs(title = "Original")

cp_ROI = spatInSituPlotPoints(xen_ROI,
                              show_image = FALSE,
                              feats = list('rna' = c("TACSTD2",
                                                     "CXCR4",
                                                     "ITGAX")),
                              feats_color_code = c("TACSTD2" = 'green',
                                                   'CXCR4' = 'blue',
                                                   'ITGAX' = 'red'),
                              point_size = 0.05,
                              show_polygon = TRUE,
                              use_overlap = FALSE,
                              background_color = "black",
                              polygon_feat_type = 'CellPose',
                              polygon_color = 'pink',
                              polygon_line_size = 0.2,
                              show_legend = FALSE,
                              coord_fix_ratio = TRUE,
                              return_plot = TRUE)+
  ggplot2::labs(title = "CellPose")

sd_ROI = spatInSituPlotPoints(xen_ROI,
                              show_image = FALSE,
                              feats = list('rna' = c("TACSTD2",
                                                     "CXCR4",
                                                     "ITGAX")),
                              feats_color_code = c("TACSTD2" = 'green',
                                                   'CXCR4' = 'blue',
                                                   'ITGAX' = 'red'),
                              point_size = 0.05,
                              show_polygon = TRUE,
                              use_overlap = FALSE,
                              background_color = "black",
                              polygon_feat_type = 'StarDist',
                              polygon_color = 'goldenrod',
                              polygon_line_size = 0.15,
                              show_legend = FALSE,
                              coord_fix_ratio = TRUE,
                              return_plot = TRUE)+
  ggplot2::labs(title = "StarDist")

br_ROI = spatInSituPlotPoints(xen_ROI,
                              show_image = FALSE,
                              feats = list('rna' = c("TACSTD2",
                                                     "CXCR4",
                                                     "ITGAX")),
                              feats_color_code = c("TACSTD2" = 'green',
                                                   'CXCR4' = 'blue',
                                                   'ITGAX' = 'red'),
                              point_size = 0.05,
                              show_polygon = TRUE,
                              use_overlap = FALSE,
                              background_color = "black",
                              polygon_feat_type = 'Baysor',
                              polygon_color = 'lightblue',
                              polygon_line_size = 0.125,
                              show_legend = FALSE,
                              coord_fix_ratio = TRUE,
                              return_plot = TRUE) +
  ggplot2::ylim(1200,1500) +
  ggplot2::labs(title = "Baysor")

all_ROI = og_ROI + cp_ROI + sd_ROI + br_ROI

ggplot2::ggsave(filename = "./s8A_og_ROI.png",
                plot = og_ROI,
                device = "png",
                dpi = "retina",
                height = 5,
                width = 5)

ggplot2::ggsave(filename = "./s8A_cp_ROI.png",
                plot = cp_ROI,
                device = "png",
                dpi = "retina",
                height = 5,
                width = 5)

ggplot2::ggsave(filename = "./s8A_sd_ROI.png",
                plot = sd_ROI,
                device = "png",
                dpi = "retina",
                height = 5,
                width = 5)


ggplot2::ggsave(filename = "./s8A_br_ROI.png",
                plot = br_ROI,
                device = "png",
                dpi = "retina",
                height = 5,
                width = 5)

ggplot2::ggsave(filename = "./ALL_ROI.png",
                plot = all_ROI,
                device = "png",
                dpi = "retina",
                height = 10,
                width = 10)




