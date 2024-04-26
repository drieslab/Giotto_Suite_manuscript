library(terra)

#----------Part 1 Process Image Registration and Load data pathes ----------
library(Giotto)
#Set up Data directory
xenium_folder = './Image_Registration/Xenium_register/outs_rep1/'
aligned_folder = './Image_Registration/Aligned_Xe_rep1/'

### Load High Resolution HE image(which is also used as a target coordinate system for image registration)
images_folder = '/Image_Registration/Xenium_register/Images/'
##Note we used a 1:10 compressed image to
HE_img_path = paste0(images_folder,'Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_compressed.png')

#Convert IF images to a point data table, then register that to H&E Image, example can be found in 1st part of fig2B_S7A_Register_Xenium_IF_Visium_to_HE.ipynb
IF_CD20_path = paste0(aligned_folder, 'CD20_aligned.csv.gz')
IF_DAPI_path = paste0(aligned_folder, 'DAPI_aligned.csv.gz')
IF_HER2_path = paste0(aligned_folder, 'HER2_aligned.csv.gz')

#Register the Xenium assay to H&E Image, example can be found in 2nd part of fig2B_S7A_Register_Xenium_IF_Visium_to_HE.ipynb
cell_bound_path = paste0(aligned_folder, 'cell_boundaries_aligned.csv.gz')
nuc_bound_path = paste0(aligned_folder, 'nucleus_boundaries_aligned.csv.gz')
tx_path = paste0(aligned_folder, 'tx_aligned.csv.gz')
feat_meta_path = paste0(xenium_folder, 'cell_feature_matrix/features.tsv.gz') # This can be found in the Xenium output bundle

##Get Xenium to Visium Register information, in this figure we just compare which registered cells can be captured by Visium assay. Example script can be found in figS7A_Xenium_to_VisiumHE.ipynb
xen_df_path = "./Image_Registration/xenium_rep1_to_visium/Xenium_Rep1_STalign_to_VisiumHE.csv"
#Visium data is raw
visium_dir = "./data/Visium/"

#----------Part 2 Read and process all segmented results in one----------

#--- Place in one Xenium Obj and visualize
# 1. set working directory
# results_folder = 'path/to/result'
results_folder = './Coregister/final'

# 2. set giotto python path
python_path = '.conda/envs/giotto/bin/python'

# 3. create giotto instructions
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = 'python_path')


#--- load HE image
#Xenium_HE_img = imager::load.image(HE_img_path)
#plot(Xenium_HE_img)
HE_gobj = createGiottoLargeImage(terra::rast(HE_img_path),
                                 negative_y = FALSE)

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


#--------------IF
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

#plot(gpoints_list$`Gene Expression`,
#     feats = c('KRT8', 'MS4A1'))
#plot(gpoints_list$IF,
#     feats = c('HER2', 'CD20'))


# Create Xenium_obj All in One
xenium_gobj = createGiottoObjectSubcellular(
  gpoints = list(rna = gpoints_list$`Gene Expression`,
                 IF = gpoints_list$IF),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  largeImages = list(HE = HE_gobj),
  instructions = instrs
)


# Process Xenium Objects
# Calculate overlapped data and
feat_type = c('IF','rna')
polygon_type = c('cell','nucleus')

for (feat in feat_type){
         cat("Generating overlapping data \n")
         cat("Start Processing Feature type: ", paste0(feat),'\n')
              for (polygon in polygon_type){
                cat("Processing Polygon type: ", paste0(feat),'\n')
                xenium_gobj = calculateOverlapRaster(xenium_gobj,
                                                     spatial_info = polygon,
                                                     feat_info = feat)
                xenium_gobj = overlapToMatrix(xenium_gobj,
                                              poly_info = polygon,
                                              feat_info = feat,
                                              name = 'raw')
                cat("Normalzie data for feature: ", paste0(feat),'with spatial unit: ',paste0(polygon))
                xenium_gobj = normalizeGiotto(xenium_gobj,
                                              spat_unit = polygon,
                                              feat_type = feat)
                cat("Add Statistics data for feature: ", paste0(feat),'with spatial unit: ',paste0(polygon),'\n\n')
                xenium_gobj = addStatistics(xenium_gobj,
                                            feat_type = feat,
                                            spat_unit = polygon)
              }
}

#showGiottoSpatialInfo(xenium_gobj)
#showGiottoSpatLocs(xenium_gobj)
#showGiottoExpression(xenium_gobj)
#showGiottoCellMetadata(xenium_gobj)
#showGiottoFeatMetadata(xenium_gobj)


HE_gobj = createGiottoLargeImage(terra::rast(HE_img_path),
                                 negative_y = FALSE)
xenium_gobj = addGiottoLargeImage(xenium_gobj,list(HE_gobj))

pl = spatInSituPlotPoints(xenium_gobj,
                     show_image = T,largeImage_name = 'image',
                     feats = list('IF' = c("CD20","HER2")),
                     feat_type = c('IF'),
                     point_size = 0.05,
                     show_polygon = FALSE,
                     feats_color_code = c('CD20' = 'green','HER2' = 'brown'),
                     return_plot = T,
                     coord_fix_ratio = T)


# Plot the overlapped xenium to visium cells
xen_df = data.table::fread(xen_df_path)
xen_df$visium_x = xen_df$aligned_x/0.092717074
xen_df$visium_y = - xen_df$aligned_y/0.092717074
library(data.table)
library(ggplot2)
library(quantreg)

## Visium Obj creation
Visium = createGiottoVisiumObject(visium_dir = visium_dir,
                                  expr_data = 'raw',
                                  png_name = 'tissue_hires_image.png',
                                  gene_column_index = 2,
                                  instructions = instrs)

# Fit quantile regressions get boundary information from projected xenium
dt = getSpatialLocations(Visium,output = 'data.table')
colnames(dt)[1:2] = c('visium_x','visium_y')
q_low <- rq(visium_y ~ visium_x, data = dt, tau = 0)
q_up <- rq(visium_y ~ visium_x, data = dt, tau = 1)
q_left <- rq(visium_x ~ visium_y, data = dt, tau = 0)
q_right <- rq(visium_x ~ visium_y, data = dt, tau = 1)

# Define conditions using the lines
condition_low_xen <- xen_df$visium_y > predict(q_low, newdata = xen_df)
condition_up_xen <- xen_df$visium_y < predict(q_up, newdata = xen_df)
condition_left_xen <- xen_df$visium_x > predict(q_left, newdata = xen_df[, .(visium_y, visium_x)])
condition_right_xen <- xen_df$visium_x < predict(q_right, newdata = xen_df[, .(visium_y, visium_x)])
# Subset xen_df based on the conditions
subset_xen_df <- xen_df[condition_low_xen & condition_up_xen & condition_left_xen & condition_right_xen,]
subset_xen_overlap = getSpatialLocations(xenium_gobj,spat_unit = 'cell',output = 'data.table')
subset_xen_overlap = subset_xen_overlap[cell_ID %in% subset_xen_df$cell_id]

#ggplot() + geom_point(data = subset_xen_overlap, aes(x = sdimx, y = sdimy), size = 0.2, alpha = 0.25,color= 'orange')

pl = spatInSituPlotPoints(xenium_gobj,
                          show_image = T,largeImage_name = 'image',
                          feats = list('IF' = c("CD20","HER2")),
                          feat_type = c('IF'),
                          point_size = 0.05,
                          show_polygon = FALSE,
                          show_legend = FALSE,
                          feats_color_code = c('CD20' = 'green','HER2' = 'red'),
                          return_plot = T,
                          coord_fix_ratio = T)

pl_final = pl + geom_point(data = subset_xen_overlap, aes(x = sdimx, y = sdimy), size = 0.2, alpha = 0.25,color= 'orange')
pl_final
