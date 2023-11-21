##%######################################################%##
#                                                          #
####      Supplementary figure 3 DBiT-seq dataset       ####
#                                                          #
##%######################################################%##


## --------------------------------------------------------------------------------------
## Download dataset
## The mouse embryo E10.5 dataset was created by [Liu, et al 2020](https://www.sciencedirect.com/science/article/pii/S0092867420313908?via%3Dihub) and downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137986

## Load data
## --------------------------------------------------------------------------------------
## RNA
rna_expression = read.table("data/GSE137986_RAW/GSM4189613_0702cL.tsv.gz",
                            sep = "\t",
                            header = TRUE)
rownames(rna_expression) = rna_expression$X

## --------------------------------------------------------------------------------------
## Protein
protein_expression = read.table("data/GSE137986_RAW/GSM4202307_0702aL.tsv.gz",
                                sep = "\t",
                                header = TRUE)
rownames(protein_expression) = protein_expression$X


## --------------------------------------------------------------------------------------
rna_expression = t(rna_expression[-1])
protein_expression = t(protein_expression[-1])


## --------------------------------------------------------------------------------------
spatial_coords = data.frame(cell_ID = colnames(rna_expression))
spatial_coords = cbind(spatial_coords,
                       tidyr::separate(spatial_coords, cell_ID, c("x","y"), sep = "x"))

spatial_coords$x = as.numeric(spatial_coords$x)
spatial_coords$y = as.numeric(spatial_coords$y)*-1


## Create Giotto object
## --------------------------------------------------------------------------------------
library(Giotto)

save_dir = 'results'
instructions = createGiottoInstructions(save_dir = save_dir,
                                        save_plot = TRUE,
                                        show_plot = TRUE)

giottoObject = createGiottoObject(expression = list(raw = rna_expression,
                                                    raw = protein_expression),
                                  expression_feat = c('rna', 'protein'),
                                  spatial_locs = spatial_coords,
                                  instructions = instructions)

## Filter
## --------------------------------------------------------------------------------------
## RNA
giottoObject = filterGiotto(gobject = giottoObject,
                            expression_threshold = 1,
                            feat_det_in_min_cells = 1,
                            min_det_feats_per_cell = 1,
                            expression_values = c('raw'),
                            verbose = TRUE)

# Number of cells removed:  0  out of  901
# Number of feats removed:  476  out of  22846

## Protein
giottoObject = filterGiotto(gobject = giottoObject,
                            spat_unit = 'cell',
                            feat_type = 'protein',
                            expression_threshold = 1,
                            feat_det_in_min_cells = 1,
                            min_det_feats_per_cell = 1,
                            expression_values = c('raw'),
                            verbose = TRUE)

# Number of cells removed:  0  out of  901
# Number of feats removed:  0  out of  22

## Normalize
## --------------------------------------------------------------------------------------
## RNA
giottoObject = normalizeGiotto(gobject = giottoObject,
                               scalefactor = 6000,
                               verbose = TRUE)

## Protein
giottoObject = normalizeGiotto(gobject = giottoObject,
                               spat_unit = 'cell',
                               feat_type = 'protein',
                               scalefactor = 6000,
                               verbose = TRUE)

## Add statistics
## --------------------------------------------------------------------------------------
## RNA
giottoObject = addStatistics(gobject = giottoObject)

spatPlot2D(giottoObject,
           spat_unit = 'cell',
           feat_type = 'rna',
           cell_color = "nr_feats",
           color_as_factor = FALSE,
           point_size = 3.5)

spatPlot2D(giottoObject,
           cell_color = "total_expr",
           color_as_factor = FALSE,
           point_size = 3.5)

## Protein
giottoObject = addStatistics(gobject = giottoObject,
                             spat_unit = 'cell',
                             feat_type = 'protein')

spatPlot2D(giottoObject,
           spat_unit = 'cell',
           feat_type = 'protein',
           cell_color = "total_expr",
           color_as_factor = FALSE,
           point_size = 3.5)

## Calculate HVF
## --------------------------------------------------------------------------------------
giottoObject = calculateHVF(gobject = giottoObject)

## PCA
## --------------------------------------------------------------------------------------
# RNA
giottoObject <- runPCA(gobject = giottoObject)

## --------------------------------------------------------------------------------------
screePlot(giottoObject, ncp = 30)

plotPCA(gobject = giottoObject)


## --------------------------------------------------------------------------------------
# Protein
giottoObject <- runPCA(gobject = giottoObject,
                       spat_unit = 'cell',
                       feat_type = 'protein')

## --------------------------------------------------------------------------------------
screePlot(giottoObject,
          spat_unit = 'cell',
          feat_type = 'protein',
          ncp = 30)

plotPCA(gobject = giottoObject,
        spat_unit = 'cell',
        feat_type = 'protein')

## UMAP
## --------------------------------------------------------------------------------------
# RNA
giottoObject <- runUMAP(giottoObject,
                        dimensions_to_use = 1:10)


## --------------------------------------------------------------------------------------
plotUMAP(gobject = giottoObject)


## --------------------------------------------------------------------------------------
# Protein
giottoObject <- runUMAP(giottoObject,
                        spat_unit = 'cell',
                        feat_type = 'protein',
                        dimensions_to_use = 1:10)


## --------------------------------------------------------------------------------------
plotUMAP(gobject = giottoObject,
         spat_unit = 'cell',
         feat_type = 'protein')

## Clustering
## --------------------------------------------------------------------------------------
# RNA
giottoObject <- createNearestNetwork(gobject = giottoObject,
                                     dimensions_to_use = 1:10,
                                     k = 30)

giottoObject <- doLeidenCluster(gobject = giottoObject,
                                resolution = 1,
                                n_iterations = 1000)

# Protein
giottoObject <- createNearestNetwork(gobject = giottoObject,
                                     spat_unit = 'cell',
                                     feat_type = 'protein',
                                     dimensions_to_use = 1:10,
                                     k = 30)

giottoObject <- doLeidenCluster(gobject = giottoObject,
                                spat_unit = 'cell',
                                feat_type = 'protein',
                                resolution = 1,
                                n_iterations = 1000)


## --------------------------------------------------------------------------------------
# RNA
plotUMAP(gobject = giottoObject,
         cell_color = 'leiden_clus',
         show_NN_network = FALSE,
         point_size = 2,
         title = "",
         axis_text = 14,
         axis_title = 18,
         legend_text = 14)

# Protein
plotUMAP(gobject = giottoObject,
         spat_unit = 'cell',
         feat_type = 'protein',
         cell_color = 'leiden_clus',
         show_NN_network = FALSE,
         point_size = 2,
         title = "",
         axis_text = 14,
         axis_title = 18,
         legend_text = 14)


## --------------------------------------------------------------------------------------
# RNA
spatPlot2D(gobject = giottoObject,
           show_image = FALSE,
           cell_color = 'leiden_clus',
           point_size = 3.5,
           title = "",
           axis_text = 14,
           axis_title = 18,
           legend_text = 14)

# Protein
spatPlot2D(gobject = giottoObject,
           spat_unit = 'cell',
           feat_type = 'protein',
           show_image = FALSE,
           cell_color = 'leiden_clus',
           point_size = 3.5,
           title = "",
           axis_text = 14,
           axis_title = 18,
           legend_text = 14)

## Multi-omics integration
## --------------------------------------------------------------------------------------
# RNA
giottoObject <- createNearestNetwork(gobject = giottoObject,
                                     type = 'kNN',
                                     dimensions_to_use = 1:10,
                                     k = 20)

# Protein
giottoObject <- createNearestNetwork(gobject = giottoObject,
                                     spat_unit = 'cell',
                                     feat_type = 'protein',
                                     type = 'kNN',
                                     dimensions_to_use = 1:10,
                                     k = 20)

## --------------------------------------------------------------------------------------
giottoObject <- runWNN(giottoObject,
                       spat_unit = "cell",
                       modality_1 = "rna",
                       modality_2 = "protein",
                       pca_name_modality_1 = "pca",
                       pca_name_modality_2 = "protein.pca",
                       k = 20,
                       verbose = TRUE)

## --------------------------------------------------------------------------------------
giottoObject <- runIntegratedUMAP(giottoObject,
                                  modality1 = "rna",
                                  modality2 = "protein",
                                  spread = 7,
                                  min_dist = 1,
                                  force = FALSE)

## --------------------------------------------------------------------------------------
giottoObject <- doLeidenCluster(gobject = giottoObject,
                                spat_unit = "cell",
                                feat_type = "rna",
                                nn_network_to_use = "kNN",
                                network_name = "integrated_kNN",
                                name = "integrated_leiden_clus",
                                resolution = 1)

## --------------------------------------------------------------------------------------
plotUMAP(gobject = giottoObject,
         spat_unit = "cell",
         feat_type = "rna",
         cell_color = 'integrated_leiden_clus',
         dim_reduction_name = "integrated.umap",
         point_size = 2,
         title = "",
         axis_text = 14,
         axis_title = 18,
         legend_text = 14)


## --------------------------------------------------------------------------------------
spatPlot2D(giottoObject,
           spat_unit = "cell",
           feat_type = "rna",
           cell_color = "integrated_leiden_clus",
           point_size = 3.5,
           show_image = FALSE,
           title = "",
           axis_text = 14,
           axis_title = 18,
           legend_text = 14)

## Deconvolution
## --------------------------------------------------------------------------------------
## Load data
gene_count_cleaned_sampled_100k <- readRDS("gene_count_cleaned_sampled_100k.RDS")

cell_annotation = read.csv("data/scrnaseq/cell_annotate.csv")
cell_annotation = cell_annotation[,c("sample", "Total_mRNAs", "num_genes_expressed", "Main_cell_type")]
colnames(cell_annotation)[1] = "cell_ID"

cell_annotation = cell_annotation[cell_annotation$cell_ID %in% colnames(gene_count_cleaned_sampled_100k),]


## --------------------------------------------------------------------------------------
## create scRNAseq Giotto object
sc_giotto = createGiottoObject(expression = gene_count_cleaned_sampled_100k)
sc_giotto = addCellMetadata(sc_giotto,
                            new_metadata = cell_annotation)


## --------------------------------------------------------------------------------------
## Normalization
sc_giotto = normalizeGiotto(sc_giotto,
                            log_norm = FALSE,
                            scale_feats = FALSE,
                            scale_cells = FALSE)

## --------------------------------------------------------------------------------------
## Find markergenes
markers_scran <- findMarkers_one_vs_all(gobject = sc_giotto,
                                        method = "scran",
                                        expression_values = "normalized",
                                        cluster_column = 'Main_cell_type',
                                        min_feats = 3)
markergenes_scran <- unique(markers_scran[, head(.SD, 30), by = "cluster"][["feats"]])


## --------------------------------------------------------------------------------------
## Create signature matrix
DWLS_matrix_direct <- makeSignMatrixDWLSfromMatrix(
  matrix = getExpression(sc_giotto,
                         values = "normalized",
                         output = "matrix"),
  cell_type = pDataDT(sc_giotto)$Main_cell_type,
  sign_gene = markergenes_scran)


## --------------------------------------------------------------------------------------
## Fix gene names
sc_gene_names = read.csv("data/scrnaseq/GSE119945_gene_annotate.csv")

ENSMUS_names = rownames(DWLS_matrix_direct)
sc_gene_names = sc_gene_names[sc_gene_names$gene_id %in% ENSMUS_names,]

rownames(DWLS_matrix_direct) = sc_gene_names$gene_short_name


## --------------------------------------------------------------------------------------
## Run DWLS using integrated leiden clusters
giottoObject <- runDWLSDeconv(gobject = giottoObject,
                              sign_matrix = DWLS_matrix_direct,
                              cluster_column = "integrated_leiden_clus")

## --------------------------------------------------------------------------------------
# Plot DWLS deconvolution result with Pie plots
spatDeconvPlot(giottoObject,
               show_image = FALSE,
               radius = 0.5,
               return_plot = TRUE,
               save_plot = TRUE,
               save_param = list(save_name = "integrated_deconvolution"),
               title = "",
               axis_text = 14,
               axis_title = 18,
               legend_text = 8)

##%######################################################%##
#                                                          #
####     Supplementary figure 3 Nanostring dataset      ####
#                                                          #
##%######################################################%##


############################## LOAD LIBRARIES ##################################

library(terra)
library(patchwork)
library(ggplot2)
library(Giotto)
# Custom color palettes from rcartocolor

# pal10 = rcartocolor::carto_pal(n = 10, name = 'Pastel')
pal10 = c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#87C55F",
          "#9EB9F3","#FE88B1","#C9DB74","#8BE0A4","#B3B3B3")

# viv10 = rcartocolor::carto_pal(n = 10, name = 'Vivid')
viv10 = c("#E58606","#5D69B1","#52BCA3","#99C945","#CC61B0",
          "#24796C","#DAA51B","#2F8AC4","#764E9F","#A5AA99")

############################ Set seed, Initialize ###############################


# set working directory to results folder for plot saving
results_folder = '/projectnb/rd-spat/HOME/mobrien2/nanostring/full_dataset/manuscript_PDFs/'
setwd(results_folder)
my_seed_num = 315
set.seed(my_seed_num)
my_python_path = NULL # alternatively, "/local/python/path/python" if desired.

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = FALSE,
                                  show_plot = TRUE,
                                  return_plot = TRUE,
                                  python_path = my_python_path)

## provide path to nanostring folder
data_path = '/projectnb/rd-spat/DATA/Public_data/Nanostring/Lung12/Lung12-Flat_files_and_images/'

########################### Create Giotto Object ###############################

## create giotto cosmx object
fov_join = createGiottoCosMxObject(cosmx_dir = data_path,
                                   data_to_use = 'subcellular', # only subcellular
                                   FOVs = c(1,2,3,4,5,6,7,8,9,10,
                                            11,12,13,14,15,16,17,18,19,20,
                                            21,22,23,24,25,26,27,28),
                                   instructions = instrs)

showGiottoFeatInfo(fov_join)
showGiottoSpatialInfo(fov_join)


id_set = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
           "14", "15", "19", "20", "21", "22", "23", "24","25", "26", "27", "28")

new_names = paste0("fov0", id_set)

## Set up vector of image names
image_names = paste0(new_names, '-image')

######################### Create Expression Matrix #############################

# Find the feature points overlapped by polygons. This overlap information is then
# returned to the relevant giottoPolygon object's overlaps slot.
fov_join = calculateOverlapRaster(fov_join, feat_info = 'rna')
fov_join = calculateOverlapRaster(fov_join, feat_info = 'neg_probe')

# Convert the overlap information into a cell by feature expression matrix which
# is then stored in the Giotto object's expression slot
fov_join = overlapToMatrix(fov_join, feat_info = 'rna')
fov_join = overlapToMatrix(fov_join, feat_info = 'neg_probe')

showGiottoExpression(fov_join)

############################# Combine Metadata #################################

# combine cell data
morphometa = combineCellData(fov_join,
                             feat_type = 'rna')

# combine feature data
featmeta = combineFeatureData(fov_join,
                              feat_type = c('rna'))

# combine overlapping feature data
featoverlapmeta = combineFeatureOverlapData(fov_join,
                                            feat_type = c('rna'))

################### Filtering, Normalization, and Stats ########################

# filter
fov_join <- filterGiotto(gobject = fov_join,
                         feat_type = 'rna',
                         expression_threshold = 1,
                         feat_det_in_min_cells = 5,
                         min_det_feats_per_cell = 5)

# normalize
# standard method of normalization (log normalization based)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'rna',
                            norm_methods = 'standard',
                            verbose = TRUE)
fov_join <- normalizeGiotto(gobject = fov_join,
                            feat_type = 'neg_probe',
                            norm_methods = 'standard',
                            library_size_norm = FALSE,
                            verbose = TRUE)



showGiottoExpression(fov_join)

# add statistics based on log normalized values for features rna and negative probes
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'rna')
fov_join = addStatistics(gobject = fov_join,
                         expression_values = 'normalized',
                         feat_type = 'neg_probe')

# View cellular data (default is feat = 'rna')
showGiottoCellMetadata(fov_join)
# View feature data
showGiottoFeatMetadata(fov_join)


############### Dimension Reductions, NN Networks, Clustering ##################

# PCA
fov_join = runPCA(fov_join,
                  scale_unit = FALSE,
                  center = FALSE,
                  expression_values = 'normalized')

# Generate UMAP from PCA
fov_join <- runUMAP(fov_join,
                    dimensions_to_use = 1:10,
                    n_threads = determine_cores(),
                    seed_number = my_seed_num,
                    set_seed = TRUE)

# Shared Nearest Neighbor Network
fov_join <- createNearestNetwork(gobject = fov_join,
                                 dimensions_to_use = 1:10,
                                 k = 10)
# Leiden Clustering
fov_join <- doLeidenCluster(gobject = fov_join,
                            resolution = 0.15,
                            n_iterations = 100,
                            seed_number = my_seed_num,
                            set_seed = TRUE)


### Optional: Visualize Dimension Reductions:
# # visualize UMAP cluster results
# plotUMAP(gobject = fov_join,
#          cell_color = 'leiden_clus',
#          cell_color_code = pal10,
#          show_NN_network = TRUE,
#          point_size = 2)
#
# # visualize UMAP and spatial results
# spatDimPlot2D(gobject = fov_join,
#               show_image = TRUE,
#               image_name = image_names,
#               cell_color = 'leiden_clus',
#               cell_color_code = pal10,
#               spat_point_size = 1)


################## Spatially Variable Gene Identification ######################

# create spatial network based on physical distance of cell centroids
fov_join = createSpatialNetwork(gobject = fov_join,
                                minimum_k = 2,
                                maximum_distance_delaunay = 50)


# Perform Binary Spatial Extraction of genes
# NOTE: Depending on your system this could take time
rank_spatialgenes = binSpect(fov_join,
                             bin_method = "rank",
                             spatial_network_name = "Delaunay_network",
                             calc_hub = T,
                             hub_min_int = 5 )

# Extract the top 500 spatially variable genes
top_spat_var_genes = rank_spatialgenes$feats[1:500]

# Find spatially correlated features
scn_DT = detectSpatialCorFeats(fov_join,
                               method = 'network',
                               spatial_network_name = 'Delaunay_network',
                               subset_feats = top_spat_var_genes)

scn_DT = clusterSpatialCorFeats(scn_DT,
                                name = 'spat_netw_clus',
                                k = 6)

# extract an equal amount of features from each cluster
cluster_genes = getBalancedSpatCoexpressionFeats(scn_DT,
                                                 maximum = 30)

# Create a metagene representation of the extracted features from above
fov_join = createMetafeats(fov_join,
                           feat_clusters = cluster_genes,
                           name = 'cluster_metagene')

svg_names = names(cluster_genes)


####################### Spatially Informed Clustering ##########################

# Provide spatially variable gene names to the feats_to_use argument
# and give this PCA a new name
fov_join = runPCA(gobject = fov_join,
                  feats_to_use = svg_names,
                  name = "custom_pca",
                  set_seed = T,
                  seed_number = my_seed_num)


fov_join = runUMAP(gobject = fov_join,
                   dim_reduction_name = "custom_pca",
                   dimensions_to_use = 1:10,
                   n_neighbors = 5,
                   name = "custom_UMAP",
                   set_seed = T,
                   seed_number = my_seed_num)

fov_join = createNearestNetwork(gobject = fov_join,
                                dim_reduction_name = "custom_pca",
                                dimensions_to_use = 1:10,
                                name = "custom_NN")

fov_join = doLeidenClusterIgraph(gobject = fov_join,
                                 network_name = "custom_NN",
                                 resolution_parameter = 0.3,
                                 n_iterations = 100,
                                 name = "custom_leiden",
                                 set_seed = T,
                                 seed_number = my_seed_num)


################## Visualize Spatially Informed Clustering #####################

pdf("./spatPlot_custom_clustering.pdf", width = 20, height = 10)
sp_pl = spatPlot2D(fov_join,
                   cell_color = "custom_leiden",
                   point_size = 4,
                   show_plot = TRUE,
                   save_plot = FALSE,
                   return_plot = TRUE)
dev.off()

leiden_colors = getDistinctColors(length(unique(getCellMetadata(fov_join)[]$custom_leiden)))
names(leiden_colors) = sort(unique(getCellMetadata(fov_join)[]$custom_leiden))

sisp_pl = spatInSituPlotPoints(fov_join,
                               show_polygon = TRUE,
                               polygon_color = 'white',
                               polygon_line_size = 0.05,
                               polygon_fill = 'custom_leiden',
                               polygon_fill_code = leiden_colors,
                               polygon_fill_as_factor = TRUE,
                               show_plot = TRUE,
                               return_plot = TRUE,
                               save_plot = FALSE)

ggplot2::ggsave(filename = "s3_nanostring_spatInSituPlot_custom_clustering.png",
                plot = sisp_pl,
                device = "png",
                dpi = "retina",
                width = 11,
                height = 11,
                units = "in")


# Subset and visualize polygons with transcripts

ROI = subsetGiottoLocs(gobject = fov_join,
                       spat_unit = "cell",
                       feat_type = "rna",
                       spat_loc_name = "raw",
                       x_min = 9000,
                       x_max = 12000,
                       y_min = -138000,
                       y_max = -120000,
                       return_gobject = TRUE)

roi_sisp_pl = spatInSituPlotPoints(ROI,
                                   spat_unit = "cell",
                                   feat_type = "rna",
                                   spat_loc_name = "raw",
                                   feats = list("RPL21",
                                                "IGHG1",
                                                "MALAT1",
                                                "KRT19",
                                                "COL1A1",
                                                "CD74"),
                                   point_size = 0.25,
                                   show_polygon = TRUE,
                                   use_overlap = FALSE,
                                   polygon_color = 'white',
                                   polygon_line_size = 0.05,
                                   polygon_alpha = 0.33,
                                   polygon_fill = 'custom_leiden',
                                   polygon_fill_code = leiden_colors,
                                   polygon_fill_as_factor = TRUE,
                                   show_legend = FALSE,
                                   show_plot = TRUE,
                                   return_plot = TRUE,
                                   save_plot = FALSE)


ggplot2::ggsave(filename = "s3_nanostring_ROI_spatInSituPlot_custom_clustering.png",
                plot = roi_sisp_pl,
                device = "png",
                dpi = "retina",
                width = 6,
                height = 6,
                units = "in")



