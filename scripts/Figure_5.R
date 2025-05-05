##%######################################################%##
#                                                          #
####          Figure 5 Multi-omics integration          ####
#                                                          #
##%######################################################%##

############################## Download data  ##################################
## Download the 10X Genomics multi-modal Visium CytAssist Human Tonsil dataset.
## https://www.10xgenomics.com/resources/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard

dir.create("data/Visium_Tonsil")

download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_raw_feature_bc_matrix.tar.gz",
              destfile = "data/Visium_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_raw_feature_bc_matrix.tar.gz")
download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FFPE_Protein_Expression_Human_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_spatial.tar.gz",
              destfile = "data/Visium_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_spatial.tar.gz")

untar(
    tarfile = "data/Visium_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_raw_feature_bc_matrix.tar.gz",
    exdir = "data/Visium_Tonsil/"
)
untar(
    tarfile = "data/Visium_Tonsil/CytAssist_FFPE_Protein_Expression_Human_Tonsil_spatial.tar.gz",
    exdir = "data/Visium_Tonsil/"
)

## Download the Mouse embryo ME13 spatial rna-atac seq was originally published by Zhang et al.
## Spatial epigenome–transcriptome co-profiling of mammalian tissues, 2023
## https://www.nature.com/articles/s41586-023-05795-1
## The data is available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205055
## For running this tutorial, download the following files:
## - GSM6799937_ME13_50um_matrix_merge.tsv.gz
## - GSM6801813_ME13_50um_fragments.tsv.gz
## - GSM6799937_ME13_50um_spatial.tar.gz

dir.create("data/ME13")

library(Giotto)

# Set instructions to automatically save figures
save_dir <- "results/figure_5/"

instructions <- createGiottoInstructions(
    save_dir = save_dir,
    save_plot = TRUE,
    show_plot = TRUE
)

# Provide path to visium folder
data_path <- "data/Visium_Tonsil"

##%######################################################%##
#                                                          #
####         Figure 5b RNA-Protein integration          ####
#                                                          #
##%######################################################%##

############################## Create the object  ##############################

# Create Giotto object
visium <- createGiottoVisiumObject(
    visium_dir = data_path,
    expr_data = "raw",
    png_name = "tissue_lowres_image.png",
    gene_column_index = 2,
    instructions = instructions
)

# show aligned image
spatPlot(
    gobject = visium,
    cell_color = "in_tissue",
    show_image = TRUE,
    point_size = 2.5,
    point_alpha = 0.7
)

## Subset on spots that were covered by tissue
metadata <- pDataDT(visium)
in_tissue_barcodes <- metadata[in_tissue == 1]$cell_ID
visium <- subsetGiotto(visium, cell_ids = in_tissue_barcodes)

#################################### Filtering  ################################

## RNA
visium <- filterGiotto(
    gobject = visium,
    expression_threshold = 1,
    feat_det_in_min_cells = 50,
    min_det_feats_per_cell = 1000,
    expression_values = "raw",
    verbose = TRUE
)

## Protein
visium <- filterGiotto(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    expression_threshold = 1,
    feat_det_in_min_cells = 50,
    min_det_feats_per_cell = 1,
    expression_values = "raw",
    verbose = TRUE
)

################################ Normalization  ################################

## RNA
visium <- normalizeGiotto(
    gobject = visium,
    scalefactor = 6000,
    verbose = TRUE
)

## Protein
visium <- normalizeGiotto(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    scalefactor = 6000,
    verbose = TRUE
)

################################# Statistics  ##################################

## RNA
visium <- addStatistics(gobject = visium)

## Protein
visium <- addStatistics(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein"
)

## RNA
spatPlot2D(
    gobject = visium,
    point_alpha = 0.7,
    cell_color = "nr_feats",
    color_as_factor = FALSE
)

## Protein
spatPlot2D(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    point_alpha = 0.7,
    cell_color = "nr_feats",
    color_as_factor = FALSE
)

############################## Dimension reduction  ############################

## RNA
visium <- calculateHVF(gobject = visium)

visium <- runPCA(gobject = visium)

screePlot(visium, ncp = 30)

plotPCA(gobject = visium)

## Protein
visium <- runPCA(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein"
)

screePlot(visium,
          spat_unit = "cell",
          feat_type = "protein",
          ncp = 30
)

plotPCA(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein"
)

################################### Clustering  ################################

## RNA
visium <- runUMAP(visium,
                  dimensions_to_use = 1:10
)

plotUMAP(gobject = visium)

## Protein
visium <- runUMAP(visium,
                  spat_unit = "cell",
                  feat_type = "protein",
                  dimensions_to_use = 1:10
)

plotUMAP(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein"
)

## RNA
visium <- runtSNE(visium,
                  dimensions_to_use = 1:10
)

plotTSNE(gobject = visium)

## Protein
visium <- runtSNE(visium,
                  spat_unit = "cell",
                  feat_type = "protein",
                  dimensions_to_use = 1:10
)

plotTSNE(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein"
)

## RNA
visium <- createNearestNetwork(
    gobject = visium,
    dimensions_to_use = 1:10,
    k = 30
)

visium <- doLeidenCluster(
    gobject = visium,
    resolution = 1,
    n_iterations = 1000
)

## Protein
visium <- createNearestNetwork(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    dimensions_to_use = 1:10,
    k = 30
)

visium <- doLeidenCluster(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    resolution = 1,
    n_iterations = 1000
)

## RNA
plotUMAP(
    gobject = visium,
    cell_color = "leiden_clus",
    show_NN_network = FALSE,
    point_size = 2,
    title = "",
    axis_text = 14,
    axis_title = 18,
    legend_text = 14
)

## Protein
plotUMAP(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    cell_color = "leiden_clus",
    show_NN_network = FALSE,
    point_size = 2,
    title = "",
    axis_text = 14,
    axis_title = 18,
    legend_text = 14
)

## RNA
plotTSNE(
    gobject = visium,
    cell_color = "leiden_clus",
    show_NN_network = FALSE,
    point_size = 2
)

## Protein
plotTSNE(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    cell_color = "leiden_clus",
    show_NN_network = FALSE,
    point_size = 2
)

## RNA
spatPlot2D(
    gobject = visium,
    show_image = TRUE,
    cell_color = "leiden_clus",
    point_size = 2.5,
    title = "",
    axis_text = 14,
    axis_title = 18,
    legend_text = 14
)

## Protein
spatPlot2D(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    show_image = TRUE,
    cell_color = "leiden_clus",
    point_size = 2.5,
    title = "",
    axis_text = 14,
    axis_title = 18,
    legend_text = 14
)

########################### Multi-omics integration  ###########################

## RNA
visium <- createNearestNetwork(
    gobject = visium,
    type = "kNN",
    dimensions_to_use = 1:10,
    k = 20
)

## Protein
visium <- createNearestNetwork(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "protein",
    type = "kNN",
    dimensions_to_use = 1:10,
    k = 20
)


visium <- runWNN(visium,
                 spat_unit = "cell",
                 feat_types = c("rna", "protein"),
                 reduction_methods = c("pca", "pca"),
                 reduction_names = c("pca", "protein.pca"),
                 k = 20,
                 verbose = TRUE
)

visium <- runIntegratedUMAP(visium,
                            feat_types = c("rna", "protein"),
                            spread = 7,
                            min_dist = 1,
                            force = FALSE
)

visium <- doLeidenCluster(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "rna",
    nn_network_to_use = "kNN",
    network_name = "integrated_kNN",
    name = "integrated_leiden_clus",
    resolution = 1
)

plotUMAP(
    gobject = visium,
    spat_unit = "cell",
    feat_type = "rna",
    cell_color = "integrated_leiden_clus",
    dim_reduction_name = "integrated.umap",
    point_size = 2,
    title = "",
    axis_text = 14,
    axis_title = 18,
    legend_text = 14
)

spatPlot2D(visium,
           spat_unit = "cell",
           feat_type = "rna",
           cell_color = "integrated_leiden_clus",
           point_size = 2.5,
           show_image = TRUE,
           title = "",
           axis_text = 14,
           axis_title = 18,
           legend_text = 14
)

############################### Deconvolution  #################################

# Load reference scRNAseq
devtools::install_github("massonix/HCATonsilData",
                         build_vignettes = TRUE
)

sce <- HCATonsilData::HCATonsilData(assayType = "RNA", cellType = "All")

# 209786 cells profiled with scRNA-seq (3P) and 53513 cells profiled with multiome
table(sce$assay)
df <- as.data.frame(SingleCellExperiment::colData(sce))

# Create Giotto object with raw counts
giotto_SC <- createGiottoObject(expression = assay(sce, "counts"))

# Add normalized expression matrix
expression_object <- createExprObj(
    expression_data = assay(sce, "logcounts"),
    name = "logcounts"
)

giotto_SC <- setExpression(giotto_SC,
                           x = expression_object,
                           name = "logcounts"
)

# Add cell metadata
giotto_SC <- addCellMetadata(giotto_SC,
                             new_metadata = as.data.frame(colData(sce))
)

# Add umap
umap_dimObj <- createDimObj(
    coordinates = SingleCellExperiment::reducedDim(sce, "UMAP"),
    name = "umap"
)

giotto_SC <- set_dimReduction(giotto_SC,
                              dimObject = umap_dimObj,
                              reduction_method = "umap",
                              name = "umap"
)

# Select only cells from scRNA-seq (3P) experiment
cells_3P <- rownames(df[df$assay == "3P", ])

giotto_SC <- subsetGiotto(giotto_SC,
                          cell_ids = cells_3P
)

markers_scran <- findMarkers_one_vs_all(
    gobject = giotto_SC,
    method = "scran",
    expression_values = "logcounts",
    cluster_column = "annotation_level_1",
    min_feats = 3
)
markergenes_scran <- unique(markers_scran[, head(.SD, 30), by = "cluster"][["feats"]])

# Create signature matrix
DWLS_matrix_direct <- makeSignMatrixDWLSfromMatrix(
    matrix = 2**(getExpression(giotto_SC,
                               values = "logcounts",
                               output = "matrix"
    )) - 1,
    cell_type = pDataDT(giotto_SC)$annotation_level_1,
    sign_gene = markergenes_scran
)

# Run DWLS using integrated leiden clusters
visium <- runDWLSDeconv(
    gobject = visium,
    sign_matrix = DWLS_matrix_direct,
    cluster_column = "integrated_leiden_clus"
)

# Plot DWLS deconvolution result with Pie plots dataset 1
spatDeconvPlot(visium,
               show_image = TRUE,
               radius = 150,
               return_plot = TRUE,
               save_plot = TRUE,
               title = "",
               axis_text = 14,
               axis_title = 18,
               legend_text = 14
)

##%######################################################%##
#                                                          #
####           Figure 5c RNA-ATAC integration           ####
#                                                          #
##%######################################################%##

library(Giotto)

##############################  Read RNA matrix  ###############################

rna_matrix <- read.delim("data/ME13/GSM6799937_ME13_50um_matrix_merge.tsv.gz")

########################  Read spatial locations file  #########################

spatial_locs <- read.csv("data/ME13/ME13_50um_spatial/tissue_positions_list.csv",
                         header = FALSE)
spatial_locs <- spatial_locs[, c("V1", "V6", "V5")]
colnames(spatial_locs) <- c("cell_ID", "sdimx", "sdimy")
rownames(spatial_locs) <- spatial_locs$cell_ID
spatial_locs <- spatial_locs[colnames(rna_matrix),]

###############################  Read ATAC data  ###############################

library(ArchR)

## Add reference genome
ArchR::addArchRGenome('mm10')

## Creating Arrow Files
arrow_file <- ArchR::createArrowFiles(
    inputFiles = "data/ME13/GSM6801813_ME13_50um_fragments.tsv.gz",
    sampleNames = "atac",
    minTSS = 0, # Minimum TSS enrichment score
    minFrags = 0, # Minimum number of fragments
    maxFrags = 1e+07, # Maximum number of fragments
    TileMatParams = list(tileSize = 5000), # Tile size for creating accessibility matrix
    offsetPlus = 0,
    offsetMinus = 0,
    force = TRUE # Overwrite existing files
)

## Create ArchR project from Arrow files
proj <- ArchR::ArchRProject(
    arrow_file,
    outputDirectory = "ArchRProject"
)

## Get Tile matrix from Arrow file
tile_matrix <- ArchR::getMatrixFromArrow(
    ArrowFile = "ArchRProject/ArrowFiles/atac.arrow",
    useMatrix = "TileMatrix",
    binarize = TRUE
)

tile_matrix_dg <- tile_matrix@assays@data$TileMatrix # Dimensions: 526765 x 2172

rm(tile_matrix)

cell_ids <- colnames(tile_matrix_dg)
cell_ids <- gsub(pattern = "atac#", replacement = "", x = cell_ids)
cell_ids <- gsub(pattern = "-1", replacement = "", x = cell_ids)
colnames(tile_matrix_dg) <- cell_ids

## keep only cells in common
atac_cells <- colnames(tile_matrix_dg)
rna_cells <- colnames(rna_matrix)

all_cells <- c(rna_cells, atac_cells)
commom_cells <- all_cells[duplicated(all_cells)]

spatial_locs <- spatial_locs[commom_cells, ]
rna_matrix <- rna_matrix[, commom_cells]
atac_matrix <- tile_matrix_dg[, commom_cells]

###########################  Create Giotto object  #############################

instructions <- createGiottoInstructions(save_plot = TRUE,
                                         save_dir = "results/figure_5",
                                         python_path = NULL)

g <- createGiottoObject(expression = list("raw" = rna_matrix,
                                          "raw" = atac_matrix),
                        expression_feat = c("rna", "atac"),
                        spatial_locs = spatial_locs,
                        instructions = instructions)

g <- flip(g,
          direction = "vertical")

# add image
g_image <- createGiottoImage(gobject = g,
                             mg_object = "data/ME13/ME13_50um_spatial/tissue_lowres_image.png")
g <- addGiottoImage(g,
                    images = list(g_image))

#################################### Filtering  ################################

## RNA
g <- filterGiotto(g,
                  min_det_feats_per_cell = 50,
                  feat_det_in_min_cells = 50)

################################ Normalization  ################################

## RNA
g <- normalizeGiotto(g,
                     feat_type = "rna")

############################## Dimension reduction  ############################

## RNA (PCA)
g <- runPCA(g,
            feat_type = "rna")

screePlot(g,
          feat_type = "rna")

## ATAC (iterative LSI)
g <- runIterativeLSI(
    gobject = g,
    feat_type = "atac",
    expression_values = "raw",
    lsi_method = 2,
    resolution = 0.2,
    sample_cells_pre = 20000,
    var_features = 30000,
    dims = 30
)

################################### Clustering  ################################

## RNA
g <- runUMAP(g,
             feat_type = "rna")

## ATAC
g <- runUMAP(g,
             feat_type = "atac",
             dim_reduction_to_use = "lsi",
             dim_reduction_name = "atac_iterative_lsi")

## RNA
g <- createNearestNetwork(g,
                          feat_type = "rna")

## ATAC
g <- createNearestNetwork(g,
                          feat_type = "atac",
                          dim_reduction_to_use = "lsi",
                          dim_reduction_name = "atac_iterative_lsi")

## RNA
g <- doLeidenCluster(g,
                     feat_type = "rna",
                     resolution = 0.8)

## ATAC
g <- doLeidenCluster(g,
                     feat_type = "atac",
                     network_name = "sNN.lsi",
                     resolution = 0.6,
                     name = "atac_leiden_clus")

# Plots

## RNA
plotUMAP(g,
         feat_type = "rna",
         cell_color = "leiden_clus",
         point_size = 3,
         cell_color_code = c("1" = "#FFDE7A",
                             "2" = "#FEFF05",
                             "3" = "#8BA0D1",
                             "4" = "#FF9C7A",
                             "5" = "#BFFF04",
                             "6" = "#FF8002",
                             "7" = "#01FFFF",
                             "8" = "#FF0000",
                             "9" = "#00BDBE",
                             "10" = "#C06BAA",
                             "11" = "#4000FF"
         ),
         title = "RNA UMAP")

spatPlot2D(g,
           feat_type = "rna",
           cell_color = "leiden_clus",
           show_image = TRUE,
           point_size = 5,
           point_alpha = 0.5,
           cell_color_code = c("1" = "#FFDE7A",
                               "2" = "#FEFF05",
                               "3" = "#8BA0D1",
                               "4" = "#FF9C7A",
                               "5" = "#BFFF04",
                               "6" = "#FF8002",
                               "7" = "#01FFFF",
                               "8" = "#FF0000",
                               "9" = "#00BDBE",
                               "10" = "#C06BAA",
                               "11" = "#4000FF"
           ),
           title = "RNA Leiden clusters")

spatFeatPlot2D(g,
               expression_values = "scaled",
               feats = "Pax6",
               point_size = 1.5)

## ATAC
plotUMAP(g,
         feat_type = "atac",
         dim_reduction_name = "atac.umap",
         cell_color = "atac_leiden_clus",
         point_size = 3,
         cell_color_code = c("1" = "#8BA0D1",
                             "2" = "#FEFF05",
                             "3" = "#FF9C7A",
                             "4" = "#4a70f7",
                             "5" = "#FF0000",
                             "6" = "#FF8002",
                             "7" = "#C06BAA"
         ),
         title = "ATAC UMAP")

spatPlot2D(g,
           feat_type = "atac",
           cell_color = "atac_leiden_clus",
           point_size = 5,
           point_alpha = 0.5,
           show_image = TRUE,
           cell_color_code = c("1" = "#8BA0D1",
                               "2" = "#FEFF05",
                               "3" = "#FF9C7A",
                               "4" = "#4a70f7",
                               "5" = "#FF0000",
                               "6" = "#FF8002",
                               "7" = "#C06BAA"),
           title = "ATAC Leiden clusters")

########################### Multi-omics integration  ###########################

g <- createNearestNetwork(g,
                          feat_type = "rna",
                          type = "kNN",
                          dimensions_to_use = 1:10,
                          k = 10)

g <- createNearestNetwork(g,
                          feat_type = "atac",
                          dim_reduction_to_use = "lsi",
                          dim_reduction_name = "atac_iterative_lsi",
                          type = "kNN",
                          dimensions_to_use = 1:10,
                          k = 10)

g <- runWNN(
    gobject = g,
    spat_unit = "cell",
    feat_types = c("rna", "atac"),
    reduction_methods = c("pca", "lsi"),
    reduction_names = c("pca", "atac_iterative_lsi"),
    k = 10,
    integrated_feat_type = NULL,
    matrix_result_name = NULL,
    w_names = c(NULL, NULL),
    verbose = FALSE
)

g <- runIntegratedUMAP(
    gobject = g,
    spat_unit = "cell",
    feat_types = c("rna", "atac"),
    integrated_feat_type = "rna_atac",
    integration_method = "WNN",
    matrix_result_name = "theta_weighted_matrix",
    k = 10,
    spread = 10,
    min_dist = 0.01,
    force = TRUE,
    seed = 1234
)

g <- doLeidenCluster(g,
                     feat_type = "rna",
                     nn_network_to_use = "kNN",
                     network_name = "integrated_kNN",
                     resolution = 1.2,
                     name = "integrated_leiden_clus")

# Plots
plotUMAP(g,
         feat_type = "rna",
         dim_reduction_name = "integrated.umap",
         cell_color = "integrated_leiden_clus",
         point_size = 3,
         cell_color_code = c("1" = "#BFFF04",
                             "2" = "#FEFF05",
                             "3" = "#FFDE7A",
                             "4" = "#FF00FF",
                             "5" = "#8BA0D1",
                             "6" = "#01FFFF",
                             "7" = "#4000FF",
                             "8" = "#FF8002",
                             "9" = "#FF0000",
                             "10" = "#f0b699",
                             "11" = "#4a70f7",
                             "12" = "#f09567",
                             "13" = "green",
                             "14" = "#1F8A42",
                             "15" = "#00BDBE",
                             "16" = "#C06BAA",
                             "17" = "#FF9C7A"
         ),
         title = "Integrated UMAP")

spatPlot2D(g,
           feat_type = "rna",
           cell_color = "integrated_leiden_clus",
           point_size = 5,
           cell_color_code = c("1" = "#BFFF04",
                               "2" = "#FEFF05",
                               "3" = "#FFDE7A",
                               "4" = "#FF00FF",
                               "5" = "#8BA0D1",
                               "6" = "#01FFFF",
                               "7" = "#4000FF",
                               "8" = "#FF8002",
                               "9" = "#FF0000",
                               "10" = "#f0b699",
                               "11" = "#4a70f7",
                               "12" = "#f09567",
                               "13" = "green",
                               "14" = "#1F8A42",
                               "15" = "#00BDBE",
                               "16" = "#C06BAA",
                               "17" = "#FF9C7A"
           )
           ,
           title = "Integrated Leiden clusters")

