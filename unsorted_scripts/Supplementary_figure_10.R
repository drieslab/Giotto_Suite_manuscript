## -----------------------------------------------------------------------------
library(Giotto)


## -----------------------------------------------------------------------------
# Set instructions to automatically save figures
save_dir = 'results_paper'
instructions = createGiottoInstructions(save_dir = save_dir,
                                        save_plot = TRUE,
                                        show_plot = TRUE)

# Provide path to visium folder
data_path = 'data'

# Create Giotto object
visium = createGiottoVisiumObject(visium_dir = data_path,
                                  expr_data = 'raw',
                                  png_name = 'tissue_lowres_image.png',
                                  gene_column_index = 2,
                                  instructions = instructions)

# show aligned image
spatPlot(gobject = visium, 
         cell_color = 'in_tissue', 
         show_image = TRUE,
         point_size = 2.5,
         point_alpha = 0.7)


## -----------------------------------------------------------------------------
## Subset on spots that were covered by tissue
metadata = pDataDT(visium)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium = subsetGiotto(visium, cell_ids = in_tissue_barcodes)


## -----------------------------------------------------------------------------
# RNA feature
visium <- filterGiotto(gobject = visium,
                       expression_threshold = 1,
                       feat_det_in_min_cells = 50,
                       min_det_feats_per_cell = 1000,
                       expression_values = c('raw'),
                       verbose = TRUE)

visium <- normalizeGiotto(gobject = visium, 
                          scalefactor = 6000, 
                          verbose = TRUE)

visium <- addStatistics(gobject = visium)

# Protein feature
visium <- filterGiotto(gobject = visium,
                       spat_unit = 'cell',
                       feat_type = 'protein',
                       expression_threshold = 1,
                       feat_det_in_min_cells = 50,
                       min_det_feats_per_cell = 1,
                       expression_values = c('raw'),
                       verbose = TRUE)

visium <- normalizeGiotto(gobject = visium, 
                          spat_unit = 'cell',
                          feat_type = 'protein',
                          scalefactor = 6000, 
                          verbose = TRUE)

visium <- addStatistics(gobject = visium,
                        spat_unit = 'cell',
                        feat_type = 'protein')


## -----------------------------------------------------------------------------
# RNA
spatPlot2D(gobject = visium, 
           point_alpha = 0.7,
           cell_color = 'nr_feats', 
           color_as_factor = FALSE)

# Protein
spatPlot2D(gobject = visium, 
           spat_unit = 'cell',
           feat_type = 'protein',
           point_alpha = 0.7,
           cell_color = 'nr_feats', 
           color_as_factor = FALSE)


## -----------------------------------------------------------------------------
visium <- calculateHVF(gobject = visium)


## -----------------------------------------------------------------------------
# RNA
visium <- runPCA(gobject = visium)

screePlot(visium, ncp = 30)

plotPCA(gobject = visium)


## -----------------------------------------------------------------------------
# Protein
visium <- runPCA(gobject = visium,
                 spat_unit = 'cell',
                 feat_type = 'protein')

screePlot(visium, 
          spat_unit = 'cell',
          feat_type = 'protein',
          ncp = 30)

plotPCA(gobject = visium,
        spat_unit = 'cell',
        feat_type = 'protein')


## -----------------------------------------------------------------------------
# RNA
visium <- runUMAP(visium, 
                  dimensions_to_use = 1:10)

plotUMAP(gobject = visium)


## -----------------------------------------------------------------------------
# Protein
visium <- runUMAP(visium,
                  spat_unit = 'cell',
                  feat_type = 'protein',
                  dimensions_to_use = 1:10)

plotUMAP(gobject = visium,
         spat_unit = 'cell',
         feat_type = 'protein')

## -----------------------------------------------------------------------------
# RNA
visium <- runtSNE(visium, 
                  dimensions_to_use = 1:10)

plotTSNE(gobject = visium)


## -----------------------------------------------------------------------------
# Protein
visium <- runtSNE(visium, 
                  spat_unit = 'cell',
                  feat_type = 'protein',
                  dimensions_to_use = 1:10)

plotTSNE(gobject = visium,
         spat_unit = 'cell',
         feat_type = 'protein')

## -----------------------------------------------------------------------------
# RNA
visium <- createNearestNetwork(gobject = visium, 
                               dimensions_to_use = 1:10, 
                               k = 30)

visium <- doLeidenCluster(gobject = visium, 
                          resolution = 1,
                          n_iterations = 1000)

# Protein
visium <- createNearestNetwork(gobject = visium, 
                               spat_unit = 'cell', 
                               feat_type = 'protein', 
                               dimensions_to_use = 1:10, 
                               k = 30)

visium <- doLeidenCluster(gobject = visium, 
                          spat_unit = 'cell', 
                          feat_type = 'protein', 
                          resolution = 1,
                          n_iterations = 1000)


## -----------------------------------------------------------------------------
# RNA
plotUMAP(gobject = visium, 
         cell_color = 'leiden_clus', 
         show_NN_network = FALSE, 
         point_size = 2,
         title = "",
         axis_text = 14,
         axis_title = 18,
         legend_text = 14)

# Protein
plotUMAP(gobject = visium, 
         spat_unit = 'cell', 
         feat_type = 'protein',
         cell_color = 'leiden_clus', 
         show_NN_network = FALSE, 
         point_size = 2,
         title = "",
         axis_text = 14,
         axis_title = 18,
         legend_text = 14)


## -----------------------------------------------------------------------------
# RNA
plotTSNE(gobject = visium, 
         cell_color = 'leiden_clus', 
         show_NN_network = FALSE, 
         point_size = 2)

# Protein
plotTSNE(gobject = visium, 
         spat_unit = 'cell', 
         feat_type = 'protein',
         cell_color = 'leiden_clus', 
         show_NN_network = FALSE, 
         point_size = 2)


## -----------------------------------------------------------------------------
# RNA
spatPlot2D(gobject = visium, 
           show_image = TRUE,
           cell_color = 'leiden_clus',
           point_size = 2.5,
           title = "",
           axis_text = 14,
           axis_title = 18,
           legend_text = 14)

# Protein
spatPlot2D(gobject = visium, 
           spat_unit = 'cell', 
           feat_type = 'protein',
           show_image = TRUE,
           cell_color = 'leiden_clus',
           point_size = 2.5,
           title = "",
           axis_text = 14,
           axis_title = 18,
           legend_text = 14)


## -----------------------------------------------------------------------------
# RNA
visium <- createNearestNetwork(gobject = visium, 
                               type = 'kNN',
                               dimensions_to_use = 1:10, 
                               k = 20)

# Protein
visium <- createNearestNetwork(gobject = visium, 
                               spat_unit = 'cell', 
                               feat_type = 'protein', 
                               type = 'kNN',
                               dimensions_to_use = 1:10, 
                               k = 20)


## -----------------------------------------------------------------------------
visium <- runWNN(visium,
                 spat_unit = "cell",
                 modality_1 = "rna",
                 modality_2 = "protein",
                 pca_name_modality_1 = "pca",
                 pca_name_modality_2 = "protein.pca",
                 k = 20,
                 verbose = TRUE)


## -----------------------------------------------------------------------------
visium <- runIntegratedUMAP(visium,
                            modality1 = "rna",
                            modality2 = "protein",
                            spread = 7,
                            min_dist = 1,
                            force = FALSE)


## -----------------------------------------------------------------------------
visium <- doLeidenCluster(gobject = visium,
                          spat_unit = "cell",
                          feat_type = "rna",
                          nn_network_to_use = "kNN",
                          network_name = "integrated_kNN",
                          name = "integrated_leiden_clus",
                          resolution = 1)


## -----------------------------------------------------------------------------
plotUMAP(gobject = visium,
         spat_unit = "cell",
         feat_type = "rna",
         cell_color = 'integrated_leiden_clus',
         dim_reduction_name = "integrated.umap",
         point_size = 2,
         title = "",
         axis_text = 14,
         axis_title = 18,
         legend_text = 14)


## -----------------------------------------------------------------------------
spatPlot2D(visium,
           spat_unit = "cell",
           feat_type = "rna",
           cell_color = "integrated_leiden_clus",
           point_size = 2.5,
           show_image = TRUE,
           title = "",
           axis_text = 14,
           axis_title = 18,
           legend_text = 14)


## -----------------------------------------------------------------------------
# Load reference scRNAseq
devtools::install_github("massonix/HCATonsilData", 
                         build_vignettes = TRUE)

sce <- HCATonsilData::HCATonsilData(assayType = "RNA", cellType = "All")

# 209786 cells profiled with scRNA-seq (3P) and 53513 cells profiled with multiome
table(sce$assay) 
df <- as.data.frame(SingleCellExperiment::colData(sce))

# Create Giotto object with raw counts
giotto_SC = createGiottoObject(expression = assay(sce, "counts"))

# Add normalized expression matrix
expression_object = createExprObj(expression_data = assay(sce, "logcounts"),
                                  name = "logcounts")

giotto_SC = setExpression(giotto_SC,
                          x = expression_object,
                          name = "logcounts")

# Add cell metadata
giotto_SC = addCellMetadata(giotto_SC,
                            new_metadata = as.data.frame(colData(sce)))

# Add umap
umap_dimObj = createDimObj(coordinates = SingleCellExperiment::reducedDim(sce, "UMAP"),
                           name = "umap")

giotto_SC = set_dimReduction(giotto_SC,
                             dimObject = umap_dimObj,
                             reduction_method = "umap",
                             name = "umap")

# Select only cells from scRNA-seq (3P) experiment
cells_3P = rownames(df[df$assay == "3P",])

giotto_SC = subsetGiotto(giotto_SC,
                         cell_ids = cells_3P)


## -----------------------------------------------------------------------------
markers_scran <- findMarkers_one_vs_all(gobject = giotto_SC,
                                        method = "scran",
                                        expression_values = "logcounts",
                                        cluster_column = 'annotation_level_1',
                                        min_feats = 3)
markergenes_scran <- unique(markers_scran[, head(.SD, 30), by = "cluster"][["feats"]])


## -----------------------------------------------------------------------------
# Create signature matrix
DWLS_matrix_direct <- makeSignMatrixDWLSfromMatrix(
    matrix = 2**(getExpression(giotto_SC,
                               values = "logcounts",
                               output = "matrix"))-1,
    cell_type = pDataDT(giotto_SC)$annotation_level_1,
    sign_gene = markergenes_scran)


## -----------------------------------------------------------------------------
# Run DWLS using integrated leiden clusters
visium <- runDWLSDeconv(gobject = visium,
                        sign_matrix = DWLS_matrix_direct,
                        cluster_column = "integrated_leiden_clus")


## -----------------------------------------------------------------------------
# Plot DWLS deconvolution result with Pie plots dataset 1
spatDeconvPlot(visium,
               show_image = TRUE,
               radius = 150,
               return_plot = TRUE,
               save_plot = TRUE,
               save_param = list(save_name = "integrated_deconvolution_annotation_level_1"),
               title = "",
               axis_text = 14,
               axis_title = 18,
               legend_text = 14)

