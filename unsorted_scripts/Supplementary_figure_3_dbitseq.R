# Supplementary figure 3 DBiT-seq dataset

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

