## %######################################################%##
#                                                          #
####   Supplementary figure 3 Spatial Cite-seq dataset  ####
#                                                          #
## %######################################################%##

############################## Download dataset  ###############################
## Download the Skin dataset from the article
## [High-plex protein and whole transcriptome co-mapping at cellular resolution with spatial CITE-seq](https://www.nature.com/articles/s41587-023-01676-0)


############################## Create the object  ##############################

library(Giotto)

instructions <- createGiottoInstructions(
    save_plot = TRUE,
    save_dir = "results",
    show_plot = TRUE
)

x <- data.table::fread("data/GSM6578065_humanskin_RNA.tsv.gz")
spatial_coords <- data.frame(cell_ID = x$X)
spatial_coords <- cbind(
    spatial_coords,
    stringr::str_split_fixed(spatial_coords$cell_ID,
        pattern = "x",
        n = 2
    )
)

colnames(spatial_coords)[2:3] <- c("sdimx", "sdimy")
spatial_coords$sdimx <- as.integer(spatial_coords$sdimx)
spatial_coords$sdimy <- as.integer(spatial_coords$sdimy)
spatial_coords$sdimy <- spatial_coords$sdimy * (-1)

rna_matrix <- data.table::fread("data/GSM6578065_humanskin_RNA.tsv.gz")
rna_matrix <- rna_matrix[rna_matrix$X %in% spatial_coords$cell_ID, ]
rna_matrix <- rna_matrix[match(spatial_coords$cell_ID, rna_matrix$X), ]
rna_matrix <- t(rna_matrix[, -1])
colnames(rna_matrix) <- spatial_coords$cell_ID

protein_matrix <- data.table::fread("data/GSM6578074_humanskin_protein.tsv.gz")
protein_matrix <- protein_matrix[protein_matrix$X %in% spatial_coords$cell_ID, ]
protein_matrix <- protein_matrix[match(spatial_coords$cell_ID, protein_matrix$X), ]
protein_matrix <- t(protein_matrix[, -1])
colnames(protein_matrix) <- spatial_coords$cell_ID

my_giotto_image <- createGiottoLargeImage(
    raster_object = "data/skin.jpg",
    scale_factor = 0.5,
    negative_y = TRUE
)

my_giotto_object <- createGiottoObject(
    expression = list(
        rna = list(raw = rna_matrix),
        protein = list(raw = protein_matrix)
    ),
    expression_feat = list("rna", "protein"),
    spatial_locs = spatial_coords,
    largeImages = list(my_giotto_image),
    instructions = instructions
)

#################################### Filtering  ################################

## RNA
my_giotto_object <- filterGiotto(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    expression_threshold = 1,
    feat_det_in_min_cells = 1,
    min_det_feats_per_cell = 1
)


## Protein
my_giotto_object <- filterGiotto(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    expression_threshold = 1,
    feat_det_in_min_cells = 1,
    min_det_feats_per_cell = 1
)

################################ Normalization  ################################

## RNA
my_giotto_object <- normalizeGiotto(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    norm_methods = "standard",
    scalefactor = 10000,
    verbose = TRUE
)

## Protein
my_giotto_object <- normalizeGiotto(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    scalefactor = 6000,
    verbose = TRUE
)

################################# Statistics  ##################################

## RNA
my_giotto_object <- addStatistics(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna"
)


## Protein
my_giotto_object <- addStatistics(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    expression_values = "normalized"
)

############################## Dimension reduction  ############################

## RNA
my_giotto_object <- runPCA(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    expression_values = "normalized",
    reduction = "cells",
    name = "rna.pca"
)

## Protein
my_giotto_object <- runPCA(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    expression_values = "normalized",
    scale_unit = TRUE,
    center = FALSE,
    method = "factominer"
)

################################### Clustering  ################################

## RNA
my_giotto_object <- runUMAP(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    expression_values = "normalized",
    reduction = "cells",
    dimensions_to_use = 1:10,
    dim_reduction_name = "rna.pca"
)


## Protein
my_giotto_object <- runUMAP(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    expression_values = "normalized",
    dimensions_to_use = 1:10
)


## RNA
my_giotto_object <- createNearestNetwork(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    type = "sNN",
    dim_reduction_to_use = "pca",
    dim_reduction_name = "rna.pca",
    dimensions_to_use = 1:10,
    k = 20
)


## Protein
my_giotto_object <- createNearestNetwork(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    type = "sNN",
    name = "protein_sNN.pca",
    dimensions_to_use = 1:10,
    k = 20
)


## RNA
my_giotto_object <- doLeidenCluster(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    nn_network_to_use = "sNN",
    name = "leiden_clus",
    resolution = 1
)


## Protein
my_giotto_object <- doLeidenCluster(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    nn_network_to_use = "sNN",
    network_name = "protein_sNN.pca",
    name = "leiden_clus",
    resolution = 1
)


############################ Multi-omics integration  ##########################

## RNA
my_giotto_object <- createNearestNetwork(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    type = "kNN",
    dim_reduction_name = "rna.pca",
    name = "rna_kNN.pca",
    dimensions_to_use = 1:10,
    k = 20
)

## Protein
my_giotto_object <- createNearestNetwork(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "protein",
    type = "kNN",
    name = "protein_kNN.pca",
    dimensions_to_use = 1:10,
    k = 20
)

## WNN
my_giotto_object <- runWNN(my_giotto_object,
    feat_types = c("rna", "protein"),
    reduction_methods = c("pca", "pca"),
    reduction_names = c("rna.pca", "protein.pca"),
    pca_name_modality_2 = "protein.pca",
    k = 20
)

my_giotto_object <- runIntegratedUMAP(my_giotto_object,
    feat_types = c("rna", "protein")
)

my_giotto_object <- doLeidenCluster(
    gobject = my_giotto_object,
    spat_unit = "cell",
    feat_type = "rna",
    nn_network_to_use = "kNN",
    network_name = "integrated_kNN",
    name = "integrated_leiden_clus",
    resolution = 0.7
)

######################### Spatially informed clusters  #########################

my_giotto_object <- createSpatialNetwork(
    gobject = my_giotto_object,
    method = "kNN",
    k = 6,
    maximum_distance_knn = 5,
    name = "spatial_network"
)

ranktest <- binSpect(my_giotto_object,
    bin_method = "rank",
    calc_hub = TRUE,
    hub_min_int = 5,
    spatial_network_name = "spatial_network"
)

spatFeatPlot2D(my_giotto_object,
    expression_values = "scaled",
    feats = ranktest$feats[1:6],
    cow_n_col = 2,
    point_size = 1.5
)


## cluster the top 1500 spatial genes into 20 clusters
ext_spatial_genes <- ranktest[1:500, ]$feats

# here we use existing detectSpatialCorGenes function to calculate pairwise distances between genes (but set network_smoothing=0 to use default clustering)
spat_cor_netw_DT <- detectSpatialCorFeats(my_giotto_object,
    method = "network",
    spatial_network_name = "spatial_network",
    subset_feats = ext_spatial_genes
)

# identify potential spatial co-expression
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT,
    name = "spat_netw_clus",
    k = 3
)

# visualize clusters
heatmSpatialCorFeats(my_giotto_object,
    spatCorObject = spat_cor_netw_DT,
    use_clus_name = "spat_netw_clus",
    heatmap_legend_param = list(title = NULL),
    save_param = list(
        base_height = 6,
        base_width = 8,
        units = "cm"
    )
)


# create metagenes/co-expression modules
cluster_genes <- getBalancedSpatCoexpressionFeats(spat_cor_netw_DT,
    maximum = 30
)

my_giotto_object <- createMetafeats(my_giotto_object,
    feat_clusters = cluster_genes,
    name = "cluster_metagene"
)

my_spatial_genes <- names(cluster_genes)

my_giotto_object <- runPCA(
    gobject = my_giotto_object,
    feats_to_use = my_spatial_genes,
    name = "custom_pca"
)

my_giotto_object <- runUMAP(my_giotto_object,
    dim_reduction_name = "custom_pca",
    dimensions_to_use = 1:20,
    name = "custom_umap"
)

my_giotto_object <- createNearestNetwork(
    gobject = my_giotto_object,
    dim_reduction_name = "custom_pca",
    dimensions_to_use = 1:20,
    k = 3,
    name = "custom_NN"
)

my_giotto_object <- doLeidenCluster(
    gobject = my_giotto_object,
    network_name = "custom_NN",
    resolution = 0.1,
    n_iterations = 1000,
    name = "custom_leiden"
)

spatPlot2D(my_giotto_object,
    show_image = FALSE,
    cell_color = "custom_leiden",
    cell_color_code = c(
        "#eb4034",
        "#5877e8",
        "#ebd834",
        "#9beb34",
        "#6fab6a",
        "#24703f",
        "#58e8cb",
        "#58d0e8",
        "#eb8f34",
        "#7f58e8",
        "#d758e8",
        "#e85892"
    ),
    point_size = 3.5,
    background_color = "black",
    title = "Spatially informed clustering"
)

################################# Session info  ################################

sessionInfo()

# R version 4.3.3 (2024-02-29)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.4.1
#
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# time zone: America/New_York
# tzcode source: internal
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] Giotto_4.0.6      GiottoClass_0.2.4
#
# loaded via a namespace (and not attached):
#   [1] colorRamp2_0.1.0            bitops_1.0-7                deldir_2.0-4
# [4] rlang_1.1.3                 magrittr_2.0.3              GiottoUtils_0.1.6
# [7] matrixStats_1.3.0           compiler_4.3.3              systemfonts_1.0.6
# [10] png_0.1-8                   vctrs_0.6.5                 stringr_1.5.1
# [13] pkgconfig_2.0.3             SpatialExperiment_1.12.0    crayon_1.5.2
# [16] fastmap_1.1.1               backports_1.4.1             magick_2.8.3
# [19] XVector_0.42.0              labeling_0.4.3              utf8_1.2.4
# [22] rmarkdown_2.26              ragg_1.3.0                  xfun_0.43
# [25] zlibbioc_1.48.2             beachmat_2.18.1             GenomeInfoDb_1.38.8
# [28] jsonlite_1.8.8              flashClust_1.01-2           DelayedArray_0.28.0
# [31] BiocParallel_1.36.0         terra_1.7-71                irlba_2.3.5.1
# [34] parallel_4.3.3              cluster_2.1.6               R6_2.5.1
# [37] stringi_1.8.3               reticulate_1.36.0           GenomicRanges_1.54.1
# [40] estimability_1.5            Rcpp_1.0.12                 SummarizedExperiment_1.32.0
# [43] knitr_1.45                  R.utils_2.12.3              FNN_1.1.4
# [46] IRanges_2.36.0              Matrix_1.6-5                igraph_2.0.3
# [49] tidyselect_1.2.1            rstudioapi_0.16.0           abind_1.4-5
# [52] yaml_2.3.8                  codetools_0.2-20            lattice_0.22-6
# [55] tibble_3.2.1                Biobase_2.62.0              withr_3.0.0
# [58] evaluate_0.23               pillar_1.9.0                MatrixGenerics_1.14.0
# [61] checkmate_2.3.1             DT_0.33                     stats4_4.3.3
# [64] generics_0.1.3              dbscan_1.1-12               sp_2.1-3
# [67] RCurl_1.98-1.14             S4Vectors_0.40.2            ggplot2_3.5.0
# [70] munsell_0.5.1               scales_1.3.0                gtools_3.9.5
# [73] xtable_1.8-4                leaps_3.1                   glue_1.7.0
# [76] emmeans_1.10.1              scatterplot3d_0.3-44        tools_4.3.3
# [79] GiottoVisuals_0.1.6         data.table_1.15.4           ScaledMatrix_1.10.0
# [82] mvtnorm_1.2-4               cowplot_1.1.3               grid_4.3.3
# [85] colorspace_2.1-0            SingleCellExperiment_1.24.0 GenomeInfoDbData_1.2.11
# [88] BiocSingular_1.18.0         cli_3.6.2                   rsvd_1.0.5
# [91] textshaping_0.3.7           fansi_1.0.6                 S4Arrays_1.2.1
# [94] dplyr_1.1.4                 uwot_0.1.16                 gtable_0.3.4
# [97] R.methodsS3_1.8.2           digest_0.6.35               BiocGenerics_0.48.1
# [100] SparseArray_1.2.4           ggrepel_0.9.5               farver_2.1.1
# [103] FactoMineR_2.10             rjson_0.2.21                htmlwidgets_1.6.4
# [106] htmltools_0.5.8.1           R.oo_1.26.0                 lifecycle_1.0.4
# [109] multcompView_0.1-10         MASS_7.3-60.0.1
