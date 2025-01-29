## %######################################################%##
#                                                           #
####   Supplementary figure 3 Slide-seq dataset          ####
#                                                           #
## %######################################################%##

############################## Download dataset  ###############################
## Download the mouse brain dataset deposited in the NeMO database
## by the Macosko lab under the grant rf1_macosko

data_path <- "data"

## Get the expression data
download.file(url = "https://data.nemoarchive.org/biccn/grant/rf1_macosko/macosko/spatial_transcriptome/cellgroup/Slide-seq/mouse/processed/counts/2020-12-19_Puck_201112_26.matched.digital_expression.mex.tar.gz",
              destfile = file.path(data_path, "2020-12-19_Puck_201112_26.matched.digital_expression.mex.tar.gz"))

## Get the spatial coordinates
download.file(url = "https://data.nemoarchive.org/biccn/grant/rf1_macosko/macosko/spatial_transcriptome/cellgroup/Slide-seq/mouse/processed/other/2020-12-19_Puck_201112_26.BeadLocationsForR.csv.tar",
              destfile = file.path(data_path, "2020-12-19_Puck_201112_26.BeadLocationsForR.csv.tar"))

## Untar the expression files
untar(tarfile = file.path(data_path, "2020-12-19_Puck_201112_26.matched.digital_expression.mex.tar.gz"),
      exdir = data_path)

############################## Create the object  ##############################
library(Giotto)

instructions <- createGiottoInstructions(save_plot = TRUE,
                                         save_dir = "results",
                                         show_plot = FALSE,
                                         return_plot = FALSE,
                                         python_path = NULL)

expression_matrix <- get10Xmatrix(file.path(data_path, "2020-12-19_Puck_201112_26.matched.digital_expression"))

spatial_locs <- data.table::fread(file.path(data_path, "2020-12-19_Puck_201112_26.BeadLocationsForR.csv.tar"))

spatial_locs <- spatial_locs[spatial_locs$barcodes %in% colnames(expression_matrix),]

giotto_object <- createGiottoObject(expression = expression_matrix,
                                    spatial_locs = spatial_locs,
                                    instructions = instructions)

#################################### Filtering  ################################

giotto_object <- filterGiotto(giotto_object,
                              min_det_feats_per_cell = 10,
                              feat_det_in_min_cells = 10)

################################ Normalization  ################################

giotto_object <- normalizeGiotto(giotto_object,
                                 scale_feats = FALSE,
                                 scale_cells = FALSE)

################################# Statistics  ##################################

giotto_object <- addStatistics(giotto_object)

############################ Find spatial genes ################################

giotto_object <- createSpatialNetwork(gobject = giotto_object,
                                      method = "kNN",
                                      k = 6,
                                      name = "spatial_network")

ranktest <- binSpect(giotto_object,
                     bin_method = "rank",
                     calc_hub = TRUE,
                     hub_min_int = 5,
                     spatial_network_name = "spatial_network")

## cluster the top 500 spatial genes into clusters
ext_spatial_genes <- ranktest[1:500, ]$feats

# calculate pairwise distances between genes
spat_cor_netw_DT <- detectSpatialCorFeats(giotto_object,
                                          method = "network",
                                          spatial_network_name = "spatial_network",
                                          subset_feats = ext_spatial_genes)

# identify potential spatial co-expression
spat_cor_netw_DT <- clusterSpatialCorFeats(spat_cor_netw_DT,
                                           name = "spat_netw_clus",
                                           k = 10)

# create metagenes/co-expression modules
cluster_genes <- getBalancedSpatCoexpressionFeats(spat_cor_netw_DT,
                                                  maximum = 30)

giotto_object <- createMetafeats(giotto_object,
                                 feat_clusters = cluster_genes,
                                 name = "cluster_metagene")

spatial_genes <- names(cluster_genes)

################### Calculate spatial informed clusters ########################

giotto_object <- runPCA(gobject = giotto_object,
                        feats_to_use = spatial_genes,
                        name = "custom_pca")

giotto_object <- runUMAP(giotto_object,
                         dim_reduction_name = "custom_pca",
                         dimensions_to_use = 1:20,
                         name = "custom_umap")

giotto_object <- createNearestNetwork(gobject = giotto_object,
                                      dim_reduction_name = "custom_pca",
                                      dimensions_to_use = 1:20,
                                      k = 30,
                                      name = "custom_NN")

giotto_object <- doLeidenCluster(gobject = giotto_object,
                                 network_name = "custom_NN",
                                 resolution = 0.5,
                                 name = "custom_leiden")

spatPlot2D(giotto_object,
           show_image = FALSE,
           cell_color = "custom_leiden",
           point_size = 1,
           background_color = "black")

################################# Session info  ################################

sessionInfo()

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sequoia 15.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
  [1] Giotto_4.2.0      GiottoClass_0.4.7

loaded via a namespace (and not attached):
  [1] colorRamp2_0.1.0            remotes_2.5.0
[3] rlang_1.1.5                 magrittr_2.0.3
[5] RcppAnnoy_0.0.22            GiottoUtils_0.2.3
[7] matrixStats_1.5.0           compiler_4.4.2
[9] systemfonts_1.2.1           png_0.1-8
[11] callr_3.7.6                 vctrs_0.6.5
[13] pkgconfig_2.0.3             SpatialExperiment_1.14.0
[15] crayon_1.5.3                fastmap_1.2.0
[17] backports_1.5.0             magick_2.8.5
[19] XVector_0.44.0              labeling_0.4.3
[21] UCSC.utils_1.0.0            ps_1.8.1
[23] ragg_1.3.3                  purrr_1.0.2
[25] zlibbioc_1.50.0             beachmat_2.20.0
[27] GenomeInfoDb_1.40.1         jsonlite_1.8.9
[29] DelayedArray_0.30.1         BiocParallel_1.38.0
[31] terra_1.8-15                irlba_2.3.5.1
[33] parallel_4.4.2              R6_2.5.1
[35] RColorBrewer_1.1-3          reticulate_1.40.0
[37] GenomicRanges_1.56.1        scattermore_1.2
[39] Rcpp_1.0.14                 SummarizedExperiment_1.34.0
[41] R.utils_2.12.3              IRanges_2.38.1
[43] Matrix_1.7-2                igraph_2.1.4
[45] tidyselect_1.2.1            rstudioapi_0.17.1
[47] abind_1.4-8                 codetools_0.2-20
[49] curl_6.2.0                  processx_3.8.5
[51] pkgbuild_1.4.6              lattice_0.22-6
[53] tibble_3.2.1                Biobase_2.64.0
[55] withr_3.0.2                 desc_1.4.3
[57] pillar_1.10.1               MatrixGenerics_1.16.0
[59] checkmate_2.3.2             stats4_4.4.2
[61] plotly_4.10.4               generics_0.1.3
[63] dbscan_1.2.2                S4Vectors_0.42.1
[65] ggplot2_3.5.1               munsell_0.5.1
[67] scales_1.3.0                GiottoData_0.2.14
[69] gtools_3.9.5                glue_1.8.0
[71] lazyeval_0.2.2              tools_4.4.2
[73] GiottoVisuals_0.2.11        data.table_1.16.4
[75] ScaledMatrix_1.12.0         cowplot_1.1.3
[77] grid_4.4.2                  tidyr_1.3.1
[79] colorspace_2.1-1            SingleCellExperiment_1.26.0
[81] GenomeInfoDbData_1.2.12     BiocSingular_1.20.0
[83] cli_3.6.3                   rsvd_1.0.5
[85] textshaping_1.0.0           S4Arrays_1.4.1
[87] viridisLite_0.4.2           dplyr_1.1.4
[89] uwot_0.2.2                  gtable_0.3.6
[91] R.methodsS3_1.8.2           digest_0.6.37
[93] BiocGenerics_0.50.0         SparseArray_1.4.8
[95] ggrepel_0.9.6               farver_2.1.2
[97] rjson_0.2.23                htmlwidgets_1.6.4
[99] htmltools_0.5.8.1           R.oo_1.27.0
[101] lifecycle_1.0.4             httr_1.4.7

