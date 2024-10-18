## %######################################################%##
#                                                           #
####   Supplementary figure 3 Visium HD dataset          ####
#                                                           #
## %######################################################%##

# Download data
## Download the Visium HD Spatial Gene Expression Library,
## Human Colorectal Cancer (FFPE) from the 10X genomics website
## https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc

dir.create("data")

## Get the data
download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_spatial.tar.gz",
              destfile = file.path("data", "Visium_HD_Human_Colon_Cancer_spatial.tar.gz"))

download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz",
              destfile = file.path("data", "Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz"))

## Untar the file
untar(tarfile = file.path("data", "Visium_HD_Human_Colon_Cancer_spatial.tar.gz"),
      exdir = data_path)

untar(tarfile = file.path("data", "Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz"),
      exdir = data_path)

# Create the object
library(Giotto)

# Set instructions
results_folder <- "results/"

python_path <- NULL

data_path <- "data/binned_outputs/square_002um/"

## create instructions
instructions <- createGiottoInstructions(save_dir = results_folder,
                                         save_plot = TRUE,
                                         show_plot = FALSE,
                                         return_plot = FALSE,
                                         python_path = python_path)

## Read directly from the visium folder
readerHD <- importVisiumHD()

readerHD$data_path <- data_path

visiumHD <- readerHD$create_gobject(
  data_path = data_path,
  png_name = "tissue_lowres_image.png",
  gene_column_index = 2)

# Filtering
visiumHD <- filterGiotto(gobject = visiumHD,
                         expression_threshold = 1,
                         feat_det_in_min_cells = 40,
                         min_det_feats_per_cell = 50,
                         expression_values = "raw",
                         feat_type = "rna",
                         verbose = TRUE)

# Normalization
visiumHD <- normalizeGiotto(gobject = visiumHD,
                            scalefactor = 6000,
                            feat_type = "rna",
                            verbose = TRUE)

# Statistics
visiumHD <- addStatistics(gobject = visiumHD,
                          feat_type = "rna")

# Dimension reduction
visiumHD <- calculateHVF(gobject = visiumHD,
                         spat_unit = "hexagon400",
                         feat_type = "rna",
                         show_plot = TRUE)

visiumHD <- runPCA(gobject = visiumHD,
                   feat_type = "rna")

# Clustering
visiumHD <- runUMAP(visiumHD,
                    dimensions_to_use = 1:10,
                    feat_type = "rna")

visiumHD <- createNearestNetwork(gobject = visiumHD,
                                 dimensions_to_use = 1:10,
                                 k = 30,
                                 feat_type = "rna")

visiumHD <- doLeidenCluster(gobject = visiumHD,
                            feat_type = "rna",
                            resolution = 1,
                            n_iterations = 10)

# Identify Spatial Genes

featData <- fDataDT(visiumHD)
hvf_genes <- featData[hvf == "yes"]$feat_ID

visiumHD <- createSpatialNetwork(visiumHD,
                                 name = "kNN_network",
                                 method = "kNN",
                                 k = 8)

ranktest <- binSpect(visiumHD,
                     subset_feats = hvf_genes,
                     bin_method = "rank",
                     calc_hub = FALSE,
                     do_fisher_test = TRUE,
                     spatial_network_name = "kNN_network")

# Plot spatial gene groups

balanced_genes <- getBalancedSpatCoexpressionFeats(spatCorObject = spat_cor_netw_DT,
                                                   maximum = 5)
selected_feats <- names(balanced_genes)

# give genes from same cluster same color
giotto_colors <- getDistinctColors(n = 20)
names(giotto_colors) <- 1:20

my_colors <- giotto_colors[balanced_genes]
names(my_colors) <- names(balanced_genes)

spatInSituPlotPoints(visiumHD,
                     show_image = FALSE,
                     feats = list("rna" = selected_feats),
                     feats_color_code = my_colors,
                     show_legend = FALSE,
                     point_size = 0.20,
                     show_polygon = FALSE,
                     use_overlap = FALSE,
                     polygon_feat_type = "hexagon400",
                     polygon_bg_color = NA,
                     polygon_color = "white",
                     polygon_line_size = 0.01,
                     jitter = c(25,25))

# Session info
sessionInfo()

R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 15.0.1

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
    [1] Giotto_4.1.0      GiottoClass_0.3.5

loaded via a namespace (and not attached):
    [1] colorRamp2_0.1.0            deldir_2.0-4
[3] rlang_1.1.4                 magrittr_2.0.3
[5] RcppAnnoy_0.0.22            GiottoUtils_0.1.12
[7] matrixStats_1.4.1           compiler_4.4.1
[9] systemfonts_1.1.0           png_0.1-8
[11] vctrs_0.6.5                 reshape2_1.4.4
[13] stringr_1.5.1               pkgconfig_2.0.3
[15] SpatialExperiment_1.14.0    crayon_1.5.3
[17] fastmap_1.2.0               backports_1.5.0
[19] magick_2.8.5                XVector_0.44.0
[21] labeling_0.4.3              utf8_1.2.4
[23] UCSC.utils_1.0.0            ragg_1.3.3
[25] purrr_1.0.2                 zlibbioc_1.50.0
[27] beachmat_2.20.0             GenomeInfoDb_1.40.1
[29] jsonlite_1.8.9              DelayedArray_0.30.1
[31] BiocParallel_1.38.0         terra_1.7-78
[33] irlba_2.3.5.1               parallel_4.4.1
[35] R6_2.5.1                    stringi_1.8.4
[37] RColorBrewer_1.1-3          reticulate_1.39.0
[39] GenomicRanges_1.56.1        scattermore_1.2
[41] Rcpp_1.0.13                 SummarizedExperiment_1.34.0
[43] R.utils_2.12.3              IRanges_2.38.1
[45] Matrix_1.7-0                igraph_2.0.3
[47] tidyselect_1.2.1            rstudioapi_0.16.0
[49] abind_1.4-8                 codetools_0.2-20
[51] lattice_0.22-6              tibble_3.2.1
[53] plyr_1.8.9                  Biobase_2.64.0
[55] withr_3.0.1                 pillar_1.9.0
[57] MatrixGenerics_1.16.0       checkmate_2.3.2
[59] stats4_4.4.1                plotly_4.10.4
[61] generics_0.1.3              dbscan_1.2-0
[63] sp_2.1-4                    S4Vectors_0.42.1
[65] ggplot2_3.5.1               munsell_0.5.1
[67] scales_1.3.0                gtools_3.9.5
[69] glue_1.8.0                  lazyeval_0.2.2
[71] tools_4.4.1                 GiottoVisuals_0.2.4
[73] data.table_1.16.0           ScaledMatrix_1.12.0
[75] cowplot_1.1.3               grid_4.4.1
[77] tidyr_1.3.1                 colorspace_2.1-1
[79] SingleCellExperiment_1.26.0 GenomeInfoDbData_1.2.12
[81] BiocSingular_1.20.0         cli_3.6.3
[83] rsvd_1.0.5                  textshaping_0.4.0
[85] fansi_1.0.6                 S4Arrays_1.4.1
[87] viridisLite_0.4.2           dplyr_1.1.4
[89] uwot_0.2.2                  gtable_0.3.5
[91] R.methodsS3_1.8.2           digest_0.6.37
[93] BiocGenerics_0.50.0         SparseArray_1.4.8
[95] ggrepel_0.9.6               rjson_0.2.23
[97] htmlwidgets_1.6.4           farver_2.1.2
[99] htmltools_0.5.8.1           R.oo_1.26.0
[101] lifecycle_1.0.4             httr_1.4.7
