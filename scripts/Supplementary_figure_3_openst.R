## %######################################################%##
#                                                           #
####   Supplementary figure 3 Open-ST dataset            ####
#                                                           #
## %######################################################%##

############################## Download dataset  ###############################
## Download the adult mouse hippocampus dataset from the paper
## [Open-ST: High-resolution spatial transcriptomics in 3D](https://www.sciencedirect.com/science/article/pii/S0092867424006366)
## available in this link https://rajewsky-lab.github.io/openst/latest/examples/datasets/
## Download the file adult_mouse_hippocampus_by_cell.h5ad.tar.gz and untar it to run this tutorial.

dir.create("data")

data_path <- "data"

## Get the data
download.file(url = "http://bimsbstatic.mdc-berlin.de/rajewsky/openst-public-data/adult_mouse_by_cell.h5ad.tar.gz",
              destfile = file.path(data_path, "adult_mouse_by_cell.h5ad.tar.gz"))

## Untar the file
untar(tarfile = file.path(data_path, "adult_mouse_by_cell.h5ad.tar.gz"),
      exdir = data_path)

############################## Create the object  ##############################

library(Giotto)

## create the object directly from the h5ad file
giotto_object <- anndataToGiotto(file.path(data_path, "openst_demo_adult_mouse_by_cell.h5ad"),
                                 python_path = NULL)

## update instructions
instructions(giotto_object, "save_plot") <- TRUE
instructions(giotto_object, "save_dir") <- "results"
instructions(giotto_object, "show_plot") <- FALSE
instructions(giotto_object, "return_plot") <- FALSE

## The cell with the identifier 0 will have all the transcripts belonging to the background. Therefore, make sure to omit it from the dataset.
cell_metadata <- pDataDT(giotto_object)
cell_metadata <- cell_metadata[cell_metadata$cell_ID != "0",]

giotto_object <- subsetGiotto(giotto_object,
                              cell_ids = cell_metadata$cell_ID)

#################################### Filtering  ################################

giotto_object <- filterGiotto(gobject = giotto_object,
                              expression_threshold = 1,
                              feat_det_in_min_cells = 25,
                              min_det_feats_per_cell = 100,
                              verbose = TRUE)

################################ Normalization  ################################

giotto_object <- normalizeGiotto(giotto_object,
                                 scale_feats = FALSE,
                                 scale_cells = FALSE)

################################# Statistics  ##################################

giotto_object <- addStatistics(giotto_object)

############################## Dimension reduction  ############################

giotto_object <- runPCA(giotto_object,
                        ncp = 50)

################################### Clustering  ################################

giotto_object <- runUMAP(giotto_object,
                         dimensions_to_use = 1:10)

giotto_object <- createNearestNetwork(giotto_object)

giotto_object <- doLeidenCluster(giotto_object,
                                 resolution = 0.5)

spatPlot2D(giotto_object,
           cell_color = "leiden_clus",
           point_size = 0.7,
           background_color = "black")

################################# Session info  ################################

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
