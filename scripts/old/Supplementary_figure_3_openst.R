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

############################ Find spatial genes ################################

giotto_object <- createSpatialNetwork(gobject = giotto_object,
                                      method = "kNN",
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
           point_size = 1,
           cell_color = "custom_leiden",
           background_color = "black")

################################# Session info  ################################

sessionInfo()

R version 4.4.2 (2024-10-31)
Platform: x86_64-apple-darwin20
Running under: macOS Sequoia 15.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
    [1] Giotto_4.2.0      GiottoClass_0.4.7

loaded via a namespace (and not attached):
    [1] tidyselect_1.2.1            viridisLite_0.4.2          
[3] farver_2.1.2                dplyr_1.1.4                
[5] GiottoVisuals_0.2.11        R.utils_2.12.3             
[7] fastmap_1.2.0               SingleCellExperiment_1.28.1
[9] lazyeval_0.2.2              digest_0.6.37              
[11] rsvd_1.0.5                  lifecycle_1.0.4            
[13] terra_1.8-15                magrittr_2.0.3             
[15] dbscan_1.2-0                compiler_4.4.2             
[17] rlang_1.1.5                 tools_4.4.2                
[19] igraph_2.1.4                data.table_1.16.4          
[21] labeling_0.4.3              S4Arrays_1.6.0             
[23] htmlwidgets_1.6.4           reticulate_1.40.0          
[25] DelayedArray_0.32.0         RColorBrewer_1.1-3         
[27] abind_1.4-8                 BiocParallel_1.40.0        
[29] withr_3.0.2                 purrr_1.0.2                
[31] R.oo_1.27.0                 BiocGenerics_0.52.0        
[33] grid_4.4.2                  stats4_4.4.2               
[35] beachmat_2.22.0             colorspace_2.1-1           
[37] ggplot2_3.5.1               scales_1.3.0               
[39] gtools_3.9.5                SummarizedExperiment_1.36.0
[41] cli_3.6.3                   crayon_1.5.3               
[43] ragg_1.3.3                  generics_0.1.3             
[45] rstudioapi_0.17.1           httr_1.4.7                 
[47] rjson_0.2.23                zlibbioc_1.52.0            
[49] parallel_4.4.2              XVector_0.46.0             
[51] matrixStats_1.5.0           vctrs_0.6.5                
[53] Matrix_1.7-2                jsonlite_1.8.9             
[55] BiocSingular_1.22.0         IRanges_2.40.1             
[57] S4Vectors_0.44.0            ggrepel_0.9.6              
[59] irlba_2.3.5.1               scattermore_1.2            
[61] systemfonts_1.2.1           magick_2.8.5               
[63] GiottoUtils_0.2.3           plotly_4.10.4              
[65] tidyr_1.3.1                 glue_1.8.0                 
[67] codetools_0.2-20            uwot_0.2.2                 
[69] cowplot_1.1.3               RcppAnnoy_0.0.22           
[71] gtable_0.3.6                GenomeInfoDb_1.42.1        
[73] GenomicRanges_1.58.0        UCSC.utils_1.2.0           
[75] ScaledMatrix_1.14.0         munsell_0.5.1              
[77] tibble_3.2.1                pillar_1.10.1              
[79] htmltools_0.5.8.1           GenomeInfoDbData_1.2.13    
[81] R6_2.5.1                    textshaping_1.0.0          
[83] lattice_0.22-6              Biobase_2.66.0             
[85] png_0.1-8                   R.methodsS3_1.8.2          
[87] backports_1.5.0             SpatialExperiment_1.16.0   
[89] Rcpp_1.0.14                 SparseArray_1.6.0          
[91] checkmate_2.3.2             colorRamp2_0.1.0           
[93] MatrixGenerics_1.18.0       pkgconfig_2.0.3  

