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

my_spatial_genes <- names(cluster_genes)

giotto_object <- runPCA(gobject = giotto_object,
                        feats_to_use = my_spatial_genes,
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
           point_size = 0.7,
           background_color = "black",
           title = "Spatially informed clustering"
)

################################# Session info  ################################

sessionInfo()

