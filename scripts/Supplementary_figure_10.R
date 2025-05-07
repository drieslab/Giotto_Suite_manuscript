######################## Download Visium Brain data #######################

# Download data from:
# https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain

########################### Load Libraries ################################

library(Giotto)
if (!requireNamespace("gstat", quietly = TRUE)) {
    install.packages("gstat")
}
if (!requireNamespace("future", quietly = TRUE)) {
    install.packages(c("future", "future.apply"))
}

# setup parallelization
future::plan(future::multisession())

################## Instantiate Data/Save Locations #######################

my_seed_num <- 315
data_directory <- "PATH TO DOWNLOADED DATA"
stardist_cell_poly_path = "scripts/CREATE/EXT_DATA/stardist_only_cell_bounds.geojson"
save_directory <- "scripts/FIGURE/S10/"
interp_directory <- "scripts/CREATE/OUTS/interp_rasters/"

########################### core kriging function ############################

# (this is now implemented more conveniently as `Giotto::interpolateFeature()`)
krige <- function(gobject,
    feats,
    ext, # region to krige across
    spat_unit = "cell",
    spat_loc_name = "raw",
    buffer = 50,
    rastersize = 500,
    overwrite = FALSE,
    savedir = tempdir()) {
    
    # Create directory if it doesn't exist
    if (!dir.exists(savedir)) dir.create(savedir, recursive = TRUE)
    
    # xy coordinates to initiate with
    sl <- getSpatialLocations(gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table",
        set_defaults = FALSE,
        copy_obj = FALSE,
        verbose = FALSE
    )
    
    # values to krige
    data <- spatValues(gobject, 
        feats = feats,
        spat_unit = spat_unit,
        expression_values = "raw", # pull raw counts
        verbose = FALSE
    )
    
    # combine xy and values into single table
    annotatedlocs <- merge(sl, data, by = "cell_ID")
    
    # create numeric representation since SpatExtent is not passable to
    # workers
    e_numeric <- ext[]
    
    future.apply::future_lapply(feats, function(feat) {
        filename <- file.path(savedir, paste0(feat, ".tif"))
        
        # create subset table with only relevant data
        data <- annotatedlocs[, c("cell_ID", feat, "sdimx", "sdimy"), with = FALSE]
        data.table::setnames(data, old = feat, new = "count")
        
        # model to use
        model <- gstat::gstat(
            id = feat,
            formula = count ~ 1,
            locations = ~ sdimx + sdimy,
            data = data,
            nmax = 7,
            set = list(
                idp = 0.5
            )
        )
        
        # regenerate extent
        interp_e <- terra::ext(e_numeric)
        
        # generate raster to interpolate on
        res <- max(range(interp_e) / rastersize)
        r <- terra::rast(x = interp_e, res = res)
        
        # perform interpolation and save to disk
        terra::interpolate(
            object = r,
            model = model,
            xyNames = c("sdimx", "sdimy"),
            debug.level = 0,
            index = 1L,
            filename = filename,
            overwrite = overwrite,
            wopt = list(
                progress = FALSE,
                filetype = "COG"
            )
        )
    },
    future.seed = TRUE,
    future.globals = list(
        e_numeric = e_numeric,
        annotatedlocs = annotatedlocs,
        savedir = savedir,
        overwrite = overwrite
    ),
    future.packages = c(
        "terra", "gstat", "data.table"
    ))
    return(invisible())
}

################## Create Instructions ##############################

instrs <- createGiottoInstructions(
    show_plot = TRUE,
    save_plot = TRUE,
    return_plot = FALSE,
    save_dir = save_directory
)

################## Create Giotto Object #############################

v_brain = createGiottoVisiumObject(data_directory, 
    gene_column_index = 2,
    instructions = instrs
)

######### Filter, normalize, ID hvf, dimension reductions, clustering ##########


# Subset to in tissue only
v_brain <- subset(v_brain, subset = in_tissue == 1)

# Filter
v_brain = filterGiotto(v_brain,
    expression_threshold = 1,
    feat_det_in_min_cells = 50,
    min_det_feats_per_cell = 1000,
    expression_values = c('raw')
)

# Normalize
v_brain = normalizeGiotto(v_brain, scalefactor = 6000, verbose = TRUE)

# Add stats
v_brain = addStatistics(gobject = v_brain)

# ID HVF
v_brain = calculateHVF(gobject = v_brain)

# Dimension Reductions
v_brain = runPCA(gobject = v_brain)

v_brain = runUMAP(v_brain,
    dimensions_to_use = 1:10,
    n_neighbors = 15,
    set_seed = TRUE,
    seed_number = my_seed_num
)

# NN Network
v_brain = createNearestNetwork(v_brain,
    dimensions_to_use = 1:10,
    k = 15
)
# Leiden Cluster
v_brain = doLeidenCluster(v_brain,
    resolution = 0.4,
    n_iterations = 200,
    seed_number = my_seed_num
)

# Spatial Network (kNN)
v_brain <- createSpatialNetwork(v_brain,
    method = 'kNN',
    k = 5,
    maximum_distance_knn = 400,
    name = 'spatial_network'
)

# Find Spatially Variable Features
ranktest = binSpect(v_brain,
    bin_method = 'rank',
    calc_hub = TRUE,
    hub_min_int = 5,
    spatial_network_name = 'spatial_network',
    seed = my_seed_num
)

# !!! Here we will proceed with these top 500 features from binspect.
top_spatial_features = head(ranktest$feats, 500)
# you can instead further refine the top features for for genes that are
# representative of the main spatial patterns

REFINE_TOP_SPATIAL_FEATURES <- FALSE
if (REFINE_TOP_SPATIAL_FEATURES) {
    # cluster the top 2000 spatial features from binSpect() into 20 clusters
    # Find pairwise distances between them
    spat_cor_netw_DT = detectSpatialCorFeats(v_brain,
        method = 'network',
        spatial_network_name = 'spatial_network',
        subset_feats = head(ranktest$feats, 2000)
    )
    
    # cluster spatially correlated features into similar spatial signatures
    spat_cor_clus = clusterSpatialCorFeats(spat_cor_netw_DT,
        name = 'spat_netw_clus',
        k = 20
    )
    
    heatmSpatialCorFeats(v_brain,
        spatCorObject = spat_cor_clus,
        use_clus_name = 'spat_netw_clus',
        heatmap_legend_param = list(title = NULL),
        show_plot = TRUE,
        save_plot = TRUE,
        return_plot = FALSE
    )
    
    # select at most 25 features from each of the 20 spatial clusters to get a
    # balanced representation of the spatial variation.
    # top 497 spatially variable features
    svgs <- getBalancedSpatCoexpressionFeats(spatCorObject = spat_cor_clus,
        seed = my_seed_num, 
        maximum = 25)
    top_spatial_features <- names(svgs)
}


################## kriging for `top_spatial_features` ######################

# Takes roughly 5 min
GENERATE_INTERP_RASTERS <- FALSE
if (GENERATE_INTERP_RASTERS) {
    krige(v_brain,
        spat_unit = "cell",
        feats = top_spatial_features,
        ext = ext(getPolygonInfo(v_brain)),
        savedir = interp_directory
    )
}

# load created rasters
img_list <- lapply(top_spatial_features, function(f) {
    createGiottoLargeImage(file.path(interp_directory, sprintf("%s.tif", f)), 
        name = f,
        use_rast_ext = TRUE, 
        verbose = FALSE)
})
v_brain <- setGiotto(v_brain, img_list)



# attach stardist segmentations
stardist <- createGiottoPolygon(stardist_cell_poly_path,
    name = "stardist_cell", 
    remove_background_polygon = TRUE,
    make_valid = TRUE
)
# this is inverted compared to the rest of the data
stardist <- flip(stardist)
v_brain <- setGiotto(v_brain, stardist)


# calculate overlaps on the interpolated data

v_brain <- calculateOverlap(v_brain,
    name_overlap = "interp_rna",
    spatial_info = "cell",
    image_names = top_spatial_features
)

v_brain <- calculateOverlap(v_brain,
    name_overlap = "interp_rna",
    spatial_info = "stardist_cell",
    image_names = top_spatial_features
)

v_brain <- overlapToMatrix(v_brain,
    poly_info = "cell",
    type = "intensity",
    feat_info = "interp_rna"
)

v_brain <- overlapToMatrix(v_brain,
    poly_info = "stardist_cell",
    type = "intensity",
    feat_info = "interp_rna"
)




############# Subcellular, Cell: Dimension reductions and Clustering ###########

activeSpatUnit(v_brain) <- "cell" # visium spots
activeFeatType(v_brain) <- "interp_rna"

v_brain = normalizeGiotto(v_brain)

# PCA
v_brain = runPCA(v_brain,
    set_seed = TRUE,
    seed_number = my_seed_num
)

# UMAP
v_brain = runUMAP(v_brain,
    dimensions_to_use = 1:10,
    n_neighbors = 25,
    min_dist = 0.5,
    spread = 1,
    seed_number = my_seed_num
)

# NN Network
v_brain = createNearestNetwork(v_brain,
    dimensions_to_use = 1:10,
    k = 25
)

# Leiden Cluster
v_brain = doLeidenCluster(v_brain,
    resolution = 0.15,
    n_iterations = 100,
    seed_number = my_seed_num
)

v_brain = doKmeans(v_brain,
    centers = 12,
    seed_number = my_seed_num
)




########### Subcellular, StarDist: Dimension reductions and Clustering #########

activeSpatUnit(v_brain) <- "stardist_cell"
activeFeatType(v_brain) <- "interp_rna"

# Filter out any polygons with expression of 0
v_brain = filterGiotto(v_brain,
    expression_values = "raw",
    expression_threshold = 1,
    feat_det_in_min_cells = 0,
    min_det_feats_per_cell = 1
)

v_brain = normalizeGiotto(v_brain)

# PCA
v_brain = runPCA(v_brain,
    feats_to_use = NULL,
    set_seed = TRUE,
    seed_number = my_seed_num
)

# UMAP
v_brain = runUMAP(v_brain,
    dimensions_to_use = 1:10,
    n_neighbors = 25,
    min_dist = 0.01,
    spread = 1,
    seed_number = my_seed_num
)

# NN Network
v_brain = createNearestNetwork(v_brain,
    dimensions_to_use = 1:10,
    k = 25
)
    
# Leiden Cluster
v_brain = doLeidenCluster(v_brain,
    resolution = 0.2,
    n_iterations = 100,
    seed_number = my_seed_num
)

# K-Means
v_brain = doKmeans(v_brain,
    centers = 12,
    iter_max = 50,
    nstart = 100,
    seed_number = my_seed_num
)



############################ Visualize and Compare #############################

# leiden clustering results #########
spatInSituPlotPoints(v_brain,
    polygon_feat_type = "cell",
    feat_type = "interp_rna",
    polygon_fill = "leiden_clus",
    polygon_fill_as_factor = TRUE,
    polygon_line_size = 0.1,
    polygon_color = "black",
    polygon_alpha = 1,
    show_image = TRUE,
    save_plot = FALSE,
    image_name = "image",
)

spatInSituPlotPoints(v_brain,
    polygon_feat_type = "stardist_cell",
    feat_type = "interp_rna",
    polygon_fill = "leiden_clus",
    polygon_fill_as_factor = TRUE,
    polygon_line_size = 0,
    polygon_alpha = 1,
    show_image = TRUE,
    image_name = "image",
    save_plot = FALSE
)


# FIG S10 C kmeans clustering results #########

# COLORS FOR kmeans=12 mapping

visium_colors <- getDistinctColors(12)
stardist_colors <- visium_colors[c(8, 7, 9, 2, 11, 10, 4, 12, 1, 3, 5, 6)]

# visium_colors = getDistinctColors(12)
# stardist_colors = visium_colors
# 
names(visium_colors) <- names(stardist_colors) <- seq(12)
# names(stardist_colors) <- seq(12)

spatInSituPlotPoints(v_brain,
    polygon_feat_type = "cell",
    feat_type = "interp_rna",
    polygon_fill = "kmeans",
    polygon_fill_as_factor = TRUE,
    polygon_alpha = 1,
    polygon_line_size = 0.1,
    polygon_color = "black",
    show_image = TRUE,
    image_name = "image",
    polygon_fill_code = visium_colors,
    save_param = list(
        save_name = "C_visium",
        save_format = "svg"
    )
)

spatInSituPlotPoints(v_brain,
    polygon_feat_type = "stardist_cell",
    feat_type = "interp_rna",
    polygon_fill = "kmeans",
    polygon_fill_as_factor = TRUE,
    polygon_alpha = 1,
    polygon_line_size = 0,
    show_image = TRUE,
    image_name = "image",
    polygon_fill_code = stardist_colors,
    save_param = list(
        save_name = "C_stardist",
        save_format = "svg"
    )
)

spatInSituPlotPoints(v_brain,
    polygon_feat_type = "cell",
    feat_type = "interp_rna",
    polygon_fill = "kmeans",
    polygon_fill_as_factor = TRUE,
    polygon_alpha = 1,
    polygon_line_size = 0.3,
    polygon_color = "black",
    show_image = TRUE,
    image_name = "image",
    polygon_fill_code = visium_colors,
    xlim = c(2500, 3500),
    ylim = c(-5000, -3000),
    save_param = list(
        save_name = "C_visium_zoom",
        save_format = "svg",
        base_width = 5.5
    )
)

spatInSituPlotPoints(v_brain,
    polygon_feat_type = "stardist_cell",
    feat_type = "interp_rna",
    polygon_fill = "kmeans",
    polygon_fill_as_factor = TRUE,
    polygon_alpha = 1,
    polygon_line_size = 0.3,
    polygon_color = "black",
    show_image = TRUE,
    image_name = "image",
    polygon_fill_code = stardist_colors,
    xlim = c(2500, 3500),
    ylim = c(-5000, -3000),
    save_param = list(
        save_name = "C_stardist_zoom",
        save_format = "svg",
        base_width = 5.5
    )
)

# FIG S10 D Feature plots ###########

spatInSituPlotPoints(v_brain,
    polygon_feat_type = "cell",
    feat_type = "interp_rna",
    polygon_fill = "Pantr1",
    expression_values = "normalized",
    polygon_fill_gradient_style = "sequential",
    polygon_alpha = 1,
    polygon_line_size = 0.1,
    polygon_color = "black",
    background_color = "black",
    save_param = list(
        save_name = "D_visium_pantr1",
        save_format = "svg",
        base_width = 8
    )
)

spatInSituPlotPoints(v_brain,
    polygon_feat_type = "stardist_cell",
    feat_type = "interp_rna",
    polygon_fill = "Pantr1",
    expression_values = "normalized",
    polygon_fill_gradient_style = "sequential",
    polygon_alpha = 1,
    polygon_line_size = 0,
    background_color = "black",
    save_param = list(
        save_name = "D_stardist_pantr1",
        save_format = "svg",
        base_width = 8
    )
)
