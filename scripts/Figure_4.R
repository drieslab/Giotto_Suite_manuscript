## %######################################################%##
#                                                          #
####            Figure 4 Multi-Segmentation             ####
#                                                          #
## %######################################################%##

## Download data
## Original Xenium data (In situ sample1, replicate 1) can be downloaded from 10X:
## https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
##
## After download, unzip it to "./Xenium/"
##
## For Julia and BAYSOR-0.5.2 installation, 10x does provide some instructions:
## https://www.10xgenomics.com/analysis-guides/using-baysor-to-perform-xenium-cell-segmentation.
## However we used some slightly modified steps as detailed below.

dir.create("data/Xenium")

download.file(
    url = "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip",
    destfile = "data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip"
)
download.file(
    url = "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.ome.tif",
    destfile = "data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.ome.tif"
)
download.file(
    url = "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_imagealignment.csv",
    destfile = "data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_imagealignment.csv"
)

untar(
    tarfile = "data/Xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip",
    exdir = "data/Xenium/"
)

############################# Segmentation Steps #############################
library(Giotto)

# ensure that a compatible Giotto segmentation python env exists
# see https://github.com/drieslab/Giotto/blob/suite/inst/python/configuration/genv_segmentation.yml
# a python env .yml

# for this script, the full environment will be installed
if (!reticulate::condaenv_exists("giotto_segmentation")) {
    installGiottoEnvironment(
        envname = "giotto_segmentation",
        python_version = "3.8",
    )
    reticulate::conda_install(
        envname = "giotto_segmentation",
        packages = c("tifffile", "deepcell", "stardist", "cellpose==3.1.0"),
        pip = TRUE
    )
}

# Multi Segmentation -------------------------------------------------- #
# Some of these steps may take significant time and computational resources.
# Pre-made outputs for these steps are available on Zenodo.
# --------------------------------------------------------------------- #


# Mesmer segmentation (below) requires a deepcell access token
# If you don't yet have a token, you can create one at
# https://users.deepcell.org
# ** This must be set before the python env is attached **
Sys.setenv(DEEPCELL_ACCESS_TOKEN = "???")

# attach segmentation python env
set_giotto_python_path("giotto_segmentation")



# filepaths
data_dir <- "data/Xenium/"
results_folder <- "results/figure_4/"
Segmentation_dir <- file.path(results_folder, "Segmentation_dir")
Xenium_dir <- file.path(data_dir, "outs") # xenium output bundle

dir.create(results_folder, recursive = TRUE)
dir.create(Segmentation_dir, recursive = TRUE)


# 1. Load needed transcripts info for Baysor
## transcripts
x <- importXenium(Xenium_dir)
x$qv <- 20 # default
tx_dt <- x$load_transcripts(
    flip_vertical = FALSE,
    # split_keyword only works with `giottoPoints` outputs
    output = "data.table" # default output would have been `giottoPoints`
)

# Prep transcript data for baysor (see 10x instructions for Baysor for details)
# * qv filter is already done
# * filter out control probes
tx_dt <- tx_dt[!grepl("BLANK|NegControlCodeword|NegControlProbe|antisense", feat_ID)]
# * set cell_ID -1 values to 0
tx_dt[cell_id == -1, cell_id := 0]
# * write to csv for baysor processing
data.table::fwrite(tx_dt, file = file.path(Segmentation_dir, "transcripts_flip.csv"))
rm(tx_dt)

# 2. Load aligned IF cell bounds image for image segmentation methods
IF_xen <- x$load_aligned_image(
    path = file.path(data_dir, "/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.ome.tif"),
    imagealignment_path = file.path(data_dir, "/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_imagealignment.csv"),
    name = "bound"
)
IF_xen <- flip(IF_xen) # flip across y to match original outputs

############## Multi_segmentation #############################################

# Set mini ROI to use
zoom <- ext(6000, 7000, 3800, 5100)
img_cropped <- crop(IF_xen, zoom)
# Xenium aligned images are aligned with affine transforms. Affines are a
# delayed transform, and are only performed according to what is needed.
# * render the delayed affine transform with `doDeferred()`
# * `size` param is the allowed number of pixels (full size)
IF_rast <- doDeferred(img_cropped, size = prod(dim(img_cropped)))
IF_subset_path <- file.path(Segmentation_dir, "IF_cropped.tif")

terra::writeRaster(
    IF_rast[], # drop to SpatRaster
    filename = IF_subset_path,
    overwrite = TRUE
)

# preserve affine image extent (this may be slightly different from
# the zoom extent due to increased margin size from the image transform)
write.csv(ext(IF_rast)[], file = file.path(Segmentation_dir, "zoomext.csv"))

doMesmerSegmentation(
    input = IF_subset_path,
    mask_output = file.path(Segmentation_dir, "mesmer.tiff"),
    nucleus_channel = 3,
    membrane_channel = 2,
    micron_scale = 0.2125,
    python_env = "giotto_segmentation"
)

doCellposeSegmentation(
    input = IF_subset_path,
    mask_output = file.path(Segmentation_dir, "cellpose.tiff"),
    model_name = "cyto3",
    channel_1 = 2,
    channel_2 = 3,
    python_env = "giotto_segmentation"
)

doStardistSegmentation(
    input = IF_subset_path,
    mask_output = file.path(Segmentation_dir, "stardist.tiff"),
    model_name = "2D_demo",
    nuclei_channel = 3,
    python_env = "giotto_segmentation"
)

# For baysor (v0.5.2):
# A high performance cluster is recommended.
#
# 1. Download the relevant asset from the Baysor github releases, then unzip it.
# https://github.com/kharchenkolab/Baysor/releases/tag/v0.5.2
# (we used `baysor_ubuntu-latest_x64_build.zip`)
#
# 2. Either load Julia as a module or install it.
# $ module load julia/1.7.3
#
# 3. Build binary using the downloaded files (this may take a while)
# $ julia ./bin/build.jl
#
# 4. Add built Baysor binary to PATH
#
# 5. Run baysor on the `transcripts.csv` that was written out earlier.
# We included the shell script (`scripts/Xenium_supporting_files/Baysor_run.sh`)
# that was run on our cluster as an example
# This shell script was moved to the `Segmentation_dir` and then run
# Will likely take a day or so of processing.
#
# 6. Place the Baysor output "segmentation_polygons.json" as "baysor.json"
# within the Segmentation_dir location.

# plotting steps ---------------------------------------------------------- #

# restart session to return to standard giotto_env python environment
.rs.restartR()

# reload values due to session restart
library(Giotto)

data_dir <- "data/Xenium/"
results_folder <- "results/figure_4/"
Segmentation_dir <- file.path(results_folder, "Segmentation_dir")
Xenium_dir <- file.path(data_dir, "outs") # xenium output bundle
IF_subset_path <- file.path(Segmentation_dir, "IF_cropped.tif")

instructions <- createGiottoInstructions(
    show_plot = FALSE,
    return_plot = FALSE,
    save_plot = TRUE,
    save_dir = results_folder
)

# Load transcripts as giottoPoints
x <- importXenium(Xenium_dir)
x$qv <- 20 # default

tx <- x$load_transcripts(
    flip_vertical = FALSE,
    split_keyword = list(
        c("BLANK"),
        c("NegControlCodeword"),
        c("NegControlProbe", "antisense")
    )
)

# original cell segmentation
cell <- x$load_polys()
cell <- flip(cell) # flip to match the data y orientation in the manuscript

# load extent info
zoomext <- ext(read.csv(file.path(Segmentation_dir, "zoomext.csv"))$x)
# load in cropped image
IF_rast <- createGiottoLargeImage(IF_subset_path, extent = zoomext)

# read in generated segmentation outputs
stardist_poly <- createGiottoPolygonsFromMask(
    file.path(Segmentation_dir, "stardist.tiff"),
    shift_vertical_step = FALSE,
    shift_horizontal_step = FALSE,
    calc_centroids = TRUE
)

mesmer_poly <- createGiottoPolygonsFromMask(
    file.path(Segmentation_dir, "mesmer.tiff"),
    shift_vertical_step = FALSE,
    shift_horizontal_step = FALSE,
    calc_centroids = TRUE
)

cellpose_poly <- createGiottoPolygonsFromMask(
    file.path(Segmentation_dir, "cellpose.tiff"),
    shift_vertical_step = FALSE,
    shift_horizontal_step = FALSE,
    calc_centroids = TRUE
)

# original
cell_ROI <- crop(cell, zoomext)

ext(stardist_poly) <- zoomext
# expand slightly for cell as opposed to nucleus
stardist_poly <- buffer(stardist_poly, 5)
ext(mesmer_poly) <- zoomext
ext(cellpose_poly) <- zoomext

# BAYSOR
# No need to register as it was generated using transcript coordinates
Baysor_JSON <- file.path(Segmentation_dir, "segmentation_polygons.json")
Baysor_poly <- createGiottoPolygonsFromGeoJSON(
    Baysor_JSON,
    name = "Baysor",
    calc_centroids = TRUE,
    make_valid = TRUE # needed for baysor 0.5.2 outputs
)
Baysor_ROI <- crop(Baysor_poly, zoomext)

# preview plots to check alignment
# original
plot(IF_rast)
plot(cell_ROI, add = TRUE, border = "magenta", lwd = 0.5)
# stardist
plot(IF_rast)
plot(stardist_poly, add = TRUE, border = "magenta", lwd = 0.5)
# mesmer
plot(IF_rast)
plot(mesmer_poly, add = TRUE, border = "magenta", lwd = 0.5)
# cellpose
plot(IF_rast)
plot(cellpose_poly, add = TRUE, border = "magenta", lwd = 0.5)
# baysor
plot(IF_rast)
plot(Baysor_ROI, add = TRUE, border = "magenta", lwd = 0.5)

# crop transcripts to match
rna <- crop(tx$rna, zoomext)

############## Figure 4a ######################################################

## Create Giotto Xenium Object
xen_ROI <- createGiottoObjectSubcellular(
    gpoints = list("rna" = rna),
    gpolygons = list(
        "cell" = cell_ROI,
        "Baysor" = Baysor_ROI,
        "CellPose" = cellpose_poly,
        "StarDist" = stardist_poly,
        "mesmer" = mesmer_poly
    ),
    instructions = instructions
)

# further cropping to make segmentation visible on figure
xen_ROI <- subsetGiottoLocs(xen_ROI,
    spat_unit = ":all:",
    x_min = 6650, x_max = 7000,
    y_min = 3900, y_max = 4350
)

for (segmentation in c("cell", "StarDist", "mesmer", "CellPose", "Baysor")) {
    spatInSituPlotPoints(xen_ROI,
        show_image = FALSE,
        feats = list("rna" = c("TACSTD2", "CXCR4", "ITGAX")),
        feats_color_code = c(
            "TACSTD2" = "magenta",
            "CXCR4" = "yellow",
            "ITGAX" = "orange"
        ),
        point_size = 0.8,
        polygon_line_size = 0.5,
        polygon_alpha = 0,
        polygon_feat_type = segmentation,
        polygon_color = "cyan",
        use_overlap = FALSE,
        save_param = list(
            save_name = sprintf("f4a_%s", segmentation),
            save_format = "pdf"
        )
    )
}

############## Data Analysis Steps ########################################

# create full dataset giotto analysis... ----------------------------------- #

# Preprocess
xen_cell <- createGiottoObjectSubcellular(
    gpoints = tx$rna,
    gpolygons = cell
) |>
    calculateOverlapRaster() |>
    overlapToMatrix()

objName(Baysor_poly) <- "cell" # change name to be easier to join
xen_baysor <- createGiottoObjectSubcellular(
    gpoints = tx$rna,
    gpolygons = Baysor_poly
) |>
    calculateOverlapRaster() |>
    overlapToMatrix()

# drop transcripts info to lower memory footprint
xen_cell@feat_info <- NULL
xen_baysor@feat_info <- NULL

# Joined process
join_xen <- joinGiottoObjects(
    gobject_list = list(xen_cell, xen_baysor),
    gobject_names = c("cell", "baysor"),
    join_method = "no_change"
)
instructions(join_xen) <- instructions



# process data ---------------------------------------------------------- #
join_xen <- filterGiotto(
    gobject = join_xen,
    expression_threshold = 1,
    feat_det_in_min_cells = 3,
    min_det_feats_per_cell = 5
)

join_xen <- normalizeGiotto(
    gobject = join_xen,
    scalefactor = 5000
)

join_xen <- addStatistics(join_xen)

subset_25_percent <- round(ncol(join_xen) * 0.25)

join_xen <- runPCAprojection(join_xen,
    feats_to_use = NULL, # no hvfs
    random_subset = subset_25_percent
)

screePlot(join_xen, dim_reduction_name = "pca.projection")

dims_to_use <- seq_len(15) # for UMAP & nearest network

join_xen <- runUMAPprojection(join_xen,
    dim_reduction_name = "pca.projection",
    name = "umap.projection",
    random_subset = subset_25_percent,
    dimensions_to_use = dims_to_use,
    n_neighbors = 15,
    min_dist = 0.15,
    n_epochs = 200,
    spread = 1
)

# clustering ---------------------------------------------------------- #

join_xen <- createNearestNetwork(join_xen,
    dim_reduction_name = "pca.projection",
    dimensions_to_use = dims_to_use,
    k = 30
)

join_xen <- doLeidenClusterIgraph(join_xen,
    resolution = 0.15,
    n_iterations = 100,
    name = "leiden_clus"
)

dimPlot2D(join_xen,
    dim_reduction_name = "umap.projection",
    point_border_stroke = 0,
    point_size = 0.1,
    cell_color = "leiden_clus",
    save_param = list(
        save_name = "f4_leiden"
    )
)

# cell typing ---------------------------------------------------------- #

# plot general known markers
dimFeatPlot2D(join_xen,
    dim_reduction_name = "umap.projection",
    feats = c(
        "CD3E", "CD8A", "CD4", # T cell
        "LILRA4", "IL3RA", # plasmacytoid Dendritic Cell (pDC)
        "MZB1", "CD79A", # plasmablast
        "MS4A1", # B
        "PDGFRB", "LUM", "POSTN", # CAF (stromal)
        "RAMP2", # endothelial
        "TYROBP", # myeloid
        "CPA3", "TPSAB1", # mast cell
        "MKI67", "KRT8", # epithelial/cancer
        "MYH11", "KRT14" # myoepithelial
    ),
    point_size = 0.05,
    point_border_stroke = 0,
    background_color = "black",
    cow_n_col = 4,
    show_legend = FALSE,
    cell_color_gradient = c("#555555", "red"),
    save_param = list(
        save_name = "f4_markers",
        base_height = 12
    )
)

ann <- c(
    "stromal",
    "DCIS_myoepithelial",
    "myeloid",
    "B_and_DC_cells",
    "plasmablast",
    "invasive_tumor",
    "endothelial",
    "T Cells",
    "mast cells"
)
names(ann) <- seq_along(ann)

join_xen <- annotateGiotto(join_xen,
    annotation_vector = ann,
    cluster_column = "leiden_clus",
    name = "major_celltype"
)

# major_celltype annotations
dimPlot2D(join_xen,
    cell_color = "major_celltype",
    dim_reduction_name = "umap.projection",
    cell_color_code = getColors("Vivid", 9),
    point_border_stroke = 0,
    point_size = 0.3,
    show_center_label = FALSE,
    save_param = list(
        save_name = "f4_major_celltype"
    )
)

# subclustering ----------------------------------------------------------- #

subclus_process <- function(x, dims_to_use = 1:10) {
    x <- filterGiotto(x,
        expression_threshold = 1,
        feat_det_in_min_cells = 1,
        min_det_feats_per_cell = 1
    )
    x <- normalizeGiotto(x)
    x <- runPCA(x, feats_to_use = NULL)
    x <- runUMAP(x, dimensions_to_use = dims_to_use)
    x
}

subclus_cluster <- function(x, dims_to_use = 1:10, k = 20, res = 0.1) {
    x |>
        createNearestNetwork(
            dimensions_to_use = dims_to_use,
            k = k
        ) |>
        doLeidenClusterIgraph(resolution = res)
}

# T Cell subclustering

subclus_t <- subset(join_xen, subset = major_celltype == "T Cells")
subclus_t <- subclus_process(subclus_t)
subclus_t <- subclus_cluster(subclus_t, res = 0.08)

plotUMAP(subclus_t, cell_color = "leiden_clus", point_border_stroke = 0)
res_g <- findMarkers_one_vs_all(subclus_t, cluster_column = "leiden_clus")
dimFeatPlot2D(subclus_t,
    point_border_stroke = 0,
    point_size = 0.3,
    gradient_style = "sequential",
    feats = c(
        "CD3", "CD4", "CD8A",
        "IL7R", # memory
        "FOXP3", # regulatory
        "KLRG1", # T effector
        "GATA3", # Th2
        "BCL6", # Tfh
        "PDCD1", # exhaustion
        "GNLY", "KLRD1", # NK
        "CCR7" # naive/central mem
    )
)

ann <- c("CD8T", "CD4T", "NK")
names(ann) <- seq_along(ann)
subclus_t <- annotateGiotto(
    subclus_t,
    annotation_vector = ann,
    cluster_column = "leiden_clus",
    name = "subclus_t"
)

plotUMAP(subclus_t, cell_color = "subclus_t", point_border_stroke = 0)
subclus_t_res <- spatValues(subclus_t, feats = "subclus_t")

join_xen <- addCellMetadata(join_xen, new_metadata = subclus_t_res, by_column = TRUE)

join_xen$major_celltype[join_xen$cell_ID %in% subclus_t_res$cell_ID] <-
    join_xen$subclus_t[join_xen$cell_ID %in% subclus_t_res$cell_ID]

# B/DC subclustering

subclus_bdc <- subset(join_xen, subset = major_celltype == "B_and_DC_cells")
subclus_bdc <- subclus_process(subclus_bdc)
subclus_bdc <- subclus_cluster(subclus_bdc, res = 0.08)

plotUMAP(subclus_bdc, cell_color = "leiden_clus", point_border_stroke = 0)
dimFeatPlot2D(
    subclus_bdc,
    feats = c(
        "MS4A1", # B cell
        "LILRA4", "IL3RA" # plasmacytoid DC (pDC)
    )
)

ann <- c("B", "pDC")
names(ann) <- seq_along(ann)
subclus_bdc <- annotateGiotto(subclus_bdc,
    annotation_vector = ann,
    cluster_column = "leiden_clus", name = "subclus_bdc"
)
subclus_bdc_res <- spatValues(subclus_bdc, feats = "subclus_bdc")
join_xen <- addCellMetadata(join_xen, new_metadata = subclus_bdc_res, by_column = TRUE)
join_xen$major_celltype[join_xen$cell_ID %in% subclus_bdc_res$cell_ID] <-
    join_xen$subclus_bdc[join_xen$cell_ID %in% subclus_bdc_res$cell_ID]

# DCIS/myo subclustering

subclus_dcis <- subset(join_xen, subset = major_celltype == "DCIS_myoepithelial")
subclus_dcis <- subclus_process(subclus_dcis)
subclus_dcis <- subclus_cluster(subclus_dcis, res = 0.08)

plotUMAP(subclus_dcis, cell_color = "leiden_clus", point_border_stroke = 0)
dimFeatPlot2D(
    subclus_dcis,
    feats = c(
        "MYH11", "KRT14", "ACTA2",
        "ERBB2", "ESR1", # DCIS
        "GATA3", "CEACAM6", "SERPINA3", "TCIM", # lum DCIS
        "KRT8", # luminal epithelial DCIS
        "FOXA1", "CCND1" # lum progenitor like
    ),
    point_border_stroke = 0,
    point_size = 0.3,
    show_legend = FALSE
)

ann <- c("DCIS3", "myoepithelial", "DCIS2", "DCIS1")
names(ann) <- seq_along(ann)
subclus_dcis <- annotateGiotto(subclus_dcis,
    annotation_vector = ann,
    cluster_column = "leiden_clus", name = "subclus_dcis"
)

subclus_dcis_res <- spatValues(subclus_dcis, feats = "subclus_dcis")
join_xen <- addCellMetadata(join_xen, new_metadata = subclus_dcis_res, by_column = TRUE)
join_xen$major_celltype[join_xen$cell_ID %in% subclus_dcis_res$cell_ID] <-
    join_xen$subclus_dcis[join_xen$cell_ID %in% subclus_dcis_res$cell_ID]

# stromal subclustering
subclus_stromal <- subset(join_xen, subset = major_celltype == "stromal")
subclus_stromal <- subclus_process(subclus_stromal)
subclus_stromal <- subclus_cluster(subclus_stromal, res = 0.08)

plotUMAP(subclus_stromal, cell_color = "leiden_clus", point_border_stroke = 0)
dimFeatPlot2D(
    subclus_stromal,
    feats = c(
        "LUM", "POSTN", "FAP", "PDGFRB", # CAF
        "MYH11", "MYLK", "ACTG2", # PVL
        "TOP2A", "TUBB2B", "HMGA1" # low quality
    ),
    point_border_stroke = 0,
    point_size = 0.5,
    show_legend = FALSE
)

ann <- c("low_quality", "stromal", "stromal")
names(ann) <- seq_along(ann)

subclus_stromal <- annotateGiotto(subclus_stromal,
    annotation_vector = ann,
    cluster_column = "leiden_clus", name = "subclus_stromal"
)

subclus_stromal_res <- spatValues(subclus_stromal, feats = "subclus_stromal")
join_xen <- addCellMetadata(join_xen,
                            new_metadata = subclus_stromal_res,
                            by_column = TRUE)

join_xen$major_celltype[join_xen$cell_ID %in% subclus_stromal_res$cell_ID] <-
    join_xen$subclus_stromal[join_xen$cell_ID %in% subclus_stromal_res$cell_ID]


# plot final celltypes

# cleanup names
join_xen$major_celltype <- gsub("_", " ", join_xen$major_celltype)

plotUMAP(join_xen,
    dim_reduction_name = "umap.projection",
    cell_color = "major_celltype",
    point_border_stroke = 0,
    point_size = 0.3,
    show_legend = FALSE
)

mycolors <- c(
    "low quality" = "#666666",
    "DCIS3" = "cyan4",
    "DCIS2" = "#8FA242",
    "DCIS1" = "#2F4F30",
    "endothelial" = "#FFCFDF",
    "myoepithelial" = "blue",
    "stromal" = "purple",
    "invasive tumor" = "#AF008D",
    "NK" = "orange",
    "CD4T" = "brown",
    "CD8T" = "red",
    "pDC" = "salmon",
    "mast cells" = "#FF4F09",
    "plasmablast" = "#00F09F",
    "B" = "cyan",
    "myeloid" = "yellow"
)

xen_list <- splitGiotto(join_xen, by = "list_ID")

############## Figure 4b ######################################################

spatInSituPlotPoints(xen_list$cell,
    polygon_fill = "major_celltype",
    polygon_fill_code = mycolors,
    polygon_fill_as_factor = TRUE,
    polygon_line_size = 0,
    polygon_alpha = 1,
    show_legend = FALSE,
    save_param = list(
        save_name = "f4b_spat_original",
        base_height = 14,
        base_width = 20,
        save_format = "png"
    )
)

spatInSituPlotPoints(xen_list$baysor,
    polygon_fill = "major_celltype",
    polygon_fill_code = mycolors,
    polygon_fill_as_factor = TRUE,
    polygon_line_size = 0,
    polygon_alpha = 1,
    show_legend = FALSE,
    save_param = list(
        save_name = "f4b_spat_baysor",
        base_height = 14,
        base_width = 20,
        save_format = "png"
    )
)

########################### ------------- Zoomed In

roi_colors <- list(
    "ERBB2" = "magenta",
    "LUM" = "cyan",
    "CEACAM6" = "yellow"
)

rna_zoom <- crop(tx$rna, ext(3950, 4250, 4650, 4950))

baysor_zoom <- subsetGiottoLocs(xen_list$baysor,
    x_min = 3950, x_max = 4250,
    y_min = 4650, y_max = 4950
)
baysor_zoom <- setGiotto(baysor_zoom, rna_zoom)

spatInSituPlotPoints(baysor_zoom,
    feats = list("rna" = c("ERBB2", "LUM", "CEACAM6")),
    feats_color_code = roi_colors,
    point_size = 0.5,
    background_color = "black",
    use_overlap = FALSE,
    show_polygon = TRUE,
    polygon_color = "grey",
    polygon_line_size = 1,
    polygon_fill_as_factor = TRUE,
    polygon_fill = "major_celltype",
    polygon_alpha = 0.5,
    polygon_fill_code = mycolors,
    plot_last = "points",
    show_legend = FALSE,
    save_param = list(
        save_name = "f4b_baysor_zoom",
        save_format = "pdf"
    )
)

cell_zoom <- subsetGiottoLocs(xen_list$cell,
    x_min = 3950, x_max = 4250,
    y_min = 4650, y_max = 4950
)
cell_zoom <- setGiotto(cell_zoom, rna_zoom)

spatInSituPlotPoints(cell_zoom,
    feats = list("rna" = c("ERBB2", "LUM", "CEACAM6")),
    feats_color_code = roi_colors,
    point_size = 0.5,
    background_color = "black",
    use_overlap = FALSE,
    show_polygon = TRUE,
    polygon_color = "grey",
    polygon_line_size = 1,
    polygon_fill_as_factor = TRUE,
    polygon_fill = "major_celltype",
    polygon_alpha = 0.5,
    polygon_fill_code = mycolors,
    plot_last = "points",
    show_legend = FALSE,
    save_param = list(
        save_name = "f4b_original_zoom",
        save_format = "pdf"
    )
)

############## Figure 4c ######################################################

plotUMAP(xen_list$cell,
    dim_reduction_name = "umap.projection",
    cell_color = "major_celltype",
    cell_color_code = mycolors,
    point_border_stroke = 0,
    point_size = 1,
    show_center_label = FALSE,
    show_legend = FALSE,
    title = "cell_umap",
    save_param = list(
        save_name = "f4c_original_umap",
        base_height = 12,
        base_width = 12,
        save_format = "png"
    )
)

plotUMAP(xen_list$baysor,
    dim_reduction_name = "umap.projection",
    cell_color = "major_celltype",
    cell_color_code = mycolors,
    point_border_stroke = 0,
    point_size = 1,
    show_center_label = FALSE,
    show_legend = FALSE,
    title = "baysor_umap",
    save_param = list(
        save_name = "f4c_baysor_umap",
        base_height = 12,
        base_width = 12,
        save_format = "png"
    )
)

############## Figure 4d ######################################################

# area values were calculated during `addStatistics()` on join_xen in the S8B figure generation
cell_area <- data.frame(
    segmentation = "Original",
    area = xen_list$cell$area
)

baysor_area <- data.frame(
    segmentation = "Baysor",
    area = xen_list$baysor$area
)

poly_area <- rbind(cell_area, baysor_area)
poly_area[, "segmentation"] <- factor(poly_area[, "segmentation"], c("Original", "Baysor"))

violin_area <- ggpubr::ggviolin(poly_area,
    x = "segmentation", y = "area", fill = "segmentation",
    palette = c("#00AFBB", "#FC4E07"),
    xlab = "Segmentation Method",
    ylab = "Polygon Area",
    add = "boxplot",
    add.params = list(size = 0.25)
) +
    ggpubr::stat_compare_means(
        comparisons = list(c("Original", "Baysor")),
        label = "p.signif"
    ) +
    ggpubr::stat_compare_means(label.x = 1.5, label.y = 125)

ggplot2::ggsave(file.path(results_folder, "f4d_area_violin.pdf"), violin_area)

############## Figure 4e ######################################################

# Calculate counts and percentages
cell_meta <- spatValues(join_xen, feats = c("major_celltype", "list_ID"))
cell_meta[list_ID == "cell", list_ID := "original"]

ct_counts <- cell_meta[, .N, by = c("list_ID", "major_celltype")]
ct_counts[, pct := (N / sum(N)) * 100, by = "list_ID"]
cell_type_order <- ct_counts[order(-pct)][, .SD[1], by = major_celltype]$major_celltype
ct_counts[, major_celltype := factor(major_celltype, levels = cell_type_order)]
data.table::setorder(ct_counts, list_ID, -pct)

library(ggplot2)

ct_comp <- ggplot(ct_counts, aes(x = major_celltype, y = pct, fill = list_ID)) +
    geom_bar(
        stat = "identity",
        position = position_dodge()
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45),
        axis.ticks.length.x = unit(1.5, "cm")
    ) +
    labs(
        title = "Cell Types by Segmentation",
        x = "Cell Type",
        y = "Percentage of Segmented Cells"
    )

ggplot2::ggsave(
    filename = file.path(results_folder, "f4e_celltype_comp.pdf"),
    plot = ct_comp,
    dpi = 300,
    units = "in",
    height = 7,
    width = 14
)

############## Figure 4f ######################################################

ct_stack <- ggplot(ct_counts, aes(
    x = factor(list_ID),
    y = pct, fill = major_celltype
)) +
    geom_bar(
        stat = "identity",
        position = position_stack()
    ) +
    theme_classic() +
    scale_fill_manual(values = mycolors) +
    theme(
        axis.text.x = element_text(size = 10),
        axis.ticks.length.x = unit(1.5, "cm")
    ) +
    labs(
        title = "Cell Types by Segmentation",
        x = "Segmentation",
        y = "Percentage of Segmented Cells"
    ) +
    guides(fill = guide_legend(title = "Cell Types"))

ggplot2::ggsave(
    filename = file.path(results_folder, "f4f_celltype_stack.pdf"),
    plot = ct_stack
)
