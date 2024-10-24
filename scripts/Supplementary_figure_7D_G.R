## %######################################################%##
#                                                          #
####      Supplementary figure 7 Xenium Co-Register        ####
#                                                          #
## %######################################################%##

############################## Download dataset  ###############################
## Original Xenium data can be downloaded from 10X: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## This example script download orginal Xenium and Visium data to './Image_Registration/Xenium_register/outs_rep1/' and "./data/Visium/"
## Registered intermediate files
## Registered files can be downloaded from Zenodo in order to reproduce the figure: 10.5281/zenodo.11075079
## These files include transcripts_Xenium_to_visium_rep1.csv
## Place the downloaded aligned files into './Image_Registration/xenium_rep1_to_visium/' to use the script

library(Giotto)

############ Process Image Registration and Load data paths  ###################

# Visium data should be raw data
visium_dir <- "./data/Visium/"

# feat meta data could also be find in the Xenium output bundle
feat_meta_path <- "./data/outs_rep1/cell_feature_matrix/features.tsv.gz"

# Note that coordinates for transcripts needs to be aligned to Visium HE images,
# example register script can be found in figS7A_Xenium_to_VisiumHE.ipynb,
# which provides an example output for tx_path:
tx_path <- "./Image_Registration/xenium_rep1_to_visium/transcripts_Xenium_to_visium_rep1.csv"

results_folder <- "./Visium_Xen_Integration"
python_path <- NULL

instructions <- createGiottoInstructions(
    save_dir = results_folder,
    save_plot = TRUE,
    show_plot = FALSE,
    python_path = python_path
)


############################ Visium Obj creation  ##############################

Visium <- createGiottoVisiumObject(
    visium_dir = visium_dir,
    expr_data = "raw",
    png_name = "tissue_hires_image.png",
    gene_column_index = 2,
    instructions = instructions
)

Visium <- normalizeGiotto(gobject = Visium,
                          scalefactor = 6000,
                          verbose = TRUE)

Visium <- addStatistics(gobject = Visium)

pl <- spatPlot2D(
    gobject = Visium,
    cell_color = "nr_feats",
    show_image = TRUE,
    point_alpha = 0.7,
    color_as_factor = FALSE,
    return_plot = TRUE,
    save_plot = FALSE
)

pl

############################## Xen Obj creation  ###############################

# featmeta
feature_dt <- data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) <- c("feat_ID", "feat_name", "feat_type")

# find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types <- names(feature_dt[, table(feat_type)])
feat_types_IDs <- lapply(
    feat_types, function(type) feature_dt[feat_type == type, unique(feat_name)]
)
names(feat_types_IDs) <- feat_types

# feat_tx
tx_dt <- data.table::fread(tx_path)
# Image dimentions/scalefactor = tx coords
tx_dt$Visium_x <- tx_dt$aligned_x / 0.092717074
tx_dt$Visium_y <- -tx_dt$aligned_y / 0.092717074

data.table::setnames(
    x = tx_dt,
    old = c("feature_name", "Visium_x", "Visium_y"),
    new = c("feat_ID", "x", "y")
)

# filter by qv (Phred score)
tx_dt_filtered <- tx_dt[qv >= 30]

# separate detections by feature type
tx_dt_types <- lapply(
    feat_types_IDs, function(types) tx_dt_filtered[feat_ID %in% types]
)

gpoints_list <- lapply(
    tx_dt_types, function(x) createGiottoPoints(x = x)
)

# Create Visium spots polygons
centers <- getSpatialLocations(Visium, output = "data.table")
center_to_center <- sqrt((15274 - 3769)^2 + (16784 - 16873)^2) / 63
radius <- center_to_center * 55 / 200
stamp_dt <- circleVertices(radius = radius, npoints = 100)
visium_poly <- createGiottoPolygonsFromDfr(polyStamp(stamp_dt, centers), calc_centroids = T)


xenium_visium <- createGiottoObjectSubcellular(
    gpoints = list(rna = gpoints_list$`Gene Expression`),
    gpolygons = list(Visium = visium_poly),
    instructions = instructions
)

xenium_visium <- calculateOverlapRaster(xenium_visium,
    spatial_info = "Visium",
    feat_info = "rna"
)
xenium_visium <- overlapToMatrix(xenium_visium,
    poly_info = "Visium",
    feat_info = "rna",
    name = "raw"
)

xenium_visium <- normalizeGiotto(xenium_visium)

xenium_visium <- addStatistics(xenium_visium,
                               spat_unit = "Visium")

dt <- pDataDT(xenium_visium)
dt <- dt[, names(dt) := lapply(.SD, function(col) ifelse(is.na(col), 0, col))]

xenium_visium@cell_metadata$Visium$rna[] <- dt

spatPlot2D(
    gobject = xenium_visium,
    spat_unit = "Visium",
    cell_color = "nr_feats",
    color_as_factor = FALSE,
    point_size = 3,
    point_shape = "no_border",
    background_color = "black",
)

######################## Filter No Xenium present area  ########################

xenium_visium <- filterGiotto(
    gobject = xenium_visium,
    spat_unit = "Visium",
    poly_info = "Visium",
    expression_threshold = 1,
    feat_det_in_min_cells = 1,
    min_det_feats_per_cell = 1
)

Visium <- subsetGiotto(Visium,
                       cell_ids = pDataDT(xenium_visium)$cell_ID)

spatPlot2D(
    gobject = xenium_visium,
    spat_unit = "Visium",
    cell_color = "nr_feats",
    color_as_factor = FALSE,
    point_size = 3,
    point_shape = "no_border",
    background_color = "black",
)

############################## Correlation plot  ###############################

xen_feats <- levels(factor(tx_dt_types$`Gene Expression`$feat_ID))

xen_visium_exprs <- getExpression(xenium_visium,
                                  output = "matrix",
                                  values = "raw")

visium_exprs <- getExpression(Visium,
                              output = "matrix",
                              values = "raw")

visium_exprs <- visium_exprs[rownames(visium_exprs) %in% xen_feats, ]

sum(colnames(xen_visium_exprs) == colnames(visium_exprs))
xen_visium_exprs_mtx <- as.matrix(xen_visium_exprs)
visium_exprs_mtx <- as.matrix(visium_exprs)
intersect_genes <- intersect(rownames(visium_exprs_mtx), rownames(xen_visium_exprs_mtx))
cor_df <- data.frame(genes = intersect_genes)

for (i in 1:length(intersect_genes)) {
    gene <- intersect_genes[i]
    x <- visium_exprs_mtx[gene, ]
    y <- xen_visium_exprs_mtx[gene, ]
    cor_value <- cor(log10(x + 1), log10(y + 1))
    # cor_value <- cor(x,y)
    cor_df$cor_value[i] <- cor_value
    cor_df$mean_exprs_visium[i] <- mean(x)
    cor_df$mean_exprs_xenium[i] <- mean(y)
}

df_sorted <- cor_df[order(-cor_df$cor_value), ]
rownames(df_sorted) <- NULL

library(ggplot2)

df_sorted$mean_exprs_visium_log <- log10(df_sorted$mean_exprs_visium + 1)
df_sorted$mean_exprs_xenium_log <- log10(df_sorted$mean_exprs_xenium + 1)

df_sorted$cor_rank <- 1:nrow(df_sorted)

ggplot(df_sorted, aes(x = reorder(genes, -cor_value), y = cor_value)) +
    geom_bar(stat = "identity", color = "white") +
    theme_void() +
    labs(
        title = "Correlation Scores Between Registered Xenium and Visium ",
        x = "Gene",
        y = "Correlation Score"
    ) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        axis.text.y = element_text(angle = 0, hjust = 2, size = 8),
        axis.title.y = element_text(angle = 90),
        axis.line = element_line(),
        axis.title = element_text(), legend.position = "none"
    )

ggplot(df_sorted) +
    geom_point(aes(x = cor_rank, y = mean_exprs_xenium_log, color = "blue"),
               shape = 5, alpha = 0.5) +
    geom_smooth(aes(x = cor_rank, y = mean_exprs_xenium_log, color = "blue"),
        method = "loess", se = TRUE
    ) +
    geom_point(aes(x = cor_rank, y = mean_exprs_visium_log, color = "red"),
               shape = 17, alpha = 0.5) +
    geom_smooth(aes(x = cor_rank, y = mean_exprs_visium_log, color = "red"),
        method = "loess", se = TRUE
    ) +
    theme_void() +
    labs(
        title = "Loess Trends for Xenium and Visium Expressions vs. Rank",
        x = "Correlation Rank",
        y = "Average Log Expression"
    ) +
    scale_color_identity(
        name = "",
        breaks = c("red", "blue"),
        labels = c("Visium expression", "Xenium expression"),
        guide = "legend"
    ) +
    theme(
        axis.text = element_text(),
        axis.line = element_line(),
        axis.title.y = element_text(angle = 90),
        axis.title = element_text(), legend.position = "top"
    )


################################### Comparison  ################################

# Get log expression for both assay

visium_exprs <- getExpression(Visium,
                              output = "matrix",
                              values = "raw")

Visium <- setExpression(Visium,
                        x = createExprObj(log10(visium_exprs + 1)),
                        name = "log")

xen_visium_exprs <- getExpression(xenium_visium,
                                  output = "matrix",
                                  values = "raw")

xenium_visium <- setExpression(xenium_visium,
                               x = createExprObj(log10(xen_visium_exprs + 1)),
                               spat_unit = "Visium",
                               name = "log")

## Rotate coordinates to match Xenium all in one Visually(Fig 2E)
xdt_raw <- getSpatialLocations(xenium_visium, output = "data.table")
xdt <- xdt_raw
xdt$newx <- -xdt_raw$sdimy
xdt$newy <- -xdt_raw$sdimx
xdt <- xdt[, .(cell_ID, newx, newy)]
data.table::setnames(xdt, old = c("newy", "newx"), new = c("sdimy", "sdimx"))

xenium_visium <- setSpatialLocations(xenium_visium,
                                     createSpatLocsObj(xdt),
                                     spat_unit = "Visium",
                                     name = "rotate")

## Rotate coordinates to match Xenium all in one Visually(Fig 2E)
vdt_raw <- getSpatialLocations(Visium,
                               output = "data.table")

vdt <- vdt_raw
vdt$newx <- -vdt_raw$sdimy
vdt$newy <- -vdt_raw$sdimx
vdt <- vdt[, .(cell_ID, newx, newy)]
data.table::setnames(vdt, old = c("newy", "newx"), new = c("sdimy", "sdimx"))

Visium <- setSpatialLocations(Visium,
                              createSpatLocsObj(vdt),
                              name = "rotate")

gene <- "FASN"

spatFeatPlot2D(
    gobject = xenium_visium,
    show_image = FALSE,
    expression_values = "log",
    spat_unit = "Visium",
    spat_loc_name = "rotate",
    feats = gene,
    point_size = 3,
    point_shape = "no_border",
    cell_color_gradient = c("blue", "white", "red"),
    background_color = "black",
    return_plot = TRUE,
    show_legend = TRUE,
    save_param = list(base_width = 10)
)

spatFeatPlot2D(
    gobject = Visium,
    show_image = FALSE,
    feats = gene,
    expression_values = "log",
    spat_loc_name = "rotate",
    cell_color_gradient = c("blue", "white", "red"),
    point_size = 3,
    point_shape = "no_border",
    background_color = "black",
    return_plot = TRUE,
    show_legend = TRUE,
    save_param = list(base_width = 10)
)

gene <- "HDC"

spatFeatPlot2D(
    gobject = xenium_visium,
    show_image = FALSE,
    expression_values = "log",
    spat_unit = "Visium",
    spat_loc_name = "rotate",
    feats = gene,
    point_size = 3,
    point_shape = "no_border",
    cell_color_gradient = c("blue", "white", "red"),
    background_color = "black",
    return_plot = TRUE,
    show_legend = TRUE,
    save_param = list(base_width = 10)
)

spatFeatPlot2D(
    gobject = Visium,
    show_image = FALSE,
    feats = gene,
    spat_loc_name = "rotate",
    expression_values = "log",
    cell_color_gradient = c("blue", "white", "red"),
    point_size = 3,
    point_shape = "no_border",
    background_color = "black",
    return_plot = TRUE,
    show_legend = TRUE,
    save_param = list(base_width = 10)
)
