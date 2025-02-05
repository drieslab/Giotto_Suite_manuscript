## %######################################################%##
#                                                           #
####      Supplementary figure 7 Xenium Co-Register      ####
#                                                           #
## %######################################################%##

## Download data
## Original Xenium data(including H&E,IF) can be downloaded from 10X: 
## https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## In situ sample1, replicate 1
## Original Visium data can be downloaded from same link,
## This example script download orginal Xenium data to './Xenium/'
## This example script download orginal Visium data to './Visium/'

############################# Preprocess and Load ##############################
## Assign the data directory
Xenium_dir <- 'Xenium/'
Visium_dir <- 'Visium/'
results_folder <- "scripts/FIGURE/S7"

library(Giotto)

####### Registering Adjacent Visium to coordinate system(Xenium) #############

instrs <- instructions(
    save_dir = results_folder,
    save_plot = TRUE,
    return_plot = FALSE,
    show_plot = FALSE
)

# load visium data
vis <- createGiottoVisiumObject(Visium_dir,
    gene_column_index = 2,
    png_name = 'tissue_hires_image.png',
    instructions = instrs
)
# extract visium HE image
img_vis <- vis[["images"]][[1]]

# load xenium run data
xen <- createGiottoXeniumObject(
    file.path(Xenium_dir, "outs"),
    feat_type = c("rna", "UnassignedCodeword", "NegControlCodeword", "NegControlProbe"),
    split_keyword = list(
        c("BLANK"),
        c("NegControlCodeword"),
        c("NegControlProbe", "antisense")
    ),
    instructions = instrs
)

# load post-xenium IF images
# * The Xenium also has a paired H&E, but the paired alignment data is not great.
#
# ometif to tif conversion
img_xen_path <- file.path(Xenium_dir, "Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.ome.tif")
for (page in 1:3) {
    GiottoClass::ometif_to_tif(img_xen_path, page = page)
}

img_xen <- lapply(
    sprintf("Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image_000%d.tif", seq_len(3)),
    function(fbasename) {
        read10xAffineImage(
            file = file.path(Xenium_dir, "tif_exports", fbasename),
            imagealignment_path = file.path(Xenium_dir,"Xenium_FFPE_Human_Breast_Cancer_Rep1_if_imagealignment.csv"),
            micron = 0.2125
        )
    }
)
objName(img_xen) <- c("CD20", "HER2", "DAPI") # set image names
# append the post-xenium images
xen <- setGiotto(xen, img_xen)

# set up affine transform to flip whole xenium object to older data orientation
aff <- affine()
ext(aff) <- ext(xen)
aff <- flip(aff)
# perform flip of whole xenium object
xen <- affine(xen, aff)

# extract xenium dapi image
# dapi (lowercase) created during xenium run and was used in segmentation
img_xen <- xen[[,"dapi"]][[1]]

# increase image sampling to help with manual alignment
options("giotto.plot_img_max_sample" = 5e6)

# Select landmarks and calculate affine transform

landmarks <- interactiveLandmarkSelection(img_vis, img_xen) # you can skip this and load the premade
options("giotto.plot_img_max_sample" = 5e5) # revert to default sampling
# landmarks <- readRDS("scripts/Xenium_supporting_files/landmarks.RDS") # premade landmarks
affine_mtx <- calculateAffineMatrixFromLandmarks(landmarks[[1]], landmarks[[2]])


# import polys to show how well the registration worked
cell_poly <- xen[[, "cell"]][[1]]

# align the visium data
vis_aligned <- affine(vis, affine_mtx, pre_multiply = TRUE)
aff_img_vis <- vis_aligned[["images"]][[1]]

# Figure 27A ###############################################################
plot(aff_img_vis)
plot(cell_poly, add = TRUE, border = "blue", lwd = 0.1)
plot(vis_aligned@spatial_info$cell, add = TRUE, col = "orange", lwd = 0.3)

# Figure 2E ################################################################
plot(aff_img_vis)
plot(xen@images$HER2, add = TRUE, alpha = 0.7, col = getMonochromeColors("cyan"))
plot(cell_poly, add = TRUE, border = "blue", lwd = 0.1)
plot(vis_aligned@spatial_info$cell, add = TRUE, col = "orange", lwd = 0.1)

# Figure S7A ###############################################################

spatInSituPlotPoints(vis,
    polygon_bg_color = "black",
    polygon_line_size = 0,
    show_image = TRUE,
    xlim = c(0, 19510),
    ylim = c(-21570, 0),
    save_param = list(
        save_name = "S7A_visium"
    )
)

spatInSituPlotPoints(xen,
    polygon_alpha = 0,
    polygon_color = "blue",
    polygon_line_size = 0.1,
    show_image = TRUE,
    image_name = "dapi",
    save_param = list(
        save_name = "S7A_xenium"
    )
)

# Figure 2F #############################################################
options("giotto.plot_img_max_sample" = 5e6)
spatInSituPlotPoints(xen,
    show_image = TRUE,
    spat_unit = "cell",
    polygon_line_size = 0.05,
    polygon_color = "cyan",
    polygon_alpha = 0,
    image_name = "HER2",
    feats = list(rna = "ERBB2"),
    feats_color_code = "red",
    feat_shape_code = ".",
    xlim = c(-300, 7822),
    ylim = c(-100, 5580),
    point_size = 0.01,
    use_overlap = FALSE,
    show_legend = FALSE,
    save_param = list(
        save_name = "2F_rna_prot_cell",
        base_height = 18,
        base_width = 26
    )
)

spatInSituPlotPoints(xen,
    show_image = TRUE,
    spat_unit = "cell",
    polygon_line_size = 0.7,
    polygon_color = "cyan",
    polygon_alpha = 0,
    image_name = "HER2",
    feats = list(rna = "ERBB2"),
    feats_color_code = "red",
    xlim = c(6150, 6350),
    ylim = c(3950, 4150),
    point_size = 2,
    use_overlap = FALSE,
    show_legend = FALSE,
    save_param = list(
        save_name = "2F_rna_prot_cell_zoom",
        base_height = 10,
        base_width = 10
    )
)
options("giotto.plot_img_max_sample" = 5e5)

# bin data (feature points and raster intensities) ##########################

# extent to cover
e <- ext(xen, prefer = "points") # extent of the points data

binsizes <- c(8, 16, 32, 64, 128, 256, 512)
bins <- lapply(binsizes, function(size) {
    tessellate(extent = e, shape_size = size, shape = "hexagon", id_prefix = "bin_", name = sprintf("hex_%d", size))
})


# intensities data overlap 
# (not fully streamlined to work in the giotto object yet)
CD20 <- xen[[, "CD20"]][[1]]
HER2 <- xen[[, "HER2"]][[1]]
names(CD20) <- "CD20"
names(HER2) <- "HER2"

bins <- lapply(bins, calculateOverlap, CD20)
bins <- lapply(bins, calculateOverlap, HER2)

CD20_cts <- lapply(bins, overlapToMatrix, type = "intensity", feat_info = "CD20")
HER2_cts <- lapply(bins, overlapToMatrix, type = "intensity", feat_info = "HER2")

protein <- lapply(seq_along(binsizes), function(i) {
    expr_mat <- rbind(CD20_cts[[i]], HER2_cts[[i]])
    createExprObj(expr_mat,
        name = "raw",
        spat_unit = objName(bins)[[i]],
        feat_type = "protein"
    )
})

xen <- setGiotto(xen, bins)
xen <- setGiotto(xen, protein)
# points data overlap
for (bin in objName(bins)) {
    xen <- calculateOverlap(xen, spatial_info = bin, feat_info = "rna")
    xen <- overlapToMatrix(xen, poly_info = bin, feat_info = "rna")
}

activeSpatUnit(xen) <- "hex_8"

# Figure S7B ##############################################################
spatFeatPlot2D(xen,
    # plot_method = "scattermore",
    expression_values = "raw",
    spat_unit = "hex_8", 
    feat_type = "protein",
    feats = "CD20", 
    point_size = 0.5,
    point_border_stroke = 0, 
    gradient_style = "s",
    background_color = "black",
    save_param = list(
        save_name = "S7B_CD20"
    )
)

spatFeatPlot2D(xen, 
    expression_values = "raw",
    spat_unit = "hex_8", 
    feat_type = "protein",
    feats = "HER2", 
    point_size = 0.5,
    point_border_stroke = 0, 
    gradient_style = "s",
    background_color = "black",
    save_param = list(
        save_name = "S7B_HER2"
    )
)


spatFeatPlot2D(xen, 
    expression_values = "raw",
    spat_unit = "hex_8", 
    feat_type = "rna",
    feats = "MS4A1",
    point_size = 0.5,
    point_border_stroke = 0, 
    gradient_style = "s",
    background_color = "black",
    save_param = list(
        save_name = "S7B_MS4A1"
    )
)

spatFeatPlot2D(xen, 
    expression_values = "raw",
    spat_unit = "hex_8", 
    feat_type = "rna",
    feats = "ERBB2", 
    point_size = 0.5,
    point_border_stroke = 0, 
    gradient_style = "s",
    background_color = "black",
    save_param = list(
        save_name = "S7B_ERBB2"
    )
)



# find multiscale correlations between feature types

# Find bin spatial neighbors then create weight matrices
for (bin_i in seq_along(bins)) {
    su <- objName(bins)[[bin_i]]
    xen <- createSpatialNetwork(xen,
        spat_unit = su,
        method = "kNN",
        k = 6,
        # all connections should be at the tessellation shape_size dist at
        # maximum to get hexagonal neighbors. We multiply by 1.1 to increase
        # slightly just in case
        maximum_distance_knn = binsizes[[bin_i]] * 1.1
    )
    
    xen <- createSpatialWeightMatrix(xen,
        spat_unit = su,
        spatial_network_to_use = "kNN_network"
    )
}

r2p_MS4A1 <- lapply(seq_along(bins), function(bin_i) {
    su <- objName(bins)[[bin_i]]
    
    x <- spatValues(xen,
        spat_unit = su,
        feat_type = "rna",
        feats = c("MS4A1"),
        expression_values = "raw"
    )
    y <- spatValues(xen,
        spat_unit = su,
        feat_type = "protein",
        feats = c("CD20"),
        expression_values = "raw"
    )
    
    # scale 0 - 1
    x[, MS4A1 := MS4A1 / max(MS4A1)]
    y[, CD20 := CD20 / max(CD20)]

    wm <- getSpatialNetwork(xen, spat_unit = su)@misc$weight_matrix$spat_weights
    id_order <- colnames(wm)
    x <- x[match(id_order, cell_ID),]
    y <- y[match(id_order, cell_ID),]
    
    callSpdep("moran_bv", x = x$MS4A1, y = y$CD20, listw = wm, nsim = 200)
})

r2p_ERBB2 <- lapply(seq_along(bins), function(bin_i) {
    su <- objName(bins)[[bin_i]]
    
    x <- spatValues(xen,
        spat_unit = su,
        feat_type = "rna",
        feats = c("ERBB2"),
        expression_values = "raw"
    )
    y <- spatValues(xen,
        spat_unit = su,
        feat_type = "protein",
        feats = c("HER2"),
        expression_values = "raw"
    )
    
    # scale 0 - 1
    x[, ERBB2 := ERBB2 / max(ERBB2)]
    y[, HER2 := HER2 / max(HER2)]
    
    wm <- getSpatialNetwork(xen, spat_unit = su)@misc$weight_matrix$spat_weights
    id_order <- colnames(wm)
    x <- x[match(id_order, cell_ID),]
    y <- y[match(id_order, cell_ID),]
    
    callSpdep("moran_bv", x = x$ERBB2, y = y$HER2, listw = wm, nsim = 200)
})

# directly find poly area from first polygon of each bin level
# you can also get polygon area from addStatistics()
bin_area <- sapply(bins, function(b) terra::expanse(b[][1]))

morans_bv_results <- data.table::data.table(
    "binarea" = bin_area,
    "MS4A1-CD20" = sapply(1:7, function(i) r2p_MS4A1[[i]]$t0),
    "ERBB2-HER2" = sapply(1:7, function(i) r2p_ERBB2[[i]]$t0)
)

# Figure S7C ##############################################################
library(ggplot2)

multiscale_mbv <- ggplot(morans_bv_results) +
    geom_point(aes(x = log2(binarea), y = `MS4A1-CD20`), color = "red") +
    geom_line(aes(x = log2(binarea), y = `MS4A1-CD20`), color = "red") +
    geom_point(aes(x = log2(binarea), y = `ERBB2-HER2`), color = "blue") +
    geom_line(aes(x = log2(binarea), y = `ERBB2-HER2`), color = "blue") +
    ylab("Bivariate Moran's I") +
    xlab("log2(bin area) (microns)") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank())

ggsave(multiscale_mbv, filename = file.path(results_folder, "S7C_bivariate_morans.png"),
       width = 5, height = 4)

# max morans scores for each comparison
morans_bv_results[which(`MS4A1-CD20` == max(`MS4A1-CD20`)),]
#    binarea MS4A1-CD20 ERBB2-HER2
#      <num>      <num>      <num>
# 1: 3547.24  0.4841736  0.6749982
morans_bv_results[which(`ERBB2-HER2` == max(`ERBB2-HER2`)),]
#     binarea MS4A1-CD20 ERBB2-HER2
#       <num>      <num>      <num>
# 1: 221.7025  0.3272328  0.7163517




############################# Visium and Pseudo Visium comparison ############################# 

# the visium object was aligned to the xenium object earlier
# we can extract the visium polygons in order to compare the visium adjacent
# slice data against the visium poly-aggregated xenium data.

vis_poly <- vis_aligned[["spatial_info"]][[1]]
plot(vis_poly)

# set a new name so that spat_unit naming does not conflict
objName(vis_poly) <- "visium_spot"

# append to xenium object
xen <- setGiotto(xen, vis_poly)

# aggregate xenium values with visium polys
xen <- calculateOverlap(xen, spatial_info = "visium_spot", feat_info = "rna")
xen <- overlapToMatrix(xen, poly_info = "visium_spot", feat_info = "rna")

# pull out only relevant data to make easier to subset
xen_comp <- sliceGiotto(xen, spat_unit = "visium_spot", feat_type = "rna")
activeSpatUnit(xen) <- "visium_spot"

# only use shared features
shared_feats <- intersect(featIDs(xen_comp), featIDs(vis_aligned))
length(shared_feats)

xen_comp <- xen_comp[shared_feats,]
vis_aligned <- vis_aligned[shared_feats,]

# lenient filtering since we are only comparing expression levels
xen_comp <- filterGiotto(xen_comp,
    expression_threshold = 1,
    feat_det_in_min_cells = 1,
    min_det_feats_per_cell = 1
)
vis_aligned <- filterGiotto(vis_aligned,
    expression_threshold = 1,
    feat_det_in_min_cells = 1,
    min_det_feats_per_cell = 1
)

# ensure same polys are used in both datasets
shared_spots <- intersect(spatIDs(xen_comp), spatIDs(vis_aligned))

xen_comp <- xen_comp[, shared_spots]
vis_aligned <- vis_aligned[, shared_spots]

# extract values
xen_expr <- xen_comp[["expression", "raw"]][[1]]
vis_expr <- vis_aligned[["expression", "raw"]][[1]]

# log transform
xen_expr[] <- log2(xen_expr[] + 1)
vis_expr[] <- log2(vis_expr[] + 1)

objName(xen_expr) <- "log"
objName(vis_expr) <- "log"
xen_comp <- setGiotto(xen_comp, xen_expr)
vis_aligned <- setGiotto(vis_aligned, vis_expr)

# ensure same spot order
xen_expr <- xen_expr[, shared_spots]
vis_expr <- vis_expr[, shared_spots]

# correlate values
cor_dt <- data.table::data.table(genes = shared_feats)
for (i in seq_len(nrow(cor_dt))) {
    gene <- cor_dt$genes[i]
    x <- vis_expr[][gene,]
    y <- xen_expr[][gene,]
    cor_value <- cor(x, y)
    cor_dt$cor_value[i] <- cor_value
    cor_dt$mean_expr_visium[i] <- mean(x)
    cor_dt$mean_expr_xenium[i] <- mean(y)
    cor_dt$perc_expressed_visium[i] <- sum(x >= max(x) * 0.05) / length(x)
    cor_dt$perc_expressed_xenium[i] <- sum(y >= max(y) * 0.05) / length(y)
}

dt_sorted <- cor_dt[order(-cor_dt$cor_value), ]
rownames(dt_sorted) = NULL
dt_sorted$cor_rank = 1:nrow(dt_sorted)
dt_sorted[, perc_logfc := log2(perc_expressed_xenium / perc_expressed_visium)]
dt_sorted[, max_perc_expressed := pmax(perc_expressed_visium, perc_expressed_xenium)]

# median cor value
cor_dt[, median(cor_value)] # 0.41433

# cor vs log2(xen % spots expressed / vis % spots expressed)
cor_perc <- ggplot(dt_sorted, aes(x = cor_value, y = perc_logfc)) +
    xlab("pearson correlation") +
    ylab("log2FC of % spots expressing (Xenium / Visium)") +
    geom_point(aes(color = max_perc_expressed)) +
    geom_hline(yintercept = 0, color = "red") +
    labs(color = "max % expressed")

ggsave(cor_perc, filename = file.path(results_folder, "cor_perc.pdf"),
       width = 8.5, height = 4)

# Figure S7D ##############################################################
rankplot <- ggplot(dt_sorted, aes(x = reorder(genes, -cor_value), y = cor_value)) + 
    geom_bar(stat = "identity", color = 'white', width = 1) +
    theme_void() +
    labs(title = "Correlation Scores Between Registered Xenium and Visium ",
         x = "Gene",
         y = "Correlation Score") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5), 
          axis.text.y = element_text(angle = 0, hjust = 2, size = 8),
          axis.title.y = element_text(angle = 90),
          axis.line = element_line(), 
          axis.title = element_text(), legend.position = "none") 

ggsave(rankplot, filename = file.path(results_folder, "S7D_rankplot.pdf"),
       width = 8.5, height = 4)

# Figure S7E ##############################################################
loessplot <- ggplot(dt_sorted) + 
    geom_point(aes(x = cor_rank, y = mean_expr_xenium, color = "blue"), shape = 5, alpha = 0.5) +
    geom_smooth(aes(x = cor_rank, y = mean_expr_xenium, color = "blue"), 
                method = "loess", se = TRUE) +
    geom_point(aes(x = cor_rank, y = mean_expr_visium, color = "red"), shape = 17, alpha = 0.5) +
    geom_smooth(aes(x = cor_rank, y = mean_expr_visium, color = "red"), 
                method = "loess", se = TRUE) +
    theme_void() +
    labs(title = "Loess Trends for Xenium and Visium Expressions vs. Rank",
         x = "Correlation Rank",
         y = "Average Log Expression") +
    scale_color_identity(name = "", 
                         breaks = c("red", "blue"),
                         labels = c("Visium expression", "Xenium expression"),
                         guide = "legend") +
    theme(axis.text = element_text(), 
          axis.line = element_line(),
          axis.title.y = element_text(angle = 90),
          axis.title = element_text(), legend.position = "top") 

ggsave(loessplot, filename = file.path(results_folder, "S7E_loess.pdf"),
       width = 5, height = 4)


# Figure S7F ################################################################## 

highcorgene = "FASN"
spatFeatPlot2D(gobject = xen_comp,
    expression_values = 'log',
    feats = highcorgene,
    point_size = 3,
    point_shape = 'no_border',
    gradient_style = "sequential",
    background_color = 'black',
    show_legend = TRUE,
    save_param = list(
        base_width = 10,
        save_name = "S7F_FASN_xen"
    )
)
spatFeatPlot2D(gobject = vis_aligned,
    feats = highcorgene,
    expression_values = 'log',
    point_size = 3,
    point_shape = 'no_border',
    gradient_style = "sequential",
    background_color = 'black',
    show_legend = TRUE,
    save_param = list(
        base_width = 10,
        save_name = "S7F_FASN_vis"
    )
)

# Figure S7G ################################################################## 

lowcorgene = "HDC"
spatFeatPlot2D(xen_comp,
    expression_values = 'log',
    feats = lowcorgene,
    point_size = 3,
    point_shape = 'no_border',
    gradient_style = "sequential",
    background_color = 'black',
    show_legend = TRUE,
    save_param = list(
        base_width = 10,
        save_name = "S7G_HDC_xen"
    )
)
spatFeatPlot2D(vis_aligned,
    feats = lowcorgene,
    expression_values = 'log',
    gradient_style = "sequential",
    point_size = 3,
    point_shape = 'no_border',
    background_color = 'black',
    show_legend = TRUE,
    save_param = list(
        base_width = 10,
        save_name = "S7G_HDC_vis"
    )
)


