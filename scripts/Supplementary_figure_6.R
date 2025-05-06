## %##########################################################%##
#                                                              #
####                 Supplementary figure 4                    #
#                                                              #
## %##########################################################%##

library(Giotto)

# create dataset from vizgen FFPE breast cancer dataset ####

# see scripts/CREATE/create_vizgen_mini7.R

############################## Load giotto object  #############################

python_path <- NULL

g <- loadGiotto(
    path_to_folder = "scripts/CREATE/OUTS/mini7/",
    python_path = python_path
)

instructions(g, "save_dir") <- "results"
instructions(g, "show_plot") <- TRUE
instructions(g, "return_plot") <- FALSE
instructions(g, "save_plot") <- TRUE

fig_dir <- file.path("scripts", "FIGURE", "S4")
dir.create(fig_dir, showWarnings = FALSE)

z_digits <- 0:6
z_levels <- paste0("z", z_digits)

# calculate get counts overlapped per cell polygon to produce expression matrix
for (z_digit in z_digits) {
    z_level <- paste0("z", z_digit)

    cat("level ", z_level, "\n")

    g <- calculateOverlap(
        g,
        spatial_info = z_level,
        feat_info = "rna",
        # feat_subset params subset the points info used when calculating
        # overlaps to only those relevant to that z layer"s polygon data
        feat_subset_column = "z",
        feat_subset_ids = z_digit
    )

    g <- overlapToMatrix(
        g,
        poly_info = z_level,
        feat_info = "rna",
        name = "raw"
    )
}

g <- aggregateStacks(
    gobject = g,
    spat_units = z_levels,
    feat_type = "rna",
    values = "raw",
    summarize_expression = "sum",
    summarize_locations = "mean",
    new_spat_unit = "aggregate"
)

# polygon differences #

# get polys
poly_list <- g@spatial_info[z_levels]

# set up rasterizations to be 1000 x 1000
r <- terra::rast(nrows = 1e3, ncols = 1e3, extent = ext(poly_list[[1]]))

# perform rasterizations on polygons
poly_rasters <- lapply(poly_list, function(p) {
    terra::rasterize(p@spatVector, r)
})

# add all rasterized polys together to show overlapped regions
poly_rasterstack <- Reduce("c", poly_rasters)
poly_rasterstack[is.na(poly_rasterstack)] <- 0
poly_rs_sum <- terra::app(poly_rasterstack, sum)

vcolors <- viridisLite::turbo(8)
vcoltab <- data.table::data.table(sum = 0:7, color = vcolors)
terra::coltab(poly_rs_sum) <- vcoltab

##### Fig S4 B ---------------------------------------------------------- ####

svg(file.path(fig_dir, "B.svg"), height = 4, width = 4.1)
plot(poly_rs_sum) # plot of multiple layers
dev.off()

## rasterized feature comparisons ##

# decide on the image resolution
featinfo <- getFeatureInfo(g, return_giottoPoints = FALSE)
rasterimage <- terra::rast(x = featinfo, nrows = 20, ncols = 20)
# this rasterizes the feature points into a 20x20 image
# The resolution is arbitrary, but requires that the resolution is high enough
# that it makes some visual sense, while also being low enough that the
# rasterization will aggregate counts

# create z specific feature information layers
z_layers <- sort(unique(featinfo$z))
z_featinfo_list <- lapply(z_layers, FUN = function(x) {
    layerinfo <- featinfo[featinfo$z == x]
})

# rasterize each gene in all z-layers
# Will print index of gene being processed / layer to show progress
genes <- unique(featinfo$feat_ID)
rast_list <- lapply(1:length(z_featinfo_list), FUN = function(i) {
    cat("for layer ", i, "\n")

    z_layer <- z_featinfo_list[[i]]

    final_z_collection <- list()

    for (gene_i in 1:length(genes)) {
        print(gene_i)
        gene <- genes[gene_i]
        genevec <- z_layer[z_layer$feat_ID == gene]
        generast <- terra::rasterize(x = genevec, y = rasterimage, fun = "sum", na.rm = TRUE)
        final_z_collection[[gene_i]] <- generast
    }

    final_z_collection
})

# combine the genes together
new_savelist <- lapply(1:length(genes), FUN = function(x) {
    print(x)

    # aggregate each gene across all rasterized layers
    comb_generast <- NULL
    for (i_layer in 1:length(rast_list)) {
        generast <- rast_list[[i_layer]][[x]]
        comb_generast <- c(comb_generast, generast)
    }

    comb_generast <- do.call("c", comb_generast)
})


# new_savelist
# each gene is converted into a spatraster with defined dimensions (rastermimage)
# each spatraster has 7 layers (each z-stack)

# example:
one_gene <- new_savelist[[2]]
plot(one_gene)


# Function to plot the rasterized feature info for specific genes
# with individual plots for each z layer.
plot_layers_gene <- function(rasterized_feats,
    feats_order = NULL,
    feat = "Gfap",
    pal = "Sunset",
    n = 100,
    rev = F) {
    feat_index <- which(feats_order == feat)
    spatRasterGene <- rasterized_feats[[feat_index]]

    min_value <- min(unlist(lapply(1:terra::nlyr(spatRasterGene), FUN = function(x) {
        terra::minmax(spatRasterGene[[x]])[[1]]
    })))

    max_value <- max(unlist(lapply(1:terra::nlyr(spatRasterGene), FUN = function(x) {
        terra::minmax(spatRasterGene[[x]])[[2]]
    })))

    mycolors <- GiottoVisuals::getColors(pal = pal, n = n, rev = rev)
    plot(spatRasterGene, range = c(min_value, max_value), col = mycolors, main = feat)
}

# examples
plot_layers_gene(new_savelist, genes, "Gfap")
plot_layers_gene(new_savelist, genes, "Htr1a")
plot_layers_gene(new_savelist, genes, "Cxcl12")

## example for all genes

# calculate the total sum
sum_spat_genes <- lapply(X = new_savelist, function(x) {
    sum_gene <- sum(x, na.rm = TRUE)
})
sum_spat_genes_vector <- do.call("c", sum_spat_genes)
sum_all <- sum(sum_spat_genes_vector, na.rm = TRUE)

mycolors <- hcl.colors(n = 200, palette = "Sunset", rev = TRUE)

# set missing values to 0
sum_all[is.na(sum_all)] <- 0
sum_range <- terra::minmax(sum_all)

cell_polys <- g@spatial_info[paste0("z", 0:6)]
cell_centroids <- lapply(cell_polys, GiottoClass::centroids)


##### Fig S4 C ---------------------------------------------------------- ####
# plot the total sum across all layers
svg(file.path(fig_dir, "C.svg"), height = 9, width = 9)
plot(sum_all, col = mycolors, range = sum_range, main = "all layers")
plot(g@spatial_info$aggregate, add = TRUE)
plot(centroids(g@spatial_info$aggregate), add = TRUE, cex = 0.5)
dev.off()

##### Fig S4 D ---------------------------------------------------------- ####

# plot the total sum for each layer (figure S4 D continued)
vals <- list()
for (i in 1:7) {
    sum_spat_genes_layer <- lapply(X = new_savelist, function(x) {
        selected_layer <- terra::subset(x = x, i)
    })
    sum_spat_genes_vector_layer <- do.call("c", sum_spat_genes_layer)
    sum_spat_genes_vector_layer[is.na(sum_spat_genes_vector_layer)] <- 0
    sum_all_layer <- sum(sum_spat_genes_vector_layer)
    vals[[i]] <- sum_all_layer |>
        terra::values() |>
        as.numeric()

    # plot
    svg(file.path(fig_dir, sprintf("D%d.svg", i)), height = 9, width = 9)
    plot(sum_all_layer, col = mycolors, main = paste0("layer = ", i))
    plot(cell_polys[[i]], add = T)
    plot(cell_centroids[[i]], add = T, cex = 0.5)
    dev.off()
}


##### Fig S4 E ---------------------------------------------------------- ####

# plot as boxplot
vals_DF <- Reduce(cbind, vals)
colnames(vals_DF) <- z_levels

svg(file.path(fig_dir, "E.svg"), height = 9, width = 9)
boxplot(vals_DF)
dev.off()

# COV
cov_spat_genes <- lapply(X = new_savelist, function(x) {
    std_gene <- terra::stdev(x, na.rm = T)
    mean_gene <- terra::mean(x, na.rm = T)
    cov_gene <- std_gene / mean_gene
    cov_gene
})

cov_spat_genes_vector <- do.call("c", cov_spat_genes)
cov_spat_genes_vector[is.na(cov_spat_genes_vector)] <- 0
mean_cov_all <- terra::mean(cov_spat_genes_vector, na.rm = TRUE)

mycolors <- hcl.colors(10, palette = "Temps")

# plot mean cov
plot(mean_cov_all, col = mycolors, main = "mean COV")
plot(g@spatial_info$aggregate, add = TRUE)
plot(g@spatial_info$aggregate, type = "centroid", add = TRUE, cex = 0.5)

sum_cov_all <- sum(cov_spat_genes_vector, na.rm = TRUE)

##### Fig S4 F ---------------------------------------------------------- ####

# plot total cov
svg(file.path(fig_dir, "F.svg"), height = 9, width = 9)
plot(sum_cov_all, col = mycolors, main = "total COV")
plot(g@spatial_info$aggregate, add = TRUE)
plot(g@spatial_info$aggregate,
    type = "centroid",
    add = TRUE, cex = 0.5
)
dev.off()

spatial_diff <- data.table::data.table(genes = genes, sum = unlist(sum_spat_genes), cov = unlist(cov_spat_genes))
spatial_diff[, sumt := sum(terra::values(sum[[1]], na.rm = T)), by = 1:nrow(spatial_diff)]
spatial_diff[, covt := sum(terra::values(cov[[1]], na.rm = T)), by = 1:nrow(spatial_diff)]

data.table::setorder(spatial_diff, -covt)
spatial_diff[1:50]

spatial_diff[sumt > 1500]

##### Fig S4 G ---------------------------------------------------------- ####

library(ggplot2)
pl <- ggplot()
pl <- pl + geom_point(data = spatial_diff, aes(x = log(sumt + 1), y = log(covt + 1)))
pl <- pl + theme_bw() + xlab("log(counts +1)") + ylab("log(COV+1)")
pl

# with Gad1 label
pl <- pl + geom_point(
    data = spatial_diff[genes == "Gad1"],
    aes(
        x = log(sumt + 1), y = log(covt + 1)
    ),
    color = "red"
) +
    geom_text(
        data = spatial_diff[genes == "Gad1"],
        aes(
            x = log(sumt + 1), y = log(covt + 1)
        ),
        color = "red",
        label = "Gad1",
        hjust = -0.2
    )
pl
ggsave(file.path(fig_dir, "G.svg"), plot = pl, height = 9, width = 9)

plot_layers_gene(new_savelist, genes, "Cspg5")

##### Fig S4 H ---------------------------------------------------------- ####

# example of gene with differences in distribution across the z slices
svg(file.path(fig_dir, "H.svg"), height = 9, width = 9)
plot_layers_gene(new_savelist, genes, "Gad1")
dev.off()
