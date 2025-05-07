##%##########################################################%##
#                                                              #
####          Figure 2 H - M                                   #
#                                                              #
##%##########################################################%##



library(Giotto)

# ensure required packages are present
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler") # v4.12.6
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt") # v2.60.1
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db") # v3.19.1
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
    BiocManager::install("enrichplot") # v1.24.4
}

# set working directory to project
setwd("PATH TO Giotto_Suite_manuscript REPO")
fig_dir <- file.path("scripts", "FIGURE", "F2")
dir.create(fig_dir, showWarnings = FALSE)


# create dataset from vizgen FFPE breast cancer dataset ####

# see scripts/CREATE/create_vizgen_hBC_mini.R

# load giotto object ####
python_path <- NULL
g = loadGiotto(path_to_folder = 'scripts/CREATE/OUTS/vizgen_ffpe_BC_mini/',
               python_path = python_path)
instructions(g, 'save_dir') = fig_dir
instructions(g, 'show_plot') = TRUE
instructions(g, 'return_plot') = FALSE
instructions(g, 'save_plot') = TRUE

set.seed(1234) # for reproducibility

# get points
gpoints = getFeatureInfo(g, return_giottoPoints = TRUE)

# get polygons
cell_poly = getPolygonInfo(g,
    polygon_name = 'cell',
    return_giottoPolygon = TRUE
)
nuc_poly = getPolygonInfo(g,
    polygon_name = 'nucleus',
    return_giottoPolygon = TRUE
)

# check data overlaps
# cell segmentations
plot(g@images$cellbound3_z0, range = c(0, 24e2))
terra::lines(cell_poly@spatVector, col = 'red', alpha = 0.3, lwd = 0.5)

# nuclear manual segmentations
plot(g@images$dapi_z0, range = c(0, 15e3))
terra::lines(nuc_poly@spatVector, col = 'red', alpha = 0.5, lwd = 0.5)


# giottoPoints assigns each transcript an unique ID (feat_ID_uniq)
force(gpoints)


# find feature overlaps ####

# GiottoClass::calculateOverlap() finds which features are overlapped by polygons.
# When no polygons overlap a feature, an NA value is returned

# example with the first 10 polys and first 100 feature points
ex <- GiottoClass::calculateOverlap(
  cell_poly[1:10], gpoints[1:100],
  verbose = FALSE, return_gpoly = FALSE
)
force(ex)
# The overlaps information is returned as the overlapping poly's ID, the feat_ID,
# and the point's unique identifier. Individual overlaps results can be converted
# into expression matrices using GiottoClass::overlapToMatrix(), but here we
# want to compare multiple overlaps results against each other, so we can
# convert the results to data.table and merge the results together by the
# feat_ID_uniq.

# wrapper function around GiottoClass::calculateOverlap() to run for each polygon 
# in a provided list and then merge the results by `feat_ID_uniq`
feat_overlap <- function(x, y) {
  if (!inherits(x, "list")) x <- list(x)
  
  ovlp <- lapply(x, function(seg) {
    res <- calculateOverlap(
      seg, y, return_gpoly = FALSE, verbose = FALSE
    ) |> # return only the overlaps info
      data.table::as.data.table() |>
      data.table::setnames(old = "poly_ID", new = objName(seg))
    return(res[, c("feat_ID_uniq", objName(seg)), with = FALSE])
  })

  if (length(ovlp) > 1L) {
    ovlp <- Reduce(function(x, y) merge(x, y, by = "feat_ID_uniq"), ovlp)
  }
  data.table::setnames(ovlp, old = "feat_ID_uniq", new = objName(y))
}

rel_DT <- feat_overlap(list(cell_poly, nuc_poly), gpoints)
# This produces a data.table of relations between the feature points and the
# two sets of polygons
force(rel_DT)


# determine compartments of interest
in_cell = rel_DT[!is.na(cell), c('rna', 'cell')]
outside_cell = rel_DT[is.na(cell), c('rna', 'cell')]

in_nuc = rel_DT[!is.na(nucleus), c('rna', 'nucleus')]
in_cyto = rel_DT[!is.na(cell) & is.na(nucleus), c('rna', 'nucleus')]

# check overlap results by plotting ####
# inside/outside
plot(gpoints@spatVector[match(in_cell$rna, gpoints$feat_ID_uniq)], background = "black", col = "cyan", alpha = 0.5, cex = 0.2)
plot(gpoints@spatVector[match(outside_cell$rna, gpoints$feat_ID_uniq)], add = TRUE, col = "magenta", alpha = 0.5, cex = 0.2)
terra::lines(cell_poly@spatVector, col = "yellow")
# nuc/cyto
plot(gpoints@spatVector[match(in_nuc$rna, gpoints$feat_ID_uniq)], background = "black", col = "cyan", alpha = 0.5, cex = 0.2)
plot(gpoints@spatVector[match(in_cyto$rna, gpoints$feat_ID_uniq)], add = TRUE, col = "magenta", alpha = 0.5, cex = 0.2)
terra::lines(nuc_poly@spatVector, col = "yellow")
terra::lines(cell_poly@spatVector, col = "white")





# convert overlaps to matrix ####
# Tally up the detections of each gene overlapped by each set of polygons so
# that we can compare the enrichments of features that fall within each
# compartment.
# 
# The rna col from in_cell, outside_cell, in_nuc, and in_cyto was generated
# from the `feat_ID_uniq` attribute in gpoints. So we match those together in
# `tally_overlaps()`
tally_overlaps <- function(gpoints, feat_idx) {
  table(gpoints[match(feat_idx, gpoints$feat_ID_uniq)]$feat_ID) |>
    as.matrix()
}

cell_feats = tally_overlaps(gpoints, in_cell$rna)
outside_feats = tally_overlaps(gpoints, outside_cell$rna)
nuc_feats = tally_overlaps(gpoints, in_nuc$rna)
cyto_feats = tally_overlaps(gpoints, in_cyto$rna)

# combine into expr matrix #
identical(rownames(cell_feats), rownames(outside_feats))
identical(rownames(nuc_feats), rownames(cyto_feats))

inside_outside_expr = cbind(cell_feats, outside_feats)
colnames(inside_outside_expr) = c('cell', 'extracellular')

nuc_cyto_expr = cbind(nuc_feats, cyto_feats)
colnames(nuc_cyto_expr) = c('nucleus', 'cytoplasm')



# return to Giotto for expr normalization ####
# generate exprObj
inside_outside = Giotto::createExprObj(
  expression_data = Matrix::Matrix(inside_outside_expr),
  name = 'raw',
  spat_unit = 'inside_outside',
  feat_type = 'rna',
  provenance = c('cell')
)

nuc_cyto = Giotto::createExprObj(
  expression_data = Matrix::Matrix(nuc_cyto_expr),
  name = 'raw',
  spat_unit = 'nuc_cyto',
  feat_type = 'rna',
  provenance = c('nucleus', 'cell')
)


g = setExpression(g, inside_outside)
g = setExpression(g, nuc_cyto)

# No filtering needed. All feats are present for compartments looked at.
# Perform only library normalization
# This is done by dividing them (per compartment) by the the total expression of
# that compartment, then multiplying by an arbitrary scalefactor (default is 6000)
# 
# Do not perform log normalization
g <- processExpression(g, normParam("library"), spat_unit = "inside_outside")
g <- processExpression(g, normParam("library"), spat_unit = "nuc_cyto")





# gsea analysis  ####


## extract libnorm expression of inside/outside ####

# Values are not log normalized.
viz_gsea_inside_outside = getExpression(g,
  spat_unit = "inside_outside", 
  values = "normalized",
  output = 'matrix'
)

viz_gsea_nuc_cyto = getExpression(g,
  spat_unit = 'nuc_cyto',
  values = 'normalized',
  output = 'matrix'
)

mat_to_dt = function(mat) {
  as.matrix(mat) |>
    data.table::as.data.table()
}

dt_compartment_ranked_lfc = function(expr_mat, a, b) {
  # create DT of dgeMatrix 
  DT = mat_to_dt(expr_mat)
  DT[, 'genes' := rownames(expr_mat)]
  
  # step 2
  # calculate log2 ratios of compartment a vs b values to get log fold change
  DT[, 'lfc_ab' := log2(get(a)/get(b))]
  
  # rank values
  data.table::setorder(DT, -lfc_ab)
}


# function to create log fold change plot
compartment_lfc_plot = function(
    x, a, b, ylab = NULL, xlab = NULL, main = NULL
) {
  op = par()
  on.exit(par(mar = op$mar))
  par(mar = c(5, 6, 4, 0.5) + 0.1)
  plot(x$lfc_ab,
       ylab = ylab,
       xlab = xlab,
       main = main)
  abline(h = 0, col = "red")
}

viz_gsea_inside_outside_dt = dt_compartment_ranked_lfc(
  expr_mat = viz_gsea_inside_outside, a = 'cell', b = 'extracellular'
)
viz_gsea_nuc_cyto_dt = dt_compartment_ranked_lfc(
  expr_mat = viz_gsea_nuc_cyto, a = 'nucleus', b = 'cytoplasm'
)


##### Fig 2 J ---------------------------------------------------------- ####
svg(file.path(fig_dir, "J.svg"), height = 5, width = 4)
compartment_lfc_plot(
    viz_gsea_nuc_cyto_dt,
    a = 'nucleus',
    b = 'cytoplasm',
    ylab = 'Normalized Expression - Log2FC(Nucleus/Cytoplasm)',
    xlab = 'Gene',
    main = 'Log2FC Nuclear/Cytoplasm'
)
dev.off()

##### Fig 2 K ---------------------------------------------------------- ####
svg(file.path(fig_dir, "K.svg"), height = 5, width = 4)
compartment_lfc_plot(
    x = viz_gsea_inside_outside_dt,
    a = "cell",
    b = "extracellular",
    ylab = "Normalized Expression - Log2FC(Inside/Outside)",
    xlab = 'Gene',
    main = "Log2FC Inside/Outside Cell"
)
dev.off()





# get ensembl IDs for the genes we have.
# returns data.table
get_ensembl_ID = function(gene_names) {
  mart = biomaRt::useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
  out = biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "external_gene_name",
                       values = rownames(viz_gsea_inside_outside),
                       mart = mart) |>
    data.table::setDT() |>
  data.table::setnames(new = c('genes', 'ensembl'))
  return(out)
}

eID = get_ensembl_ID(featIDs(gpoints))
# you will note that there are more rows in eID than there are unique features
# this is due to multiple ensembl ID matches.


# BM is the lookup for full name vs ensembl from biomart generated in prev step
# DT is a data.table with gene symbols to add an ensemble ID col to
append_ensembl_id <- function(BM, DT, gene_ID_col = 'genes') {
  # filter for first unique ensembl ID
  BM = BM[!duplicated(genes)]
  DT = merge(DT, BM, by.x = gene_ID_col, by.y = 'genes')
}

get_ranked_ensembl_id_dt <- function(BM, DT, gene_ID_col = 'genes', order_by = 'lfc_ab') {
  res = append_ensembl_id(BM = BM, DT = DT, gene_ID_col = gene_ID_col)
  data.table::setorderv(res, cols = order_by, order = '-1')
  res
}

ranked_inside_outside = get_ranked_ensembl_id_dt(BM = eID, viz_gsea_inside_outside_dt)
ranked_nuc_cyto = get_ranked_ensembl_id_dt(BM = eID, viz_gsea_nuc_cyto_dt)

# no NA values are present from the log2FC
# already sorted for decreasing on inside / outside

# GSEA of gene ontology ####
# Required input is a list of values named by ensembl ID that are ranked according from enriched
# to un-enriched
GO <- function(DT, lfc_col = 'lfc_ab', eID_col = 'ensembl', lfc_cutoff = NULL) {
    
    if (!is.null(lfc_cutoff)) {
        DT <- DT[abs(get(lfc_col)) >= lfc_cutoff]
    }
    
    io_gene_list = DT[[lfc_col]]
    names(io_gene_list) = DT[[eID_col]]
    
    clusterProfiler::gseGO(
        geneList = io_gene_list,
        ont = "CC",
        keyType = "ENSEMBL",
        minGSSize = 3,
        maxGSSize = 500,
        pvalueCutoff = 0.2,
        verbose = TRUE,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        pAdjustMethod = "BH",
        seed = TRUE
    )
}

# these take a while
gsea_inside_outside = GO(ranked_inside_outside)
gsea_nuc_cyto = GO(ranked_nuc_cyto)




# plot GSEA ####

##### Fig 2 L ---------------------------------------------------------- ####
eplot_l <- enrichplot::dotplot(gsea_nuc_cyto, showCategory = 6, split = ".sign") + 
    enrichplot::facet_grid(.~.sign)
ggplot2::ggsave(
    plot = eplot_l,
    filename = file.path(fig_dir, "L.svg"),
    width = 8, 
    height = 8
)

##### Fig 2 M ---------------------------------------------------------- ####

eplot_m <- enrichplot::dotplot(gsea_inside_outside, showCategory = 6, split=".sign") + 
  enrichplot::facet_grid(.~.sign)
ggplot2::ggsave(
    plot = eplot_m,
    filename = file.path(fig_dir, "M.svg"),
    width = 8, 
    height = 8
)




# Plot enriched GO results ####

GO_genes <- function(GO, desc) {
  res <- data.table::setDT(GO@result)
  GO_ID <- res[match(desc, Description), ID]
  
  GO@geneSets[[GO_ID]]
}

in_genes <- GO_genes(gsea_inside_outside, "chromosome, centromeric region") # inside top
out_genes <- GO_genes(gsea_inside_outside, "collagen-containing extracellular matrix") # outside top

nuc_genes <- GO_genes(gsea_nuc_cyto, "nuclear lumen") # nuc top
cyto_genes <- GO_genes(gsea_nuc_cyto, "COPII-coated ER to Golgi transport vesicle") # cyto top



# find symbol names
my_in_genes = ranked_inside_outside[ensembl %in% in_genes, genes]
my_out_genes = ranked_inside_outside[ensembl %in% out_genes, genes]

my_nuc_genes = ranked_nuc_cyto[ensembl %in% nuc_genes, genes]
my_cyto_genes = ranked_nuc_cyto[ensembl %in% cyto_genes, genes]



##### Fig 2H ---------------------------------------------------------- ####
svg(file.path(fig_dir, "H.svg"), height = 7, width = 9)
plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints, feats = my_nuc_genes, col = 'magenta', add = TRUE, alpha = 0.4, raster = FALSE)
plot(gpoints, feats = my_cyto_genes, col = 'cyan', add = TRUE, alpha = 0.5, raster = FALSE)
terra::lines(nuc_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.5)
dev.off()

##### Fig 2I ---------------------------------------------------------- ####
svg(file.path(fig_dir, "I.svg"))
plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints, feats = my_in_genes, col = 'magenta', add = TRUE, alpha = 1, raster = FALSE)
plot(gpoints, feats = my_out_genes, col = 'cyan', add = TRUE, alpha = 0.4, raster = FALSE)
terra::lines(cell_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.7)
dev.off()






# plot the above independently

in_cell_bool <- (gpoints$feat_ID_uniq %in% in_cell$rna) & (gpoints$feat_ID %in% my_in_genes)
out_cell_bool <- (gpoints$feat_ID_uniq %in% in_cell$rna) & (gpoints$feat_ID %in% my_out_genes)

nuc_cell_bool <- (gpoints$feat_ID_uniq %in% in_nuc$rna) & (gpoints$feat_ID %in% my_nuc_genes)
cyto_cell_bool <- (gpoints$feat_ID_uniq %in% in_cyto$rna) & (gpoints$feat_ID %in% my_cyto_genes)


# cell vs extracellular
plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints[in_cell_bool], feats = my_in_genes, col = 'magenta', add = TRUE, alpha = 1, raster = FALSE, cex = 0.15)
plot(gpoints[!in_cell_bool], feats = my_in_genes, col = 'cyan', add = TRUE, alpha = 1, raster = FALSE, cex = 0.15)
terra::lines(cell_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.4)

plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints[!out_cell_bool], feats = my_out_genes, col = 'magenta', add = TRUE, alpha = 1, raster = FALSE)
plot(gpoints[out_cell_bool], feats = my_out_genes, col = 'cyan', add = TRUE, alpha = 1, raster = FALSE)
terra::lines(cell_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.4)


# nuclear vs cyto
plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints[nuc_cell_bool], feats = my_nuc_genes, col = 'magenta', add = TRUE, alpha = 1, raster = FALSE)
plot(gpoints[!nuc_cell_bool], feats = my_nuc_genes, col = 'cyan', add = TRUE, alpha = 1, raster = FALSE)
terra::lines(nuc_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.5)


plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints[cyto_cell_bool], feats = my_cyto_genes, col = 'magenta', add = TRUE, alpha = 1, raster = FALSE)
plot(gpoints[!cyto_cell_bool], feats = my_cyto_genes, col = 'cyan', add = TRUE, alpha = 1, raster = FALSE)
terra::lines(nuc_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.5)




