# Looking at transcript enrichment within the Vizgen FFPE Breast Cancer data release



data_dir = 'DATA/vizgen_ffpe_BC_mini/'

# load data ####
viz_ffpe_brca_sub = loadGiotto(paste0(data_dir))


# get points
gpoints = getFeatureInfo(viz_ffpe_brca_sub, return_giottoPoints = TRUE)

# get polygons
cell_poly = getPolygonInfo(viz_ffpe_brca_sub,
                           polygon_name = 'cell',
                           return_giottoPolygon = TRUE)
nuc_poly = getPolygonInfo(viz_ffpe_brca_sub,
                          polygon_name = 'nucleus',
                          return_giottoPolygon = TRUE)



# check data overlaps
# range param has an issue
plot(viz_ffpe_brca_sub@largeImages$cellbound3_z0, range = c(0, 24e2))
terra::lines(cell_poly@spatVector, col = 'red', alpha = 0.3)

plot(viz_ffpe_brca_sub@largeImages$dapi_z0, range = c(0, 15e3))
terra::lines(nuc_poly@spatVector, col = 'red', alpha = 0.5)




# find feature overlaps ####
# TODO clean this up and use rasterization approach to avoid multiple polygon to single feature point relationships
feat_overlap = function(x, y) {
  if(!inherits(x, 'list')) x = list(x)
  
  overlap_DT = lapply(x, function(sel) {
    out = terra::extract(sel@spatVector, y@spatVector)
    data.table::setDT(out)
    data.table::setnames(out, old = 'poly_ID', new = objName(sel))
    out
  })
  if(length(overlap_DT) > 1L) {
    overlap_DT = Reduce(function(x,y) merge(x, y, by = 'id.y'), overlap_DT)
  }
  data.table::setnames(overlap_DT, old = 'id.y', new = objName(y))
  overlap_DT
}

rel_DT = feat_overlap(list(cell_poly, nuc_poly), gpoints)





# determine compartments of interest
in_cell = rel_DT[!is.na(cell), c('rna', 'cell')]
outside_cell = rel_DT[is.na(cell), c('rna', 'cell')]

in_nuc = rel_DT[!is.na(nucleus), c('rna', 'nucleus')]
in_cyto = rel_DT[!is.na(cell) & is.na(nucleus), c('rna', 'nucleus')]


# convert gpoint overlaps to matrix ####
gpoints_to_count_matrix = function(gpoints, feat_idx) {
  table(gpoints@spatVector[feat_idx]$feat_ID) %>%
    as.matrix()
}

cell_feats = gpoints_to_count_matrix(gpoints, in_cell$rna)
outside_feats = gpoints_to_count_matrix(gpoints, outside_cell$rna)

nuc_feats = gpoints_to_count_matrix(gpoints, in_nuc$rna)
cyto_feats = gpoints_to_count_matrix(gpoints, in_cyto$rna)

# combine into expr matrix #
inside_outside_expr = cbind(cell_feats, outside_feats)
colnames(inside_outside_expr) = c('cell', 'extracellular')

nuc_cyto_expr = cbind(nuc_feats, cyto_feats)
colnames(nuc_cyto_expr) = c('nucleus', 'cytoplasm')
# TODO
# rowname order IS identical for both right now so cbind is possible

# generate exprObj ####
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


# return to Giotto for normalization ####
viz_ffpe_brca_sub = setExpression(viz_ffpe_brca_sub, inside_outside)
viz_ffpe_brca_sub = setExpression(viz_ffpe_brca_sub, nuc_cyto)

# no filtering needed. All feats are present for compartments looked at
viz_ffpe_brca_sub = normalizeGiotto(viz_ffpe_brca_sub, spat_unit = 'inside_outside', log_norm = FALSE)
viz_ffpe_brca_sub = normalizeGiotto(viz_ffpe_brca_sub, spat_unit = 'nuc_cyto', log_norm = FALSE)






# TODO move more of this analysis and especially the plotting step into Giotto

# gsea analysis ####

library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)

## extract expression of inside/outside ####

# values have only been scaled. Not log normalized
viz_gsea_inside_outside = getExpression(viz_ffpe_brca_sub, 
                                        spat_unit = "inside_outside", 
                                        values = "normalized",
                                        output = 'matrix')

viz_gsea_nuc_cyto = getExpression(viz_ffpe_brca_sub,
                                  spat_unit = 'nuc_cyto',
                                  values = 'normalized',
                                  output = 'matrix')

mat_to_dt = function(mat) {
  as.matrix(mat) |>
    data.table::as.data.table()
}

dt_compartment_ranked_lfc = function(expr_mat, a, b) {
  # create DT of dgeMatrix 
  DT = mat_to_dt(expr_mat)
  DT[, 'genes' := rownames(expr_mat)]
  
  # step 2
  # calculate log2 ratios of inside and outside (normalized) values
  DT[, 'lfc_ab' := log2(get(a)/get(b))]
  # DT[, outside_inside := log2(b/a), with = FALSE]
  
  # rank values
  data.table::setorder(DT, -lfc_ab)
}

compartment_lfc_plot = function(x, a, b, ylab = NULL, xlab = NULL, main = NULL) {
  op = par()
  par(mar = c(5, 6, 4, 0.5) + 0.1)
  plot(x$lfc_ab,
       ylab = ylab,
       xlab = xlab,
       main = main)
  abline(h = 0, col = "red")
  par(mar = op$mar)
}

viz_gsea_inside_outside_dt = dt_compartment_ranked_lfc(expr_mat = viz_gsea_inside_outside, a = 'cell', b = 'extracellular')
compartment_lfc_plot(x = viz_gsea_inside_outside_dt,
                     ylab = "Normalized Expression - Log2FC(Inside/Outside)",
                     xlab = 'Gene',
                     main = "Log2FC Inside/Outside Cell")
viz_gsea_nuc_cyto_dt = dt_compartment_ranked_lfc(expr_mat = viz_gsea_nuc_cyto, a = 'nucleus', b = 'cytoplasm')
compartment_lfc_plot(viz_gsea_nuc_cyto_dt,
                     a = 'nucleus',
                     b = 'cytoplasm',
                     ylab = 'Normalized Expression - Log2FC(Nucleus/Cytoplasm)',
                     xlab = 'Gene',
                     main = 'Log2FC Nuclear/Cytoplasm')




# get ensembl IDs for the genes we have.
# returns data.table
get_ensembl_ID = function(gene_names) {
  mart = biomaRt::useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
  out = biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                       filters = "external_gene_name",
                       values = rownames(viz_gsea_inside_outside),
                       mart = mart) |>
    data.table::setDT()
  data.table::setnames(out, new = c('genes', 'ensembl')) 
  out
}

eID = get_ensembl_ID(rownames(viz_gsea_inside_outside))



# BM needs the dictionary for full name vs ensembl from biomart generated in prev
# step
append_ensembl_id = function(BM, DT, gene_ID_col = 'genes') {
  # filter for first unique ensembl ID
  BM = BM[!duplicated(genes)]
  DT = merge(DT, BM, by.x = gene_ID_col, by.y = 'genes')
}

get_ranked_ensembl_id_dt = function(BM, DT, gene_ID_col = 'genes', order_by = 'lfc_ab') {
  res = append_ensembl_id(BM = BM, DT = DT, gene_ID_col = gene_ID_col)
  data.table::setorderv(res, cols = order_by, order = '-1')
  res
}

ranked_inside_outside = get_ranked_ensembl_id_dt(BM = eID, viz_gsea_inside_outside_dt)
ranked_nuc_cyto = get_ranked_ensembl_id_dt(BM = eID, viz_gsea_nuc_cyto_dt)


# no NA values are present from the log2FC
# already sorted for decreasing on inside / outside

# run gene ontology
# Required input is a list of values named by ensembl ID that are ranked according from enriched
# to un-enriched
get_gsea = function(DT, lfc_col = 'lfc_ab', eID_col = 'ensembl') {
  io_gene_list = DT[[lfc_col]]
  names(io_gene_list) = DT[[eID_col]]
  
  clusterProfiler::gseGO(
    geneList = io_gene_list,
    ont = "ALL",
    keyType = "ENSEMBL",
    minGSSize = 3,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    OrgDb = org.Hs.eg.db, 
    pAdjustMethod = "none"
  )
}

gsea_inside_outside = get_gsea(ranked_inside_outside)
gsea_nuc_cyto = get_gsea(ranked_nuc_cyto)




# plot GSEA ####

# library(DOSE)
enrichplot::dotplot(gsea_inside_outside, showCategory = 6, split=".sign") + facet_grid(.~.sign)

enrichplot::dotplot(gsea_nuc_cyto, showCategory = 6, split = ".sign") + facet_grid(.~.sign)

# enrichplot::emapplot(enrichplot::pairwise_termsim(gsea_inside_outside),
#                      cex_label_category = 0.7,
#                      label_format = 1)
# enrichplot::emapplot(enrichplot::pairwise_termsim(gsea_nuc_cyto),
#                      cex_label_category = 0.7,
#                      label_format = 1)


# find major driving genes for GO

gsea_io_res = gsea_inside_outside@result |>
  data.table::setDT() |>
  data.table::setorder(-rank)
gsea_nc_res = gsea_nuc_cyto@result |>
  data.table::setDT() |>
  data.table::setorder(-rank)

# get definitions
GO_defs = gsea_io_res[, c('ONTOLOGY', 'ID', 'Description')]
GO_defs[match(c('ribonucleotide binding', # inside top
                'connective tissue development'), # outside top
              Description)]

GO_defs = gsea_nc_res[, c('ONTOLOGY', 'ID', 'Description')]
GO_defs[match(c('programmed cell death', # nuc top
                'COPII-coated ER to Golgi transport vesicle'), # cyto top
              Description)]

ribo_genes = gsea_inside_outside@geneSets$`GO:0032553`
conn_genes = gsea_inside_outside@geneSets$`GO:0061448`

pcd_genes = gsea_nuc_cyto@geneSets$`GO:0012501`
cop_genes = gsea_nuc_cyto@geneSets$`GO:0030134`

my_ribo_genes = ranked_inside_outside[ensembl %in% ribo_genes]
my_conn_genes = ranked_inside_outside[ensembl %in% conn_genes]

my_pcd_genes = ranked_nuc_cyto[ensembl %in% pcd_genes]
my_cop_genes = ranked_nuc_cyto[ensembl %in% cop_genes]




plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints, feats = my_ribo_genes[, genes], col = 'magenta', raster = FALSE, add = TRUE, alpha = 0.5)
plot(gpoints, feats = my_conn_genes[, genes], col = 'cyan', raster = FALSE, add = TRUE, alpha = 0.5)
terra::lines(cell_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.7)

plot(cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(gpoints, feats = my_pcd_genes[, genes], col = 'magenta', raster = FALSE, add = TRUE, alpha = 0.5)
plot(gpoints, feats = my_cop_genes[, genes], col = 'cyan', raster = FALSE, add = TRUE, alpha = 0.5)
# terra::lines(cell_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.7)
terra::lines(nuc_poly@spatVector, col = 'yellow')


ultra_mini_crop1 = terra::ext(5700, 5800, 3800, 3900)
ultra_mini_crop1 = terra::ext(5700, 5800, 3800, 3900)

um_gpoints = terra::crop(gpoints@spatVector, ultra_mini_crop)
um_gpoints = giottoPoints(spatVector = um_gpoints,
                          unique_ID_cache = unique(um_gpoints$feat_ID))
um_nuc_poly = terra::crop(nuc_poly@spatVector, ultra_mini_crop)
um_cell_poly = terra::crop(cell_poly@spatVector, ultra_mini_crop)

plot(um_cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(um_gpoints, feats = my_ribo_genes[, genes], col = 'magenta', raster = FALSE, add = TRUE, alpha = 0.7, cex = 0.2)
plot(um_gpoints, feats = my_conn_genes[, genes], col = 'cyan', raster = FALSE, add = TRUE, alpha = 0.7, cex = 0.2)
terra::lines(um_cell_poly, col = 'yellow', lwd = 1.4, alpha = 0.7)

plot(um_cell_poly, background = 'black', col = 'darkgreen', alpha = 0.7,)
plot(um_gpoints, feats = my_pcd_genes[, genes], col = 'magenta', raster = FALSE, add = TRUE, alpha = 0.7, cex = 0.2)
plot(um_gpoints, feats = my_cop_genes[, genes], col = 'cyan', raster = FALSE, add = TRUE, alpha = 0.7, cex = 0.2)
# terra::lines(cell_poly@spatVector, col = 'yellow', lwd = 0.8, alpha = 0.7)
terra::lines(um_nuc_poly, col = 'yellow', lwd = 1.4)








