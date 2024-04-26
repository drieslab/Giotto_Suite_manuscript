# ----- Dataset -----
# This scirpt take the Xenium Human Breast Cancer Data as an example.
# The download link is https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast

# ----- Creat Python Environment -----
# #run in terminal
# conda create -n giotto python=3.10.2 r=4.2 -y
# conda activate giotto
# pip install pandas==1.5.1 python-igraph==0.10.2 networkx==2.8.8 # python-louvain==0.16 leidenalg==0.9.0 scikit-learn==1.1.3 smfishHmrf
# pip install git+https://github.com/wwang-chcn/bento-tools.git@giotto_install


# ----- Processing 10x Xenium Human Breast Cancer Data -----
library(Giotto)
library(reticulate)

results_folder = paste0(getwd(),'/Xenium_results')
python_path = '.conda/envs/giotto/bin/python'  # change .conda to your conda path, which can be found by running 'conda env list' in terminal

# 3. Create Giotto instructions
# Directly saving plots to the working directory without rendering them in the editor saves time.
my_instrs = createGiottoInstructions(python_path = python_path,
                                     save_dir = results_folder,
                                     save_plot = TRUE,
                                     show_plot = TRUE,
                                     return_plot = FALSE)

# ** SET PATH TO FOLDER CONTAINING XENIUM DATA **
xenium_folder = paste0(getwd(),'/Xenium/')
list.files(xenium_folder)
# general files (some are supplemental files)
settings_path = paste0(xenium_folder, 'experiment.xenium')
he_img_path = paste0(xenium_folder, 'Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif')
if_img_path = paste0(xenium_folder, 'Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.tif')
panel_meta_path = paste0(xenium_folder, 'Xenium_FFPE_Human_Breast_Cancer_Rep1_panel.tsv') # (optional)

xenium_out_folder = paste0(xenium_folder,'outs/')
# files (SUBCELLULAR): (tutorial focuses on working with these files)
cell_bound_path = paste0(xenium_out_folder, 'cell_boundaries.csv.gz')
nuc_bound_path = paste0(xenium_out_folder, 'nucleus_boundaries.csv.gz')
tx_path = paste0(xenium_out_folder, 'transcripts.csv.gz')
feat_meta_path = paste0(xenium_out_folder, 'cell_feature_matrix/features.tsv.gz') # (also used in aggregate)

# files (AGGREGATE):
expr_mat_path = paste0(xenium_out_folder, 'cell_feature_matrix')
cell_meta_path = paste0(xenium_out_folder, 'cells.csv.gz') # contains spatlocs

# load features metadata
# (make sure cell_feature_matrix folder is unpacked)
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('feat_ID','feat_name','feat_type')

# find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])

feat_types_IDs = lapply(
  feat_types, function(type) feature_dt[feat_type == type, unique(feat_name)] # change feat_ID to feat_name
)
names(feat_types_IDs) = feat_types

tx_dt = data.table::fread(tx_path)
data.table::setnames(x = tx_dt,
                     old = c('feature_name', 'x_location', 'y_location'),
                     new = c('feat_ID', 'x', 'y'))
cat('Transcripts info available:\n ', paste0('"', colnames(tx_dt), '"'), '\n',
    'with', tx_dt[,.N], 'unfiltered detections\n')

# filter by qv (Phred score)
tx_dt_filtered = tx_dt[qv >= 20]
cat('and', tx_dt_filtered[,.N], 'filtered detections\n\n')

# separate detections by feature type
tx_dt_types = lapply(
  feat_types_IDs, function(types) tx_dt_filtered[feat_ID %in% types]
)

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], 'detections: ', tx_dt_types[[x]][,.N], '\n')
}))

gpoints_list = lapply(
  tx_dt_types, function(x) createGiottoPoints(x = x)
) # 208.499 sec elapsed

# preview QC probe detections
plot(gpoints_list$`Blank Codeword`,
     point_size = 0.3,
     main = 'Blank Codeword')
plot(gpoints_list$`Negative Control Codeword`,
     point_size = 0.3,
     main = 'Negative Control Codeword')
plot(gpoints_list$`Negative Control Probe`,
     point_size = 0.3,
     main = 'Negative Control Probe')

# preview two genes (slower)
plot(gpoints_list$`Gene Expression`,  # 77.843 sec elapsed
     feats = c('KRT8', 'MS4A1'))
tx_dt_types$`Gene Expression`[feat_ID %in% c('KRT8', 'MS4A1'), table(feat_ID)]


cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(cellPoly_dt,
                     old = c('cell_id', 'vertex_x', 'vertex_y'),
                     new = c('poly_ID', 'x', 'y'))
data.table::setnames(nucPoly_dt,
                     old = c('cell_id', 'vertex_x', 'vertex_y'),
                     new = c('poly_ID', 'x', 'y'))

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)

plot(x = gpoly_nucs, point_size = 0.1, type = 'centroid')

xenium_gobj = createGiottoObjectSubcellular(
  gpoints = list(rna = gpoints_list$`Gene Expression`,
                 blank_code = gpoints_list$`Blank Codeword`,
                 neg_code = gpoints_list$`Negative Control Codeword`,
                 neg_probe = gpoints_list$`Negative Control Probe`),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  instructions = my_instrs
)


# ----- subset gobj -----
subset_xenium_gobj <- subsetGiottoLocs(xenium_gobj, spat_unit='cell', feat_type='rna',
                                       x_max=200,x_min=0,y_max=200,y_min=0)


png("Supple14_D_region.png", width = 3000, height = 3000, res=600)
plot(x = subset_xenium_gobj@feat_info$rna@spatVector, point_size = 0.1, type = 'centroid', xlim=c(0,200), ylim=c(0,200), col='gray', pch=20, cex=0.2)
plot(x = subset_xenium_gobj@spatial_info$cell@spatVector, point_size = 0.1, type = 'centroid', xlim=c(0,200), ylim=c(0,200), add=TRUE)
plot(x = subset_xenium_gobj@spatial_info$nucleus@spatVector, point_size = 0.1, type = 'centroid', xlim=c(0,200), ylim=c(0,200), add=TRUE)
dev.off()

# ----- create bento -----
bento_adata <- createBentoAdata(gobject = subset_xenium_gobj, env_to_use='/sc/arion/work/wangw32/conda-env/envs/giotto')


# ----- analysis -----
bento_analysis_path <- system.file("python","python_bento_analysis.py",package="Giotto")
reticulate::source_python(bento_analysis_path)


analysis_shape_features(adata=bento_adata)
plot_shape_features_analysis_results(adata=bento_adata, fname='test_shape_features.pdf')

analysis_points_features(adata=bento_adata)
plot_points_features_analysis_results(adata=bento_adata, fname='test_points_features.pdf')

analysis_rna_forest(adata=bento_adata)
plot_rna_forest_analysis_results(adata=bento_adata, fname1='Supple14_D_radvis.pdf', fname2='Supple14_E.pdf')

analysis_colocalization(adata=bento_adata, fname='test_colocalization_knee_pos.pdf', ranks=seq(10))
plot_colocalization_analysis_results(adata=bento_adata, rank=5, fname='Supple14_F.pdf')
