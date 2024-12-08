## %######################################################%##
#                                                          #
####      Supplementary figure 7 Xeinium Co-Register        ####
#                                                          #
## %######################################################%##

## Download data
## Original Xenium data(including H&E,IF) can be downloaded from 10X: https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
## In situ sample1, replicate 1
## Original Visium data can be downloaded from same link,
## This example script download orginal Xenium data to './Xenium/'
## This example script download orginal Visium data to './Visium/'
## Download the landmarks.RDS 

############################# Preprocess and Load ##############################
## Assign the data directory
Xenium_dir <- paste0(data_dir,'/Xenium/')
Visium_dir <- paste0(data_dir,'/Visium/')
download_dir <- getwd()
library(Giotto)
library(ggplot2)

############################# Create Xenium object as target coordinate system ############################# 
#use the ggplot of centroids to pin landmarks
xen_cell_df <- read.csv(paste0(Xenium_dir,"/outs/cells.csv.gz"))
xen_cell_pl <- ggplot2::ggplot() + ggplot2::geom_point(data = xen_cell_df, ggplot2::aes(x = x_centroid , y = y_centroid),size = 1e-150,,color = 'orange') + ggplot2::theme_classic()


#Load transcripts
x <- importXenium(paste0(Xenium_dir,'/outs/'))
x$qv <- 20 # default
tx <- x$load_transcripts(split_keyword =  list(
    c("BLANK"),
    c("NegControlCodeword"),
    c("NegControlProbe", "antisense")
))
cell <- x$load_polys()

#Flip vertically
aff <- affine() |> flip(direction = 'vertical')
cell <- affine(cell,aff)
tx_pts <- affine(tx[[1]]$rna,aff)


############################################# Figure 7 B and C ############################################# 
############################# Registering IF to coordinate system(Xenium) ##############################
IF_xen <- read10xAffineImage(file = paste0(Xenium_dir, "/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_image.ome.tif"),
                             imagealignment_path = paste0(Xenium_dir,"/Xenium_FFPE_Human_Breast_Cancer_Rep1_if_imagealignment.csv"),
                             micron = 0.2125)

IF_xen <- IF_xen |> flip(direction = "vertical")
gimg_rast <- IF_xen@funs$realize_magick(size = round(dim(IF_xen)[1] *dim(IF_xen)[2]))
plot(gimg_rast)

### Figure B MS4A1, CD20
tx_vect <- tx_pts['MS4A1'][] 
tx_ext <- ext(gimg_rast)

temp_rast <- terra::rast(ext = tx_ext, nrow = 500, ncol = 700)
tx_rast <- terra::rasterize(tx_vect, temp_rast, fun = sum)
ab_rast <- terra::resample(gimg_rast@raster_object[[1]], temp_rast, method = "bilinear")
combined_rast = c(tx_rast, ab_rast)
names(combined_rast) = c("MS4A1", "CD20")
combined_rast[is.na(combined_rast)] <- 0
pearson_R = terra::layerCor(combined_rast, fun = "pearson")


par(mfrow = c(1,2))
plot(combined_rast$CD20,
     col = GiottoVisuals::getColors(pal = "viridis"),
     main = 'CD20',
     background = "black",
     plg = list(loc = "bottom",
                title = "Intensity",
                title.cex = 2.5,
                inset = -0.1,
                cex = 2,
                at = c(0, 50, 150, 200, 250)))

plot(combined_rast$MS4A1,
     col = GiottoVisuals::getColors(pal = "viridis"),
     main = 'MS4A1',
     background = "black",
     plg = list(loc = "bottom",
                title = "Counts",
                title.cex = 2.5,
                inset = -0.1,
                cex = 2))

mtext(sprintf("MS4A1 and CD20 (r = %s)",
              round(pearson_R$correlation[1,2],2)),
      side = 3,
      line = -3,
      cex = 3,
      outer = TRUE)



subset_plot_ext = ext(1500, 2500, 2200, 3700)


par(mfrow = c(1,2))
plot(combined_rast$CD20,
     ext = subset_plot_ext,
     col = GiottoVisuals::getColors(pal = "viridis"),
     background = "black")

plot(combined_rast$MS4A1,
     ext = subset_plot_ext,
     cex = 1.5,
     col = GiottoVisuals::getColors(pal = "viridis"))



### Figure C, ERBB2, HER2
tx_vect <- tx_pts['ERBB2'][] 
tx_ext <- ext(gimg_rast)
temp_rast <- terra::rast(ext = tx_ext, nrow = 500, ncol = 700)
tx_rast <- terra::rasterize(tx_vect, temp_rast, fun = sum)
ab_rast <- terra::resample(gimg_rast@raster_object[[2]], temp_rast, method = "bilinear")
combined_rast = c(tx_rast, ab_rast)
names(combined_rast) = c("ERBB2", "HER2")
combined_rast[is.na(combined_rast)] <- 0
pearson_R = terra::layerCor(combined_rast, fun = "pearson")

par(mfrow = c(1,2))
plot(combined_rast$HER2,
     col = GiottoVisuals::getColors(pal = "viridis"),
     main = 'HER2',
     background = "black",
     plg = list(loc = "bottom",
                title = "Intensity",
                title.cex = 2.5,
                inset = -0.1,
                cex = 2,
                at = c(0, 50, 150, 200, 250)))

plot(combined_rast$ERBB2,
     col = GiottoVisuals::getColors(pal = "viridis"),
     main = 'ERBB2',
     background = "black",
     plg = list(loc = "bottom",
                title = "Counts",
                title.cex = 2.5,
                inset = -0.1,
                cex = 2))

mtext(sprintf("HER2 and ERBB2 (r = %s)",
              round(pearson_R$correlation[1,2],2)),
      side = 3,
      line = -3,
      cex = 3,
      outer = TRUE)


par(mfrow = c(1,1))








############################# Registering Adjacent Visium to coordinate system(Xenium) ##############################
G_visium <- createGiottoVisiumObject(visium_dir = Visium_dir,
                                     gene_column_index = 2,
                                     png_name = 'tissue_hires_image.png',
                                     instructions = NULL)
# In the meantime, calculate statistics for easier plot showing
G_visium <- normalizeGiotto(G_visium)
G_visium <- addStatistics(G_visium)
V_origin <- spatPlot2D(G_visium,show_image = T,point_size = 0,return_plot = T)

# create affine2d
aff <- affine(diag(c(1,1)))
aff <- aff |> 
    spin(90) |>
    flip(direction = "horizontal")
force(aff)

# Apply the transform to better pin landmarks
V_tansformed <- affine(G_visium,aff)
spatplot_to_register <- spatPlot2D(V_tansformed,show_image = T,point_size = 0,return_plot = T)

# Select landmarks and calculate affine transform
landmarks <- interactiveLandmarkSelection(spatplot_to_register, xen_cell_pl)
#landmarks<- readRDS(paste0(download_dir,'/Visium_to_Xen_Landmarks.rds'))
affine_mtx <- calculateAffineMatrixFromLandmarks(landmarks[[1]],landmarks[[2]])
V_final <- affine(G_visium,affine_mtx %*% aff@affine)
spatplot_final <- spatPlot2D(V_final,show_image = T,point_size = 0,show_plot = F) 

############################# Visium and Pseudo Visium comparison ############################# 
# Create PseudoVisium Polygons from Registered Visium using the k = 1 spatial network
V_final <- createSpatialNetwork(V_final, k = 1)
spat_network <- getSpatialNetwork(V_final,output = 'networkDT')
center_to_center <- min(spat_network$distance)
radius <- center_to_center*55/200
Visium_centroid <- getSpatialLocations(V_final,output = 'data.table')
stamp_dt <- circleVertices(radius = radius, npoints = 100)
pseudo_visium_dt <- polyStamp(stamp_dt, Visium_centroid)
pseudo_visium_poly <- createGiottoPolygonsFromDfr(pseudo_visium_dt,calc_centroids = T,name = 'visium')
overlap_extent <- ext(tx_pts@spatVector)
pseudo_visium_poly <- crop(pseudo_visium_poly,overlap_extent)
plot(pseudo_visium_poly)



# Now, modify Xen object with new spat unit
Xen = createGiottoObjectSubcellular(gpoints = list('rna' = tx_pts),
                                    gpolygons = list('visium' = pseudo_visium_poly))
Xen <- calculateOverlap(Xen,
                            feat_info = 'rna',
                            spatial_info = 'visium')
Xen <- overlapToMatrix(x = Xen,
                           type = "point", 
                           poly_info = "visium", 
                           feat_info = "rna",
                           aggr_function = "sum")

tmp_exprs <- getExpression(Xen,
                           feat_type = 'rna',
                           spat_unit = 'visium',
                           output = 'matrix')
Xen <- setExpression(Xen,
                         x = createExprObj(log(tmp_exprs+1)),
                         feat_type = 'rna',
                         spat_unit = 'visium',
                         name = 'log')

# Now, modify Visium object to get the subsetted polygons
sub_visium <- subsetGiotto(V_final,spat_unit = 'cell',cell_ids = pseudo_visium_poly$poly_ID)
tmp_exprs <- getExpression(sub_visium,
                           feat_type = 'rna',
                           output = 'matrix')
sub_visium <- setExpression(sub_visium,
                     x = createExprObj(log(tmp_exprs+1)),
                     name = 'log')

############################################# Figure 7 D and E ############################################# 
### correlation plot -------
xen_visium_exprs = getExpression(Xen,spat_unit = 'visium',output = 'matrix',values = 'raw')
visium_exprs = getExpression(sub_visium,output = 'matrix',values = 'raw')
intersect_genes = intersect(fDataDT(Xen,spat_unit = 'visium')$feat_ID,rownames(visium_exprs))
xen_visium_exprs = xen_visium_exprs[intersect_genes ,]
visium_exprs = visium_exprs[intersect_genes,]


sum(colnames(xen_visium_exprs) != colnames(visium_exprs))
xen_visium_exprs_mtx = as.matrix(xen_visium_exprs)
visium_exprs_mtx = as.matrix(visium_exprs)
cor_df = data.frame(genes = intersect_genes)

for (i in 1:length(intersect_genes)){
    gene = intersect_genes[i]
    x = visium_exprs_mtx[gene,]
    y = xen_visium_exprs_mtx[gene,]
    cor_value <- cor(log10(x+1), log10(y+1))
    #cor_value <- cor(x,y)
    cor_df$cor_value[i] = cor_value
    cor_df$mean_exprs_visium[i] = mean(x)
    cor_df$mean_exprs_xenium[i] = mean(y)
}

df_sorted <- cor_df[order(-cor_df$cor_value), ]
rownames(df_sorted) = NULL

df_sorted$mean_exprs_visium_log = log10(df_sorted$mean_exprs_visium+1)
df_sorted$mean_exprs_xenium_log = log10(df_sorted$mean_exprs_xenium+1)

df_sorted$cor_rank = 1:nrow(df_sorted)

ggplot(df_sorted, aes(x=reorder(genes, -cor_value), y=cor_value)) + 
    geom_bar(stat="identity",color = 'white')+
    theme_void() +
    labs(title = "Correlation Scores Between Registered Xenium and Visium ",
         x = "Gene",
         y = "Correlation Score") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 5), 
          axis.text.y = element_text(angle = 0, hjust = 2,size = 8),
          axis.title.y = element_text(angle = 90),
          axis.line = element_line(), 
          axis.title = element_text(), legend.position = "none") 

ggplot(df_sorted) + 
    geom_point(aes(x = cor_rank, y = mean_exprs_xenium_log, color = "blue"),shape = 5, alpha = 0.5) +
    geom_smooth(aes(x = cor_rank, y = mean_exprs_xenium_log, color = "blue"), 
                method = "loess", se = T) +
    geom_point(aes(x = cor_rank, y = mean_exprs_visium_log, color = "red"),shape = 17, alpha = 0.5) +
    geom_smooth(aes(x = cor_rank, y = mean_exprs_visium_log, color = "red"), 
                method = "loess", se = T) +
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


############################################# Figure 7 F and G ############################################# 
gene = 'FASN'
p1 <- spatFeatPlot2D(gobject = Xen,
               show_image = F,
               expression_values = 'log',
               spat_unit = 'visium',
               feats = gene,
               point_size = 3,
               point_shape = 'no_border',
               cell_color_gradient = c("blue", "white", "red"),
               background_color = 'black',
               return_plot = T,
               show_legend = TRUE)
p2 <- spatFeatPlot2D(gobject = sub_visium,
               show_image = F,
               feats = gene,
               expression_values = 'log',
               cell_color_gradient = c("blue", 'white', "red"),
               point_size = 3,
               point_shape = 'no_border',
               background_color = 'black',
               return_plot = T,
               show_legend = TRUE)

gene = 'HDC'
p3 <- spatFeatPlot2D(gobject = Xen,
               show_image = F,
               expression_values = 'log',
               spat_unit = 'visium',
               feats = gene,
               point_size = 3,
               point_shape = 'no_border',
               cell_color_gradient = c("blue", "white", "red"),
               background_color = 'black',
               return_plot = T,
               show_legend = TRUE,
               save_param = list(base_width = 10))
p4 <- spatFeatPlot2D(gobject = sub_visium,
               show_image = F,
               feats = gene,
               expression_values = 'log',
               cell_color_gradient = c("blue", 'white', "red"),
               point_size = 3,
               point_shape = 'no_border',
               background_color = 'black',
               return_plot = T,
               show_legend = TRUE,
               save_param = list(base_width = 10))

gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)









