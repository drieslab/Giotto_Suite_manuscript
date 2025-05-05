##%######################################################%##
#                                                          #
####      Supplementary figure 8 Interactive tools      ####
#                                                          #
##%######################################################%##

####################### Supplementary figure 8a ################################
# Pre-run the Parts 1-4 from the [Mouse Visium Brain tutorial](https://drieslab.github.io/Giotto_website/articles/visium_mouse_brain.html)

library(Giotto)

# Plot spatial cells/spots without tissue on the background
my_spatPlot <- spatPlot2D(gobject = visium_brain,
                          point_size = 1.5,
                          point_alpha = 0.1,
                          show_image = FALSE)


# Run the Shiny app and save the spots coordinates
my_polygon_coordinates <- plotInteractivePolygons(my_spatPlot,
                                                  color = 'black')


####################### Supplementary figure 8b ################################

## Create a Giotto polygon object
my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coordinates,
                                                  name = 'selections')

## Add the polygons to the Giotto object
visium_brain <- addGiottoPolygons(gobject = visium_brain,
                                  gpolygons = list(my_giotto_polygons))

## Add selected cells to the Giotto object
visium_brain <- addPolygonCells(visium_brain,
                                polygon_name = 'selections')

## Print the modified cell_metadata
pDataDT(visium_brain)

## Compare gene expression
comparePolygonExpression(visium_brain,
                         selected_feats = c('Mbp', 'Mobp', 'Tcf712', 'Prkcd'))

## Compare the cell type abundances
compareCellAbundance(visium_brain, cell_type_column = 'leiden_clus')

####################### Supplementary figure 8c ################################

# Pre-run the Parts 1-7 from the [merFISH Mouse Hypothalmic Preoptic Region tutorial](https://drieslab.github.io/Giotto_website/articles/merfish_mouse_hypothalamic.html)

## Create a color-code vector
mycolorcode = c('red', 'lightblue', 'yellowgreen','purple', 'darkred',
                'magenta', 'mediumblue', 'yellow', 'gray')
names(mycolorcode) = c('Inhibitory', 'Excitatory','OD Mature', 'OD Immature',
                       'Astrocyte', 'Microglia', 'Ependymal','Endothelial', 'Ambiguous')

## Run the 3D interactive tool,
## Cell-types selection will include only Astrocyte and Endothelial cells, between z-coords 80 and 145.
my_slice_coordinates <- plotInteractive3D(merFISH_gobject,
                                          cell_color = 'cell_types',
                                          cell_color_code = mycolorcode,
                                          width = "100%", height = "600px")

## Subset object
merFISH_gobject_110 <- subsetGiotto(merFISH_gobject,
                                    cell_ids = my_slice_coordinates$cell_ID)

## Compare expression between cell types
scran_markers = findMarkers_one_vs_all(merFISH_gobject_110,
                                       method = "scran",
                                       cluster_column = 'cell_types',
                                       min_feats = 5,
                                       logFC = 1.5)
top_genes <- scran_markers[, head(.SD, 2), by = 'cluster']$feats

## Plot heatmap
plotMetaDataHeatmap(merFISH_gobject_110,
                    selected_feats = top_genes,
                    metadata_cols = 'cell_types',
                    show_values = 'zscores')

