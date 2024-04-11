# **************************************************************************** #
# This script contains a set of functions to ingest Vizgen MERSCOPE data outputs
# in a way that matches with the original HDF5 boundary polygons outputs.
#
# The `\cell_boundaries\` folder of these MERSCOPE dataset outputs contains an
# HDF5 file for each of the FOVs in the dataset. Each HDF5 contains the polygons
# segmented for that FOV, separated by z layer in the internal hierarchical
# structure.
#
# These functions index through the HDF5 files and return a manifest of all the
# polygons contained within, which file_ID/FOV they are associated with, and
# the spatial extent covered by the FOV. All of this information is then
# organized within an S4 `H5TileProxy` that is defined in this script.
#
# Queries with spatial extents can then check for overlaps with the extents of
# the FOVs that have been indexed, allowing fast access to subsets of the full
# dataset.
#
#
#
# Moving forward, these functions are not necessary and a simpler workflow
# can be used now that Vizgen has switched to using parquet for their outputs.
#
# Other scripts can source() this script in order to use its functionalities.
# **************************************************************************** #





## Retrieve values in extent ####

# Takes the crop extent and extracts the desired spatial region from the HDF5
# polygon values using the H5TileProxy generated from the vizgen dataset and
# returns a terra SpatVector object.
#
# NOTE: SpatRaster input is required (the image created with createMerscopeLargeImage)
# This is so that spatial flipping can take the entire dataset's spatial bounds
# into account.
H5TPqueryPolys = function(vizproxy, SpatRaster, poly_crop_ext) {
  e = terra::ext(SpatRaster)
  # find y midline
  y_range = c(e$ymin, e$ymax)
  names(y_range) = NULL
  y_midline = mean(y_range)

  # flip crop extent about the y midline
  poly_crop_ext_flip = flip_extent(poly_crop_ext, direction = 'vertical', y0 = y_midline)

  out = get_values_extent(vizproxy, poly_crop_ext_flip)

  if(nrow(out) == 0L) stop('no geometries found within specified extent')
  out_p = Giotto::createGiottoPolygonsFromDfr(out)
  out_p@spatVector = terra::makeValid(out_p@spatVector)
  out_p@spatVector = terra::crop(out_p@spatVector, poly_crop_ext_flip)
  # perform remaining needed steps to flip
  # terra flips spatvectors over the minimum y extent value
  # this means that the extent of the resultant spatvector points (specifically) are needed for downstream positioning
  out_p@spatVector = terra::flip(out_p@spatVector)
  out_p@spatVector = terra::shift(out_p@spatVector, dy = (poly_crop_ext$ymax - terra::ext(out_p@spatVector)$ymax))
}


# Using the crop extent, find the matching location in the Vizgen transcript
# detections output (should be provided as a data.table). Returns the matching
# crop region in data.table format.
# This requires flipping of the crop extent to match the coordinates orientation
# of the points data. The values are selected and then the y values are inverted
# so that the cropped points information are in the same spatial orientation as
# the images and the polygons information.
ext_query_points = function(xy, SpatRaster, crop_ext) {
  e = terra::ext(SpatRaster)
  # find y midline
  y_range = c(e$ymin, e$ymax)
  names(y_range) = NULL
  y_midline = mean(y_range)

  # flip crop extent about the y midline
  crop_ext_flip = flip_extent(crop_ext, direction = 'vertical', y0 = y_midline)

  out = xy[x >= crop_ext_flip$xmin & x <= crop_ext_flip$xmax & y >= crop_ext_flip$ymin & y <= crop_ext_flip$ymax]
  out[, y := -y + (crop_ext$ymin - flip_extent(crop_ext_flip)$ymin)]
}


# A helper function to return SpatExtent that has been flipped across the specified
# y or x value. This has been incorporated into Giotto as a method for SpatExtent
# objects.
flip_extent = function(x, direction = c('vertical', 'horizontal'), x0 = 0, y0 = 0) {
  direction = match.arg(direction, choices = c('vertical', 'horizontal'))
  x_range = c(x$xmin, x$xmax) ; names(x_range) = NULL
  y_range = c(x$ymin, x$ymax) ; names(y_range) = NULL
  switch(direction,
         'vertical' = {
           y_range = -(y_range - y0) + y0
           terra::ext(c(x_range, y_range[2], y_range[1]))
         },
         'horizontal' = {
           x_range = -(x_range - x0) + x0
           terra::ext(c(x_range[2], x_range[1], y_range))
         })
}




# Function that retrieves the values into memory from an H5TileProxy object based
# on the spatial extent that is desired.
# If any filters for which set of data should be retrieved then they are also
# applied.
# Filters should be applied directly to the H5TileProxy prior to this step.
get_values_extent = function(h5tp, extent) {
  parser = h5tp@parser # get parser function
  root = h5tp@root # get root search path
  file_list = h5tp@files
  x_range = c(extent$xmin, extent$xmax) ; names(x_range) = NULL
  y_range = c(extent$ymin, extent$ymax) ; names(y_range) = NULL

  # select relevant FOV based on extent
  fov_ext = h5tp@extents
  fov_sel = fov_ext[xmin <= x_range[2] &
                      xmax >= x_range[1] &
                      ymin <= y_range[2] &
                      ymax >= y_range[1],
                    fov]
  manifest = h5tp@manifest[name == 'coordinates' & fov %in% fov_sel]

  # additional filtering based on tokens
  filters = h5tp@filter
  if(!is.null(filters)) {
    for(token in names(filters)) {
      manifest = manifest[eval(as.name(token)) %in% filters[[token]]]
    }
  }

  # find specific fovs after filtering
  poly_fovs = manifest[, unique(fov)]
  # function to get vertices for each polygon within an fov dataset
  .fov_polys = function(h5id, dfns, pids) {
    poly_xy = lapply(seq_along(dfns), function(poly_i) {
      d = rhdf5::h5read(h5id, dfns[[poly_i]]) # read hdf5 to get array
      d = data.table::as.data.table(parser(d)) # parse array into xy values
      d[, poly_ID := pids[[poly_i]]] # append poly_ID information
      d
    })
    data.table::rbindlist(poly_xy) # combine into single data.table and return
    # DT contains xy and poly_ID of every polygon within FOV
  }


  xy = future.apply::future_lapply(
    poly_fovs,
    future.packages = c('rhdf5', 'data.table'),
    future.seed = TRUE,
    function(fov_i) {

      file_i = manifest[fov == fov_i, unique(file_ID)]
      Dfullnames = manifest[fov == fov_i, fullname]
      poly_IDs = manifest[fov == fov_i, V1] # hardcoded
      fid = rhdf5::H5Fopen(file_list[[file_i]], flags = 'H5F_ACC_RDONLY')
      on.exit(rhdf5::H5Fclose(fid))
      if(!is.na(root)) {
        gid = rhdf5::H5Gopen(fid, name = root)
        on.exit(rhdf5::H5Gclose(gid), add = TRUE, after = FALSE)
        fov_poly_DT = .fov_polys(h5id = gid, dfns = Dfullnames, pids = poly_IDs)
      } else {
        fov_poly_DT = .fov_polys(h5id = fid, dfns = Dfullnames, pids = poly_IDs)
      }
      return(fov_poly_DT)
    }
  )
  # combine poly information across all selected FOV
  data.table::rbindlist(xy)
}













## manifest generation ####

# Generate manifest of available information within the HDF5 file including all
# of the H5 nesting for each dataset. This nesting information can be selected
# for by using the filters.
# The depth at which spatial information can be set using the root param to
# account for HDF5 datasets where the spatial information is only within one
# subdirectory.
# Scanning is parallelized using future_lapply (see get_manifest_single) across
# sets of multiple HDF5 files.
# Returns data.table


#' @param file_list list of files
#' @param root root hierarchy name in .h5 from which to begin cataloging
#' @param token_names names to assign the token/colnames
#' @param fov_tokens combination of tokens/colnames in output that uniquely identifies
#' a FOV. Defaults to selecting each file as its own fov by passing 'file_ID'.
#' Select specific tokens
#' @param verbose be verbose
get_manifest = function(file_list, root = NA_character_, token_names = NULL, fov_tokens = c('file_ID'), verbose = TRUE) {
  if(verbose) {
    s_time = proc.time()
    GiottoUtils::vmsg(.v = "log", 'cataloguing...\n [start] ')
  }

  out = future.apply::future_lapply(
    seq_along(file_list),
    future.packages = c('data.table', 'rhdf5'),
    future.seed = TRUE,
    function(root, file_list, file_i) {
      get_manifest_single(i = file_i, file_list = file_list, root = root)
    },
    root = root,
    file_list = file_list
  )

  manifest = data.table::rbindlist(out)

  manifest[, group := gsub('^/', '', group)] # remove leading '/' character in name
  manifest = cbind(manifest[, data.table::tstrsplit(group, split = '/')], manifest[, c('name', 'group', 'file_ID')])

  # create fullname column
  manifest[, group := do.call(paste, c(.SD, sep = '/')), .SDcols = c('group', 'name')]
  data.table::setnames(manifest, 'group', 'fullname')

  # add fov identifier
  manifest[, fov := .GRP, by = fov_tokens]

  # apply token_names if available
  if(!is.null(token_names)) {
    token_oldnames = paste0('V', seq(length(token_names)))
    data.table::setnames(manifest, old = token_oldnames, new = token_names)
  }

  # format for return
  data.table::setcolorder(manifest, c('file_ID', 'fov'))
  key_names = unique(c('file_ID', 'fov', names(manifest)[names(manifest) != 'fullname']))
  data.table::setkeyv(manifest, key_names)
  if(verbose) {
    GiottoUtils::vmsg(.v = "log", ' [finish] ', data.table::timetaken(s_time))
  }
  return(manifest)
}

# HDF5 manifest scanning operation for a single file.
get_manifest_single = function(i, file_list, root = NA_character_) {
  fid = rhdf5::H5Fopen(file_list[[i]], flags = 'H5F_ACC_RDONLY')
  on.exit(rhdf5::H5Fclose(fid))
  if(!is.na(root)) {
    gid = rhdf5::H5Gopen(fid, name = root)
    on.exit(rhdf5::H5Gclose(gid), add = TRUE, after = FALSE)
    manifest = data.table::setDT(h5ls(gid))[otype == 'H5I_DATASET', .SD, .SDcols = c('group', 'name')]
  } else {
    manifest = data.table::setDT(h5ls(fid))[otype == 'H5I_DATASET', .SD, .SDcols = c('group', 'name')]
  }

  manifest[, c('file_ID', 'fov') := i]
  # debug
  print(basename(file_list[[i]]))

  manifest
}







## parser ####
# Callback functions that are used to read and/or format the values within the
# HDF5 into expected formats when they are requested.

# parse functions should be attached to the H5TileProxy object

# Spline polynomial to smooth polygon inputs
spline_poly <- function(xy, vertices = 20, k = 3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.

  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }

  # Spline the x and y coordinates.
  data.spline <- stats::spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- stats::spline(1:(n+2*k), data[,2], n=vertices, ...)$y

  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}

# parser function for Vizgen HDF5 arrays
parse_array = function(a) {
  xy = data.table::data.table(
    x = a[1,,],
    y = a[2,,]
  )
  xy = as.matrix(xy)
  xy = spline_poly(xy, vertices = 60)
  colnames(xy) = c('x', 'y')
  return(xy)
}

# terra geom matrix formatting
# geom_format = function(xy) {
#   xy = cbind(1,1, xy, 0)
#   colnames(xy) = c('geom', 'part', 'x', 'y', 'hole')
#   xy
# }








## extent scanning ####

# extent scanner determines spatial extent for each fov (combination of manifest
# item and any indexing or chunking possible that might exist. Chunk detection
# might be use case specific for the moment.)
#
# Extent scanning is performed using he manifest of datasets.

# extent scan function for Vizgen HDF5 arrays
scan_extent_vizgen = function(manifest, file_list,
                              root = NA_character_,
                              verbose = TRUE) {

  manifest = manifest[name == 'coordinates']

  if(verbose) {
    s_time = proc.time()
    GiottoUtils::vmsg(.v = "log", 'scanning extents...\n [start] ')
  }

  out = future.apply::future_lapply(
    unique(manifest$fov), # fov_i
    future.packages = c('data.table', 'rhdf5'),
    future.seed = TRUE,
    function(fov_i) {
      scan_extent_single_vizgen(manifest = manifest, fov_i = fov_i, file_list = file_list,
                                parser = parse_array, root = root)
    }
  )

  fov_ext = do.call(rbind, out)
  fov_ext = data.table::as.data.table(fov_ext)
  data.table::setnames(fov_ext, c('fov', 'xmin', 'xmax', 'ymin', 'ymax'))
  GiottoUtils::vmsg(.v = "log", ' [finish] ', data.table::timetaken(s_time))
  if(verbose) {
    message("[finish]", data.table::timetaken(s_time))
  }
  
  return(fov_ext)
}




# method of getting x and y values are from vizgen specific parser
# filtering of manifest info is vizgen specific
# parser for array data is also vizgen specific
scan_extent_single_vizgen = function(manifest, parser = parse_array, ...) {

  scan_extent_single(manifest = manifest, parser = parser, ...)
}


# Function to find the extent of a whole fov
# Framework function from which method-specific scanners can be made
scan_extent_single = function(manifest,
                              fov_i,
                              file_list,
                              parser = function(a) {
                                xy = data.table::data.table(
                                  x = a[1,,],
                                  y = a[2,,]
                                )
                                as.matrix(xy)
                              },
                              root = NA_character_) {

  .find_fov_ext = function(h5id, dfns) {
    poly_ext = lapply(dfns, function(dfn) {
      d = rhdf5::h5read(h5id, dfn)
      d = parser(d)
      c(min(d[,'x']), max(d[,'x']), min(d[,'y']), max(d[,'y']))
    })
    ext_array = do.call(rbind, poly_ext)
    c(min(ext_array[,1]), max(ext_array[,2]),
      min(ext_array[,3]), max(ext_array[,4]))
  }

  file_i = manifest[fov == fov_i, unique(file_ID)]
  Dfullnames = manifest[fov == fov_i, fullname]
  fid = rhdf5::H5Fopen(file_list[[file_i]], flags = 'H5F_ACC_RDONLY')
  on.exit(rhdf5::H5Fclose(fid))
  if(!is.na(root)) {
    gid = rhdf5::H5Gopen(fid, name = root)
    on.exit(rhdf5::H5Gclose(gid), add = TRUE, after = FALSE)
    fov_ext = .find_fov_ext(h5id = gid, dfns = Dfullnames)
  } else {
    fov_ext = .find_fov_ext(h5id = fid, dfns = Dfullnames)
  }

  # return FOV extent
  c(fov_i, fov_ext)
}



# HDF5 Vizgen table

## H5TileProxy class ####

# Generalized framework for H5

# Scan through and compile a DT filterable list of value locations
# FOV, tokens --- 1 row / cell, key on FOV then cell_ID then all other tokens
# store filepaths separately as vector (more space efficient)
# store EXT separately as table w/ FOV key (more space efficient)

# generic [] accept terra::ext()

# pkgs: rhdf5, data.table, future, future.callr

# root_ID - a base query. Each may or may not be a separate file. Could also be
# different base queries within a single file

# build a retrieval statement:
# file | root | manifest (fullname) | chunk idx | value idx

h5TileProxy = setClass(
  'h5TileProxy',
  slots = c(
    files = 'character', # map file to file_ID
    root = 'character', # root search path as ordered character vector which manifest builds upon
    manifest = 'data.table', # map entries to fov and file_ID
    extents = 'data.table', # map fov to extents
    filter = 'list', # list character vectors of selected token values
    parser = 'function',
    init = 'logical'
  ),
  prototype = list(
    files = NA_character_,
    root = NA_character_,
    manifest = data.table::data.table(),
    extents = data.table::data.table(),
    filter = list(),
    parser = parse_array,
    init = FALSE
  )
)



## H5TileProxy methods ####

methods::setMethod('initialize', signature('h5TileProxy'), function(.Object, verbose = TRUE, ...) {
  .Object = methods::callNextMethod(.Object, ...)
  if(!.Object@init) {
    if(!is.na(.Object@files)) {

      .Object@files = list.files(.Object@files, full.names = TRUE)
      .Object@manifest = get_manifest(file_list = .Object@files, root = .Object@root, verbose = verbose)
      .Object@extents = scan_extent_vizgen(manifest = .Object@manifest[name == 'coordinates'],
                                           file_list = .Object@files,
                                           root = .Object@root, verbose = verbose)
      .Object@init = TRUE
    }
  }
  methods::validObject(.Object)
  .Object
})

methods::setMethod('names', signature('h5TileProxy'), function(x) {
  names(x@manifest)
})
methods::setReplaceMethod('names', signature(x = 'h5TileProxy', value = 'character'), function(x, value) {
  op = names(x@manifest)
  names(x@manifest) = value
  for(i in seq_along(value)) {
    names(x@select)[which(names(x@select) == op[[i]])] = value[[i]]
  }
})
methods::setMethod('show', signature('h5TileProxy'), function(object) {
  e = object@extents
  cat('Object of class ', class(object), '\n')
  cat('dir       : ', gsub(basename(object@files)[1], '', object@files[1]), '\n')
  cat('files     : ', length(object@files), '\n')
  cat('fovs      : ', nrow(e), '\n')
  cat('extent    : ', min(e$xmin), max(e$xmax), min(e$ymin), max(e$ymax), '(xmin, xmax, ymin, ymax)\n')
  cat('datasets  : ', nrow(object@manifest), '\n')
  cat('tokens    : \n')
  for(token in names(object@filter)) {
    cat('            ', token, '~', object@filter[[token]], '\n')
  }
})
library(terra)
methods::setMethod('ext', signature('h5TileProxy'), function(x, ...) {
  e = x@extents
  terra::ext(min(e$xmin), max(e$xmax), min(e$ymin), max(e$ymax))
})



# change the filepaths recorded within the H5TileProxy @files slot so that
# they point to a new location
H5TPmigratePaths <- function(x, new_dir) {
  path <- x@files[[1L]]
  old_dir <- gsub(basename(path), "", path)
  base_paths <- gsub(old_dir, "", x@files)
  x@files <- file.path(new_dir, base_paths)
  return(x)
}




# TileProxy Generation Workflow ####

# Generate an H5TileProxy that is provides a spatial index of where the polygonal
# information is within the Vizgen HDF5 files folder


# poly_dir directory containing the Vizgen HDF5 boundary files
# out_dir is where to save the FOV index
fovIndexVizgenHDF5 <- function(poly_dir, out_dir) {

  # polygons ####
  # create H5TileProxy
  mb_vizproxy = h5TileProxy(root = 'featuredata',
                            files = poly_dir,
                            parser = parse_array)

  token_filter = list(V3 = 'p_0',
                      V2 = 'zIndex_0')
  mb_vizproxy@filter = token_filter

  qs::qsave(mb_vizproxy, file = file.path(out_dir, 'vizH5Proxy.qs'))
}



