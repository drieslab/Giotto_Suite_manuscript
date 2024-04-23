if (!requireNamespace("GiottoDB", quietly = TRUE)) {
  remotes::install_github("drieslab/GiottoDB")
}
if (!requireNamespace("duckdb", quietly = TRUE)) {
  BiocManager::install("duckdb")
}


library(GiottoDB)


# The location to place the duckdb instance that represents the values in this
# dataset
database_path <- 'scripts/CREATE/OUTS/stereoseq_E16.5_E2S6/'
dir.create(database_path)

b = GiottoDB::dbBackend(
  drv = duckdb::duckdb(),
  dbdir = database_path
)


# bin1 matrix GiottoDB ingestion ####
# This file is downloadable from:
# https://db.cngb.org/stomics/mosta/download/
f <- "[PATH TO E16.5_E2S6_GEM_bin1.tsv.gz FILE]"

# take a quick look at the data (uncomment to run)
# data.table::fread(f, nrows = 0, header = TRUE)
# expected colnames: "geneID", "x", "y", "MIDCounts"


# select a read function ####
# chunk reader functions should exposes the params:
# `x` : filepath
# `n` : n lines to process/chunk
# `i` : ith chunk
# 
# GiottoDB provides chunked readers that work through data.table and arrow
# This reader function returns a read-in data.table chunk of the file
reader <- function(x, n, i) {
  GiottoDB::stream_reader_fread(x = x, n = n, i = i)
}

# callback ####
# a callback function can be used to edit the read-in chunk before it is saved 
# to the database.
# Here we define one that applies the colnames that GiottoDB expects with
# dbSpatProxy points data.
callback <- function(chunk) {
  cnames <- c("feat_ID", "x", "y", "count")
  data.table::setnames(chunk, new = cnames)
  return(chunk)
}

# run the chunked read and write into the database as a dbSpatProxy points
# ofbject. We also set the geom param as x and y to tell it which cols have xy 
# info.
s <- GiottoDB::dbvect(
  x = f,
  db = b, 
  remote_name = 'rna',
  type = 'points',
  geom = c('x', 'y'),
  read_fun = reader,
  n = 1e7,
  overwrite = TRUE,
  callback = callback,
  report_n_chunks = 1L, # report every chunk processed
  stop_cond = function(x) nrow(x) == 0L # stop when number of rows = 0
)

GiottoDB::gdbReadLog() # opens and readlines the last logfile made

# plot rasterization of stereoseq data
GiottoDB::plot(s)

# check extent of stereoseq data
e <- GiottoDB::ext(s)
e

bInfo <- GiottoDB::getBackendInfo(b)

# save object
saveRDS(s, file = 'scripts/CREATE/OUTS/stereoseq_E16.5_E2S6/ss_points.rds')
saveRDS(bInfo, file = 'scripts/CREATE/OUTS/stereoseq_E16.5_E2S6/backend_info.rds')

