if (!requireNamespace("GiottoDB", quietly = TRUE)) {
  remotes::install_github("drieslab/GiottoDB")
}
if (!requireNamespace("rhdf5", quietly = TRUE)) {
  BiocManager::install("rhdf5")
}
if (!requireNamespace("duckdb", quietly = TRUE)) {
  BiocManager::install("duckdb")
}


library(Giotto)


# .gef GiottoDB ingestion

gef_file <- "[PATH TO E16.5_E2S6.gef FILE]"
gef_manifest = rhdf5::h5ls(gef_file)


# The location to place the duckdb instance that represents the values in this
# dataset
database_path <- 'scripts/CREATE/OUTS/stereoseq_E16.5_E2S6/'
dir.create(database_path)

b = GiottoDB::dbBackend(
  drv = duckdb::duckdb(),
  dbdir = database_path
)


# declare reader function
read_fun <- function(x, n, i) {
  GiottoDB::stream_reader_gef_tx(
    x = x,
    n = n,
    i = i,
    bin_size = "bin1",
    output = "DT")
}

s <- GiottoDB::dbvect(
  x = gef_file,
  db = b, 
  remote_name = 'rna',
  type = 'points',
  geom = c('x', 'y'),
  read_fun = read_fun,
  n = 1000L,
  overwrite = TRUE,
  callback = NULL,
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

