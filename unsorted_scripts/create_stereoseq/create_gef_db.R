
# .gef GiottoDB ingestion testing
# 292300971 bin1
# gef_file = '/projectnb/rd-spat/HOME/ecruiz/StereoSeq/GiottoPaper/data/gef_files/E16.5_E2S6.bgef'
# 359630461 bin1
gef_file = '/projectnb/rd-spat/HOME/ecruiz/StereoSeq/GiottoPaper/data/gef_files/E16.5_E2S5.bgef'

gef_manifest = rhdf5::h5ls(gef_file)

library(Giotto)
library(GiottoDB)

# database_path <- 'DATA/stereoseq/'
database_path <- 'DATA/stereoseq/E16.5_E2S5/'

b = createBackend(
  drv = duckdb::duckdb(),
  dbdir = database_path
)


# declare reader function
read_fun <- function(x, n, i) {
  stream_reader_gef_tx(
    x = x,
    n = n,
    i = i,
    bin_size = "bin1",
    output = "DT")
}

s <- dbvect(
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

gdbReadLog() # opens and readlines the last logfile made

# plot rasterization of stereoseq data
plot(s)

# check extent of stereoseq data
e <- ext(s)
e

bInfo <- getBackendInfo(b)

# save object
saveRDS(s, file = 'DATA/stereoseq/ss_points.rds')
saveRDS(bInfo, file = 'DATA/stereoseq/backend_info.rds')

