##%##########################################################%##
#                                                              #
####      Supplementary figure 4 Multi-scale analysis       ####
#                                                              #
##%##########################################################%##

library(Giotto)

# set working directory to project
setwd("PATH TO Giotto_Suite_manuscript REPO")

# 1. create dataset from vizgen FFPE breast cancer dataset ####
# ---------------------------------------------------------------- #

source("aux_scripts/source/vizgen_hdf5.R")

fovIndexVizgenHDF5()
