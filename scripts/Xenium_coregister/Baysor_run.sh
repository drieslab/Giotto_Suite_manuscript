#!/bin/bash -l

#$ -l h_rt=72:00:00   # Specify the hard time limit for the job
#$ -P rd-spat
#$ -pe omp 28
#$ -l mem_per_core=16G

#JULIA_NUM_THREADS=28 Baysor preview -x x_location -y y_location -g feature_name Xenium_gene.csv 

JULIA_NUM_THREADS=28 Baysor run -x x_location -y y_location -g feature_name --plot --save-polygons=geojson -m 30 --prior-segmentation-confidence 0.5 Xenium_gene.csv :cell_id

