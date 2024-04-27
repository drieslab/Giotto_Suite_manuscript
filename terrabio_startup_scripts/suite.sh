#!/usr/bin/env bash

pip install --upgrade pip
pip install libgsl-dev
pip install --upgrade libgsl-dev
pip install --upgrade libglpk-dev
pip install --upgrade libgdal-dev
pip install --upgrade imagemagick
pip install --upgrade libmagick++-dev
pip install --upgrade igraph
pip install --upgrade leidenalg
pip install --upgrade networkx
pip install --upgrade community


RUN R -e "BiocManager::install('beachmat', dependencies=TRUE)"
RUN R -e "BiocManager::install('HDF5Array', dependencies=TRUE)"

R -e "install.packages('pak', dependencies = TRUE)"


R -e "pak::pkg_install('drieslab/GiottoUtils')"
R -e "pak::pkg_install('drieslab/GiottoClass')"
R -e "pak::pkg_install('drieslab/GiottoVisuals')"
R -e "pak::pkg_install('drieslab/GiottoData')"

R -e "pak::pkg_install('drieslab/Giotto')"


