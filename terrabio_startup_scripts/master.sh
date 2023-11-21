#!/usr/bin/env bash

pip install --upgrade libgsl-dev
pip install --upgrade libglpk-dev
pip install --upgrade libgdal-dev
pip install --upgrade imagemagick
pip install --upgrade libmagick++-dev
pip install --upgrade igraph
pip install --upgrade leidenalg
pip install --upgrade networkx
pip install --upgrade community

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.18')"
RUN R -e "BiocManager::install('beachmat', dependencies=TRUE)"

R -e "install.packages('remotes', dependencies=TRUE)"
R -e "remotes::install_github('drieslab/Giotto@master')"


