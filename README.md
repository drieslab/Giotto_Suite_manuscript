# Giotto_Suite_manuscript
Repository containing all files related to the Giotto Suite manuscript.
The [scripts](scripts/) repository contains all the code used to generate the figures presented.

## important links
- Giotto Suite [repo](http://giottosuite.com)
- Giotto Suite [preprint](https://www.biorxiv.org/content/10.1101/2023.11.26.568752v1)

## folder structure
- [scripts](scripts/) - main directory containing reproducible scripts and associated subdirectories
  - [CREATE](scripts/CREATE/) - intermediate and additional files
    - [EXT_DATA](scripts/CREATE/EXTDATA/) - pre-made additional files
    - [OUTS](scripts/CREATE/OUTS/) - intermediate files, some of which are provided with the repo
  - [SOURCE](scripts/SOURCE/) - scripts to source for additional functionality
  - [Xenium_supporting_files](scripts/Xenium_supporting_files/) - files used with xenium co register and segmentation
 
- [terrabio_startup_scripts](terrabio_startup_scripts/) - directory containing .sh scripts that can be uploaded to [Terra cloud platform](https://app.terra.bio/) for creating a cloud environment with Giotto installed. There is a master.sh file for creating the cloud environment with the previous version of Giotto and suite.sh for creating the cloud environment with the most recent version of Giotto suite.

## data used
- The following spatial datasets were used in this manuscript: 
- Spatial Genomics dataset. The mouse kidney fresh frozen dataset was downloaded from the Spatial Genomics website at https://db.cngb.org/stomics/mosta/download/.
- DBiT-seq dataset. The mouse embryo E10.5 dataset was downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137986.
- Nanostring CosMx dataset. The CosMx FFPE Non-Small Cell Lung Cancer dataset for lung sample 12 was downloaded from the Nanostring website at  https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/ 
- Seq-Scope dataset. The Seq-Scope liver dataset was downloaded directly from the following link, https://deepblue.lib.umich.edu/data/concern/data_sets/9c67wn05f.
- Vizgen dataset. The MERSCOPE/MERFISH FF mouse brain (data release v1.0, May 2021) and FFPE breast cancer (May 2022) datasets were downloaded directly from the Vizgen website at https://info.vizgen.com/mouse-brain-data and https://info.vizgen.com/ffpe-showcase, respectively. 
- 10X Genomics Xenium dataset. Xenium, corresponding images and Visium datasets for human breast cancer were downloaded directly from the 10X website at https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast. 
- 10X Genomics multi-modal Visium CytAssist Human Tonsil dataset. The multi-modal Visium CytAssist Human Tonsil dataset was downloaded from the 10X genomics website at https://www.10xgenomics.com/resources/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard. 
- 10X Genomics Visium Mouse Brain Section (Coronal) dataset. The Visium Mouse Brain Section (Coronal) dataset was downloaded from the 10X Genomics website at https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Adult_Mouse_Brain.
- 10X Genomics Visium HD Human Colorectal Cancer. The Visium HD Human Colorectal Cancer (FFPE) dataset was downloaded from the 10X genomics website at https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc.
- OpenST adult mouse hippocampus dataset. The OpenST adult mouse hippocampus dataset was downloaded from https://rajewsky-lab.github.io/openst/latest/examples/datasets/.
- Slide-seq mouse brain. The processed counts of the Slide-seq mouse brain were downloaded from the NeMO database using the following link https://data.nemoarchive.org/biccn/grant/rf1_macosko/macosko/spatial_transcriptome/cellgroup/Slide-seq/mouse/.
- Stereo-seq dataset. The bin1 matrix files (i.e. *_GEM_bin1.tsv.gz) were downloaded from the CNGB portal at the following link (https://db.cngb.org/stomics/mosta/download/).  
- Imaging Mass Cytometry dataset. Intensity images of human lymph node FFPE tissue were downloaded from a repository created by Bost et. al https://data.mendeley.com/datasets/ncfgz5xxyb/1. 
- Single Cell Mouse Brain Dataset. A single-cell reference dataset published by Manno et al. 2021  was used to identify developmental mouse brain cell types for spatial DWLS deconvolution with Stereo-seq data in this study. This data can be downloaded in the form of a .loom file from the Mousebrain.org website (http://mousebrain.org/development/downloads.html). 
- Single-cell Human Tonsil Dataset. The Atlas of Cells in the Human Tonsil published by Massoni-Badosa et al 2022 containing the annotation of >357,000 cells was used for spatial DWLS deconvolution of Visum CytAssist data in this study. The annotated SpatialExperiment object was downloaded from https://github.com/massonix/HCATonsilData.
- Spatial RNA-ATAC seq mouse embryo ME13. The Mouse embryo ME13 spatial RNA-ATAC seq was originally published by Zhang et al. 2023 and downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205055.   
