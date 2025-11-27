# Skrupskelyte_etal_SI_code





<a href="https://doi.org/10.5281/zenodo.17734456"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17734456.svg" alt="DOI"></a>



## Overview

This repository contains R code used to produce analyses and figures for the supplementary information accompanying Skrupskelyte et al. The scripts perform data processing, statistical analyses and generate plots used in the manuscript's supplementary materials.

## Copynumber

### Install dependencies

Run the following in R to install the common dependencies:

```r
# Install CRAN packages (grid is part of base R)
cran_pkgs <- c("png", "gridExtra", "gtools", "ggplot2", "logging")
install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c("QDNAseq", "Biobase", "VariantAnnotation")
BiocManager::install(bioc_pkgs, ask = FALSE)
```


### Running the R script

Example shell script for the main off-target copy-number script:

```bash
BED_FILE=baitset.bed
SEGMENT_GAMMA=10
SEG_KMIN=3
BINS="qdnaseq/precalculated_windows/QNDAseq_bins100.RData"
RUN_NAME="SEG_GAMMA_10"
SEX="NULL"

Rscript qdnaseq_call_offtarget_read_cn.R \
    ${ID} \
    ${TUMOUR} \
    ${NORMAL} \
    ${BED_FILE} \
    ${SEX} \
    ${SEGMENT_GAMMA} \
    ${SEG_KMIN} \
    ${BINS} \
    ${RUN_NAME}
```
QDNAseq bins and the bed file used are provided at the Figshare link here : ................

A BED file for the targeted sequencing regions is provided in the `References` folder using GRCM38 coordinates.

