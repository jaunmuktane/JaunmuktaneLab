# Analysing GeoMx DSP dataset with standR

## Overview

This is a repository for the R Code used to analyse Nanostring GeoMx Data from the Nanostring GeoMx DSP platform.

Data QC was performed using the GeoMx-NGS RNA Expression Data with GeomxTools. Please see - https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

Data normalisation and batch correction using methods in the `standR` package as well as differential expression using 'voom' and 'standR' were performed as per - 
https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html

## Pre-requisites 

This code is for researchers who are interested in analysing GeoMx transcriptomics data. 

## _R_ packages used

The following key R packages were used: 

* `NanoStringNCTools`
* `GeomxTools`
* `GeoMxWorkflows`
* `gridExtra`
* `standR`
* `edgeR`
* `limma`
* `msigdb`
* `GSEABase`
* `igraph`
* `vissE`
* `SpatialExperiment`
* `Scater`
* `htmlwidgets`
* `svDialogs`
* `msigdb`
* `DT`
* `EnhancedVolcano`
*  `xlsx`
* `tidyverse`

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Authors

* **Hemanth Nelvagal** - *Initial work*
* **Toby Curless** - *Initial work* 
* **Zane Jaunmuktane** - *Initial work* 

## Acknowledgments

* https://bioconductor.org/packages/devel/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html
* https://davislaboratory.github.io/GeoMXAnalysisWorkflow/articles/GeoMXAnalysisWorkflow.html

