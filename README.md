# DEx Multi-Omics in Parkinson Disease (PD)

## Description

In this repository is presented demo files and R scripts to reproduce one of the bioinformatic pipelines (framework A) and some of the visualization available in [insert DOI].

Particularly, here you can find functions to ellaborate differential expression analyses of RNA-seq data (total and small RNA) using the Bioconductor packages [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/manuals/TCGAbiolinks/man/TCGAbiolinks.pdf).

## Usage 

Demo files are provided to conduct the differential expression analysis. This can be found in the `Data/` folder. 
Nevertheless, please check the manuscript to access the multi-omics data used for the study.

This project was conducted in [R software](https://www.r-project.org). 
As mention previously, here we present one of the bioinformatic pipelines available in our manuscript (framework A), in the `R Scripts/` folder. However, for the other bioinformatic methodology (framework B), please revert to this [GitHub page](https://github.com/gauravj49/BulkRnaseqDE).

Bear in mind, here is only provided an adaptation of the original source code. 

All the necessary R package dependencies are

* ggplot2
* magrittr
* ggpubr
* readxl
* readr
* dplyr
* nortest
* tidyverse
* plyr
* ashr
* plm

These packages and dependencies should be installed a priori using the install.packages() function and complemented by the library() function to be ready to use. 

Furthermore, we also leveraged packages available in the Bioconductor:

* edgeR
* limma
* biomaRt
* DESeq2
* apeglm
* vsn 
* TCGAbiolinks

Similarly, to employ these packages, first install the BiocManager package using install.packages("BiocManager"), and later the packages above using BiocManager::install(). 

For the RNA decomposition, we used the immunedeconv R package. This can be installed using the query:

```
install.packages("remotes")
remotes::install_github("icbi-lab/immunedeconv")
```

## Contact

For any inquiries related to this work, please contact me via e-mail ana.galhoz@helmholtz-muenchen.de
