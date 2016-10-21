---
output: html_document
---

# Technical requirements

This course is based on the popular programming language [R](https://www.r-project.org/). However, one of the methods that we describe ([SNN-Cliq](http://bioinfo.uncc.edu/SNNCliq/)) is only partly R-based. It makes a simple _python_ call from R and requires a user to have write permissions to the working directory.

To be able to run all code chunks of the course one needs to clone or download the [course GitHub repository](https://github.com/hemberg-lab/scRNA.seq.course) and start an R session in the cloned folder. 

One also needs to install the following R packages (ordered by purposes):

## General

[devtools](https://cran.r-project.org/web/packages/devtools/index.html) for installing packages from GitHub:

```r
install.packages("devtools")
```

[BiocInstaller](http://bioconductor.org/packages/release/bioc/html/BiocInstaller.html) for installing packages from BioConductor:

```r
source('https://bioconductor.org/biocLite.R')
biocLite('BiocInstaller')
```

[scRNA.seq.funcs](https://github.com/hemberg-lab/scRNA.seq.funcs) - R package containing some special functions used in this course:

```r
devtools::install_github("hemberg-lab/scRNA.seq.funcs")
```

## Plotting

[ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) for plotting general plots:

```r
install.packages("ggplot2")
```

[pheatmap](https://cran.r-project.org/web/packages/pheatmap/index.html) for plotting heatmaps:

```r
install.packages("pheatmap")
```

[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) for plotting Venn diagrams:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
```

## QC and normalisation

[scater](http://bioconductor.org/packages/release/bioc/html/scater.html) is a single-cell analysis toolkit for expression:

```r
source('https://bioconductor.org/biocLite.R')
biocLite('scater')
```

[mvoutlier](https://cran.r-project.org/web/packages/mvoutlier/index.html) - for an automatic outlier detection used by the [scater](https://github.com/davismcc/scater) package.

```r
install.packages("mvoutlier")
```

[statmod](https://cran.r-project.org/web/packages/statmod/index.html) - a dependency for [mvoutlier](https://cran.r-project.org/web/packages/mvoutlier/index.html).

```r
install.packages("statmod")
```

[Rtsne](https://cran.r-project.org/web/packages/pheatmap/index.html) for Rtsne data embedding:

```r
install.packages("Rtsne")
```

[scran](http://bioconductor.org/packages/release/bioc/html/scran.html) for a new single cell normalisation method (LSF):

```r
source('https://bioconductor.org/biocLite.R')
biocLite('scran')
```

[RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html) for normalization using ERCC controls:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("RUVSeq")
```

## Clustering

[pcaReduce](https://github.com/JustinaZ/pcaReduce) for unsupervised clustering of scRNA-seq data ([bioRxiv](http://biorxiv.org/content/early/2015/09/08/026385)):

```r
devtools::install_github("JustinaZ/pcaReduce")
```

[pcaMethods](http://bioconductor.org/packages/release/bioc/html/pcaMethods.html) is a [pcaReduce](https://github.com/JustinaZ/pcaReduce) dependency:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")
```

[SC3](http://bioconductor.org/packages/devel/bioc/html/SC3.html) for unsupervised clustering of scRNA-seq data:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("SC3")
```

* Before running __SC3__ for the first time __only__, please start R and enter:

```r
RSelenium::checkForServer()
```

[SEURAT](https://github.com/Puriney/seurat) for density clustering of scRNA-seq data:

```r
devtools::install_github('Puriney/seurat')
```

## Dropouts

[M3Drop](https://github.com/tallulandrews/M3D) for identification of important and DE genes:

```r
devtools::install_github("tallulandrews/M3D")
```

## Pseudotime

[TSCAN](https://bioconductor.org/packages/release/bioc/html/TSCAN.html) for pseudotime analysis:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("TSCAN")
```

[destiny](https://bioconductor.org/packages/release/bioc/html/destiny.html) for pseudotime analysis:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("destiny")
```


## Differential Expression

[ROCR](https://cran.r-project.org/web/packages/ROCR/index.html) for performance estimations:

```r
install.packages("ROCR")
```

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for identification of differentially expressed genes:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

[scde](http://hms-dbmi.github.io/scde/) for identification of differentially expressed genes:

```r
devtools::install_github("hms-dbmi/scde", build_vignettes = FALSE)
```

* Installation on Mac OS X may require [this additional gfortran library](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/):
```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

* See the [help page](http://hms-dbmi.github.io/scde/help.html) for additional support.



