---
knit: bookdown::preview_chapter
---

# Remove confounders using controls

## Introduction

In the previous chapter we saw that the library size normalisation is able to remove major confounders. Let us see whether further use of the controls (ERCC and MT) can provide us with more insights into the data.

Several methods (eg. [BASiCS](https://github.com/catavallejos/BASiCS), [scLVM](https://github.com/PMBio/scLVM), [RUV](http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)) have been developed for normalisation using ERCCs using different noise models and different fitting procedures.

We will demonstrate some of the methods starting from the simplest one proposed by [Brennecke et al.](http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2645.html), which identifies genes with significant variation above technical noise (ERCCs).




```r
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("blischak/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
```

## Brennecke method

To use the method, we first normalize for library size then calculate
the mean and the square coefficient of variation (variation divided by
the squared mean expression). A curve is fit for the relationship
between these two variables for the ERCC spike-in (subject to just
technical variation) then a chi-square test is used to find genes
significantly above the curve. This has been provided for you as the
Brenneck_getVariableGenes(counts, spikes) function.


```r
umi.qc <- 
    scater::normaliseExprs(umi.qc,
                           method = "RLE",
                           feature_set = endog_genes,
                           lib.size = rep(1, ncol(umi.qc)))
erccs <- grep("ERCC-", rownames(assayData(umi.qc)$norm_counts))
highly.var.genes <- scRNA.seq.funcs::Brennecke_getVariableGenes(
            assayData(umi.qc)$norm_counts,
            erccs
)
```

```
## Loading required package: statmod
```

```
## Warning in scRNA.seq.funcs::Brennecke_getVariableGenes(assayData(umi.qc)
## $norm_counts, : Only 21 spike-ins to be used in fitting, may result in poor
## fit.
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-remove-conf_files/figure-latex/rm-conf-brennecke-1} 

}

\caption{Results of using the Brennecke method on the Blischak dataset}(\#fig:rm-conf-brennecke)
\end{figure}

In the figure above blue points are the ERCC spike-ins. The red curve
is the fitted technical noise model and the dashed line is the 95%
CI. Pink dots are the genes with significant biological variability
after multiple-testing correction. Since our dataset is relatively
homogeneous only 148 genes are identified as significantly
variable.

## Remove Unwanted Variation

Factors contributing to technical noise frequently appear as "batch
effects" where cells processed on different days or by different
technicians systematically vary from one another. Removing technical
noise and correcting for batch effects can frequently be performed
using the same tool or slight variants on it. We will be considering
the [Remove Unwanted Variation (RUV)](http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)
method which uses [singular value decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition)
(similar to PCA) on subsets of the dataset which should be invariant
(e.g. [ERCC spike-ins](https://www.thermofisher.com/order/catalog/product/4456740)). Then the method removes the identified unwanted factors.


```r
assayData(umi.qc)$ruv_counts <- RUVSeq::RUVg(
    round(assayData(umi.qc)$norm_counts),
    erccs,
    k = 1)$normalizedCounts
```

## Effectiveness 1

We evaluate the effectiveness of the normalization by inspecting the
PCA plot where shape corresponds the technical replicate and colour
corresponds to different biological samples (individuals from whom the
iPSC lines where derived). Separation of biological samples and
interspersed batches indicates that technical variation has been
removed. 


```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-remove-conf_files/figure-latex/rm-conf-pca-rle-1} 

}

\caption{PCA plot of the blischak data after RLE normalisation}(\#fig:rm-conf-pca-rle)
\end{figure}


```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "ruv_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-remove-conf_files/figure-latex/rm-conf-pca-rle-ruv-1} 

}

\caption{PCA plot of the blischak data after RLE and RUV normalisations}(\#fig:rm-conf-pca-rle-ruv)
\end{figure}

## Effectiveness 2

We can also examine the relative log expression (RLE) across cells to
confirm technical noise has been removed from the dataset.


```r
boxplot(list(scRNA.seq.funcs::calc_cell_RLE(assayData(umi.qc)$norm_counts),
             scRNA.seq.funcs::calc_cell_RLE(assayData(umi.qc)$ruv_counts)))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-remove-conf_files/figure-latex/rm-conf-rle-comp-1} 

}

\caption{Comparison of the relative log expression of the blischak data before and after the RUV normalisation}(\#fig:rm-conf-rle-comp)
\end{figure}

## Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter). Additionally, experiment with other combinations of normalizations and compare the results.
