---
knit: bookdown::preview_chapter
---

# Expression overview

## Introduction

Here we will continue to work with the filtered __blischak__ dataset produced in the previous chapter. We will look at what happened to the expression matrix after the quality control step.




```r
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("blischak/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
```

## Top 500 genes

With scater one can look at the proportion of reads accounted by the top 500 expressed genes.

### Before QC


```r
scater::plot(umi[endog_genes, ],
             block1 = "individual",
             block2 = "replicate",
             colour_by = "use",
             exprs_values = "counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-top500-before-qc-1} 

}

\caption{Proportion of reads accounted by the top 500}(\#fig:expr-overview-top500-before-qc)
\end{figure}

### After QC


```r
scater::plot(umi.qc[endog_genes, ],
             block1 = "individual",
             block2 = "replicate",
             colour_by = "use",
             exprs_values = "counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-top500-after-qc-1} 

}

\caption{Proportion of reads accounted by the top 500}(\#fig:expr-overview-top500-after-qc)
\end{figure}

## PCA plot

The easiest thing to overview the data is to transform it using the principal component analysis and then visualize the first two principal components. Again the [scater](https://github.com/davismcc/scater) package provides several very useful functions to make this analysis straightforward.

### Before QC


```r
scater::plotPCA(umi[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                exprs_values = "counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-pca-before-qc-1} 

}

\caption{PCA plot of the blischak data}(\#fig:expr-overview-pca-before-qc)
\end{figure}

### After QC


```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                exprs_values = "counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-pca-after-qc-1} 

}

\caption{PCA plot of the blischak data}(\#fig:expr-overview-pca-after-qc)
\end{figure}

## Diffusion map

Another way of representing the data is a diffusion map. It is very useful if the cells represent a continuous process (e.g. development).

### Before QC


```r
scater::plotDiffusionMap(umi[endog_genes, ],
                         colour_by = "batch",
                         size_by = "total_features",
                         exprs_values = "counts",
                         rand_seed = 123456)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-diff-before-qc-1} 

}

\caption{Diffusion map of the blischak data}(\#fig:expr-overview-diff-before-qc)
\end{figure}

### After QC


```r
scater::plotDiffusionMap(umi.qc[endog_genes, ],
                         colour_by = "batch",
                         size_by = "total_features",
                         exprs_values = "counts",
                         rand_seed = 123456)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-diff-after-qc-1} 

}

\caption{Diffusion map of the blischak data}(\#fig:expr-overview-diff-after-qc)
\end{figure}

## tSNE map

Another way of representing the data is a tSNE map.

### Before QC


```r
scater::plotTSNE(umi[endog_genes, ],
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "counts",
                 rand_seed = 123456)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-tsne-before-qc-1} 

}

\caption{tSNE map of the blischak data}(\#fig:expr-overview-tsne-before-qc)
\end{figure}

### After QC


```r
scater::plotTSNE(umi.qc[endog_genes, ],
                 colour_by = "batch",
                 size_by = "total_features",
                 exprs_values = "counts",
                 rand_seed = 123456)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-overview_files/figure-latex/expr-overview-tsne-after-qc-1} 

}

\caption{tSNE map of the blischak data}(\#fig:expr-overview-tsne-after-qc)
\end{figure}

## Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter).
