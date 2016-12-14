---
output: html_document
---

# Dealing with confounders

## Introduction

In the previous chapter we normalized for library size, effectively removing it as a confounder. Now we will consider removing other less well defined confounders from our data. Technical confounders (aka batch effects) can arise from difference in reagents, isolation methods, the lab/experimenter who performed the experiment, even which day/time the experiment was performed. Accounting for technical confounders, and batch effects particularly, is a large topic that also involves principles of experimental design. Here we address approaches that can be taken to account for confounders when the experimental design is appropriate.

Fundamentally, accounting for technical confounders involves identifying and, ideally, removing sources of variation in the expression data that are not related to (i.e. are confounding) the biological signal of interest. Various approaches exist, some of which use spike-in or housekeeping genes, and some of which use endogenous genes.

The use of spike-ins as control genes is appealing, since the same amount of ERCC (or other) spike-in was added to each cell in our experiment. In principel, all the variablity we observe for these genes is due to technical noise; whereas endogenous genes are affected by both technical noise and biological variability. Technical noise can be removed by fitting a model to the spike-ins and "substracting" this from the endogenous genes. There are several methods available based on this premise (eg. [BASiCS](https://github.com/catavallejos/BASiCS), [scLVM](https://github.com/PMBio/scLVM), [RUVg](http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)); each using different noise models and different fitting procedures. Alternatively, one can identify genes which exhibit significant variation beyond technical noise (eg. Distance to median, [Highly variable genes](http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2645.html)). However, there are issues with the use of spike-ins for normalisation (particularly ERCCs, derived from bacterial sequences), including that their variability can, for various reasons, actually be *higher* than that of endogenous genes.

Given the issues with using spike-ins, better results can often be obtained by using endogenous genes instead. Where we have a large number of endogenous genes that, on average, do not vary systematically between cells and where we expect technical effects to affect a large number of genes (a very common and reasonable assumption), then such methods (for example, the RUVs method) can perform well. 

We explore both general approaches below.





```r
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater, quietly = TRUE)
library(scran)
library(edgeR)
options(stringsAsFactors = FALSE)
umi <- readRDS("blischak/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
erccs <- fData(umi.qc)$is_feature_control
```

## Remove Unwanted Variation

Factors contributing to technical noise frequently appear as "batch
effects" where cells processed on different days or by different
technicians systematically vary from one another. Removing technical
noise and correcting for batch effects can frequently be performed
using the same tool or slight variants on it. We will be considering
the [Remove Unwanted Variation (RUVSeq)](http://bioconductor.org/packages/RUVSeq). Briefly, RUVSeq works as follows. For $n$ samples and $J$ genes, consider the following generalized linear model (GLM), where the RNA-Seq read counts are regressed on both the known covariates of interest and unknown factors of unwanted variation:
\[\log E[Y|W,X,O] = W\alpha + X\beta + O\]
Here, $Y$ is the $n \times J$ matrix of observed gene-level read counts, $W$ is an $n \times k$ matrix corresponding to the factors of “unwanted variation” and $O$ is an $n \times J$ matrix of offsets that can either be set to zero or estimated with some other normalization procedure (such as upper-quartile normalization). The simultaneous estimation of $W$, $\alpha$, $\beta$, and $k$ is infeasible. For a given $k$, instead the following three
approaches to estimate the factors of unwanted variation $W$ are used:

* _RUVg_ uses negative control genes (e.g. ERCCs), assumed to have constant expression across samples;
* _RUVs_ uses centered (technical) replicate/negative control samples for which the covariates of interest are
constant;
* _RUVr_ uses residuals, e.g., from a first-pass GLM regression of the counts on the covariates of interest.

We will concentrate on the first two approaches.

### RUVg


```r
ruvg <- RUVg(counts(umi.qc), erccs, k = 1)
set_exprs(umi.qc, "ruvg1") <- ruvg$normalizedCounts
ruvg <- RUVg(counts(umi.qc), erccs, k = 2)
set_exprs(umi.qc, "ruvg2") <- ruvg$normalizedCounts
set_exprs(umi.qc, "ruvg2_logcpm") <- log2(t(t(ruvg$normalizedCounts) / 
                                           colSums(ruvg$normalizedCounts)) + 1)
```

### RUVs


```r
scIdx <- matrix(-1, ncol = max(table(umi.qc$individual)), nrow = 3)
tmp <- which(umi.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(umi.qc)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
set_exprs(umi.qc, "ruvs1") <- ruvs$normalizedCounts
ruvs <- RUVs(counts(umi.qc), cIdx, k = 2, scIdx = scIdx, isLog = FALSE)
set_exprs(umi.qc, "ruvs2") <- ruvs$normalizedCounts
set_exprs(umi.qc, "ruvs2_logcpm") <- log2(t(t(ruvs$normalizedCounts) / 
                                           colSums(ruvs$normalizedCounts)) + 1)
```

## Effectiveness 1

We evaluate the effectiveness of the normalization by inspecting the
PCA plot where colour corresponds the technical replicates and shape
corresponds to different biological samples (individuals). Separation of biological samples and
interspersed batches indicates that technical variation has been
removed. 


```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg1") +
    ggtitle("PCA - RUVg normalisation: k = 1")
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-5-1} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg2") +
    ggtitle("PCA - RUVg normalisation: k = 2")
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-5-2} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs1") +
    ggtitle("PCA - RUVs normalisation: k = 1")
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-5-3} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2") +
    ggtitle("PCA - RUVs normalisation: k = 2")
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-5-4} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2_logcpm") +
    ggtitle("PCA - RUVs normalisation log2-cpm: k = 2")
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-5-5} \end{center}

Plotting log2-normalized CPM from RUVs with k = 2 looks to give the best separation of cells by individual.

## Effectiveness 2

We can also examine the effectiveness of correction using the relative log expression (RLE) across cells to confirm technical noise has been removed from the dataset.


```r
boxplot(
    list(
        "Raw counts" = calc_cell_RLE(counts(umi.qc), erccs),
        "RUVg (k = 1)" = calc_cell_RLE(assayData(umi.qc)$ruvg1, erccs),
        "RUVg (k = 2)" = calc_cell_RLE(assayData(umi.qc)$ruvg2, erccs),
        "RUVs (k = 1)" = calc_cell_RLE(assayData(umi.qc)$ruvs1, erccs),
        "RUVs (k = 2)" = calc_cell_RLE(assayData(umi.qc)$ruvs2, erccs)
    )
)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-6-1} \end{center}

## Effectiveness 3

Another way of evaluating the effectiveness of correction is to look at the differentially expressed (DE) genes among the batches of the same individual Theoretically, these batches should not differ from each other. Let's take the most promising individual (__NA19101__, whose batches are the closest to each other) and check whether it is true.

For demonstration purposes we will only use a subset of cells. You should not do that with your real dataset, though.

```r
keep <- c(
    sample(which(umi.qc$batch == "NA19101.r1"), 20), 
    sample(which(umi.qc$batch == "NA19101.r2"), 20),
    sample(which(umi.qc$batch == "NA19101.r3"), 20)
)
design <- model.matrix(~umi.qc[, keep]$batch)
```

We will use the [edgeR](http://bioconductor.org/packages/edgeR) package to calculate DE genes between plates for this particular individual. Recall that the input data for edgeR (and similar methods like DESeq2) must always be raw counts.

The particular coefficient that we test for DE in each case below tests to for genes that show a difference in expression between replicate plate 3 and replicate plate 1.

### DE (raw counts)

```r
dge1 <- DGEList(
    counts = counts(umi.qc[, keep]), 
    norm.factors = rep(1, length(keep)),
    group = umi.qc[, keep]$batch
)
dge1 <- estimateDisp(dge1, design = design, trend.method = "none")
plotBCV(dge1)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-8-1} \end{center}

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1)
topTags(res1)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC   logCPM       LR       PValue          FDR
## ENSG00000185885 -1.4046076 8.833223 84.29007 4.272558e-20 6.008499e-16
## ENSG00000008311 -1.1922888 7.601355 48.76077 2.891662e-12 2.033272e-08
## ENSG00000115875  0.5873815 9.214540 47.13951 6.610958e-12 3.098997e-08
## ENSG00000131969 -1.2144672 7.343158 42.01740 9.046496e-11 3.180522e-07
## ENSG00000164265 -2.3287998 5.622384 39.97597 2.571066e-10 7.231380e-07
## ENSG00000110931 -1.2820663 6.602960 32.02143 1.524810e-08 3.573901e-05
## ENSG00000182718 -1.0350975 6.875030 30.94252 2.657839e-08 5.339599e-05
## ENSG00000163106 -1.0810252 7.106678 30.11108 4.079929e-08 7.172005e-05
## ENSG00000214265 -0.4331774 9.206615 27.26517 1.773779e-07 2.771629e-04
## ENSG00000125144 -2.0551927 5.891827 26.71781 2.354416e-07 3.311016e-04
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   105
## 0  13889
## 1     69
```

```r
plotSmear(
    res1, lowess = TRUE,
    de.tags = rownames(topTags(res1, n = sum(abs(decideTestsDGE(res1))))$table)
)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-8-2} \end{center}

### DE (RUVg, k = 2)

```r
design_ruvg <- model.matrix(~ruvg$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvg)
```

```
##   (Intercept) ruvg$W[keep, ]W_1 ruvg$W[keep, ]W_2
## 1           1        0.03199250       0.057216310
## 2           1        0.08126848       0.012635639
## 3           1        0.07081795      -0.005430952
## 4           1        0.02317790       0.055762151
## 5           1        0.06273906       0.046101318
## 6           1       -0.03599555       0.066168939
##   umi.qc[, keep]$batchNA19101.r2 umi.qc[, keep]$batchNA19101.r3
## 1                              0                              0
## 2                              0                              0
## 3                              0                              0
## 4                              0                              0
## 5                              0                              0
## 6                              0                              0
```

```r
dge_ruvg <- estimateDisp(dge1, design = design_ruvg, trend.method = "none")
plotBCV(dge_ruvg)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-9-1} \end{center}

```r
fit2 <- glmFit(dge_ruvg, design_ruvg)
res2 <- glmLRT(fit2)
topTags(res2)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000185885 -1.3998098  8.833058 83.69968 5.759463e-20 8.099532e-16
## ENSG00000008311 -1.1701306  7.600703 46.36870 9.796781e-12 6.888606e-08
## ENSG00000115875  0.5724239  9.215062 44.19051 2.979214e-11 1.396556e-07
## ENSG00000131969 -1.1728415  7.342296 40.36088 2.111264e-10 7.422676e-07
## ENSG00000164265 -2.3298101  5.622166 39.48067 3.313235e-10 9.318804e-07
## ENSG00000182718 -1.1521898  6.875872 38.18088 6.448115e-10 1.511331e-06
## ENSG00000125144 -2.7033767  5.892039 36.27606 1.712548e-09 3.440509e-06
## ENSG00000181163  0.1887665 11.483769 32.76984 1.037421e-08 1.776276e-05
## ENSG00000187193 -1.3350901  7.848223 32.59207 1.136776e-08 1.776276e-05
## ENSG00000110931 -1.2551852  6.602633 30.14713 4.004802e-08 5.631953e-05
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   102
## 0  13894
## 1     67
```

```r
plotSmear(
    res2, lowess = TRUE,
    de.tags = rownames(topTags(res2, n = sum(abs(decideTestsDGE(res2))))$table)
)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-9-2} \end{center}

### DE (RUVs, k = 2)

```r
design_ruvs <- model.matrix(~ruvs$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvs)
```

```
##   (Intercept) ruvs$W[keep, ]W_1 ruvs$W[keep, ]W_2
## 1           1         0.2351730       -0.07046587
## 2           1         0.2146049       -0.08745950
## 3           1         0.2267828       -0.09645385
## 4           1         0.2811607       -0.09199253
## 5           1         0.2216114       -0.07681271
## 6           1         0.3409596       -0.07033071
##   umi.qc[, keep]$batchNA19101.r2 umi.qc[, keep]$batchNA19101.r3
## 1                              0                              0
## 2                              0                              0
## 3                              0                              0
## 4                              0                              0
## 5                              0                              0
## 6                              0                              0
```

```r
dge_ruvs <- estimateDisp(dge1, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-10-1} \end{center}

```r
fit3 <- glmFit(dge_ruvs, design_ruvs)
res3 <- glmLRT(fit3)
topTags(res3)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000144713  0.5661064 10.588987 58.58420 1.947748e-14 2.739118e-10
## ENSG00000174748  0.5911278 10.633599 52.87993 3.545736e-13 2.493184e-09
## ENSG00000186468  0.3665129 11.071795 50.36227 1.278284e-12 5.992167e-09
## ENSG00000105372  0.3749875 11.571434 43.54286 4.147747e-11 1.458244e-07
## ENSG00000198918  0.6127225  9.556033 35.96502 2.008917e-09 5.650279e-06
## ENSG00000108298  0.3734569 10.700224 35.51623 2.529338e-09 5.928347e-06
## ENSG00000177954  0.5822918 11.016116 34.29092 4.745921e-09 9.534554e-06
## ENSG00000135486  0.4674854 10.012538 33.24167 8.138751e-09 1.378789e-05
## ENSG00000185885 -1.2286874  8.832730 33.08450 8.823937e-09 1.378789e-05
## ENSG00000111640  0.2875083 11.708010 31.58388 1.910089e-08 2.686158e-05
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1    65
## 0  13930
## 1     68
```

```r
plotSmear(
    res3, lowess = TRUE,
    de.tags = rownames(topTags(res3, n = sum(abs(decideTestsDGE(res3))))$table)
)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-10-2} \end{center}

In the above analyses, we have ignored size factors between cells. A typical edgeR analysis would always include these.


```r
umi.qc <- scran::computeSumFactors(umi.qc, sizes = 15)
dge_ruvs$samples$norm.factors <- sizeFactors(umi.qc)[keep]
dge_ruvs_sf <- estimateDisp(dge_ruvs, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs_sf)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-11-1} \end{center}

```r
fit4 <- glmFit(dge_ruvs_sf, design_ruvs)
res4 <- glmLRT(fit4)
topTags(res4)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000185885 -1.2391082  8.939311 33.37624 7.594632e-09 7.242602e-05
## ENSG00000144713  0.5588894 10.729895 32.56669 1.151715e-08 7.242602e-05
## ENSG00000130772  2.8829353  5.780897 31.65362 1.842712e-08 7.242602e-05
## ENSG00000177954  0.5684906 11.196333 31.43712 2.060045e-08 7.242602e-05
## ENSG00000174748  0.5833530 10.773288 30.32243 3.658708e-08 1.029048e-04
## ENSG00000198918  0.6090866  9.672585 28.61142 8.845912e-08 2.049680e-04
## ENSG00000187193 -1.6550068  7.866237 28.33518 1.020249e-07 2.049680e-04
## ENSG00000131969 -1.2676671  7.386214 26.33872 2.864867e-07 5.036079e-04
## ENSG00000115875  0.6284983  9.304699 25.65181 4.089109e-07 6.389461e-04
## ENSG00000205609  1.5836661  6.866997 24.76294 6.483206e-07 9.117333e-04
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1    47
## 0  13974
## 1     42
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{14-remove-conf_files/figure-latex/unnamed-chunk-11-2} \end{center}


## Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter). Additionally, experiment with other combinations of normalizations and compare the results.
