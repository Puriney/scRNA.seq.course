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

<img src="14-remove-conf_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg2") +
    ggtitle("PCA - RUVg normalisation: k = 2")
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-5-2.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs1") +
    ggtitle("PCA - RUVs normalisation: k = 1")
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-5-3.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2") +
    ggtitle("PCA - RUVs normalisation: k = 2")
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-5-4.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2_logcpm") +
    ggtitle("PCA - RUVs normalisation log2-cpm: k = 2")
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-5-5.png" width="672" style="display: block; margin: auto;" />

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

<img src="14-remove-conf_files/figure-html/unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

## Effectiveness 3

Another way of evaluating the effectiveness of correction is to look at the differentially expressed (DE) genes among the batches of the same individual Theoretically, these batches should not differ from each other. Let's take the most promising individual (__NA19101__, whose batches are the closest to each other) and check whether it is true.


```r
keep <- umi.qc$individual == "NA19101"
design <- model.matrix(~umi.qc[, keep]$batch)
```

We will use the [edgeR](http://bioconductor.org/packages/edgeR) package to calculate DE genes between plates for this particular individual. Recall that the input data for edgeR (and similar methods like DESeq2) must always be raw counts.

The particular coefficient that we test for DE in each case below tests to for genes that show a difference in expression between replicate plate 3 and replicate plate 1.

### DE (raw counts)

```r
dge1 <- DGEList(
    counts = counts(umi.qc[, keep]), 
    norm.factors = rep(1, sum(keep)),
    group = umi.qc[, keep]$batch
)
dge1 <- estimateDisp(dge1, design = design, trend.method = "none")
plotBCV(dge1)
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1, coef = 2)
topTags(res1)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r2 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000145423  1.0437082  8.930360 354.3682 4.741536e-79 6.668021e-75
## ENSG00000185885 -1.4704261  8.871862 312.6361 5.820571e-70 4.092734e-66
## ENSG00000104332  0.7942726  9.732365 284.9588 6.237725e-64 2.924037e-60
## ENSG00000214265 -0.7737881  9.237431 282.8562 1.791432e-63 6.298226e-60
## ENSG00000106554 -0.6836208  9.262132 241.4285 1.919657e-54 5.399226e-51
## ENSG00000167283 -0.4769965 10.083626 227.0334 2.644427e-51 6.198096e-48
## ENSG00000165502 -0.5248241 10.028994 224.9492 7.531484e-51 1.331865e-47
## ENSG00000087460  0.5510932  9.392753 224.9374 7.576560e-51 1.331865e-47
## ENSG00000121207  1.9353719  6.736511 211.5649 6.258142e-48 9.778694e-45
## ENSG00000112695 -0.5505989  9.637643 210.3241 1.167217e-47 1.641458e-44
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1]
## -1 2109
## 0  9355
## 1  2599
```

```r
plotSmear(
    res1, lowess = TRUE,
    de.tags = rownames(topTags(res1, n = sum(abs(decideTestsDGE(res1))))$table)
)
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-8-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVg, k = 2)

```r
design_ruvg <- model.matrix(~ruvg$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvg)
```

```
##   (Intercept) ruvg$W[keep, ]W_1 ruvg$W[keep, ]W_2
## 1           1        0.01010943       0.001956012
## 2           1        0.01777997       0.009123186
## 3           1        0.01968374       0.006157921
## 4           1        0.01275973       0.017130060
## 5           1        0.07081795      -0.005430952
## 6           1        0.02911818       0.021727805
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

<img src="14-remove-conf_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

```r
fit2 <- glmFit(dge_ruvg, design_ruvg)
res2 <- glmLRT(fit2)
topTags(res2)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC   logCPM        LR       PValue          FDR
## ENSG00000185885 -1.3174364 8.871753 217.14670 3.791218e-49 5.331591e-45
## ENSG00000163106 -1.3722714 7.130682 128.37173 9.307179e-30 6.544343e-26
## ENSG00000131969 -1.2168148 7.359236 118.81702 1.148503e-27 5.383799e-24
## ENSG00000008311 -1.1845327 7.534209 111.42347 4.779029e-26 1.680187e-22
## ENSG00000187193 -1.4168162 7.660033  86.83180 1.181542e-20 3.323205e-17
## ENSG00000145423  0.5749357 8.929498  85.46301 2.360764e-20 5.533238e-17
## ENSG00000110931 -1.2224083 6.553645  80.05106 3.648582e-19 7.330001e-16
## ENSG00000196683 -0.4017256 9.555568  77.23020 1.521490e-18 2.674589e-15
## ENSG00000214265 -0.4006635 9.237080  71.71226 2.489812e-17 3.890469e-14
## ENSG00000150459  0.5101314 8.387922  68.60831 1.200971e-16 1.688925e-13
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1  1039
## 0  12237
## 1    787
```

```r
plotSmear(
    res2, lowess = TRUE,
    de.tags = rownames(topTags(res2, n = sum(abs(decideTestsDGE(res2))))$table)
)
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-9-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVs, k = 2)

```r
design_ruvs <- model.matrix(~ruvs$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvs)
```

```
##   (Intercept) ruvs$W[keep, ]W_1 ruvs$W[keep, ]W_2
## 1           1         0.2619953       -0.08220987
## 2           1         0.2695799       -0.09943721
## 3           1         0.2152265       -0.09724690
## 4           1         0.2642097       -0.10269732
## 5           1         0.2267828       -0.09645385
## 6           1         0.2910819       -0.09367850
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

<img src="14-remove-conf_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

```r
fit3 <- glmFit(dge_ruvs, design_ruvs)
res3 <- glmLRT(fit3)
topTags(res3)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM        LR       PValue          FDR
## ENSG00000144713  0.4128927 10.580211 126.45399 2.446079e-29 3.439921e-25
## ENSG00000137818  0.2964018 11.606043 122.06746 2.231214e-28 1.568878e-24
## ENSG00000105372  0.3266439 11.555782 116.13967 4.429746e-27 2.076517e-23
## ENSG00000185885 -1.0832910  8.871543 107.26639 3.891781e-25 1.368253e-21
## ENSG00000162244  0.4920933  9.846306  97.28089 6.015816e-23 1.692008e-19
## ENSG00000181163  0.2865645 11.478855  94.62337 2.302830e-22 5.397449e-19
## ENSG00000174748  0.3857506 10.667147  89.43624 3.166816e-21 6.362134e-18
## ENSG00000150459  0.6417026  8.391754  85.37323 2.470423e-20 4.342695e-17
## ENSG00000240972  0.5080671 10.390759  84.80951 3.285399e-20 5.133619e-17
## ENSG00000177954  0.4557659 11.015973  84.12838 4.636678e-20 6.520560e-17
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   876
## 0  12609
## 1    578
```

```r
plotSmear(
    res3, lowess = TRUE,
    de.tags = rownames(topTags(res3, n = sum(abs(decideTestsDGE(res3))))$table)
)
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-10-2.png" width="672" style="display: block; margin: auto;" />

In the above analyses, we have ignored size factors between cells. A typical edgeR analysis would always include these.


```r
dge_ruvs$samples$norm.factors <- sizeFactors(umi.qc)[keep]
```

```
## Warning in .local(object, ...): 'sizeFactors' have not been set
```

```r
dge_ruvs_sf <- estimateDisp(dge_ruvs, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs_sf)
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

```r
fit4 <- glmFit(dge_ruvs_sf, design_ruvs)
res4 <- glmLRT(fit4)
topTags(res4)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM        LR       PValue          FDR
## ENSG00000144713  0.4128927 10.580211 126.45399 2.446079e-29 3.439921e-25
## ENSG00000137818  0.2964018 11.606043 122.06746 2.231214e-28 1.568878e-24
## ENSG00000105372  0.3266439 11.555782 116.13967 4.429746e-27 2.076517e-23
## ENSG00000185885 -1.0832910  8.871543 107.26639 3.891781e-25 1.368253e-21
## ENSG00000162244  0.4920933  9.846306  97.28089 6.015816e-23 1.692008e-19
## ENSG00000181163  0.2865645 11.478855  94.62337 2.302830e-22 5.397449e-19
## ENSG00000174748  0.3857506 10.667147  89.43624 3.166816e-21 6.362134e-18
## ENSG00000150459  0.6417026  8.391754  85.37323 2.470423e-20 4.342695e-17
## ENSG00000240972  0.5080671 10.390759  84.80951 3.285399e-20 5.133619e-17
## ENSG00000177954  0.4557659 11.015973  84.12838 4.636678e-20 6.520560e-17
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   876
## 0  12609
## 1    578
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="14-remove-conf_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />


## Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter). Additionally, experiment with other combinations of normalizations and compare the results.
