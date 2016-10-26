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

<img src="14-remove-conf_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1)
topTags(res1)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC   logCPM       LR       PValue          FDR
## ENSG00000185885 -1.3628258 8.854939 90.96117 1.465191e-21 2.060498e-17
## ENSG00000187193 -1.7554950 7.435729 73.36714 1.076432e-17 7.568932e-14
## ENSG00000131969 -1.5630889 7.325064 70.76752 4.018986e-17 1.883967e-13
## ENSG00000008311 -1.4596139 7.525723 67.83418 1.778397e-16 6.252401e-13
## ENSG00000163106 -1.6514771 7.073590 66.85171 2.927179e-16 8.232983e-13
## ENSG00000125144 -3.1737542 5.511178 52.73695 3.813490e-13 8.938184e-10
## ENSG00000198417 -3.0510227 5.393091 46.10135 1.122912e-11 2.255930e-08
## ENSG00000182463 -3.5761474 5.271596 42.14259 8.485470e-11 1.491639e-07
## ENSG00000062716 -0.8674567 8.705770 36.25206 1.733770e-09 2.709112e-06
## ENSG00000164265 -1.8660121 5.730975 33.59121 6.799932e-09 9.562744e-06
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   146
## 0  13826
## 1     91
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
## 1           1       0.027831228        0.04542965
## 2           1       0.008323015        0.03133616
## 3           1       0.018680097        0.03507238
## 4           1       0.045487344       -0.01318167
## 5           1       0.058780915        0.02744029
## 6           1       0.010699708       -0.01163183
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
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000185885 -1.4985087  8.855193 92.68745 6.124077e-22 8.612289e-18
## ENSG00000163106 -1.7669349  7.072825 58.66466 1.869696e-14 1.314677e-10
## ENSG00000131969 -1.4980658  7.323562 53.41436 2.701145e-13 1.266207e-09
## ENSG00000008311 -1.4064186  7.523867 51.39596 7.549359e-13 2.654166e-09
## ENSG00000187193 -1.4770093  7.434192 46.27074 1.029907e-11 2.896717e-08
## ENSG00000214517 -1.9782163  6.094151 35.30279 2.822276e-09 6.614945e-06
## ENSG00000198417 -2.8935373  5.392688 30.79625 2.865905e-08 5.757603e-05
## ENSG00000125144 -2.4727649  5.510436 28.73281 8.308400e-08 1.460513e-04
## ENSG00000214756 -0.8880076  7.494998 27.10534 1.926654e-07 3.010503e-04
## ENSG00000177954  0.3692837 11.054561 25.98758 3.436204e-07 4.693717e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1    98
## 0  13904
## 1     61
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
## 1           1         0.2447665       -0.10344432
## 2           1         0.2825106       -0.09287973
## 3           1         0.2134979       -0.07575268
## 4           1         0.2516191       -0.08994524
## 5           1         0.2228644       -0.09824831
## 6           1         0.2688981       -0.08130043
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
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000185966 5.9898740  5.046981 49.87011 1.642679e-12 2.310099e-08
## ENSG00000137818 0.3747621 11.606224 44.67791 2.322632e-11 1.633159e-07
## ENSG00000162244 0.6630424  9.816618 37.84494 7.659661e-10 3.590594e-06
## ENSG00000089157 0.3447326 11.186372 34.91426 3.445481e-09 1.211345e-05
## ENSG00000171863 0.4162921 10.566553 34.11631 5.191425e-09 1.395146e-05
## ENSG00000137154 0.3967160 11.179666 33.85016 5.952411e-09 1.395146e-05
## ENSG00000144713 0.4953033 10.558270 32.81459 1.013808e-08 2.036741e-05
## ENSG00000114391 0.5633371 10.161511 32.01211 1.532148e-08 2.693325e-05
## ENSG00000105372 0.3565580 11.538046 30.24876 3.800347e-08 5.458360e-05
## ENSG00000186468 0.3495589 11.075280 30.20785 3.881363e-08 5.458360e-05
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1    69
## 0  13930
## 1     64
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
umi.qc <- scran::computeSumFactors(umi.qc, sizes = 15)
dge_ruvs$samples$norm.factors <- sizeFactors(umi.qc)[keep]
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
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000185966  6.3505721  4.909325 50.20968 1.381644e-12 1.943006e-08
## ENSG00000162244  0.7268032  9.915924 32.47344 1.208326e-08 8.496344e-05
## ENSG00000114391  0.6377807 10.268325 27.85192 1.309638e-07 6.139146e-04
## ENSG00000137818  0.4259350 11.783427 26.71753 2.354749e-07 7.963257e-04
## ENSG00000124193  1.2464691  7.457052 26.06622 3.299059e-07 7.963257e-04
## ENSG00000240972  0.5872975 10.533875 26.00943 3.397536e-07 7.963257e-04
## ENSG00000185340  4.4725489  4.820466 24.06529 9.312368e-07 1.870855e-03
## ENSG00000163584  0.7273963  9.145515 23.27753 1.402289e-06 2.382867e-03
## ENSG00000198417 -3.5047887  5.274102 23.11623 1.524980e-06 2.382867e-03
## ENSG00000050405 -1.5737371  6.650556 22.88325 1.721452e-06 2.420877e-03
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1    49
## 0  13965
## 1     49
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
