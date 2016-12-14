---
output: html_document
---

# Dealing with confounders (Reads)




```r
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater, quietly = TRUE)
library(scran)
library(edgeR)
options(stringsAsFactors = FALSE)
reads <- readRDS("blischak/reads.rds")
reads.qc <- reads[fData(reads)$use, pData(reads)$use]
endog_genes <- !fData(reads.qc)$is_feature_control
erccs <- fData(reads.qc)$is_feature_control
```

## Remove Unwanted Variation

### RUVg


```r
ruvg <- RUVg(counts(reads.qc), erccs, k = 1)
set_exprs(reads.qc, "ruvg1") <- ruvg$normalizedCounts
ruvg <- RUVg(counts(reads.qc), erccs, k = 2)
set_exprs(reads.qc, "ruvg2") <- ruvg$normalizedCounts
set_exprs(reads.qc, "ruvg2_logcpm") <- log2(t(t(ruvg$normalizedCounts) / 
                                           colSums(ruvg$normalizedCounts)) + 1)
```

### RUVs


```r
scIdx <- matrix(-1, ncol = max(table(reads.qc$individual)), nrow = 3)
tmp <- which(reads.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(reads.qc)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
set_exprs(reads.qc, "ruvs1") <- ruvs$normalizedCounts
ruvs <- RUVs(counts(reads.qc), cIdx, k = 2, scIdx = scIdx, isLog = FALSE)
set_exprs(reads.qc, "ruvs2") <- ruvs$normalizedCounts
set_exprs(reads.qc, "ruvs2_logcpm") <- log2(t(t(ruvs$normalizedCounts) / 
                                           colSums(ruvs$normalizedCounts)) + 1)
```

## Effectiveness 1


```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg1") +
    ggtitle("PCA - RUVg normalisation: k = 1")
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-5-1} \end{center}

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg2") +
    ggtitle("PCA - RUVg normalisation: k = 2")
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-5-2} \end{center}

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs1") +
    ggtitle("PCA - RUVs normalisation: k = 1")
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-5-3} \end{center}

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2") +
    ggtitle("PCA - RUVs normalisation: k = 2")
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-5-4} \end{center}

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2_logcpm") +
    ggtitle("PCA - RUVs normalisation log2-cpm: k = 2")
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-5-5} \end{center}

## Effectiveness 2


```r
boxplot(
    list(
        "Raw counts" = calc_cell_RLE(counts(reads.qc), erccs),
        "RUVg (k = 1)" = calc_cell_RLE(assayData(reads.qc)$ruvg1, erccs),
        "RUVg (k = 2)" = calc_cell_RLE(assayData(reads.qc)$ruvg2, erccs),
        "RUVs (k = 1)" = calc_cell_RLE(assayData(reads.qc)$ruvs1, erccs),
        "RUVs (k = 2)" = calc_cell_RLE(assayData(reads.qc)$ruvs2, erccs)
    )
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-6-1} \end{center}

## Effectiveness 3


```r
keep <- c(
    sample(which(reads.qc$batch == "NA19101.r1"), 20), 
    sample(which(reads.qc$batch == "NA19101.r2"), 20),
    sample(which(reads.qc$batch == "NA19101.r3"), 20)
)
design <- model.matrix(~reads.qc[, keep]$batch)
```

### DE (raw counts)

```r
dge1 <- DGEList(
    counts = counts(reads.qc[, keep]), 
    norm.factors = rep(1, length(keep)),
    group = reads.qc[, keep]$batch
)
dge1 <- estimateDisp(dge1, design = design, trend.method = "none")
plotBCV(dge1)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-8-1} \end{center}

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1)
topTags(res1)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000205358 -8.361776 3.2162232 30.69522 3.019083e-08 0.0004848949
## ENSG00000078967  7.121784 2.1408893 26.53017 2.594542e-07 0.0014555353
## ENSG00000117408  6.169546 3.3722602 26.43983 2.718763e-07 0.0014555353
## ENSG00000146250  7.760828 2.5044911 24.75743 6.501769e-07 0.0016948305
## ENSG00000180938  7.373991 2.3899136 24.72537 6.610812e-07 0.0016948305
## ENSG00000188816  6.879539 1.6848100 24.70188 6.691874e-07 0.0016948305
## ENSG00000182253 -6.800201 1.5366985 24.51151 7.386722e-07 0.0016948305
## ENSG00000169247  6.937716 1.4366019 23.53462 1.226865e-06 0.0018231751
## ENSG00000064655  5.812560 0.7133747 23.45633 1.277818e-06 0.0018231751
## ENSG00000188133  6.829378 1.3563128 23.43198 1.294094e-06 0.0018231751
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   409
## 0  14855
## 1    797
```

```r
plotSmear(
    res1, lowess = TRUE,
    de.tags = rownames(topTags(res1, n = sum(abs(decideTestsDGE(res1))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-8-2} \end{center}

### DE (RUVg, k = 2)

```r
design_ruvg <- model.matrix(~ruvg$W[keep,] + reads.qc[, keep]$batch)
head(design_ruvg)
```

```
##   (Intercept) ruvg$W[keep, ]W_1 ruvg$W[keep, ]W_2
## 1           1      0.0277947090       -0.03228469
## 2           1     -0.0004362023       -0.06051295
## 3           1     -0.0154461036        0.04957296
## 4           1     -0.0105911292       -0.01936104
## 5           1      0.0793184051        0.00186142
## 6           1      0.0130598159        0.02360856
##   reads.qc[, keep]$batchNA19101.r2 reads.qc[, keep]$batchNA19101.r3
## 1                                0                                0
## 2                                0                                0
## 3                                0                                0
## 4                                0                                0
## 5                                0                                0
## 6                                0                                0
```

```r
dge_ruvg <- estimateDisp(dge1, design = design_ruvg, trend.method = "none")
plotBCV(dge_ruvg)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-9-1} \end{center}

```r
fit2 <- glmFit(dge_ruvg, design_ruvg)
res2 <- glmLRT(fit2)
topTags(res2)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000176236 -5.975226 1.0323168 43.78815 3.659133e-11 4.728708e-07
## ENSG00000165309 -5.896065 0.8655683 42.30761 7.798848e-11 4.728708e-07
## ENSG00000103888 -5.858395 0.8421441 42.06418 8.832654e-11 4.728708e-07
## ENSG00000251201 -5.589593 0.7651133 41.11874 1.432564e-10 5.752103e-07
## ENSG00000196735 -5.571323 0.6718151 40.20068 2.291685e-10 7.361350e-07
## ENSG00000116819 -4.525598 0.2020634 33.67791 6.503513e-09 1.740882e-05
## ENSG00000184786 -5.397420 0.4916194 31.38905 2.111689e-08 4.845119e-05
## ENSG00000156219 -6.600050 1.4008709 30.52734 3.291940e-08 6.608981e-05
## ENSG00000159182 -6.573696 1.0884067 30.12060 4.059953e-08 7.245211e-05
## ENSG00000078967  7.755221 2.1411162 29.77576 4.850178e-08 7.789871e-05
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   383
## 0  15291
## 1    387
```

```r
plotSmear(
    res2, lowess = TRUE,
    de.tags = rownames(topTags(res2, n = sum(abs(decideTestsDGE(res2))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-9-2} \end{center}

### DE (RUVs, k = 2)

```r
design_ruvs <- model.matrix(~ruvs$W[keep,] + reads.qc[, keep]$batch)
head(design_ruvs)
```

```
##   (Intercept) ruvs$W[keep, ]W_1 ruvs$W[keep, ]W_2
## 1           1         0.2769722         0.1985588
## 2           1         0.3111904         0.2082694
## 3           1         0.3093128         0.1938905
## 4           1         0.3176870         0.2029780
## 5           1         0.3954771         0.1603333
## 6           1         0.2873899         0.1971467
##   reads.qc[, keep]$batchNA19101.r2 reads.qc[, keep]$batchNA19101.r3
## 1                                0                                0
## 2                                0                                0
## 3                                0                                0
## 4                                0                                0
## 5                                0                                0
## 6                                0                                0
```

```r
dge_ruvs <- estimateDisp(dge1, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-10-1} \end{center}

```r
fit3 <- glmFit(dge_ruvs, design_ruvs)
res3 <- glmLRT(fit3)
topTags(res3)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000205358 -9.629375 3.2163112 39.74163 2.898796e-10 2.018330e-06
## ENSG00000127533 -6.073391 0.5559786 39.35760 3.528771e-10 2.018330e-06
## ENSG00000165309 -6.567050 0.8656534 38.86540 4.540577e-10 2.018330e-06
## ENSG00000103888 -6.527453 0.8422293 38.66688 5.026661e-10 2.018330e-06
## ENSG00000196735 -6.224409 0.6719003 37.13424 1.102704e-09 3.260125e-06
## ENSG00000131126 -6.738550 0.9402994 36.94049 1.217904e-09 3.260125e-06
## ENSG00000102554 -5.981232 0.5055222 34.74106 3.766037e-09 8.450030e-06
## ENSG00000156219 -7.268864 1.4008801 34.52461 4.208968e-09 8.450030e-06
## ENSG00000088386 -6.675796 0.9136309 33.48094 7.196572e-09 1.284268e-05
## ENSG00000104432 -5.914827 0.4962049 32.94322 9.489026e-09 1.524032e-05
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   387
## 0  15324
## 1    350
```

```r
plotSmear(
    res3, lowess = TRUE,
    de.tags = rownames(topTags(res3, n = sum(abs(decideTestsDGE(res3))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-10-2} \end{center}


```r
reads.qc <- scran::computeSumFactors(reads.qc, sizes = 15)
dge_ruvs$samples$norm.factors <- sizeFactors(reads.qc)[keep]
dge_ruvs_sf <- estimateDisp(dge_ruvs, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs_sf)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-11-1} \end{center}

```r
fit4 <- glmFit(dge_ruvs_sf, design_ruvs)
res4 <- glmLRT(fit4)
topTags(res4)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000205358 -9.575893 3.1948851 38.37163 5.847636e-10 4.953775e-06
## ENSG00000127533 -5.974917 0.4179778 37.55042 8.908033e-10 4.953775e-06
## ENSG00000165309 -6.479654 0.7214685 37.11317 1.114681e-09 4.953775e-06
## ENSG00000103888 -6.439968 0.6984075 36.91530 1.233740e-09 4.953775e-06
## ENSG00000196735 -6.136619 0.5307670 35.38790 2.701603e-09 8.629408e-06
## ENSG00000131126 -6.637499 0.7991928 35.04379 3.223738e-09 8.629408e-06
## ENSG00000156219 -7.302817 1.3229713 33.83764 5.990853e-09 1.374558e-05
## ENSG00000102554 -5.789419 0.3014686 32.46313 1.214756e-08 2.438775e-05
## ENSG00000088386 -6.635257 0.7725839 31.86489 1.652778e-08 2.949473e-05
## ENSG00000216937  7.221752 1.8483497 31.61750 1.877304e-08 3.015138e-05
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   335
## 0  15413
## 1    313
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-11-2} \end{center}
