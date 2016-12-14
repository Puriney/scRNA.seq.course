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
## ENSG00000138650 -7.441348 2.1921251 29.20068 6.525634e-08 0.0006672509
## ENSG00000160307  7.418664 1.7999768 28.02863 1.195337e-07 0.0006672509
## ENSG00000167555 -7.432707 1.8114463 27.94777 1.246344e-07 0.0006672509
## ENSG00000110195  7.344885 1.7616132 27.03319 1.999917e-07 0.0006891312
## ENSG00000186439 -6.329585 1.6129275 26.89751 2.145356e-07 0.0006891312
## ENSG00000105877 -7.260839 1.6735360 26.46016 2.690294e-07 0.0007032305
## ENSG00000196547 -7.631320 2.6455310 26.02277 3.374147e-07 0.0007032305
## ENSG00000151014 -5.953574 0.8058841 25.95053 3.502798e-07 0.0007032305
## ENSG00000163106 -3.097488 4.2586338 25.54573 4.320204e-07 0.0007709645
## ENSG00000196109 -7.557884 2.2977156 25.29799 4.912203e-07 0.0007889490
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   559
## 0  14868
## 1    634
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
## 1           1      -0.009896783        0.01987111
## 2           1       0.011533972       -0.06066909
## 3           1       0.018149399       -0.06327037
## 4           1      -0.002228917        0.03877421
## 5           1       0.034174399       -0.01684021
## 6           1      -0.010043808        0.07328756
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
## ENSG00000165309 -6.289280 0.8695936 45.88979 1.250965e-11 8.855225e-08
## ENSG00000235531 -6.521231 1.2393856 45.60397 1.447477e-11 8.855225e-08
## ENSG00000203392 -6.208378 0.8224942 45.34268 1.654049e-11 8.855225e-08
## ENSG00000023892 -6.879285 1.4435444 42.48057 7.138799e-11 2.738491e-07
## ENSG00000137252 -5.732205 0.5662601 42.09698 8.685704e-11 2.738491e-07
## ENSG00000127533 -5.684170 0.5431679 41.77693 1.023034e-10 2.738491e-07
## ENSG00000104432 -5.573076 0.4898047 41.01552 1.510255e-10 3.465172e-07
## ENSG00000054179 -7.324881 1.5139814 31.91086 1.614125e-08 3.240557e-05
## ENSG00000186439 -6.864318 1.6130818 31.58453 1.909446e-08 3.407512e-05
## ENSG00000088386 -6.256865 0.9016214 30.71564 2.987481e-08 4.798192e-05
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   284
## 0  15558
## 1    219
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
## 1           1         0.3176779         0.2226902
## 2           1         0.3775191         0.1868876
## 3           1         0.2702873         0.2169784
## 4           1         0.3444123         0.2143000
## 5           1         0.2864326         0.1671973
## 6           1         0.3290784         0.1965167
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
## ENSG00000162552  5.919575 0.8798908 31.52689 1.966972e-08 0.0001996355
## ENSG00000112541  7.388383 1.5779887 29.68909 5.071948e-08 0.0001996355
## ENSG00000196547 -8.345312 2.6457262 29.07009 6.980646e-08 0.0001996355
## ENSG00000183570 -6.690838 1.4109883 29.04321 7.078177e-08 0.0001996355
## ENSG00000151834  6.027799 0.7249389 28.99626 7.251808e-08 0.0001996355
## ENSG00000186439 -6.543542 1.6131270 28.90302 7.609428e-08 0.0001996355
## ENSG00000082497  7.652492 1.5279419 28.37037 1.001870e-07 0.0001996355
## ENSG00000038295  6.235536 0.8592432 28.29882 1.039592e-07 0.0001996355
## ENSG00000160307  7.112246 1.7996327 28.15689 1.118685e-07 0.0001996355
## ENSG00000164100  5.008404 0.3684580 27.34545 1.701645e-07 0.0002733012
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   268
## 0  15503
## 1    290
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
## ENSG00000162552  6.103844 0.9299905 31.51496 1.979097e-08 0.0003178628
## ENSG00000186439 -6.671444 1.6315203 28.97744 7.322596e-08 0.0003699401
## ENSG00000196547 -8.254273 2.2825598 28.16871 1.111877e-07 0.0003699401
## ENSG00000183570 -6.710140 1.2137842 27.99022 1.219298e-07 0.0003699401
## ENSG00000082497  7.809497 1.5300307 27.71198 1.407882e-07 0.0003699401
## ENSG00000038295  6.279808 0.7446956 27.58062 1.506812e-07 0.0003699401
## ENSG00000164100  5.190206 0.3546801 27.44970 1.612341e-07 0.0003699401
## ENSG00000160307  7.364708 2.0307282 26.81481 2.239145e-07 0.0004241947
## ENSG00000112541  7.123755 1.3614677 26.69934 2.377033e-07 0.0004241947
## ENSG00000151834  5.702388 0.4254371 26.36902 2.820276e-07 0.0004529646
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   237
## 0  15557
## 1    267
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-11-2} \end{center}
