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
## ENSG00000157014  6.663541 2.8913192 29.51408 5.551189e-08 0.0006367981
## ENSG00000145990  6.055787 2.8277740 28.57047 9.034971e-08 0.0006367981
## ENSG00000136237  6.935771 3.0702270 28.03817 1.189462e-07 0.0006367981
## ENSG00000152977 -8.282698 2.4359500 24.26421 8.398518e-07 0.0023264609
## ENSG00000131620  7.695002 2.3660074 23.92329 1.002518e-06 0.0023264609
## ENSG00000138642  6.683492 2.1868406 23.31680 1.373950e-06 0.0023264609
## ENSG00000168505 -6.964629 1.3731528 23.31095 1.378134e-06 0.0023264609
## ENSG00000186470 -6.128681 0.8351887 23.21172 1.451104e-06 0.0023264609
## ENSG00000122390 -6.538248 1.0818092 23.13374 1.511160e-06 0.0023264609
## ENSG00000148219 -5.544937 2.9610611 22.96420 1.650468e-06 0.0023264609
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   557
## 0  14875
## 1    629
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
## 1           1        0.02696928       0.014887328
## 2           1       -0.01059113      -0.019361045
## 3           1        0.05102047      -0.008166728
## 4           1        0.04848977      -0.007775879
## 5           1        0.02887529      -0.021003497
## 6           1        0.05401540      -0.031853433
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
##                     logFC      logCPM       LR       PValue          FDR
## ENSG00000105383 -3.897266  0.24192772 34.87068 3.523477e-09 2.829528e-05
## ENSG00000106384 -3.897266  0.24192772 34.87068 3.523477e-09 2.829528e-05
## ENSG00000157014  7.745171  2.89137782 33.07371 8.873040e-09 3.621373e-05
## ENSG00000010438 -3.658437  0.12149131 33.04198 9.019048e-09 3.621373e-05
## ENSG00000145990  6.336950  2.82773426 31.02203 2.551157e-08 8.194827e-05
## ENSG00000023445 -3.270633 -0.03953745 30.05149 4.207256e-08 1.126212e-04
## ENSG00000102174 -4.346068  0.46991113 28.25493 1.063437e-07 2.217577e-04
## ENSG00000158869 -3.023714 -0.12193224 28.12880 1.135039e-07 2.217577e-04
## ENSG00000257242 -2.988503 -0.13256581 27.85206 1.309544e-07 2.217577e-04
## ENSG00000102466 -2.951876 -0.14327721 27.56530 1.518793e-07 2.217577e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   520
## 0  15303
## 1    238
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
## 1           1         0.3482541         0.1820908
## 2           1         0.3176870         0.2029780
## 3           1         0.3904875         0.1400878
## 4           1         0.3147846         0.1785769
## 5           1         0.3689274         0.2143278
## 6           1         0.4046564         0.1473477
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
## ENSG00000157014  7.181485 2.8913972 37.64377 8.491768e-10 1.363863e-05
## ENSG00000145990  6.256429 2.8277211 30.58064 3.202735e-08 2.571956e-04
## ENSG00000196735 -5.495057 0.6391009 26.97599 2.059980e-07 1.102845e-03
## ENSG00000041982 -7.797254 2.0526162 26.34455 2.856238e-07 1.146851e-03
## ENSG00000026103  8.629767 2.4411156 24.89215 6.062871e-07 1.409198e-03
## ENSG00000136237  6.701656 3.0701058 24.84897 6.200218e-07 1.409198e-03
## ENSG00000154027  7.324371 1.6272007 24.84207 6.222459e-07 1.409198e-03
## ENSG00000127533 -5.236535 0.4966899 24.29101 8.282489e-07 1.409198e-03
## ENSG00000148219 -5.676093 2.9609983 24.12254 9.039610e-07 1.409198e-03
## ENSG00000131126 -5.960863 0.9058634 23.91884 1.004835e-06 1.409198e-03
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   432
## 0  15240
## 1    389
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
## ENSG00000157014  7.350360 2.7507839 37.91091 7.404989e-10 1.189315e-05
## ENSG00000145990  6.453002 2.8425562 30.90213 2.713733e-08 2.179263e-04
## ENSG00000196735 -5.475502 0.4440058 26.28267 2.949238e-07 1.054656e-03
## ENSG00000136237  6.997524 3.1793777 26.07404 3.285716e-07 1.054656e-03
## ENSG00000041982 -7.776020 1.7571825 25.86451 3.662404e-07 1.054656e-03
## ENSG00000026103  9.015648 2.3826895 25.52151 4.374779e-07 1.054656e-03
## ENSG00000154027  7.859342 1.5041142 25.42609 4.596595e-07 1.054656e-03
## ENSG00000122863  6.292541 0.7141398 24.42574 7.722981e-07 1.345959e-03
## ENSG00000138642  6.873626 2.2087605 24.31151 8.194783e-07 1.345959e-03
## ENSG00000131620  7.887993 2.1830325 24.26839 8.380291e-07 1.345959e-03
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   371
## 0  15272
## 1    418
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-11-2} \end{center}
