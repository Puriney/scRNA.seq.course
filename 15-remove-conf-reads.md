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
##                     logFC   logCPM       LR       PValue          FDR
## ENSG00000138650 -7.678806 2.003747 33.84554 5.966578e-09 9.582922e-05
## ENSG00000109743  7.667826 1.989344 31.97871 1.558716e-08 1.251727e-04
## ENSG00000044446  6.087456 3.657979 29.89003 4.572581e-08 2.448007e-04
## ENSG00000160307  6.768762 1.307775 27.95346 1.242684e-07 4.978484e-04
## ENSG00000146021  5.784734 3.135279 27.52613 1.549867e-07 4.978484e-04
## ENSG00000112137  6.132768 1.292363 27.14609 1.886459e-07 5.049737e-04
## ENSG00000214188  6.988412 1.485779 25.73225 3.922183e-07 8.999170e-04
## ENSG00000171551  6.742559 1.289356 25.12943 5.360833e-07 9.988565e-04
## ENSG00000171208  6.307105 3.643464 25.04622 5.597228e-07 9.988565e-04
## ENSG00000131849  6.618593 1.204518 23.98708 9.698418e-07 1.415977e-03
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   441
## 0  14859
## 1    761
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
## 1           1       0.016561929      -0.008540817
## 2           1      -0.025856039       0.027531992
## 3           1      -0.009896783       0.019871106
## 4           1       0.125910579       0.007593017
## 5           1      -0.006315637       0.035765872
## 6           1      -0.049085371      -0.010073119
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
## ENSG00000138650 -7.710654 2.0040978 33.39202 7.533247e-09 0.0001209915
## ENSG00000146021  6.307282 3.1352154 31.77292 1.732924e-08 0.0001391625
## ENSG00000044446  6.303099 3.6579550 30.61585 3.145133e-08 0.0001683799
## ENSG00000163219 -4.822317 0.8585435 26.69957 2.376748e-07 0.0007734088
## ENSG00000160307  6.512451 1.3076718 26.48642 2.653971e-07 0.0007734088
## ENSG00000112137  6.072095 1.2925307 26.32234 2.889267e-07 0.0007734088
## ENSG00000189134 -4.366128 0.5390822 25.74219 3.902025e-07 0.0008764743
## ENSG00000145431  7.026457 2.1802636 25.52551 4.365727e-07 0.0008764743
## ENSG00000169403  8.557565 2.2185715 25.05101 5.583337e-07 0.0009963775
## ENSG00000109743  6.501223 1.9894723 24.72481 6.612748e-07 0.0010620735
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   371
## 0  15345
## 1    345
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
## 1           1         0.3393820         0.2006604
## 2           1         0.2879353         0.1807737
## 3           1         0.3176779         0.2226902
## 4           1         0.3775897         0.1547174
## 5           1         0.3144554         0.1941110
## 6           1         0.2172195         0.1480125
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
## ENSG00000138650 -7.822534 2.0041784 35.62401 2.393199e-09 3.843717e-05
## ENSG00000144962 -7.757015 1.8291101 33.01236 9.157473e-09 7.353908e-05
## ENSG00000044446  6.258706 3.6579495 31.79771 1.710940e-08 9.159800e-05
## ENSG00000169071 -6.816668 0.8677284 30.56531 3.228138e-08 9.302654e-05
## ENSG00000128683 -5.977243 0.6514647 30.37054 3.569081e-08 9.302654e-05
## ENSG00000108231  6.633075 1.4514755 30.12366 4.053558e-08 9.302654e-05
## ENSG00000177483 -8.142323 1.9894202 30.00973 4.298845e-08 9.302654e-05
## ENSG00000085276 -7.113673 1.8199658 29.86430 4.633661e-08 9.302654e-05
## ENSG00000112619 -5.232082 0.6578571 28.80205 8.016642e-08 1.345195e-04
## ENSG00000090554 -4.896455 0.8363723 28.71723 8.375540e-08 1.345195e-04
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   331
## 0  15387
## 1    343
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
##                     logFC     logCPM       LR       PValue          FDR
## ENSG00000138650 -7.659682 1.67575067 33.87807 5.867651e-09 0.0000580736
## ENSG00000044446  6.540458 3.67746787 33.47149 7.231629e-09 0.0000580736
## ENSG00000144962 -7.676915 1.61779567 31.08434 2.470564e-08 0.0001322658
## ENSG00000108231  6.591068 1.46349823 29.94053 4.455014e-08 0.0001788800
## ENSG00000169071 -6.886334 0.72057266 29.36268 6.002242e-08 0.0001927710
## ENSG00000128683 -5.900157 0.46520459 28.96348 7.375594e-08 0.0001927710
## ENSG00000085276 -7.020925 1.55238975 28.62206 8.797410e-08 0.0001927710
## ENSG00000177483 -8.071504 1.83803493 28.45263 9.601941e-08 0.0001927710
## ENSG00000146021  5.446115 3.17203186 27.10201 1.929969e-07 0.0003256494
## ENSG00000162989 -5.075773 0.09068354 26.70524 2.369776e-07 0.0003256494
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   264
## 0  15452
## 1    345
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-11-2} \end{center}
