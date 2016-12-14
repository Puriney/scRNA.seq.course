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
## ENSG00000186439 -8.380392 3.103152 30.22830 3.840650e-08 0.0006168468
## ENSG00000205358 -9.004225 3.461390 28.60485 8.875960e-08 0.0007024891
## ENSG00000076770  7.399246 2.676850 27.25351 1.784510e-07 0.0007024891
## ENSG00000123610 -7.783018 2.086625 26.54849 2.570058e-07 0.0007024891
## ENSG00000186448  7.044079 3.610900 26.52971 2.595158e-07 0.0007024891
## ENSG00000119509  7.014250 1.522664 26.07907 3.277180e-07 0.0007024891
## ENSG00000166535 -7.818781 2.113461 25.96177 3.482456e-07 0.0007024891
## ENSG00000105974 -7.763365 2.817838 25.95256 3.499105e-07 0.0007024891
## ENSG00000204314 -7.257694 1.668619 24.79167 6.387304e-07 0.0009940588
## ENSG00000169851  6.486070 2.713638 24.79147 6.387950e-07 0.0009940588
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   479
## 0  14779
## 1    803
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
## 1           1      -0.010591129       -0.01936104
## 2           1      -0.006531984        0.05755855
## 3           1      -0.051262695        0.02268730
## 4           1       0.030955368        0.05499800
## 5           1      -0.026220790        0.03641728
## 6           1       0.005369460        0.04289753
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
## ENSG00000101951 -6.034568 1.0363284 39.30067 3.633166e-10 3.712802e-06
## ENSG00000187726 -4.865434 0.4616214 38.65121 5.067175e-10 3.712802e-06
## ENSG00000186439 -8.922191 3.1031318 38.03883 6.935064e-10 3.712802e-06
## ENSG00000185155 -5.131980 1.1809121 33.63346 6.653825e-09 2.671677e-05
## ENSG00000076770  7.824960 2.6768467 33.06597 8.908417e-09 2.861562e-05
## ENSG00000071909 -4.728837 0.2933059 32.09735 1.466374e-08 3.925239e-05
## ENSG00000105974 -7.788893 2.8178616 31.43662 2.060580e-08 4.727853e-05
## ENSG00000127329 -6.003254 1.1761452 28.91627 7.557527e-08 1.290168e-04
## ENSG00000177807 -7.211029 1.5771758 28.80617 7.999568e-08 1.290168e-04
## ENSG00000103534  8.886471 2.6852119 28.79811 8.032926e-08 1.290168e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   435
## 0  15290
## 1    336
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
## 1           1         0.3176870         0.2029780
## 2           1         0.2895892         0.1938409
## 3           1         0.3411501         0.1898753
## 4           1         0.2701635         0.1806164
## 5           1         0.3092167         0.2235698
## 6           1         0.2630351         0.1353374
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
## ENSG00000174327 -6.632293 1.0325464 43.33373 4.615570e-11 7.413066e-07
## ENSG00000205426 -6.208561 0.7649563 41.21258 1.365410e-10 7.604659e-07
## ENSG00000258466 -6.199806 0.7602025 41.13532 1.420458e-10 7.604659e-07
## ENSG00000169306 -4.829915 1.2470712 35.93229 2.042952e-09 8.202964e-06
## ENSG00000186439 -9.079850 3.1031256 34.26924 4.799066e-09 1.541556e-05
## ENSG00000109819 -5.717629 0.5843696 32.72779 1.060106e-08 2.741886e-05
## ENSG00000205358 -9.995314 3.4614250 32.33965 1.294448e-08 2.741886e-05
## ENSG00000151952 -8.583591 2.2871449 32.23548 1.365736e-08 2.741886e-05
## ENSG00000161610 -4.723280 0.1073482 31.34033 2.165344e-08 3.864177e-05
## ENSG00000205464 -6.671561 0.8499503 31.12373 2.420932e-08 3.888259e-05
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   505
## 0  15227
## 1    329
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
##                     logFC      logCPM       LR       PValue          FDR
## ENSG00000174327 -6.719538  0.92294374 42.84978 5.910877e-11 7.945056e-07
## ENSG00000205426 -6.300933  0.65267129 41.10561 1.442214e-10 7.945056e-07
## ENSG00000258466 -6.292519  0.64786229 41.04974 1.484040e-10 7.945056e-07
## ENSG00000169306 -4.617450  1.34179512 33.58038 6.837889e-09 2.745583e-05
## ENSG00000186439 -9.104446  3.21568818 33.08667 8.814097e-09 2.831264e-05
## ENSG00000109819 -5.766238  0.47758135 31.97946 1.558115e-08 3.733502e-05
## ENSG00000151952 -8.601991  2.22545510 31.89518 1.627203e-08 3.733502e-05
## ENSG00000205358 -9.824377  3.31014368 31.05555 2.507484e-08 4.527153e-05
## ENSG00000161610 -4.797050 -0.01613118 31.03295 2.536852e-08 4.527153e-05
## ENSG00000156804 -5.732683  0.62491475 29.52302 5.525637e-08 8.050215e-05
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   441
## 0  15305
## 1    315
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-11-2} \end{center}
