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

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg2") +
    ggtitle("PCA - RUVg normalisation: k = 2")
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-5-2.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs1") +
    ggtitle("PCA - RUVs normalisation: k = 1")
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-5-3.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2") +
    ggtitle("PCA - RUVs normalisation: k = 2")
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-5-4.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2_logcpm") +
    ggtitle("PCA - RUVs normalisation log2-cpm: k = 2")
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-5-5.png" width="672" style="display: block; margin: auto;" />

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

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

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

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1)
topTags(res1)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC   logCPM       LR       PValue         FDR
## ENSG00000158169  7.197222 2.247014 28.53859 9.184970e-08 0.001098626
## ENSG00000182759 -6.414136 1.058092 27.21032 1.824822e-07 0.001098626
## ENSG00000157014  6.901329 3.506767 26.69703 2.379867e-07 0.001098626
## ENSG00000118855  6.701824 3.399124 26.19298 3.089443e-07 0.001098626
## ENSG00000168778  6.614154 1.206185 25.99661 3.420167e-07 0.001098626
## ENSG00000196167 -7.427657 1.791652 25.36512 4.744192e-07 0.001167795
## ENSG00000165300  7.245209 2.450908 25.03266 5.636726e-07 0.001167795
## ENSG00000186952 -7.765114 2.775445 24.97203 5.816799e-07 0.001167795
## ENSG00000081923 -7.318848 1.698339 24.24409 8.486726e-07 0.001372503
## ENSG00000134698  6.038990 2.480129 23.86768 1.031895e-06 0.001372503
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   616
## 0  14721
## 1    724
```

```r
plotSmear(
    res1, lowess = TRUE,
    de.tags = rownames(topTags(res1, n = sum(abs(decideTestsDGE(res1))))$table)
)
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-8-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVg, k = 2)

```r
design_ruvg <- model.matrix(~ruvg$W[keep,] + reads.qc[, keep]$batch)
head(design_ruvg)
```

```
##   (Intercept) ruvg$W[keep, ]W_1 ruvg$W[keep, ]W_2
## 1           1       0.011884519        0.06237138
## 2           1      -0.010043808        0.07328756
## 3           1      -0.041163098        0.06647058
## 4           1      -0.004933132       -0.04423065
## 5           1      -0.010591129       -0.01936104
## 6           1       0.023211720       -0.02820286
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

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

```r
fit2 <- glmFit(dge_ruvg, design_ruvg)
res2 <- glmLRT(fit2)
topTags(res2)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000196503 -6.080249 0.8691154 38.47381 5.549354e-10 8.912817e-06
## ENSG00000087303 -7.950008 1.3060670 35.62274 2.394756e-09 1.923109e-05
## ENSG00000158481 -5.752882 0.6194376 31.66258 1.834224e-08 7.646286e-05
## ENSG00000205678 -4.836044 0.2170978 31.58976 1.904311e-08 7.646286e-05
## ENSG00000081923 -7.460578 1.6986486 28.57893 8.995571e-08 2.889557e-04
## ENSG00000143869 -5.155483 0.5813877 27.70982 1.409454e-07 3.315218e-04
## ENSG00000147174 -5.824378 0.6331817 27.57141 1.514004e-07 3.315218e-04
## ENSG00000215014 -5.296638 0.9623415 27.40351 1.651313e-07 3.315218e-04
## ENSG00000158169  6.882079 2.2470452 26.98753 2.047727e-07 3.624874e-04
## ENSG00000175426 -6.979055 1.3647334 26.79951 2.256942e-07 3.624874e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   407
## 0  15353
## 1    301
```

```r
plotSmear(
    res2, lowess = TRUE,
    de.tags = rownames(topTags(res2, n = sum(abs(decideTestsDGE(res2))))$table)
)
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-9-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVs, k = 2)

```r
design_ruvs <- model.matrix(~ruvs$W[keep,] + reads.qc[, keep]$batch)
head(design_ruvs)
```

```
##   (Intercept) ruvs$W[keep, ]W_1 ruvs$W[keep, ]W_2
## 1           1         0.3544075         0.2425527
## 2           1         0.3290784         0.1965167
## 3           1         0.3360008         0.1705224
## 4           1         0.2375417         0.1452930
## 5           1         0.3176870         0.2029780
## 6           1         0.3641549         0.1829880
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

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

```r
fit3 <- glmFit(dge_ruvs, design_ruvs)
res3 <- glmLRT(fit3)
topTags(res3)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000203710  5.123287 1.0375701 35.88773 2.090212e-09 3.357089e-05
## ENSG00000157014  7.514598 3.5068098 32.10915 1.457490e-08 6.101958e-05
## ENSG00000185686  4.480325 0.5505938 31.76658 1.738587e-08 6.101958e-05
## ENSG00000087303 -6.952885 1.3061282 31.49909 1.995334e-08 6.101958e-05
## ENSG00000182759 -6.575687 1.0585106 31.36179 2.141540e-08 6.101958e-05
## ENSG00000005884 -7.407967 2.6678679 31.24055 2.279543e-08 6.101958e-05
## ENSG00000196878  4.291226 0.4317795 30.55968 3.237527e-08 7.428274e-05
## ENSG00000203872  4.073434 0.3767709 29.96864 4.390915e-08 8.815311e-05
## ENSG00000169900 -7.355347 1.8172005 29.36503 5.994980e-08 1.069837e-04
## ENSG00000158169  7.248781 2.2470553 28.89720 7.632297e-08 1.225823e-04
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   390
## 0  15190
## 1    481
```

```r
plotSmear(
    res3, lowess = TRUE,
    de.tags = rownames(topTags(res3, n = sum(abs(decideTestsDGE(res3))))$table)
)
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-10-2.png" width="672" style="display: block; margin: auto;" />


```r
reads.qc <- scran::computeSumFactors(reads.qc, sizes = 15)
dge_ruvs$samples$norm.factors <- sizeFactors(reads.qc)[keep]
dge_ruvs_sf <- estimateDisp(dge_ruvs, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs_sf)
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

```r
fit4 <- glmFit(dge_ruvs_sf, design_ruvs)
res4 <- glmLRT(fit4)
topTags(res4)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000203710  5.455013 1.3048741 37.21201 1.059595e-09 1.701815e-05
## ENSG00000185686  4.784438 0.7075399 33.09979 8.754839e-09 6.096704e-05
## ENSG00000157014  7.677421 3.3696408 32.07787 1.481150e-08 6.096704e-05
## ENSG00000196878  4.590897 0.5547440 31.89324 1.628830e-08 6.096704e-05
## ENSG00000203872  4.329318 0.4843484 31.26623 2.249583e-08 6.096704e-05
## ENSG00000005884 -7.365194 2.5480047 30.95379 2.642457e-08 6.096704e-05
## ENSG00000087303 -7.045420 1.1278744 30.94301 2.657177e-08 6.096704e-05
## ENSG00000182759 -6.510981 0.8226632 30.17516 3.947322e-08 7.924743e-05
## ENSG00000158169  7.330271 2.0045313 28.96524 7.368881e-08 1.315018e-04
## ENSG00000169900 -7.331762 1.6392305 28.35914 1.007697e-07 1.618462e-04
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   348
## 0  15240
## 1    473
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />
