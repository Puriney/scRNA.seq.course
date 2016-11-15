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
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000105928  6.634751 1.2199904 32.51418 1.183263e-08 0.0001900438
## ENSG00000128591 -7.694811 2.0077190 28.68684 8.508015e-08 0.0004712027
## ENSG00000173706  5.957047 2.7767207 28.62116 8.801495e-08 0.0004712027
## ENSG00000111490  7.420324 1.7940434 27.14825 1.884354e-07 0.0005967399
## ENSG00000204103  5.546966 0.5589962 26.90375 2.138439e-07 0.0005967399
## ENSG00000196597  6.466540 1.0930981 26.82335 2.229276e-07 0.0005967399
## ENSG00000157680  7.446276 2.5413064 25.64815 4.096889e-07 0.0009400019
## ENSG00000107954  6.021619 0.8231169 25.24272 5.055002e-07 0.0010148548
## ENSG00000136010  7.631070 2.6348322 23.77237 1.084280e-06 0.0014356842
## ENSG00000114378  6.527095 1.1408013 23.59658 1.187985e-06 0.0014356842
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   357
## 0  14964
## 1    740
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
## 1           1      -0.007980392        0.01889913
## 2           1      -0.014831372        0.05209915
## 3           1       0.009836556       -0.06301022
## 4           1      -0.002009755        0.04140351
## 5           1       0.023338833       -0.03946856
## 6           1      -0.051262695        0.02268730
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
##                     logFC      logCPM       LR       PValue          FDR
## ENSG00000105383 -4.767389  0.28438872 41.48555 1.187451e-10 9.535828e-07
## ENSG00000106384 -4.767389  0.28438872 41.48555 1.187451e-10 9.535828e-07
## ENSG00000149201 -4.216091  0.07717356 37.00641 1.177418e-09 6.303502e-06
## ENSG00000160097 -5.487192  0.57681125 36.25632 1.729983e-09 6.946312e-06
## ENSG00000158869 -3.698707 -0.06921461 32.72079 1.063932e-08 2.511071e-05
## ENSG00000257242 -3.655886 -0.07952024 32.36012 1.280882e-08 2.511071e-05
## ENSG00000102466 -3.612822 -0.08989963 31.98636 1.552593e-08 2.511071e-05
## ENSG00000100427 -3.612822 -0.08989963 31.98636 1.552593e-08 2.511071e-05
## ENSG00000118307 -5.282330  0.52063787 31.98167 1.556342e-08 2.511071e-05
## ENSG00000166984 -5.379685  0.68516811 31.97281 1.563459e-08 2.511071e-05
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   440
## 0  15340
## 1    281
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
## 1           1         0.3521279         0.2034467
## 2           1         0.2947330         0.1690802
## 3           1         0.3330821         0.1983118
## 4           1         0.3527278         0.2031540
## 5           1         0.3140148         0.1752826
## 6           1         0.3411501         0.1898753
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
## ENSG00000125740  7.087729 1.8031905 30.01877 4.278834e-08 0.0001995498
## ENSG00000105928  6.466170 1.2200511 29.87280 4.613408e-08 0.0001995498
## ENSG00000206073 -6.799016 2.1261105 29.62308 5.247625e-08 0.0001995498
## ENSG00000143199  6.241256 0.7685883 29.32320 6.125802e-08 0.0001995498
## ENSG00000215504  6.514500 0.9573878 29.06568 6.996547e-08 0.0001995498
## ENSG00000128591 -7.716045 2.0079360 28.84924 7.823668e-08 0.0001995498
## ENSG00000101445 -8.430144 2.1325419 28.43582 9.685642e-08 0.0001995498
## ENSG00000115232 -7.355491 2.3203502 28.28138 1.049004e-07 0.0001995498
## ENSG00000173706  5.985095 2.7767217 28.15773 1.118204e-07 0.0001995498
## ENSG00000173261 -4.963095 0.2252402 27.72171 1.400821e-07 0.0002218025
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   241
## 0  15420
## 1    400
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
## ENSG00000105928  6.488146 1.1229689 29.39846 5.892448e-08 0.0002524942
## ENSG00000125740  7.616454 1.5354269 28.74396 8.260731e-08 0.0002524942
## ENSG00000173706  6.027704 2.6720886 28.32838 1.023837e-07 0.0002524942
## ENSG00000206073 -6.876650 2.2816134 28.18759 1.101081e-07 0.0002524942
## ENSG00000172794  5.472679 0.6197632 27.76846 1.367373e-07 0.0002524942
## ENSG00000215504  6.356196 0.6514204 27.70779 1.410933e-07 0.0002524942
## ENSG00000143199  6.098565 0.4774130 27.45275 1.609794e-07 0.0002524942
## ENSG00000115232 -7.474551 2.9836377 27.39530 1.658342e-07 0.0002524942
## ENSG00000111490  7.379566 1.6867893 27.35062 1.697104e-07 0.0002524942
## ENSG00000128591 -7.768335 1.8690044 27.34713 1.700163e-07 0.0002524942
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   198
## 0  15478
## 1    385
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />
