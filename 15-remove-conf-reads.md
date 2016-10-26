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
##                     logFC   logCPM       LR       PValue          FDR
## ENSG00000204344  7.004452 1.525096 32.47966 1.204466e-08 0.0001934492
## ENSG00000182463 -9.511724 3.664212 27.51124 1.561841e-07 0.0012542364
## ENSG00000173040 -6.844402 1.831306 26.31308 2.903151e-07 0.0015542503
## ENSG00000196167 -7.057363 1.563424 24.87986 6.101659e-07 0.0021836512
## ENSG00000185019  6.714889 1.314413 24.17922 8.777403e-07 0.0021836512
## ENSG00000091664 -6.837299 1.400167 23.72683 1.110238e-06 0.0021836512
## ENSG00000169851  7.105891 3.788751 23.46643 1.271128e-06 0.0021836512
## ENSG00000135480  6.026659 2.393240 23.39354 1.320213e-06 0.0021836512
## ENSG00000182759 -6.450023 1.139566 23.16665 1.485513e-06 0.0021836512
## ENSG00000146350  5.949527 2.241650 23.01052 1.611170e-06 0.0021836512
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   453
## 0  14823
## 1    785
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
## 1           1       0.018295663       0.063547371
## 2           1       0.030955368       0.054997999
## 3           1       0.018149399      -0.063270368
## 4           1       0.051020472      -0.008166728
## 5           1      -0.004933132      -0.044230649
## 6           1       0.012839421       0.062802508
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
## ENSG00000023892 -5.966359  1.36314159 39.97449 2.573015e-10 4.132520e-06
## ENSG00000154102 -5.061293  0.54373068 36.53241 1.501487e-09 1.189716e-05
## ENSG00000164318 -4.531673  0.34109088 35.76839 2.222245e-09 1.189716e-05
## ENSG00000185105 -3.861706  0.07823428 31.03787 2.530431e-08 8.820229e-05
## ENSG00000198844 -5.237723  0.82877346 30.87930 2.745853e-08 8.820229e-05
## ENSG00000023445 -3.771145  0.05005500 30.39778 3.519298e-08 9.420574e-05
## ENSG00000164794 -3.639360  0.08104765 29.49906 5.594374e-08 1.283589e-04
## ENSG00000158869 -3.494913 -0.02790036 28.40071 9.862893e-08 1.785974e-04
## ENSG00000257242 -3.455308 -0.03794812 28.11339 1.144115e-07 1.785974e-04
## ENSG00000115718 -4.624344  0.49554583 28.00202 1.211891e-07 1.785974e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   399
## 0  15345
## 1    317
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
## 1           1         0.2871249         0.2263144
## 2           1         0.2701635         0.1806164
## 3           1         0.2702873         0.2169784
## 4           1         0.3904875         0.1400878
## 5           1         0.2375417         0.1452930
## 6           1         0.3843119         0.1985025
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
##                     logFC   logCPM       LR       PValue          FDR
## ENSG00000204344  7.203041 1.524853 33.83950 5.985134e-09 9.612724e-05
## ENSG00000064692 -6.487226 1.140961 29.52980 5.506368e-08 4.421889e-04
## ENSG00000041982 -6.552649 1.396027 27.92519 1.260977e-07 6.570285e-04
## ENSG00000146278 -8.522043 3.828421 27.42113 1.636333e-07 6.570285e-04
## ENSG00000132436  5.628197 3.830149 25.90000 3.595683e-07 1.125507e-03
## ENSG00000051009 -8.190908 2.684092 25.35439 4.770663e-07 1.125507e-03
## ENSG00000216937  7.163452 1.833696 25.30067 4.905389e-07 1.125507e-03
## ENSG00000130513 -7.621189 1.668144 24.66667 6.815271e-07 1.368251e-03
## ENSG00000105371 -5.328311 1.315246 23.61043 1.179466e-06 1.436663e-03
## ENSG00000179862  6.387533 1.732257 23.57405 1.201977e-06 1.436663e-03
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   318
## 0  15321
## 1    422
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
## ENSG00000204344  7.424735 1.5238121 33.53966 6.982547e-09 0.0001121467
## ENSG00000064692 -6.463162 0.8971568 27.49956 1.571302e-07 0.0011319319
## ENSG00000041982 -6.671702 1.0420603 26.92568 2.114312e-07 0.0011319319
## ENSG00000132436  5.728165 3.7358726 26.17394 3.120051e-07 0.0012527783
## ENSG00000216937  7.228695 1.7957246 24.62774 6.954357e-07 0.0021989374
## ENSG00000146278 -8.342066 3.7269643 24.30684 8.214697e-07 0.0021989374
## ENSG00000135480  6.266938 2.4107527 23.63883 1.162187e-06 0.0022014596
## ENSG00000179862  6.542878 1.5642680 23.43350 1.293069e-06 0.0022014596
## ENSG00000130513 -7.706912 1.5815459 23.28053 1.400099e-06 0.0022014596
## ENSG00000042286  7.611205 1.7627526 23.19451 1.464149e-06 0.0022014596
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   244
## 0  15426
## 1    391
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="15-remove-conf-reads_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />
