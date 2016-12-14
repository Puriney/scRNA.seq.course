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
## ENSG00000113369  8.100379 2.949502 41.60426 1.117494e-10 1.794807e-06
## ENSG00000172731  7.082647 2.954827 31.71628 1.784198e-08 1.432800e-04
## ENSG00000182463 -9.541607 3.683570 27.81035 1.338081e-07 6.859104e-04
## ENSG00000166987  6.861502 3.163053 27.04573 1.986987e-07 6.859104e-04
## ENSG00000095015  5.145575 4.138913 26.59132 2.513706e-07 6.859104e-04
## ENSG00000091428 -7.299255 1.719740 26.29932 2.923914e-07 6.859104e-04
## ENSG00000174485  7.177151 2.821271 26.08082 3.274197e-07 6.859104e-04
## ENSG00000139269  6.454930 2.759370 25.78706 3.812359e-07 6.859104e-04
## ENSG00000010319  6.485168 1.470522 25.77131 3.843592e-07 6.859104e-04
## ENSG00000186815  7.231468 2.619838 25.31317 4.873709e-07 7.175933e-04
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   434
## 0  14724
## 1    903
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
## 1           1       0.011258160      -0.058359117
## 2           1       0.053689592      -0.035473558
## 3           1       0.018149399      -0.063270368
## 4           1       0.023547154       0.032317902
## 5           1      -0.006531984       0.057558548
## 6           1       0.049484081      -0.004065018
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
## ENSG00000101605 -5.735918  0.89707647 46.51267 9.102797e-12 1.462000e-07
## ENSG00000113369  8.166276  2.94955696 40.49438 1.971819e-10 1.583470e-06
## ENSG00000182600 -4.380632  0.17211666 39.11599 3.993588e-10 2.138034e-06
## ENSG00000205838 -5.401921  0.79446903 34.24512 4.858927e-09 1.950981e-05
## ENSG00000172731  7.232770  2.95477322 31.93624 1.593170e-08 5.117581e-05
## ENSG00000115718 -3.593633 -0.02353303 29.29469 6.216580e-08 1.664075e-04
## ENSG00000151577 -5.013203  0.78802233 26.99392 2.040959e-07 3.952267e-04
## ENSG00000139269  6.700945  2.75939996 26.90440 2.137721e-07 3.952267e-04
## ENSG00000166987  6.858631  3.16296112 26.70710 2.367503e-07 3.952267e-04
## ENSG00000095015  5.292069  4.13887548 26.46736 2.680286e-07 3.952267e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   314
## 0  15357
## 1    390
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
## 1           1         0.3764673         0.1991044
## 2           1         0.3837730         0.2115007
## 3           1         0.2702873         0.2169784
## 4           1         0.3470296         0.2160446
## 5           1         0.2895892         0.1938409
## 6           1         0.3781664         0.1922841
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
## ENSG00000131126 -6.568281 0.9349129 54.03695 1.967544e-13 2.444343e-09
## ENSG00000165309 -6.454869 0.8715194 53.17974 3.043824e-13 2.444343e-09
## ENSG00000184185 -5.492349 0.3584413 51.81512 6.098043e-13 3.264689e-09
## ENSG00000196735 -6.091716 0.6785122 50.38929 1.260807e-12 5.062457e-09
## ENSG00000127533 -5.817543 0.5453366 48.26163 3.729803e-12 1.198087e-08
## ENSG00000225556 -5.632087 0.4614987 46.81063 7.818799e-12 2.092962e-08
## ENSG00000203780 -5.560035 0.4303208 46.24460 1.043738e-11 2.394783e-08
## ENSG00000167131 -5.658000 0.6265610 43.25578 4.803156e-11 9.642935e-08
## ENSG00000113369  7.998970 2.9495742 40.29936 2.178799e-10 3.888188e-07
## ENSG00000137252 -5.809642 0.5795658 38.95602 4.334616e-10 6.961827e-07
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   348
## 0  15258
## 1    455
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
## ENSG00000131126 -6.419037 0.7867827 51.47495 7.251665e-13 5.976155e-09
## ENSG00000184185 -5.420805 0.2737568 51.20820 8.307089e-13 5.976155e-09
## ENSG00000165309 -6.308255 0.7242035 50.62824 1.116273e-12 5.976155e-09
## ENSG00000196735 -5.951366 0.5336947 47.87303 4.547331e-12 1.825867e-08
## ENSG00000127533 -5.680596 0.4022356 45.77365 1.327370e-11 4.263776e-08
## ENSG00000225556 -5.498108 0.3194622 44.34286 2.756152e-11 7.377759e-08
## ENSG00000203780 -5.427236 0.2886753 43.78495 3.665126e-11 8.409369e-08
## ENSG00000167131 -5.447891 0.3677534 40.99746 1.524274e-10 3.060171e-07
## ENSG00000113369  8.032599 2.8655308 40.28856 2.190880e-10 3.909747e-07
## ENSG00000168702 -7.087317 1.2505385 37.19061 1.071284e-09 1.720590e-06
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   302
## 0  15318
## 1    441
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf-reads_files/figure-latex/unnamed-chunk-11-2} \end{center}
