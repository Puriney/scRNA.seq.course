---
# knit: bookdown::preview_chapter
output: html_document
---

# DE in a real dataset



```r
library(scRNA.seq.funcs)
library(DESeq2)
library(scde)
library(ROCR)
library(limma)
set.seed(1)
```

## Introduction

The main advantage of using synthetic data is that we have full
control over all aspects of the data, and this facilitates the
interpretation of the results. However, the transcriptional bursting
model is unable to capture the full complexity of a real scRNA-seq
dataset. Next, we are going to analyze the difference between the
transcriptomes of the two stages of the mouse embryo development (_mid2cell_ and _16cell_) again from the `deng`. However, this time we need to use the raw counts instead of the normalized counts. This is a major requirement for the DE analysis tools like `DESeq2`, `edgeR`, `scde` etc.

```r
deng <- readRDS("deng/deng_raw.rds")
# select cells at different stages
deng <- deng[ , c(which(colnames(deng) == "mid2cell"),
                  which(colnames(deng) == "16cell")[1:12])]
# keep those genes that are expressed in at least 6 cells
deng <- deng[rowSums(deng > 0) > 5, ]
pheatmap::pheatmap(
    log2(deng + 1),
    cutree_cols = 2,
    kmeans_k = 100,
    show_rownames = FALSE
)
```

<img src="22-de-real_files/figure-html/de-real-deng-fpkm-1.png" width="672" style="display: block; margin: auto;" />

As you can see, the cells cluster well by their developmental stage.

We can now use the same methods as before to obtain a list of
differentially expressed genes.

Because `scde` is very slow here we will only use a subset of genes. You should not do that with your real dataset, though. Here we do it just for demostration purposes:

```r
deng <- deng[sample(1:nrow(deng), 1000), ]
```

## KS-test


```r
pVals <- rep(1, nrow(deng))
for (i in 1:nrow(deng)) {
    res <- ks.test(
        deng[i, 1:12],
        deng[i, 13:24]
    )
    pVals[i] <- res$p.value
}
# Bonferroni correction
pVals <- p.adjust(pVals, method = "bonferroni")
```

## DESeq2


```r
cond <- factor(
    c(
        rep("mid2cell", 12),
        rep("16cell", 12)
    )
)
colnames(deng) <- 1:ncol(deng)
dds <- DESeq2::DESeqDataSetFromMatrix(
    deng,
    colData = DataFrame(cond),
    design = ~ cond
)
dds <- DESeq2::DESeq(dds)
resDESeq <- DESeq2::results(dds)
pValsDESeq <- resDESeq$padj
```

## SCDE


```r
cnts <- apply(
    deng,
    2,
    function(x) {
        storage.mode(x) <- 'integer'
        return(x)
    }
)
names(cond) <- 1:length(cond)
colnames(cnts) <- 1:length(cond) 
o.ifm <- scde::scde.error.models(
    counts = cnts,
    groups = cond,
    n.cores = 1,
    threshold.segmentation = TRUE,
    save.crossfit.plots = FALSE,
    save.model.plots = FALSE,
    verbose = 0,
    min.size.entries = 2
)
priors <- scde::scde.expression.prior(
    models = o.ifm,
    counts = cnts,
    length.out = 400,
    show.plot = FALSE
)
resSCDE <- scde::scde.expression.difference(
    o.ifm,
    cnts,
    priors,
    groups = cond,
    n.randomizations = 100,
    n.cores = 1,
    verbose = 0
)
# Convert Z-scores into 2-tailed p-values
pValsSCDE <- pnorm(abs(resSCDE$cZ), lower.tail = FALSE) * 2 
pValsSCDE <- p.adjust(pValsSCDE, method = "bonferroni")
```

## Comparison of the methods

```r
vennDiagram(
    vennCounts(
        cbind(
            pVals < 0.05,
            pValsDESeq < 0.05,
            pValsSCDE < 0.05
        )
    ),
    names = c("KS-test", "DESeq2", "SCDE"),
    circle.col = c("red", "blue", "green")
)
```

<img src="22-de-real_files/figure-html/de-real-comparison-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 1:__ How does this Venn diagram correspond to what you would expect based on the synthetic data? 

## Visualisation

To further characterize the list of genes, we can calculate the
average fold-changes and compare the ones that were called as
differentially expressed to the ones that were not. 


```r
group1 <- deng[, 1:12] # mid2cell stage
group2 <- deng[, 13:24] # 16cell stage
ksGenesChangedInds <- which(pVals<.05)
deSeqGenesChangedInds <- which(pValsDESeq<.05)
scdeGenesChangedInds <- which(pValsSCDE<.05)
ksGenesNotChangedInds <- which(pVals>=.05)
deSeqGenesNotChangedInds <- which(pValsDESeq>=.05)
scdeGenesNotChangedInds <- which(pValsSCDE>=.05)
meanFoldChange <- rowSums(group1)/rowSums(group2)

par(mfrow=c(2,1))
hist(log2(meanFoldChange[ksGenesChangedInds]),
     freq = FALSE,
     xlab = "# fold-change",
     col = rgb(1, 0, 0, 1/4))
hist(log2(meanFoldChange[deSeqGenesChangedInds]),
     freq = FALSE,
     xlab = "# fold-change",
     col = rgb(0, 0, 1, 1/4))
```

<img src="22-de-real_files/figure-html/de-real-deng-hist-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 2:__ Create the histogram of fold-changes for SCDE. Compare
the estimated fold-changes between the different methods? What do the
genes where they differ have in common?

### Volcano plots

A popular method for illustrating the difference between two
conditions is the volcano plot which compares the magnitude of the
change and the significance.


```r
par(mfrow=c(2,1))
plot(log2(meanFoldChange),
     -log10(pVals),
     xlab = "mean expression change",
     ylab = "-log10(P-value), KS-test") 
points(log2(meanFoldChange[ksGenesChangedInds]),
       -log10(pVals[ksGenesChangedInds]),
       col = "red")

plot(log2(meanFoldChange), 
     -log10(pValsDESeq),
     xlab = "mean expression change",
     ylab = "-log10(P-value), DESeq2")
points(log2(meanFoldChange[deSeqGenesChangedInds]),
       -log10(pValsDESeq[deSeqGenesChangedInds]),
       col = "blue")
```

<img src="22-de-real_files/figure-html/de-real-deng-volcano-1.png" width="672" style="display: block; margin: auto;" />

### MA-plot

Another popular method for illustrating the difference between two
conditions is the MA-plot, in which the data has been transformed onto the M (log ratios) and A (mean average) scale:

```r
par(mfrow=c(2,1))
plot(log2(rowMeans(group1)),
     log2(meanFoldChange),
     ylab = "mean fold change",
     xlab = "mean expression") 
points(log2(rowMeans(group1[ksGenesChangedInds,])),
       log2(meanFoldChange[ksGenesChangedInds]),
       col = "red")

plot(log2(rowMeans(group1)),
     log2(meanFoldChange),
     ylab = "mean fold change",
     xlab = "mean expression")
points(log2(rowMeans(group1[deSeqGenesChangedInds,])),
       log2(meanFoldChange[deSeqGenesChangedInds]),
       col = "blue")
```

<img src="22-de-real_files/figure-html/de-real-deng-ma-plot-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 3:__ The volcano and MA-plots for the SCDE are missing - can
you generate them? Compare to the synthetic data, what do they tell
you about the properties of the genes that have changed?

### Heatmap of DE genes

Finally, we can plot heatmaps of the genes that were called as DE by
the intersection of the three.


```r
allChangedInds <- intersect(which(pValsDESeq<.05),
                            intersect(which(pValsSCDE<.05),
                                      which(pVals<.05)))
pheatmap::pheatmap(log2(1 + deng[allChangedInds,]),
                   cutree_cols = 2,
                   show_rownames = FALSE)
```

<img src="22-de-real_files/figure-html/de-real-deng-heatmap-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 4:__ Create heatmaps for the genes that were detected by at least 2/3 methods.

