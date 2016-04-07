---
knit: bookdown::preview_chapter
---

# Normalization for library size

## Introduction



In the previous chapter we identified important confounding factors and explanatory variables. scater allows one to account for these variables in subsequent statistical models or to condition them out using `normaliseExprs()`, if so desired. This can be done by providing a design matrix to `normaliseExprs()`. We are not covering this topic here, but you can try to do it yourself as an exercise.

Instead we will explore how simple size-factor normalisations correcting for library size can remove the effects of some of the confounders and explanatory variables.

## Library size

Library sizes vary because scRNA-seq data is often sequenced on highly multiplexed platforms the total reads which are derived from each cell may differ substantially. Some quantification methods
(eg. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/), [RSEM](http://deweylab.github.io/RSEM/)) incorporated library size when determining gene expression estimates thus do not require this normalization.

However, if another quantification method was used then library size must be corrected for by multiplying or dividing each column of the expression matrix by a "normalization factor" which is an estimate of the library size relative to the other cells. Many methods to correct for library size have been developped for bulk RNA-seq and can be equally applied to scRNA-seq (eg. UQ, SF, CPM, RPKM, FPKM, TPM). 


## Normalisations

The simplest way to normalize this data is to convert it to counts per
million (__CPM__) by dividing each column by its total then multiplying by
1,000,000. Note that spike-ins should be excluded from the
calculation of total expression in order to correct for total cell RNA
content, therefore we will only use endogenous genes. 


```r
calc_cpm <-
function (expr_mat, spikes = NULL) 
{
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
}
```

One potential drawback of __CPM__ is if your sample contains genes that are both very highly expressed and differentially expressed across the cells. In this case, the total molecules in the cell may depend of whether such genes are on/off in the cell and normalizing by total molecules may hide the differential expression of those genes and/or falsely create differential expression for the remaining genes. 

__Note__: RPKM, FPKM and TPM are variants on CPM which further adjust counts by the length of the respective gene/transcript.

To deal with this potentiality several other measures were devised:

The __size factor (SF)__ was proposed and popularized by DESeq ([Anders and Huber (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)). First the geometric mean of each gene across all cells is calculated. The size factor for each cell is the median across genes of the ratio of the expression to the gene's geometric mean. A draw back to this method is that since it uses the geometric mean only genes with non-zero expression values across all cells can be used in its calculation, making it unadvisable for large low-depth scRNASeq experiments. edgeR & scater call this method __RLE__ for "relative log expression".


```r
calc_sf <-
function (expr_mat, spikes = NULL) 
{
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
        median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
            0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
}
```

The __upperquartile (UQ)__ was proposed by [Bullard et al (2010)](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94). Here each column is divided by the 75% quantile of the counts for each library. Often the calculated quantile is scaled by the median across cells to keep the absolute level of expression relatively consistent. A draw back to this method is that for low-depth scRNASeq experiments the large number of undetected genes may result in the 75% quantile being zero (or close to it). This limitation can be overcome by generalizing the idea and using a higher quantile (eg. the 99% quantile is the default in scater) or by excluding zeros prior to calculating the 75% quantile.


```r
calc_uq <-
function (expr_mat, spikes = NULL) 
{
    UQ <- function(x) {
        quantile(x[x > 0], 0.75)
    }
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
}
```

Another method is called __TMM__ is the weighted trimmed mean of M-values (to the reference) proposed by [Robinson and Oshlack (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25). The M-values in question are the gene-wise log2 fold changes between individual cells. One cell is used as the reference then the M-values for each other cell is calculated compared  to this reference. These values are then trimmed by removing the top and bottom ~30%, and the average of the remaining values is calculated by weighting them to account for the effect of the log scale on variance. Each non-reference cell is multiplied by the calculated factor. Two potential issues with this method are insufficient non-zero genes left after trimming, and the assumption that most genes are not differentially expressed.

We will use visual inspection of PCA plots and calculation of cell-wise relative log expression (`calc_cell_RLE()`) to compare the efficiency of different normalization methods. Cells with many[few] reads have higher[lower] than median expression for most genes resulting in a positive[negative] RLE across the cell, whereas normalized cells have an RLE close to zero.


```r
calc_cell_RLE <-
function (expr_mat, spikes = NULL) 
{
    RLE_gene <- function(x) {
        if (median(unlist(x)) > 0) {
            log((x + 1)/(median(unlist(x)) + 1))/log(2)
        }
        else {
            rep(NA, times = length(x))
        }
    }
    if (!is.null(spikes)) {
        RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
        RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(cell_RLE)
}
```

_Tallulah: I need to revamp this to (1) include RLE, (2) change what is/isn't an exercise, and (3) possibly add comparison of scater/edgeR with actual functions for DESeq_

_Tallulah: Vague plan: show RLE & compare to DESeq 'truth' then have UQ, TMM, and CPM as exercises. Show everything for FPKMs_

 __scater__ acts as a wrapper for the `calcNormFactors` function from [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) which implements several library size normalization methods making it easy to apply any of these methods to our data. 

__Note:__ edgeR makes extra adjustments to some of the normalization methods which will result in somewhat different results than if the original methods are followed exactly, e.g. edgeR's and scater's "RLE" method which is based on the "size factor" used by [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) may give different results to the `estimateSizeFactorsForMatrix` method in the DESeq/DESeq2 packages. In addition, some versions of edgeR will not calculate the normalization factors correctly unless `lib.size` is set at 1 for all cells.

We will continue to work with the Blischak data that was used in the previous chapter.


```r
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("blischak/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
```

### Raw

```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "counts")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-pca-raw-1.png" alt="(\#fig:norm-pca-raw)PCA plot of the blischak data" width="90%" />
<p class="caption">(\#fig:norm-pca-raw)PCA plot of the blischak data</p>
</div>


```r
boxplot(calc_cell_RLE(counts(umi.qc)),
        col = "grey50",
        ylab = "RLE",
        main = "")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-ours-rle-raw-1.png" alt="(\#fig:norm-ours-rle-raw)Cell-wise RLE of the blischak data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-raw)Cell-wise RLE of the blischak data</p>
</div>

### CPM
scater performs this normalisation by default, you can control it by changing `exprs_values` parameter to `"exprs"`.

```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "cpm")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-pca-cpm-1.png" alt="(\#fig:norm-pca-cpm)PCA plot of the blischak data after CPM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-cpm)PCA plot of the blischak data after CPM normalisation</p>
</div>

```r
boxplot(calc_cell_RLE(cpm(umi.qc)),
        col = "grey50",
        ylab = "RLE",
        main = "")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-ours-rle-cpm-1.png" alt="(\#fig:norm-ours-rle-cpm)Cell-wise RLE of the blischak data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-cpm)Cell-wise RLE of the blischak data</p>
</div>


### TMM

```r
umi.qc <- 
    scater::normaliseExprs(umi.qc,
                           method = "TMM",
                           feature_set = endog_genes,
                           lib.size = rep(1, ncol(umi.qc)))
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-pca-tmm-1.png" alt="(\#fig:norm-pca-tmm)PCA plot of the blischak data after TMM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-tmm)PCA plot of the blischak data after TMM normalisation</p>
</div>

```r
boxplot(calc_cell_RLE(norm_counts(umi.qc)),
        col = "grey50",
        ylab = "RLE",
        main = "")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-ours-rle-tmm-1.png" alt="(\#fig:norm-ours-rle-tmm)Cell-wise RLE of the blischak data" width="90%" />
<p class="caption">(\#fig:norm-ours-rle-tmm)Cell-wise RLE of the blischak data</p>
</div>

__Exercise:__ 

Use `method = "RLE"` and `method = "upperquartile"` to perform size-factor and UQ normalizations and compare to the results above.

### RLE 

```r
umi.qc <- 
    scater::normaliseExprs(umi.qc,
                           method = "RLE",
                           feature_set = endog_genes,
                           lib.size = rep(1, ncol(umi.qc)))
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-pca-rle-1.png" alt="(\#fig:norm-pca-rle)PCA plot of the blischak data after RLE normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-rle)PCA plot of the blischak data after RLE normalisation</p>
</div>

### Upperquantile

```r
umi.qc <- 
    scater::normaliseExprs(umi.qc,
                           method = "upperquartile", 
                           feature_set = endog_genes,
                           p = 0.99,
                           lib.size = rep(1, ncol(umi.qc)))
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-pca-uq-1.png" alt="(\#fig:norm-pca-uq)PCA plot of the blischak data after UQ normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-uq)PCA plot of the blischak data after UQ normalisation</p>
</div>


## Other methods

Some methods combine library size and fragment/gene length normalization such as:

* __RPKM__ - Reads Per Kilobase Million (for single-end sequencing)
* __FPKM__ - Fragments Per Kilobase Million (same as __RPKM__ but for paired-end sequencing, makes sure that paired ends mapped to the same fragment are not counted twice)
* __TPM__ - Transcripts Per Kilobase Million (same as __RPKM__, but the order of normalizations is reversed - length first and sequencing depth second)

These methods are not applicable to our dataset since the end
of the transcript which contains the UMI was preferentially
sequenced. Furthermore in general these should only be calculated
using appropriate quantification software from aligned BAM files not
from read counts since often only a portion of the entire
gene/transcript is sequenced, not the entire length.

However, here we show how these normalisations can be calculated using scater. First, we need to find the effective transcript length in Kilobases. However, our dataset containes only gene IDs, therefore we will be using the gene lengths instead of transcripts. scater uses the [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) package, which allows one to annotate genes by other attributes:

```r
umi.qc <-
    scater::getBMFeatureAnnos(umi.qc,
                              filters = "ensembl_gene_id", 
                              attributes = c("ensembl_gene_id",
                                             "hgnc_symbol",
                                             "chromosome_name",
                                             "start_position",
                                             "end_position"), 
                              feature_symbol = "hgnc_symbol",
                              feature_id = "ensembl_gene_id",
                              biomart = "ENSEMBL_MART_ENSEMBL", 
                              dataset = "hsapiens_gene_ensembl",
                              host = "www.ensembl.org")

# If you have mouse data, change the arguments based on this example:
# scater::getBMFeatureAnnos(object,
#                           filters = "ensembl_transcript_id", 
#                           attributes = c("ensembl_transcript_id", 
#                                        "ensembl_gene_id", "mgi_symbol", 
#                                        "chromosome_name",
#                                        "transcript_biotype",
#                                        "transcript_start",
#                                        "transcript_end", 
#                                        "transcript_count"), 
#                           feature_symbol = "mgi_symbol",
#                           feature_id = "ensembl_gene_id",
#                           biomart = "ENSEMBL_MART_ENSEMBL", 
#                           dataset = "mmusculus_gene_ensembl",
#                           host = "www.ensembl.org") 
```

Some of the genes were not annotated, therefore we filter them out:

```r
umi.qc.ann <-
    umi.qc[!is.na(fData(umi.qc)$ensembl_gene_id), ]
```

Now we compute the total gene length in Kilobases by using the `end_position` and `start_position` fields:

```r
eff_length <- abs(fData(umi.qc.ann)$end_position -
                      fData(umi.qc.ann)$start_position)/1000
```

Now we are ready to perform the normalisations:

```r
tpm(umi.qc.ann) <-
    calculateTPM(
        umi.qc.ann,
        eff_length
    )
fpkm(umi.qc.ann) <-
    calculateFPKM(
        umi.qc.ann,
        eff_length
    )
```

Plot the results as a PCA plot:

```r
scater::plotPCA(umi.qc.ann,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "fpkm")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-pca-fpkm-1.png" alt="(\#fig:norm-pca-fpkm)PCA plot of the blischak data after FPKM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-fpkm)PCA plot of the blischak data after FPKM normalisation</p>
</div>


```r
scater::plotPCA(umi.qc.ann,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "tpm")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-pca-tpm-1.png" alt="(\#fig:norm-pca-tpm)PCA plot of the blischak data after TPM normalisation" width="90%" />
<p class="caption">(\#fig:norm-pca-tpm)PCA plot of the blischak data after TPM normalisation</p>
</div>

__Note:__ The PCA primarily looks for difference between cells. Since gene length is the same across cells for each gene FPKM is almost identical to the CPM plot (it is just rotated). However, if we plot the mean expression vs gene length it is apparent that normalizing for gene length should not be performed for this dataset.

```r
plot(eff_length, rowMeans(counts(umi.qc.ann)))
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/length-vs-mean-1.png" alt="(\#fig:length-vs-mean)Gene length vs Mean Expression for the raw data" width="90%" />
<p class="caption">(\#fig:length-vs-mean)Gene length vs Mean Expression for the raw data</p>
</div>

## Visualize genes

Now after the normalisation we are ready to visualise the gene expression:

### Raw

```r
scater::plotExpression(umi.qc.ann,
                       rownames(umi.qc.ann)[1:6],
                       x = "individual",
                       exprs_values = "counts",
                       colour = "batch")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-genes-raw-1.png" alt="(\#fig:norm-genes-raw)Expression of the first 6 genes of the blischak data" width="90%" />
<p class="caption">(\#fig:norm-genes-raw)Expression of the first 6 genes of the blischak data</p>
</div>

### CPM

```r
scater::plotExpression(umi.qc.ann,
                       rownames(umi.qc.ann)[1:6],
                       x = "individual",
                       exprs_values = "cpm",
                       colour = "batch")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-genes-cpm-1.png" alt="(\#fig:norm-genes-cpm)Expression of the first 6 genes of the blischak data after the CPM normalisation" width="90%" />
<p class="caption">(\#fig:norm-genes-cpm)Expression of the first 6 genes of the blischak data after the CPM normalisation</p>
</div>

### log2(CPM)

```r
scater::plotExpression(umi.qc.ann,
                       rownames(umi.qc.ann)[1:6],
                       x = "individual",
                       exprs_values = "exprs",
                       colour = "batch")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-genes-log2-cpm-1.png" alt="(\#fig:norm-genes-log2-cpm)Expression of the first 6 genes of the blischak data after the log2(CPM) normalisation" width="90%" />
<p class="caption">(\#fig:norm-genes-log2-cpm)Expression of the first 6 genes of the blischak data after the log2(CPM) normalisation</p>
</div>

### Upperquantile

```r
scater::plotExpression(umi.qc.ann,
                       rownames(umi.qc.ann)[1:6],
                       x = "individual",
                       exprs_values = "norm_counts",
                       colour = "batch")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-genes-UQ-1.png" alt="(\#fig:norm-genes-UQ)Expression of the first 6 genes of the blischak data after the UQ normalisation" width="90%" />
<p class="caption">(\#fig:norm-genes-UQ)Expression of the first 6 genes of the blischak data after the UQ normalisation</p>
</div>

### FPKM

```r
scater::plotExpression(umi.qc.ann,
                       rownames(umi.qc.ann)[1:6],
                       x = "individual",
                       exprs_values = "fpkm",
                       colour = "batch")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-genes-fpkm-1.png" alt="(\#fig:norm-genes-fpkm)Expression of the first 6 genes of the blischak data after the FPKM normalisation" width="90%" />
<p class="caption">(\#fig:norm-genes-fpkm)Expression of the first 6 genes of the blischak data after the FPKM normalisation</p>
</div>

### TPM

```r
scater::plotExpression(umi.qc.ann,
                       rownames(umi.qc.ann)[1:6],
                       x = "individual",
                       exprs_values = "tpm",
                       colour = "batch")
```

<div class="figure" style="text-align: center">
<img src="11-exprs-norm_files/figure-html/norm-genes-tpm-1.png" alt="(\#fig:norm-genes-tpm)Expression of the first 6 genes of the blischak data after the TPM normalisation" width="90%" />
<p class="caption">(\#fig:norm-genes-tpm)Expression of the first 6 genes of the blischak data after the TPM normalisation</p>
</div>

## Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter).
