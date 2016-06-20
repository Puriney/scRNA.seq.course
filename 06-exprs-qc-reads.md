---
knit: bookdown::preview_chapter
---

# Expression QC (Reads)

This chapter contains the summary plots and tables for the QC exercise based on the reads for the Bischak data discussed in the previous chapter.






Table: (\#tab:unnamed-chunk-3)A table of the first 6 rows and 3 columns of the molecules table.

                   NA19098.r1.A01   NA19098.r1.A02   NA19098.r1.A03
----------------  ---------------  ---------------  ---------------
ENSG00000237683                 0                0                0
ENSG00000187634                 0                0                0
ENSG00000188976                57              140                1
ENSG00000187961                 0                0                0
ENSG00000187583                 0                0                0
ENSG00000187642                 0                0                0



Table: (\#tab:unnamed-chunk-3)A table of the first 6 rows of the anno table.

individual   replicate   well   batch        sample_id      
-----------  ----------  -----  -----------  ---------------
NA19098      r1          A01    NA19098.r1   NA19098.r1.A01 
NA19098      r1          A02    NA19098.r1   NA19098.r1.A02 
NA19098      r1          A03    NA19098.r1   NA19098.r1.A03 
NA19098      r1          A04    NA19098.r1   NA19098.r1.A04 
NA19098      r1          A05    NA19098.r1   NA19098.r1.A05 
NA19098      r1          A06    NA19098.r1   NA19098.r1.A06 









<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/total-counts-hist-reads-1.png" alt="(\#fig:total-counts-hist-reads)Histogram of library sizes for all cells" width="90%" />
<p class="caption">(\#fig:total-counts-hist-reads)Histogram of library sizes for all cells</p>
</div>




Table: (\#tab:unnamed-chunk-9)The number of cells removed by total counts filter (FALSE)

filter_by_total_counts    Freq
-----------------------  -----
FALSE                      180
TRUE                       684

<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/total-features-hist-reads-1.png" alt="(\#fig:total-features-hist-reads)Histogram of the number of detected genes in all cells" width="90%" />
<p class="caption">(\#fig:total-features-hist-reads)Histogram of the number of detected genes in all cells</p>
</div>




Table: (\#tab:unnamed-chunk-11)The number of cells removed by total features filter (FALSE)

filter_by_expr_features    Freq
------------------------  -----
FALSE                       120
TRUE                        744

<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/total-features-vs-counts-reads-1.png" alt="(\#fig:total-features-vs-counts-reads)Library size vs number of detected genes" width="90%" />
<p class="caption">(\#fig:total-features-vs-counts-reads)Library size vs number of detected genes</p>
</div>

<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/mt-vs-counts-reads-1.png" alt="(\#fig:mt-vs-counts-reads)Percentage of counts in MT genes" width="90%" />
<p class="caption">(\#fig:mt-vs-counts-reads)Percentage of counts in MT genes</p>
</div>

<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/ercc-vs-counts-reads-1.png" alt="(\#fig:ercc-vs-counts-reads)Percentage of counts in ERCCs" width="90%" />
<p class="caption">(\#fig:ercc-vs-counts-reads)Percentage of counts in ERCCs</p>
</div>




Table: (\#tab:unnamed-chunk-13)The number of cells removed by ERCC filter (FALSE)

filter_by_ERCC    Freq
---------------  -----
FALSE              103
TRUE               761




Table: (\#tab:unnamed-chunk-15)The number of cells removed by MT filter (FALSE)

filter_by_MT    Freq
-------------  -----
FALSE             18
TRUE             846




Table: (\#tab:unnamed-chunk-17)The number of cells removed by default filter (FALSE)

Var1     Freq
------  -----
FALSE      37
TRUE      827


```
## The following cells/samples are detected as outliers:
## NA19098.r1.B10
## NA19098.r1.D07
## NA19098.r1.E04
## NA19098.r1.F06
## NA19098.r1.H08
## NA19098.r1.H09
## NA19098.r2.A01
## NA19098.r2.A06
## NA19098.r2.A09
## NA19098.r2.A12
## NA19098.r2.B01
## NA19098.r2.B11
## NA19098.r2.B12
## NA19098.r2.C04
## NA19098.r2.C09
## NA19098.r2.D02
## NA19098.r2.D03
## NA19098.r2.D09
## NA19098.r2.E04
## NA19098.r2.E07
## NA19098.r2.F01
## NA19098.r2.F11
## NA19098.r2.G01
## NA19098.r2.G05
## NA19098.r2.G10
## NA19098.r2.H01
## NA19098.r2.H07
## NA19098.r2.H08
## NA19098.r2.H12
## NA19098.r3.A05
## NA19098.r3.A07
## NA19098.r3.B02
## NA19098.r3.C07
## NA19098.r3.E05
## NA19098.r3.E08
## NA19098.r3.E09
## NA19098.r3.F11
## NA19098.r3.F12
## NA19098.r3.G02
## NA19098.r3.G03
## NA19098.r3.G04
## NA19098.r3.G11
## NA19098.r3.G12
## NA19098.r3.H08
## NA19101.r1.A01
## NA19101.r1.A12
## NA19101.r1.B01
## NA19101.r1.B06
## NA19101.r1.E09
## NA19101.r1.E11
## NA19101.r1.F05
## NA19101.r1.F10
## NA19101.r1.G01
## NA19101.r1.G06
## NA19101.r1.H04
## NA19101.r1.H09
## NA19101.r2.A03
## NA19101.r2.C10
## NA19101.r2.E05
## NA19101.r2.F02
## NA19101.r2.H04
## NA19101.r2.H10
## NA19101.r3.A02
## NA19101.r3.A03
## NA19101.r3.A05
## NA19101.r3.A09
## NA19101.r3.B05
## NA19101.r3.C01
## NA19101.r3.C09
## NA19101.r3.C12
## NA19101.r3.D01
## NA19101.r3.D04
## NA19101.r3.D07
## NA19101.r3.D09
## NA19101.r3.E08
## NA19101.r3.F09
## NA19101.r3.G09
## NA19101.r3.H01
## NA19101.r3.H03
## NA19101.r3.H07
## NA19101.r3.H09
## NA19239.r1.F05
## NA19239.r1.G05
## NA19239.r2.B01
## NA19239.r2.B03
## NA19239.r2.B10
## NA19239.r2.B11
## NA19239.r2.C03
## NA19239.r2.C06
## NA19239.r2.C08
## NA19239.r2.D07
## NA19239.r2.D09
## NA19239.r2.E09
## NA19239.r2.F04
## NA19239.r2.F06
## NA19239.r2.F07
## NA19239.r2.F12
## NA19239.r2.G03
## NA19239.r2.G08
## NA19239.r2.H02
## NA19239.r2.H03
## NA19239.r2.H07
## NA19239.r3.A01
## NA19239.r3.B09
## NA19239.r3.C04
## NA19239.r3.C07
## NA19239.r3.E01
## NA19239.r3.E03
## NA19239.r3.E12
## NA19239.r3.H02
## NA19239.r3.H10
## Variables with highest loadings for PC1 and PC2:
## 
##                                            PC1         PC2
## ---------------------------------  -----------  ----------
## pct_counts_feature_controls          0.5057646   0.2473134
## pct_counts_top_100_features          0.4888852   0.2277068
## n_detected_feature_controls          0.0231277   0.6235516
## log10_counts_feature_controls       -0.1226860   0.6576822
## total_features                      -0.4655518   0.2219694
## log10_counts_endogenous_features    -0.5223679   0.1278782
```

<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/auto-cell-filt-reads-1.png" alt="(\#fig:auto-cell-filt-reads)PCA plot used for automatic detection of cell outliers" width="90%" />
<p class="caption">(\#fig:auto-cell-filt-reads)PCA plot used for automatic detection of cell outliers</p>
</div>


Table: (\#tab:unnamed-chunk-18)The number of cells removed by automatic filter (FALSE)

Var1     Freq
------  -----
FALSE     753
TRUE      111




Table: (\#tab:unnamed-chunk-20)The number of cells removed by manual filter (FALSE)

Var1     Freq
------  -----
FALSE     259
TRUE      605

<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/cell-filt-comp-reads-1.png" alt="(\#fig:cell-filt-comp-reads)Comparison of the default, automatic and manual cell filters" width="90%" />
<p class="caption">(\#fig:cell-filt-comp-reads)Comparison of the default, automatic and manual cell filters</p>
</div>

<div class="figure" style="text-align: center">
<img src="06-exprs-qc-reads_files/figure-html/top50-gene-expr-reads-1.png" alt="(\#fig:top50-gene-expr-reads)Number of total counts consumed by the top 50 expressed genes" width="90%" />
<p class="caption">(\#fig:top50-gene-expr-reads)Number of total counts consumed by the top 50 expressed genes</p>
</div>




Table: (\#tab:unnamed-chunk-22)The number of genes removed by gene filter (FALSE)

filter_genes     Freq
-------------  ------
FALSE            2665
TRUE            16061


```
## Features  Samples 
##    16061      605
```



If you want to further check yourself you can download our [`reads`](http://hemberg-lab.github.io/scRNA.seq.course/blischak/reads.rds) object. If you followed the steps above it should be exactly the same as yours.

By comparing Figure 7.7 and Figure 6.6, it is clear that the reads based filtering removed 49 more cells than the UMI based analysis. If you go back and compare the results you should be able to conclude that the ERCC and MT filters are more strict for the reads-based analysis.
