---
knit: bookdown::preview_chapter
---

# Remove confounders using controls (Reads)




```
## Warning in scRNA.seq.funcs::Brennecke_getVariableGenes(assayData(reads.qc)
## $norm_counts, : Only 15 spike-ins to be used in fitting, may result in poor
## fit.
```

<div class="figure" style="text-align: center">
<img src="14-remove-conf-reads_files/figure-html/rm-conf-brennecke-reads-1.png" alt="(\#fig:rm-conf-brennecke-reads)Results of using the Brennecke method on the Blischak dataset" width="90%" />
<p class="caption">(\#fig:rm-conf-brennecke-reads)Results of using the Brennecke method on the Blischak dataset</p>
</div>



<div class="figure" style="text-align: center">
<img src="14-remove-conf-reads_files/figure-html/rm-conf-pca-rle-reads-1.png" alt="(\#fig:rm-conf-pca-rle-reads)PCA plot of the blischak data after RLE normalisation" width="90%" />
<p class="caption">(\#fig:rm-conf-pca-rle-reads)PCA plot of the blischak data after RLE normalisation</p>
</div>

<div class="figure" style="text-align: center">
<img src="14-remove-conf-reads_files/figure-html/rm-conf-pca-rle-ruv-reads-1.png" alt="(\#fig:rm-conf-pca-rle-ruv-reads)PCA plot of the blischak data after RLE and RUV normalisations" width="90%" />
<p class="caption">(\#fig:rm-conf-pca-rle-ruv-reads)PCA plot of the blischak data after RLE and RUV normalisations</p>
</div>

<div class="figure" style="text-align: center">
<img src="14-remove-conf-reads_files/figure-html/rm-conf-rle-comp-reads-1.png" alt="(\#fig:rm-conf-rle-comp-reads)Comparison of the relative log expression of the blischak data before and after the RUV normalisation" width="90%" />
<p class="caption">(\#fig:rm-conf-rle-comp-reads)Comparison of the relative log expression of the blischak data before and after the RUV normalisation</p>
</div>
