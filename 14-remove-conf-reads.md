---
knit: bookdown::preview_chapter
---

# Remove confounders using controls (Reads)




```
## Loading required package: statmod
```

```
## Warning in scRNA.seq.funcs::Brennecke_getVariableGenes(assayData(reads.qc)
## $norm_counts, : Only 17 spike-ins to be used in fitting, may result in poor
## fit.
```

```
## Warning in xy.coords(x, y, xlabel, ylabel, log): 7 x values <= 0 omitted
## from logarithmic plot
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-remove-conf-reads_files/figure-latex/rm-conf-brennecke-reads-1} 

}

\caption{Results of using the Brennecke method on the Blischak dataset}(\#fig:rm-conf-brennecke-reads)
\end{figure}



\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-remove-conf-reads_files/figure-latex/rm-conf-pca-rle-reads-1} 

}

\caption{PCA plot of the blischak data after RLE normalisation}(\#fig:rm-conf-pca-rle-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-remove-conf-reads_files/figure-latex/rm-conf-pca-rle-ruv-reads-1} 

}

\caption{PCA plot of the blischak data after RLE and RUV normalisations}(\#fig:rm-conf-pca-rle-ruv-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{14-remove-conf-reads_files/figure-latex/rm-conf-rle-comp-reads-1} 

}

\caption{Comparison of the relative log expression of the blischak data before and after the RUV normalisation}(\#fig:rm-conf-rle-comp-reads)
\end{figure}
