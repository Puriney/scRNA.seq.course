---
knit: bookdown::preview_chapter
---

# Expression QC (Reads)





\begin{table}

\caption{(\#tab:unnamed-chunk-3)A table of the first 6 rows and 3 columns of the molecules table.}
\centering
\begin{tabular}[t]{lrrr}
\toprule
  & NA19098.r1.A01 & NA19098.r1.A02 & NA19098.r1.A03\\
\midrule
ENSG00000237683 & 0 & 0 & 0\\
ENSG00000187634 & 0 & 0 & 0\\
ENSG00000188976 & 57 & 140 & 1\\
ENSG00000187961 & 0 & 0 & 0\\
ENSG00000187583 & 0 & 0 & 0\\
ENSG00000187642 & 0 & 0 & 0\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}

\caption{(\#tab:unnamed-chunk-3)A table of the first 6 rows of the anno table.}
\centering
\begin{tabular}[t]{lllll}
\toprule
individual & replicate & well & batch & sample\_id\\
\midrule
NA19098 & r1 & A01 & NA19098.r1 & NA19098.r1.A01\\
NA19098 & r1 & A02 & NA19098.r1 & NA19098.r1.A02\\
NA19098 & r1 & A03 & NA19098.r1 & NA19098.r1.A03\\
NA19098 & r1 & A04 & NA19098.r1 & NA19098.r1.A04\\
NA19098 & r1 & A05 & NA19098.r1 & NA19098.r1.A05\\
NA19098 & r1 & A06 & NA19098.r1 & NA19098.r1.A06\\
\bottomrule
\end{tabular}
\end{table}









\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/dropout-overview-reads-1} 

}

\caption{Dropout rate vs mean expression}(\#fig:dropout-overview-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/top50-gene-expr-reads-1} 

}

\caption{Number of total counts consumed by the top 50 expressed genes}(\#fig:top50-gene-expr-reads)
\end{figure}



\begin{table}

\caption{(\#tab:unnamed-chunk-9)The number of genes removed by gene filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_genes & Freq\\
\midrule
FALSE & 2445\\
TRUE & 16281\\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/total-counts-hist-reads-1} 

}

\caption{Histogram of library sizes for all cells}(\#fig:total-counts-hist-reads)
\end{figure}



\begin{table}

\caption{(\#tab:unnamed-chunk-11)The number of cells removed by total counts filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_total\_counts & Freq\\
\midrule
FALSE & 180\\
TRUE & 684\\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/total-features-hist-reads-1} 

}

\caption{Histogram of the number of detected genes in all cells}(\#fig:total-features-hist-reads)
\end{figure}



\begin{table}

\caption{(\#tab:unnamed-chunk-13)The number of cells removed by total features filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_expr\_features & Freq\\
\midrule
FALSE & 16\\
TRUE & 848\\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/total-features-vs-counts-reads-1} 

}

\caption{Library size vs number of detected genes}(\#fig:total-features-vs-counts-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/mt-vs-counts-reads-1} 

}

\caption{Percentage of counts in MT genes}(\#fig:mt-vs-counts-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/ercc-vs-counts-reads-1} 

}

\caption{Percentage of counts in ERCCs}(\#fig:ercc-vs-counts-reads)
\end{figure}



\begin{table}

\caption{(\#tab:unnamed-chunk-15)The number of cells removed by ERCC filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_ERCC & Freq\\
\midrule
FALSE & 103\\
TRUE & 761\\
\bottomrule
\end{tabular}
\end{table}



\begin{table}

\caption{(\#tab:unnamed-chunk-17)The number of cells removed by MT filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_MT & Freq\\
\midrule
FALSE & 18\\
TRUE & 846\\
\bottomrule
\end{tabular}
\end{table}



\begin{table}

\caption{(\#tab:unnamed-chunk-19)The number of cells removed by default filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 37\\
TRUE & 827\\
\bottomrule
\end{tabular}
\end{table}


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
## \begin{tabular}{l|r|r}
## \hline
##   & PC1 & PC2\\
## \hline
## pct\_counts\_feature\_controls & 0.5057646 & 0.2473134\\
## \hline
## pct\_counts\_top\_100\_features & 0.4888852 & 0.2277068\\
## \hline
## n\_detected\_feature\_controls & 0.0231277 & 0.6235516\\
## \hline
## log10\_counts\_feature\_controls & -0.1226860 & 0.6576822\\
## \hline
## total\_features & -0.4655518 & 0.2219694\\
## \hline
## log10\_counts\_endogenous\_features & -0.5223679 & 0.1278782\\
## \hline
## \end{tabular}
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/auto-cell-filt-reads-1} 

}

\caption{PCA plot used for automatic detection of cell outliers}(\#fig:auto-cell-filt-reads)
\end{figure}

\begin{table}

\caption{(\#tab:unnamed-chunk-20)The number of cells removed by automatic filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 753\\
TRUE & 111\\
\bottomrule
\end{tabular}
\end{table}



\begin{table}

\caption{(\#tab:unnamed-chunk-22)The number of cells removed by manual filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 254\\
TRUE & 610\\
\bottomrule
\end{tabular}
\end{table}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{06-exprs-qc-reads_files/figure-latex/cell-filt-comp-reads-1} 

}

\caption{Comparison of the default, automatic and manual cell filters}(\#fig:cell-filt-comp-reads)
\end{figure}


```
## Features  Samples 
##    16281      610
```


