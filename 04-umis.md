---
# knit: bookdown::preview_chapter
output: html_document
---



# Unique Molecular Identifiers (UMIs)

Thanks to Andreas Buness from EMBL Monterotondo for collaboration on this section.

## Introduction

Unique Molecular Identifiers are short (4-10bp) random barcodes added to transcripts during reverse-transcription. They enable sequencing reads to be assigned to individual transcript molecules and thus the removal of amplification noise and biases from scRNASeq data. 

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-protocol.png" alt="UMI sequencing protocol" width="90%" />
<p class="caption">(\#fig:intro-umi-protocol)UMI sequencing protocol</p>
</div>

When sequencing UMI containing data, techniques are used to specifically sequence only the end of the transcript containing the UMI (usually the 3' end).

## Mapping Barcodes

There is a limited number of unique barcodes (4^N, N=length of UMI) which is much smaller than the total number of molecules per cell (~10^6) hence both barcode and mapping location must be used to identify unique molecules. Thus, the first step is to map UMI reads, we recommend using STAR since it is fast and outputs good quality BAM-alignments and mapping locations can be useful for eg. identifying poorly-annotated 3' UTRs of transcripts.

UMI-sequencing generally consists of paired-end reads where one read from each pair captures the cell and UMI barcodes while the other read consists of exonic sequence from the transcript (Figure \@ref(fig:intro-umi-reads)). Note that trimming and/or filtering to remove reads containing poly-A sequence is recommended to avoid erors due to these read mapping to genes/transcripts with internal poly-A/poly-T sequences.

__Conventions__

1. UMI added to the read name of the other paired read. 

2. Reads separated into separate files by cell barcode
	+ For extremely large, shallow datasets, cell barcode may be added to the read name as well.

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-reads.png" alt="UMI sequencing reads, red lightning bolts represent different fragmentation locations" width="90%" />
<p class="caption">(\#fig:intro-umi-reads)UMI sequencing reads, red lightning bolts represent different fragmentation locations</p>
</div>

## Counting Barcodes

In theory, every unique UMI-transcript pair should represent all reads originating from a single transcript. However, in practice this is frequently not the case. 

1. __Different UMI doesn't necessarily mean different molecule__
	+ Base-pair substitutions due to PCR or sequencing errors can create new UMI sequences which originate from the same transcript. Longer UMIs give more opportunity for errors to arise, eg. based on cell barcodes we expect 7-10% of 10bp UMIs to contain at least one error.

2. __Different transcript doesn't necessarily mean different molecule__
	+ Mapping errors and/or multimapping reads may result in some UMIs being assigned to the wrong gene/transcript.

3. __Same UMI doesn't necessarily mean same molecule__
	+ Biases in UMI frequency and short UMIs can result in the same UMI being attached to different mRNA molecules from the same gene.

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-errors.png" alt="Potential Errors in UMIs" width="90%" />
<p class="caption">(\#fig:intro-umi-errors)Potential Errors in UMIs</p>
</div>

## Correcting for Errors

1. [UMI-tools'](https://github.com/CGATOxford/UMI-tools) directional-adjacency method implements a procedure which considers both the number of mismatches and the relative frequency of similar UMIs to identify likely PCR/sequencing errors.

2. Currently an open question. May be mitigated by removing UMIs with few reads to support their association with a particular gene/transcript, or by removing all multi-mapping reads.

3. Simple saturation (aka "collision probability") correction proposed by [Grun, Kester and van Oudenaarden (2014)](http://www.nature.com/nmeth/journal/v11/n6/full/nmeth.2930.html#methods) :

$$True \approx -N*ln(1 - \frac{Obs}{N})$$ 
where N = total unique UMI barcodes.
	+ Assumes all UMIs are equally frequent which in most cases is incorrect, baises w.r.t GC content have been observed. 

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-amp.png" alt="Per gene amplification rate" width="60%" />
<p class="caption">(\#fig:intro-umi-amp)Per gene amplification rate</p>
</div>

## Exercise 1
You have been provided with unique molecule counts from Blischak et al., they used 6bp UMI barcodes. Load this data using the command below.


```r
molecules <- read.table("blischak/molecules.txt", sep = "\t")
```

Correct this data for collisions and sequencing errors assuming a 1% per base-pair sequencing error rate.

## Downstream Analysis

Current UMI platforms (DropSeq, InDrop, ICell8) exhibit low and highly variable capture efficiency (Figure \@ref(fig:intro-umi-capture)). This can introduce strong biases in downstream analysis. Recent analyses often pool cells/genes together based on cell-type or biological pathway to get meaningful results. Robust statistical analyses of this data is still an open research question.

<div class="figure" style="text-align: center">
<img src="figures/UMI-Seq-capture.png" alt="Variability in Capture Efficiency" width="70%" />
<p class="caption">(\#fig:intro-umi-capture)Variability in Capture Efficiency</p>
</div>
## Exercise 2
You have been provided with the read counts from the Blischak data, it can be loaded using the command below.


```r
reads <- read.table("blischak/reads.txt", sep = "\t")
```
Using this data and the unadjusted molecule counts from above:

1. Plot the variability in capture efficiency

2. Determine the amplification rate: average number of reads per UMI.


