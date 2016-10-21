---
output: html_document
---

# "Ideal" scRNAseq pipeline



## Experimental Design

* Avoid confounding biological and batch effects (Figure \@ref(fig:pipeline-batches))
    * Multiple conditions should be captured on the same chip if possible
    * Perform multiple replicates of each condition where replicates of different conditions should performed together if possible
    * Statistics cannot correct a completely confounded experiment!
* Unique molecular identifiers
    * Use longer UMIs (~10bp)
    * Greatly reduce noise in data
    * May reduce gene detection rates (unclear if it is UMIs or other protocol differences)
    * Correct for sequencing errors in umis using UMI-tools
* Spike-ins
    * Useful for quality control
    * Can be used to approximate cell-size/RNA content (if relevant to biological question)
    * Often exhibit higher noise than endogenous genes (pipetting errors, mixture quality)
* Cell number vs Read depth
    * Gene detection plateaus starting from 1 million reads per cell
    * Transcription factor detection (regulatory networks) require high read depth and most sensitive protocols (i.e. Fluidigm C1)
    * Cell clustering & cell-type identification benefits from large number of cells and doesn't requireas high sequencing depth (~100,000 reads per cell).

<div class="figure" style="text-align: center">
<img src="figures/Pipeline-batches.png" alt="Appropriate approaches to batch effects in scRNASeq. Red arrows indicate batche effects which are (pale) or are not (vibrant) correctable through batch-correction." width="90%" />
<p class="caption">(\#fig:pipeline-batches)Appropriate approaches to batch effects in scRNASeq. Red arrows indicate batche effects which are (pale) or are not (vibrant) correctable through batch-correction.</p>
</div>

## Read QC & Trimming

FASTQC, cutadapt

## Mapping

* Small datasets or UMI datasets: align to genome/transcriptome using STAR
* Large datasets: pseudo-alignment with Salmon

## Quantification

* Small dataset, no UMIs : Cufflinks
* Large datasets, no UMIs: Salmon
* UMI dataset : UMI-Tools + featureCounts

## Cell QC

scater, consider mtRNA, rRNA, spike-ins (if available), number of detected genes per cell, total reads/molecules per cell 

## Library Size Normalization

scran

## Batch correction (if appropriate, see below)

RUVs

## Feature Selection

M3Drop

## Clustering and Marker genes

SC3

## Pseudotime

* Few or very distinct timepoints: Tscan
* Very continuous data : dpt

## Differential Expression between groups

* Fewer than 400 cells : SCDE
* Replicates with batch effects : mixture models
* Large datasets : Kruskal-Wallis test (many groups, 1 test), or Kolmogorov-Smirnov (KS)-test (compare 2 groups at a time)
