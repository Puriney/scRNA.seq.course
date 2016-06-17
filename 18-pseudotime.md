---
# knit: bookdown::preview_chapter
output: html_document
---

# Ordering of cells according to pseudotime



```r
library(TSCAN)
library(M3Drop)
library(monocle)
set.seed(1)
```

In many situations, one is studying a process where cells change
continuously. This includes for example many differentiation processes
taking place during development, where following a stimulus, cells
will change from one cell-type to another. Ideally, we would like to
monitor the expression levels of an individual cell over
time. Unfortunately, such monitoring is not possible with scRNA-seq
since the cell is lysed (destroyed) when the RNA is extracted.

Instead, we must sample at multiple time-points and obtain snapshots of
the gene expression profiles. Since some of the cells will
proceed faster along the differentiation than others, each snapshot may contain
cells at varying points along the developmental progression. We use statistical methods
 to order the cells along one or more trajectories which represent the underlying
 developmental trajectories, this ordering is referred to a "pseudotime".

In this chapter we will consider two different tools, Monocle and
TSCAN for ordering cells according to their pseudotime development. To
illustrate the methods we will be using a dataset on mouse embryonic
development that was collected by Deng et al. The dataset consists of
268 cells from 10 different time-points of early mouse development.

## TSCAN

TSCAN combines clustering with pseudotime analysis. First it clustering the cells using `mclust`, 
which is based on a mixture of normal distributions. Then it builds a minimum spanning tree to connect the clusters together. The branch of this tree which connects the largest number of clusters is the main branch which is used to determine pseudotime.

First we will try to use all genes to order the cells.

```r
deng <- readRDS("deng/deng.rds")
cellLabels <- colnames(deng)
procdeng <- TSCAN::preprocess(deng, clusternum=10)
colnames(procdeng) <- 1:268
dengclust <- TSCAN::exprmclust(procdeng, clusternum=10)
TSCAN::plotmclust(dengclust)
```

<img src="18-pseudotime_files/figure-html/tscan-all-genes-1.png" width="672" style="display: block; margin: auto;" />

```r
dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly=F)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
pseudotime_order_tscan
```

```
##   [1] "9"   "12"  "3"   "2"   "8"   "4"   "1"   "6"   "10"  "11"  "7"  
##  [12] "5"   "16"  "15"  "18"  "14"  "17"  "19"  "23"  "24"  "21"  "20" 
##  [23] "22"  "29"  "28"  "30"  "32"  "13"  "31"  "27"  "34"  "25"  "33" 
##  [34] "26"  "40"  "47"  "45"  "46"  "48"  "37"  "43"  "38"  "44"  "39" 
##  [45] "41"  "81"  "85"  "111" "97"  "104" "96"  "128" "87"  "93"  "71" 
##  [56] "99"  "114" "116" "60"  "64"  "105" "117" "100" "119" "51"  "115"
##  [67] "108" "50"  "106" "89"  "132" "112" "58"  "75"  "63"  "66"  "65" 
##  [78] "61"  "133" "120" "59"  "94"  "68"  "52"  "80"  "82"  "79"  "67" 
##  [89] "91"  "88"  "124" "56"  "122" "62"  "121" "126" "125" "69"  "113"
## [100] "123" "84"  "36"  "35"  "42"  "57"  "90"  "86"  "107" "92"  "70" 
## [111] "135" "98"  "78"  "73"  "76"  "53"  "95"  "74"  "129" "55"  "49" 
## [122] "130" "77"  "127" "54"  "134" "83"  "208" "228" "173" "205" "153"
## [133] "142" "72"  "164" "207" "152" "179" "227" "165" "144" "169" "219"
## [144] "171" "262" "170" "199" "253" "198" "200" "155" "212" "158" "197"
## [155] "258" "157" "209" "156" "194" "188" "189" "252" "160" "191" "265"
## [166] "187" "190" "195" "220" "263" "161" "213" "214" "222" "148" "196"
## [177] "225" "259" "162" "183" "251" "193" "143" "182" "140" "206" "141"
## [188] "201" "249" "261" "181" "239" "223" "255" "159" "260" "216" "151"
## [199] "145" "235" "215" "211" "240" "268" "267" "257" "245" "256" "250"
## [210] "217" "244" "237" "254" "233" "186" "234" "236" "246" "203" "180"
## [221] "177" "149" "136" "185" "184" "231" "264" "178" "243" "204" "230"
## [232] "176" "147" "202" "229" "154" "247" "172" "146" "232" "167" "238"
## [243] "218" "266" "174" "168" "248" "241" "210" "242" "163" "224" "226"
## [254] "221" "139" "150" "166"
```

We can also examine which cells have been assigned to each state:


```r
cellLabels[dengclust$clusterid == 1]
```

```
##  [1] "1" "1" "1" "1" "2" "2" "2" "2" "2" "2" "2" "2"
```
__Exercise__ Compare results for different numbers of clusters (`clusternum`).

## monocle

Monocle skips the clustering stage of TSCAN and directly builds a minimum spanning tree to connect all cells. Monocle then identifies the longest path in this tree to determine pseudotime. If the data contains diverging trajectories (i.e. one cell type differentiation into two different cell-types), monocle can identify alternative long paths in the tree using the argument `num_paths`. Each of the resulting forked paths is defined as a separate cell state, thus `num_paths=2` will identify three different cell states.

Unfortunately, Monocle does not work when all the genes are used, so
we must carry out feature selection. First, we use M3Drop:

```r
m3dGenes <- as.character(
    M3Drop::M3Drop_Differential_Expression(deng, suppress.plot=T)$Gene
)
d <- deng[which(m3dGenes %in% rownames(deng)),]
d <- d[!duplicated(rownames(d)),]
```


```r
pd <- as.data.frame(colnames(d))
names(pd) <- "timepoint"
pd <- new("AnnotatedDataFrame", data=pd)
fd <- as.data.frame(rownames(d))
names(fd) <- "gene"
fd <- new("AnnotatedDataFrame", data=fd)
colnames(d) <- 1:dim(d)[2]
rownames(d) <- 1:dim(d)[1]
dCellData <- monocle::newCellDataSet(d, phenoData=pd, featureData=fd)
dCellData <- monocle::setOrderingFilter(dCellData, 1:length(m3dGenes))
dCellDataSet <- monocle::reduceDimension(dCellData, pseudo_expr=1)
```

```
## Reducing to independent components
```

```r
dCellDataSet <- monocle::orderCells(dCellDataSet, reverse=T, num_paths=2)
monocle::plot_spanning_tree(dCellDataSet)
```

<img src="18-pseudotime_files/figure-html/monocle-all-genes-1.png" width="672" style="display: block; margin: auto;" />

```r
pseudotime_monocle <- data.frame(Cell = cellLabels, 
                                time = phenoData(dCellDataSet)$Pseudotime, 
                                State=phenoData(dCellDataSet)$State)
pseudotime_order_monocle <- rownames(pseudotime_monocle[order(pseudotime_monocle$time), ])
pseudotime_order_monocle
```

```
##   [1] "192" "263" "139" "153" "173" "137" "221" "266" "222" "259" "242"
##  [12] "246" "163" "260" "224" "166" "188" "238" "227" "248" "264" "234"
##  [23] "231" "261" "258" "243" "241" "268" "233" "247" "226" "252" "262"
##  [34] "267" "237" "150" "230" "175" "138" "256" "235" "257" "168" "236"
##  [45] "195" "250" "232" "265" "254" "255" "253" "178" "229" "223" "190"
##  [56] "196" "148" "152" "164" "147" "170" "155" "156" "140" "154" "186"
##  [67] "228" "167" "203" "245" "219" "174" "240" "244" "143" "180" "185"
##  [78] "218" "209" "191" "184" "239" "193" "199" "172" "151" "145" "179"
##  [89] "136" "189" "204" "187" "181" "197" "220" "225" "176" "157" "144"
## [100] "198" "182" "183" "200" "160" "146" "165" "159" "207" "161" "194"
## [111] "158" "142" "177" "171" "217" "249" "215" "251" "216" "211" "201"
## [122] "202" "162" "149" "210" "213" "141" "214" "206" "169" "212" "208"
## [133] "205" "131" "72"  "118" "101" "125" "135" "102" "110" "117" "124"
## [144] "133" "128" "109" "127" "120" "113" "134" "122" "129" "121" "126"
## [155] "119" "130" "132" "108" "105" "115" "114" "112" "100" "99"  "88" 
## [166] "123" "97"  "91"  "111" "66"  "70"  "96"  "71"  "94"  "63"  "98" 
## [177] "95"  "51"  "73"  "74"  "116" "92"  "61"  "55"  "54"  "56"  "93" 
## [188] "53"  "52"  "86"  "49"  "87"  "50"  "89"  "75"  "64"  "60"  "62" 
## [199] "76"  "79"  "67"  "57"  "85"  "104" "83"  "58"  "90"  "107" "106"
## [210] "68"  "65"  "59"  "77"  "84"  "82"  "69"  "78"  "80"  "81"  "38" 
## [221] "39"  "35"  "41"  "42"  "44"  "43"  "37"  "40"  "36"  "45"  "46" 
## [232] "47"  "48"  "26"  "30"  "33"  "34"  "25"  "27"  "29"  "28"  "31" 
## [243] "32"  "17"  "18"  "13"  "20"  "24"  "19"  "14"  "22"  "23"  "21" 
## [254] "15"  "16"  "5"   "8"   "2"   "1"   "103" "4"   "3"   "9"   "10" 
## [265] "11"  "12"  "7"   "6"
```

## Comparison of the methods

How do the trajectories inferred by TSCAN and Monocle compare?

```r
matched_ordering <- match(pseudotime_order_tscan, pseudotime_order_monocle)
plot(matched_ordering, xlab = "Monocle Order", ylab = "TSCAN Order")
```

<img src="18-pseudotime_files/figure-html/tscan-monocle-compare-1.png" width="672" style="display: block; margin: auto;" />

__Exercise__: Repeat the exercise using a subset of the genes, e.g. the set of highly variable genes that can be obtained using M3Drop::Brennecke_getVariableGenes

