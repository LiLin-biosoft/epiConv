# epiConv tutorial
Li Lin<br>

EpiConv is a novel algorithm to cluster scATAC-seq data and detect differentially accessible (DE) peaks. We will demonstrate how to analyze scATAC-seq data by epiConv in this tutorial.

## Installation
EpiConv is developed under the following environments:
1. CentOS release 6.9 with perl 5.10.1
2. bedtools v2.27.1
3. MACS2
4. R 3.5.1<br>
  4.1 Matrix package<br>
  4.2 bigmemory package<br>
  4.3 biganalytics package<br>
  4.4 umap or uwot packages<br>

All source codes of epiConv can be found in the folder `source/`. Here assume we put them into the folder `~/epiConv/`. Two source files need to be compiled.<br>
```
g++ ~/epiConv/matrix_sum.c -o ~/epiConv/matrix_sum
g++ ~/epiConv/matrix_sampl.c -o ~/epiConv/matrix_sampl
```
## Usage
We will use one dataset of [PBMCs](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_pbmc_5k_v1) from 10X Genomics as example. Two files are required: `Peak by cell matrix (filtered)` and `Fragments (TSV)`. Here we put them into folder `pbmc5k/`, extract the first file and rename the second file as `pbmck5k_frag.bed.gz`.

There are two versions of epiConv: epiConv-full and epiConv-simp. EpiConv-full calculates the similarities between cells from raw Tn5 insertion profiles and epiConv-simp calculates the similarities from binary matrix. We first show the analysis pipeline for epiConv-simp. It is an implemention in R. First we read the source file of epiConv and the data:
```
source("~/epiConv/epiConv_functions.R")
mat<-readMM(file="~/pbmc5k/matrix.mtx")
barcodes<-read.table(file="~/pbmc5k/barcodes.tsv",colClass="character")[,1]
peaks<-read.table(file="~/pbmc5k/peaks.bed")

```








First we use `peak_calling.sh` to call high density regions of Tn5 insertions:
```
peak_calling.sh <prefix> <extsize> <fraction of data retained>
```
`<prefix>`: the prefix of data.<br>
`<extsize>`: the same parameter in `MACS2 pileup`. For example, set `<extsize>` to 100 will make MACS2 extend each insertion from both directions by 100bp before pileup insertions.<br>
`<fraction of data retained>`: specify the fraction of data you want to use in downstream analysis. Based on our preliminary analysis, this parameter need not to be accurately specified but should be close to the fraction of fragments from nucleosome-free regions. For example, you can check the histogram of insertion length from 10X Genomics reporting summary to learn this parameter (the fraction of fragments in first peak in the histogram).<br>
For the PBMC dataset, we use the following command:
```
~/epiConv/peak_calling.sh pbmc5k/pbmc5k 100 0.7
```
When it is finished (~4 hours), there will be a new file `pbmc5k/pbmc5k_peak.bed` in bed format containing the peaks. In order to run epiConv in parallel, we need to split the job using the linux command `split`:
```
split -a2 -d -l127000 pbmc5k/pbmc5k_peak.bed pbmc5k/pbmc5k_peak.run
```
Here we split the peak file into 10 jobs, each containing 127000 peaks. Then we use `convolution.sh` to calculate the similarites between single cells:
```

```
