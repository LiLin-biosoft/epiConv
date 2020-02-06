# epiConv tutorial
Li Lin<br>

EpiConv is a novel algorithm to cluster scATAC-seq data and detect differentially accessible (DE) peaks. We will demonstrate how to analyze scATAC-seq data by epiConv in this tutorial.

## Installation
EpiConv is developed under the following environments:
1. CentOS release 6.9 with GNU 4.9.1 and perl 5.10.1
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
