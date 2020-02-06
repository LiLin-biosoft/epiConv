# epiConv tutorial
Li Lin<br>

EpiConv is a novel algorithm to cluster scATAC-seq data and detect differentially accessible (DE) peaks. We will demonstrate how to analyze scATAC-seq data by epiConv in this tutorial.

## Installation
EpiConv is developed under the following environments:
1. CentOS release 6.9 with GNU 4.9.1 and perl 5.10.1
2. bedtools v2.27.1
3. MACS2
5. R 3.5.1
  5.1 umap or uwot packages
  5.2 Matrix package
  5.3 bigmemory package
  5.4 biganalytics package
  
```{r }
install.packages("dbscan")
```
