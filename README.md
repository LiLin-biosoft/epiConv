# epiConv tutorial
Li Lin<br>

EpiConv is an algorithm to cluster scATAC-seq data and correct batch effects.

## Requirements
EpiConv is developed in R 3.5.1, the following packages are required:
`Matrix`,`SingleCellExperiment`,`Seurat`,`PRIMME`,`umap`,`bigmemory`,`biganalytics`<br>
They can be installed through CRAN or bioconductor.

We will align one co-assay data of postnatal day 0 mouse cortex onto scATAC-seq data to show the usage of epiConv. The data after basic quality control is stored in `sce.rds`.
```
source("epiConv_funs.R")
sce<-readRDS(file="sce.rds")
```
First, we need to combine RNA-seq and ATAC-seq data to resolve the relationships of single cells in co-assay data. If there are no RNA-seq profiles, this step can be skipped.
```
index_coassay<-which(colData(sce)$ident=="co-assay")
Smat<-run.epiConv(mat=assays(sce)$counts[,index_coassay],
                  lib_size=sce$lib_size[index_coassay],
                  nbootstrap=15,
                  nsample=floor(nrow(sce)*0.2),
                  bin=5000,
                  inf_replace=(-8))

feature_coassay<-cal.feature(Smat=Smat,
                                 pcs=sce@reducedDims$PCA_RNA[index_coassay,],
                                 neigs=50)
```
In the script above, we combined 50 PCs from RNA-seq `sce@reducedDims$PCA_RNA` and 50 Eigen vectors from ATAC-seq. `feature_coassay` contains 100 features of each cell of co-assay data.
+ `run.epiConv`: calculate the similarites between single cells.
  - `mat`: Matrix object constains the peak by cell matrix.
  - `lib_size`: library size of single cells.
  - `nbootstrap`: number of bootstraps performed.
  - `nsample`: number of peaks sampled in each bootstrap.
  - `bin`: many functions are with `bin` parameter. Given R does not support very long vectors, epiConv divides the matrix into several small matrix (e.g. `bin=10000` means 10,000 cells per matrix) and performs the calculation to avoid errors.
  - `inf_replace`: sometimes the similarity is -Inf, we use a small value to replace it.
  
+ `cal.feature`: calculate features from similarity matrix and combine them with PCs from RNA-seq.
  - `Smat`: the similarity matrix.
  - `pcs`: principal components from RNA-seq.
  - `neigs`: number of features to calculate.
  - the output contains the features from `pcs` and Eigen vectors from similarity matrix. Eigen vectors are scaled to make summed variance of Eigen vectors equal to that of PCs.

Next we align the co-assay data onto scATAC-seq reference.
```
Smat<-run.epiConv(mat=assays(sce)$counts,
                  lib_size=sce$lib_size,
                  nbootstrap=15,
                  nsample=floor(nrow(sce)*0.2),
                  bin=10000,
                  inf_replace=(-8))
sce@reducedDims$ATAC_umap<-run.umap_louvain(Smat=Smat,knn=20,
                                            umap_settings=NULL,resolution=NULL)
```
In the script above, we first calculated the simiarlity matrix from ATAC-seq data and perform umap. We can see obvious batch effects from the embeddings later. 
+ `run.umap_louvain`: perform umap and louvain clustering.
  - `Smat`: the similarity matrix.
  - `knn`: number of nearest neighbors for umap and louvain clustering.
  - `umap_settings`: custom umap settings by editing `umap::umap.defaults`.
  - `resolution`: the resolution parameter of louvain clustering. This parameter can be a numeric vector. When set to null, the function does not perform clustering.
  
Next we correct batch effects.
```
Smat<-Smat-median(Smat)
eigs<-dim.reduce(Smat,neigs=30)
residual_mat<-Smat-tcrossprod(eigs$vectors,t(t(eigs$vectors)*eigs$values))

batch<-factor(colData(sce)$ident)

guide_features<-sapply(levels(batch),function(x){
  index<-which(batch==x)
  eigs$vectors[index,]
})
names(guide_features)<-levels(batch)
```
We use Eigen Value Decomposition to deconvolute the similarity matrix into 30 eigen vectors and store the residuals in `residual_mat`. Batch correction is performed on Eigen vectors. These Eigen vectors are also used as guiding features.
+ `dim.reduce`: perform Eigen Value Decompostion.
  - `Smat`: the similarity matrix.
  - `neigs`: the number of Eigen vectors to calculate.
If we have co-assay data, PCs from RNA-seq are also used as guiding features.
```
guide_features[["co-assay"]]<-cbind(guide_features[["co-assay"]],
                                    sce@reducedDims$PCA_RNA[index_coassay,])  ##skip if there are no RNA-seq profiles.
```
Next we use the similarity matrix and guiding features to calculate the knn matrix.
```
knn_update<-eigs.knn(Smat=Smat,
                     features=guide_features,
                     batch=batch,
                     reference="ATAC",
                     knn_target=50,
                     knn_reference=10,
                     threshold=2)
```
+ `eigs.knn`: learn neighbors across multiple datasets.
  - `Smat`: the similarity matrix.
  - `features`: the guiding features.
  - `batch`: factor that contains batch information.
  - `reference`: which dataset in `batch` is used as reference. Reference datasets should contain all cell types. If there are more than one reference datasets (A, B and C), we can use `c("A","B","C")` to align B to A, and align C to A and B.
  - `knn_target`: number of neighbors for cells in target datasets to learn in each reference dataset.
  - `knn_reference`: number of neighbors for cells in reference datasets to learn in each target dataset.
  - `threshold`: the Z-score threshold used to filter false neighbors.
  - In output, `knn_update[["A"]][["A]]` contains the snn matrix of dataset A and `knn_update[["A"]][["B"]]` contains the knn matrix that cells from A pick their nearest neighbors in B.<br><br>

In the script above, the snn matrix is calculated from ATAC-seq. If the dataset contrains RNA-seq profiles, it is better to calculate the snn matrix based on both RNA-seq and ATAC-seq.
```
knn_update[["co-assay"]][["co-assay"]]<-cal.snn(Smat=as.matrix(dist(feature_coassay))*(-1),
                                                knn=floor(nrow(feature_coassay)*0.01))  ##skip if there are no RNA-seq profiles.
```
+ `cal.snn`: calculate shared nearest neighbors from similarity matrix.
  - `Smat`: the similarity matrix.
  - `knn`: number of neighbors. We generally set it to 1% of total cells.
Next we correct Eigen vectors and perform umap and louvain clustering.
```
eigs_corrected<-eigs.correct(eigs=eigs$vectors,
                             knn_mat=knn_update,
                             batch=batch,
                             reference="ATAC",
                             knn_transfer_correction=10)

Smat_corrected<-tcrossprod(eigs_corrected,t(t(eigs_corrected)*eigs$values))+residual_mat

temp<-run.umap_louvain(Smat=Smat_corrected,knn=20,umap_settings=NULL,resolution=0.8)
sce@reducedDims$ATACcorrected_umap<-as.matrix(temp[,1:2])
colData(sce)<-cbind(colData(sce),cluster=temp[,3])
```
+ `eigs.correct`: correct Eigen vectors.
  - `eigs`: the Eigen vectors.
  - `knn_mat`: the knn matrix.
  - `batch`: the batch information.
  - `reference`: the datasets used as reference, see above.
  - `knn_transfer_correction`: non-anchor cells need to learn the correction parameters from anchor cells. Here specify the number of nearest anchors used.

Now we see the results before and after batch correction.
```
plot(sce@reducedDims$ATAC_umap,pch="+",cex=0.5,col=batch)
plot(sce@reducedDims$ATACcorrected_umap,pch="+",cex=0.5,col=batch)
plot(sce@reducedDims$ATACcorrected_umap,pch="+",cex=0.5,
     col=rainbow(length(unique(sce$cluster)))[sce$cluster])
text(x=tapply(sce@reducedDims$ATACcorrected_umap[,1],list(sce$cluster),median),
     y=tapply(sce@reducedDims$ATACcorrected_umap[,2],list(sce$cluster),median),
     labels=levels(sce$cluster))
```

