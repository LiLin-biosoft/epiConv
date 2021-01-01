## ----echo=F,knitr-options,message=FALSE, warning=FALSE------------------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')
options(warn=-1)

## ----eval=FALSE---------------------------------------------------------------
#  mat<-readMM(file="bmmc_matrix.mtx")
#  cell_info<-read.table(file="bmmc_ident.tsv")
#  colnames(cell_info)<-c("barcode","ident")
#  peak<-read.table(file="peaks.bed")
#  colnames(peak)<-c("seqnames","start","end")
#  rownames(mat)<-paste(peak$seqnames,":",
#                       peak$start,"-",
#                       peak$end,sep="")
#  colnames(mat)<-cell_info[,1]
#  sce<-SingleCellExperiment(assays=SimpleList(counts=mat),
#                            rowRanges=GRanges(peak),colData=cell_info)
#  colData(sce)$lib_size<-lib.estimate(assays(sce)$counts)
#  sce<-subset(sce,select=colData(sce)$lib_size>500)
#  sce<-subset(sce,subset=freq.estimate(assays(sce)$counts)>0)

## ----eval=FALSE---------------------------------------------------------------
#  Smat<-run.epiConv(mat=assays(sce)$counts,
#                    lib_size=sce$lib_size,
#                    nbootstrap=15,
#                    nsample=floor(nrow(sce)*0.2),
#                    bin=5000,
#                    inf_replace=(-8))

## ----eval=FALSE---------------------------------------------------------------
#  reducedDims(sce)$umap<-run.umap_louvain(Smat=Smat,
#                                          knn=50,
#                                          umap_settings=umap::umap.defaults,
#                                          resolution=NULL)

## ----eval=FALSE---------------------------------------------------------------
#  Smat<-Smat-median(Smat)
#  eigs<-dim.reduce(Smat,neigs=30)
#  residual_mat<-Smat-tcrossprod(eigs$vectors,t(t(eigs$vectors)*eigs$values))

## ----eval=FALSE---------------------------------------------------------------
#  batch<-factor(sce$ident)
#  guide_features<-sapply(levels(batch),function(x){
#    index<-which(batch==x)
#    eigs$vectors[index,]
#  })
#  knn_update<-eigs.knn(Smat=Smat,
#                       features=guide_features,
#                       batch=batch,
#                       reference="Resting",
#                       knn_target=50,
#                       knn_reference=10,
#                       threshold=2)

## ----eval=FALSE---------------------------------------------------------------
#  eigs_corrected<-eigs.correct(eigs=eigs$vectors,
#                               knn_mat=knn_update,
#                               batch=batch,
#                               reference="Resting",
#                               knn_transfer_correction=10)

## ----eval=FALSE---------------------------------------------------------------
#  Smat_corrected<-tcrossprod(eigs_corrected,t(t(eigs_corrected)*eigs$values))+residual_mat
#  temp<-run.umap_louvain(Smat=Smat_corrected,
#                         knn=50,
#                         umap_settings=umap::umap.defaults,
#                         resolution=0.8)
#  reducedDims(sce)$umap_corrected<-as.matrix(temp[,1:2])
#  colData(sce)<-cbind(colData(sce),cluster=factor(temp[,3]))
#  
#  plot(reducedDims(sce)$umap,pch="+",cex=0.5,col=factor(sce$ident))
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,col=factor(sce$ident))
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,col=rainbow(nlevels(sce$cluster))[sce$cluster])

## ----eval=FALSE---------------------------------------------------------------
#  sce<-readRDS(file="sce.rds")
#  index_coassay<-which(colData(sce)$ident=="co-assay")
#  Smat<-run.epiConv(mat=assays(sce)$counts[,index_coassay],
#                    lib_size=sce$lib_size[index_coassay],
#                    nbootstrap=15,
#                    nsample=floor(nrow(sce)*0.2),
#                    bin=5000,
#                    inf_replace=(-8))
#  
#  feature_coassay<-cal.feature(Smat=Smat,
#                               pcs=reducedDims(sce)$PCA_RNA[index_coassay,],
#                               neigs=50)

## ----eval=FALSE---------------------------------------------------------------
#  Smat<-run.epiConv(mat=assays(sce)$counts,
#                    lib_size=sce$lib_size,
#                    nbootstrap=15,
#                    nsample=floor(nrow(sce)*0.2),
#                    bin=10000,
#                    inf_replace=(-8))
#  reducedDims(sce)$umap<-run.umap_louvain(Smat=Smat,
#                                          knn=NULL,
#                                          umap_settings=NULL,
#                                          resolution=NULL)
#  Smat<-Smat-median(Smat)
#  eigs<-dim.reduce(Smat,neigs=30)
#  residual_mat<-Smat-tcrossprod(eigs$vectors,t(t(eigs$vectors)*eigs$values))
#  
#  batch<-factor(colData(sce)$ident)
#  guide_features<-sapply(levels(batch),function(x){
#    index<-which(batch==x)
#    eigs$vectors[index,]
#  })

## ----eval=FALSE---------------------------------------------------------------
#  guide_features[["co-assay"]]<-cbind(guide_features[["co-assay"]],
#                                      reducedDims(sce)$PCA_RNA[index_coassay,])

## ----eval=FALSE---------------------------------------------------------------
#  knn_update<-eigs.knn(Smat=Smat,
#                       features=guide_features,
#                       batch=batch,
#                       reference="ATAC",
#                       knn_target=50,
#                       knn_reference=10,
#                       threshold=2)

## ----eval=FALSE---------------------------------------------------------------
#  knn_update[["co-assay"]][["co-assay"]]<-cal.snn(Smat=as.matrix(dist(feature_coassay))*(-1),
#                                                  knn=floor(nrow(feature_coassay)*0.01))

## ----eval=FALSE---------------------------------------------------------------
#  eigs_corrected<-eigs.correct(eigs=eigs$vectors,
#                               knn_mat=knn_update,
#                               batch=batch,
#                               reference="ATAC",
#                               knn_transfer_correction=10)
#  
#  Smat_corrected<-tcrossprod(eigs_corrected,t(t(eigs_corrected)*eigs$values))+residual_mat
#  
#  temp<-run.umap_louvain(Smat=Smat_corrected,
#                         knn=NULL,
#                         umap_settings=NULL,
#                         resolution=c(0.4,0.8)) ##Here we try resolution=0.4 and 0.8.
#  reducedDims(sce)$umap_corrected<-as.matrix(temp[,1:2])
#  colData(sce)<-cbind(colData(sce),cluster=factor(temp[,3]))
#  plot(reducedDims(sce)$umap,pch="+",cex=0.5,col=batch)
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,col=batch)
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,
#       col=rainbow(nlevels(sce$cluster))[sce$cluster])

