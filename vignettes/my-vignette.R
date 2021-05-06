## ----echo=F,knitr-options,message=FALSE, warning=FALSE------------------------
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')
options(warn=-1)

## ----eval=FALSE---------------------------------------------------------------
#  suppressMessages(library(epiConv))
#  suppressMessages(library(SingleCellExperiment))
#  
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
#                    inf=(-8),
#                    backingfile="bmmc.backup",
#                    descriptorfile="bmmc.backup.descriptor")

## ----eval=FALSE---------------------------------------------------------------
#  reducedDims(sce)$umap<-run.umap_louvain(Smat=Smat,
#                                          knn=50,
#                                          resolution=NULL)

## ----eval=FALSE---------------------------------------------------------------
#  batch<-factor(sce$ident)
#  Smat_corrected<-deepcopy(Smat,
#                           backingfile="bmmc2.backup",
#                           descriptorfile="bmmc2.backup.descriptor")
#  res_joint<-epiConv.joint(Smat=Smat_corrected,
#                           batch=batch,
#                           reference="Resting",
#                           neigs=30,
#                           features=NULL,
#                           knn_mat=NULL,
#                           knn_target=50,
#                           knn_reference=10,
#                           threshold=2)
#  

## ----eval=FALSE---------------------------------------------------------------
#  temp<-run.umap_louvain(Smat=Smat_corrected,
#                         knn=50,
#                         resolution=0.8)
#  reducedDims(sce)$umap_corrected<-as.matrix(temp[,1:2])
#  colData(sce)<-cbind(colData(sce),cluster=factor(temp[,3]))
#  
#  plot(reducedDims(sce)$umap,pch="+",cex=0.5,col=factor(sce$ident))
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,col=factor(sce$ident))
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,col=rainbow(nlevels(sce$cluster))[sce$cluster])

## ----eval=FALSE---------------------------------------------------------------
#  suppressMessages(library(epiConv))
#  suppressMessages(library(SingleCellExperiment))
#  
#  sce<-readRDS(file="sce.rds")
#  index_coassay<-which(colData(sce)$ident=="co-assay")
#  Smat<-run.epiConv(mat=assays(sce)$counts[,index_coassay],
#                    lib_size=sce$lib_size[index_coassay],
#                    nbootstrap=15,
#                    nsample=floor(nrow(sce)*0.2),
#                    bin=5000,
#                    inf=(-8))
#  
#  feature_coassay<-cal.feature(Smat=Smat,
#                               pcs=reducedDims(sce)$PCA_RNA[index_coassay,],
#                               neigs=50)
#  knn_mat<-list()
#  knn_mat[["co-assay"]]<-find.knn(as.matrix(dist(feature_coassay))*(-1),
#                         knn=floor(nrow(feature_coassay)*0.01))
#  guide_features<-list()
#  guide_features[["co-assay"]]<-reducedDims(sce)$PCA_RNA[index_coassay,]

## ----eval=FALSE---------------------------------------------------------------
#  Smat<-run.epiConv(mat=assays(sce)$counts,
#                    lib_size=sce$lib_size,
#                    nbootstrap=15,
#                    nsample=floor(nrow(sce)*0.2),
#                    bin=10000,
#                    inf=(-8),
#                    backingfile="brain.backup",
#                    descriptorfile="brain.backup.descriptor")
#  reducedDims(sce)$umap<-run.umap_louvain(Smat=Smat,
#                                          knn=20,
#                                          resolution=NULL)
#  Smat_corrected<-deepcopy(Smat,
#                           backingfile="brain2.backup",
#                           descriptorfile="brain2.backup.descriptor")
#  batch<-factor(colData(sce)$ident)
#  res_joint<-epiConv.joint(Smat=Smat_corrected,
#                           batch=batch,
#                           reference="ATAC",
#                           neigs=30,
#                           features=guide_features,
#                           knn_mat=knn_mat,
#                           knn_target=50,
#                           knn_reference=10,
#                           threshold=2)

## ----eval=FALSE---------------------------------------------------------------
#  temp<-run.umap_louvain(Smat=Smat_corrected,
#                         knn=20,
#                         resolution=c(0.4,0.8)) ##Here we try resolution=0.4 and 0.8.
#  reducedDims(sce)$umap_corrected<-as.matrix(temp[,1:2])
#  colData(sce)<-cbind(colData(sce),cluster=factor(temp[,3]))
#  plot(reducedDims(sce)$umap,pch="+",cex=0.5,col=batch)
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,col=batch)
#  plot(reducedDims(sce)$umap_corrected,pch="+",cex=0.5,
#       col=rainbow(nlevels(sce$cluster))[sce$cluster])

