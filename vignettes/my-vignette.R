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
#  sce<-SingleCellExperiment(assays=SimpleList(counts=as(mat,"dgCMatrix")),
#                            rowRanges=GRanges(peak),colData=cell_info)
#  colData(sce)$lib_size<-lib.estimate(assays(sce)$counts)
#  sce<-subset(sce,select=colData(sce)$lib_size>500)

## ----eval=FALSE---------------------------------------------------------------
#  
#  row_sample<-1:ncol(sce)
#  Smat<-run.epiConv(mat=assays(sce)$counts,
#                    row_sample=row_sample,
#                    lib_size=sce$lib_size,
#                    nbootstrap=3,
#                    nsample=floor(nrow(sce)*0.5),
#                    bin=1000,
#                    inf=(-10),
#                    backingfile="bmmc.backup",
#                    descriptorfile="bmmc.backup.descriptor")

## ----eval=FALSE---------------------------------------------------------------
#  snn_mat<-Smat2snn(Smat,knn=20)
#  dis<-snn_mat
#  dis@x<-max(dis@x)-dis@x+1e-9
#  umap_res<-uwot::umap(X=dis,n_neighbors=20)
#  

## ----eval=FALSE---------------------------------------------------------------
#  batch<-factor(sce$ident)
#  Smat_corrected<-deepcopy(Smat,
#                           backingfile="bmmc2.backup",
#                           descriptorfile="bmmc2.backup.descriptor")
#  res_anchor<-epiConv.anchor(Smat=Smat_corrected,
#                             row_sample=row_sample,
#                             batch=batch,
#                             reference="Resting",
#                             neigs=30,
#                             features=NULL,
#                             knn_target=50,
#                             knn_reference=20,
#                             threshold=2)

## ----eval=FALSE---------------------------------------------------------------
#  res_eigs<-epiConv.correct(Smat=Smat_corrected,
#                            row_sample=row_sample,
#                            batch=batch,
#                            reference="Resting",
#                            neigs=30,
#                            knn_update=res_anchor$knn_update)
#  

## ----eval=FALSE---------------------------------------------------------------
#  snn_mat_corrected<-Smat2snn(Smat_corrected,knn=20)
#  dis<-snn_mat_corrected
#  dis@x<-max(dis@x)-dis@x+1e-9
#  umap_res_corrected<-uwot::umap(X=dis,n_neighbors=20)
#  
#  plot(umap_res,pch="+",cex=0.5,col=factor(sce$ident))
#  plot(umap_res_corrected,pch="+",cex=0.5,col=factor(sce$ident))

## ----eval=FALSE---------------------------------------------------------------
#  clust<-epiConv.louvain(snn=snn_mat_corrected,resolution=c(0.8,0.6,0.4,0.2))
#  head(clust)

