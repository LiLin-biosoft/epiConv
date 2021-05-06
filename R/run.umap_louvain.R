#' perform umap and louvain clustering
#'
#' perform umap and louvain clustering.
#'
#' @param Smat the similarity matrix.
#' @param knn number of neighbors to calculate shared nearest neighbors matrix. If set to NULL, it will be total cells*0.01.
#' @param umap_settings umap settings passed to umap.
#' @param method method can be umap or uwot that specifies the package used.
#' @param resolution the resolution of louvain clustering. It can be numeric vectors. If set to NULL, the function does not perform clustering.
#' @param ... Arguments passed to big.matrix when creating the distance matrix.
#'
#'
#' @examples
#' \dontrun{
#'res<-run.umap_louvain(Smat=Smat,
#'                      knn=NULL,
#'                      knn_umap=20,
#'                      resolution=NULL)
#'}
#'
#' @return Return a data frame. The first two columns are coordinates of umap embeddings and other columns are results of clustering.

run.umap_louvain<-function(Smat,knn=NULL,umap_settings=NULL,method="umap",resolution=NULL,...){
  idx_valid<-which(!is.na(Smat[,1]))
  nvalid<-length(idx_valid)
  if(is.null(knn)){
    knn<-max(min(floor(nrow(Smat)*0.01),200),20)
  }
  knn_mat<-find.knn(Smat,knn=knn)
  knn_mat<-knn_mat[idx_valid,idx_valid]
  bin<-10000
  dis<-big.matrix(nrow=nvalid,ncol=nvalid,init=0,...)
  for(i in 1:ceiling(nvalid/bin)-1){
    row_index<-(i*bin+1):min((i*bin+bin),nvalid)
    for(j in (i+1):ceiling(nvalid/bin)-1){
      col_index<-(j*bin+1):min((j*bin+bin),nvalid)
      Smat_small<-0.01-Smat[idx_valid[row_index],idx_valid[col_index]]/10000
      Smat_small[is.na(Smat_small)]<-0.01
      temp<-knn-as.matrix(Matrix::crossprod(knn_mat[,row_index],knn_mat[,col_index]))+
        Smat_small
      dis[row_index,col_index]<-temp
      if(i!=j){
        dis[col_index,row_index]<-t(temp)
      }
    }
  }
  if(is.null(umap_settings)){
    umap_settings<-umap::umap.defaults
  }
  umap_settings$input<-"dist"

  if(!is.null(resolution)){
    graph<-new("dgCMatrix")
    graph@Dim<-as.integer(dim(dis))
    graph@p<-as.integer(0)
    nelement<-as.integer(0)
    for(i in 1:ncol(dis)){
      temp<-knn+0.01-dis[,i]
      temp[i]<-0
      temp<-graph.norm(temp,knn=umap_settings$n_neighbors)
      index<-sort(which(temp!=0))
      graph@i<-c(graph@i,as.integer(index-1))
      nelement<-nelement+length(index)
      graph@p<-c(graph@p,nelement)
      graph@x<-c(graph@x,temp[index])

    }
    graph<-graph+Matrix::t(graph)-graph*Matrix::t(graph)
    rownames(graph)<-as.character(1:nvalid)
    colnames(graph)<-as.character(1:nvalid)

    graph<-Seurat::as.Graph(graph)
  }


  if(method=="umap"){
    temp<-dis[1:nvalid,1:nvalid]
    rm(dis)
    umap_res<-umap::umap(temp,config=umap_settings)$layout
    rm(temp)
    gc()
  }else{
    temp<-biganalytics::apply(dis,2,function(x){
      aa<-order(x)[1:umap_settings$n_neighbors]
      bb<-x[aa]
      return(c(aa,bb))
    })
    rm(dis)
    temp<-list(idx=t(temp[1:umap_settings$n_neighbors,]),
               dist=t(temp[(umap_settings$n_neighbors+1):(2*umap_settings$n_neighbors),]))
    umap_res<-uwot::umap(X=NULL,nn_method=temp)
    rm(temp)
    gc()
  }

  if(is.null(resolution)){
    return(umap_res[match(1:nrow(Smat),idx_valid),])
  }

  clust<-matrix("",nvalid,length(resolution))
  for(i in 1:length(resolution)){
    temp<-Seurat::FindClusters(graph,resolution=resolution[i])[,1]
    clust[,i]<-sprintf(paste0("%0",ceiling(log10(nlevels(temp))),"d"),temp)
  }
  output<-data.frame(umap_res,clust)
  colnames(output)<-c("umap_1","umap_2",paste0("resolution_",resolution))
  return(output[match(1:nrow(Smat),idx_valid),])
}
