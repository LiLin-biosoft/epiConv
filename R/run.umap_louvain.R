#' perform umap and louvain clustering
#'
#' perform umap and louvain clustering.
#'
#' @param Smat the similarity matrix.
#' @param knn number of neighbors to calculate shared nearest neighbors matrix. If set to NULL, it will be total cells*0.01.
#' @param umap_settings umap settings. If set to NULL, it will be defaults of umap.
#' @param resolution the resolution of louvain clustering. It can be numeric vectors. If set to NULL, the function does not perform clustering.
#'
#' @examples
#' \dontrun{
#'res<-run.umap_louvain(Smat=Smat,
#'                      knn=NULL,
#'                      umap_settings=NULL,
#'                      resolution=NULL)
#'}
#'
#' @return Return a data frame. The first two columns are coordinates of umap embeddings and other columns are results of clustering.

run.umap_louvain<-function(Smat,knn=NULL,umap_settings=NULL,resolution=NULL){
  if(is.null(knn)){
    knn<-max(min(floor(nrow(Smat)*0.01),200),20)
  }
  snn_mat<-cal.snn(Smat,knn=knn)+Smat/1e4
  if(is.null(umap_settings)){
    umap_settings<-umap::umap.defaults
  }
  umap_settings$input<-"dist"
  umap_res<-umap::umap(max(snn_mat)-snn_mat,config=umap_settings)$layout

  if(is.null(resolution)){
    return(umap_res)
  }

  diag(snn_mat)<-0
  temp<-apply(snn_mat,2,graph.norm,knn=knn)
  temp<-temp+Matrix::t(temp)-temp*Matrix::t(temp)
  rownames(temp)<-as.character(1:nrow(Smat))
  colnames(temp)<-as.character(1:nrow(Smat))

  graph<-as.Graph(Matrix(temp,sparse=T))
  clust<-matrix("",nrow(Smat),length(resolution))
  for(i in 1:length(resolution)){
    temp<-FindClusters(graph,resolution=resolution[i])[,1]
    clust[,i]<-sprintf(paste0("%0",ceiling(log10(nlevels(temp))),"d"),temp)
  }
  output<-data.frame(umap_res,clust)
  colnames(output)<-c("umap_1","umap_2",paste0(resolution,resolution))
  return(output)
}
