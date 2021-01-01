#' Calculate shared nearest neighbors matrix from similarity matrix
#'
#' Calculate shared nearest neighbors matrix from similarity matrix
#'
#' @param Smat the similarity matrix.
#' @param knn number of neighbors. We generally set it to total cells*0.01.
#' @param bin Matrix are divided by bin to avoid very long vectors, It does not affect the results.
#'
#' @examples
#' \dontrun{
#' snn<-cal.snn(Smat=Smat,knn=50)
#'}
#'
#' @return Return a shared nearest neighbor matrix.

cal.snn<-function(Smat,knn=50,bin=5000){
  diag(Smat)<-apply(Smat,2,max)
  ncell<-nrow(Smat)
  neighbor_mat<-Matrix(apply(Smat,2,function(x,knn){
    output<-rep(0,length(x))
    output[order(x,decreasing=T)[1:knn]]<-1
    return(output)
  },knn=knn))
  output<-symm.crossprod(mat=neighbor_mat,bin=bin)
  return(output)
}
