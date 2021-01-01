#' Calculate meta features from RNA-seq and ATAC-seq
#'
#' Calculate meta features from RNA-seq and ATAC-seq
#'
#' @param Smat the similarity matrix.
#' @param pcs principal components from RNA-seq.
#' @param neigs number of Eigen vectors to calculate.
#'
#' @examples
#' \dontrun{
#'feature_coassay<-cal.feature(Smat=Smat,
#'                             pcs=prcomp(t(expr))$x,
#'                             neigs=50)
#'}
#'
#' @return Return meta features.

cal.feature<-function(Smat,pcs,neigs=50){
  eigs<-dim.reduce(Smat-median(Smat),neigs=neigs)
  ratio<-sqrt(sum(apply(pcs,2,function(x) sum(x^2)))/sum(abs(eigs$values))*abs(eigs$values))
  return(cbind(pcs,Matrix::t(Matrix::t(eigs$vectors)*ratio)))
}
