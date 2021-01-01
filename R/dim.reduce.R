#' Perform Eigenvalue decompostion.
#'
#' Perform Eigenvalue decompostion.
#'
#' @param Smat the similarity matrix.
#' @param neigs the number of Eigen vectors to calculate.
#'
#' @examples
#' \dontrun{
#' eigs<-dim.reduce(Smat,neigs=30)
#'}
#'
#' @return Return a list with Eigenvalues and Eigen vectors.

dim.reduce<-function(Smat,neigs=50){
  eigs1<-eigs_sym(Smat,NEig=neigs,which="LA")
  eigs2<-eigs_sym(Smat,NEig=neigs,which="SA")
  values<-c(eigs1$values,eigs2$values)
  vectors<-cbind(eigs1$vectors,eigs2$vectors)
  odr<-order(abs(values),decreasing=T)[1:neigs]
  return(list(vectors=vectors[,odr],values=values[odr]))
}
