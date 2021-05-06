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
  if(is.big.matrix(Smat)){
    temp<-Smat[1:nrow(Smat),1:ncol(Smat)]
    gc()
    eigs1<-PRIMME::eigs_sym(temp,NEig=neigs,which="LA")
    eigs2<-PRIMME::eigs_sym(temp,NEig=neigs,which="SA")
    rm(temp)
    gc()
    values<-c(eigs1$values,eigs2$values)
    vectors<-cbind(eigs1$vectors,eigs2$vectors)
    odr<-order(abs(values),decreasing=T)[1:neigs]
    eigs<-list(vectors=vectors[,odr],values=values[odr])
  }else{
    eigs1<-PRIMME::eigs_sym(Smat,NEig=neigs,which="LA")
    eigs2<-PRIMME::eigs_sym(Smat,NEig=neigs,which="SA")
    values<-c(eigs1$values,eigs2$values)
    vectors<-cbind(eigs1$vectors,eigs2$vectors)
    odr<-order(abs(values),decreasing=T)[1:neigs]
    eigs<-list(vectors=vectors[,odr],values=values[odr])
  }
  ratio<-sqrt(sum(apply(pcs,2,function(x) sum(x^2)))/sum(abs(eigs$values))*abs(eigs$values))
  return(cbind(pcs,t(t(eigs$vectors)*ratio)))
}
