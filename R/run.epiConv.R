#' Calculate similarities between cells
#'
#' Calculate similarities between cells
#'
#' @param mat Matrix object constains the peak by cell matrix.
#' @param lib_size library size of single cells.
#' @param nbootstrap number of bootstraps performed.
#' @param nsample number of peaks sampled in each bootstrap.
#' @param bin Matrix are divided by bin to avoid very long vectors, It does not affect the results.
#' @param inf_replace a small value used to replace -Inf.
#'
#' @examples
#' \dontrun{
#' Smat<-run.epiConv(mat=mat,
#'                   lib_size=colSums(mat),
#'                   nbootstrap=15,
#'                   nsample=floor(nrow(mat)*0.2),
#'                   bin=5000,
#'                   inf_replace=(-8))
#'}
#'
#' @return Return a similarity matrix.


run.epiConv<-function(mat,lib_size,nbootstrap=15,nsample=floor(nrow(mat)*0.2),bin=10000,inf_replace=(-8)){
  mat<-tfidf.norm(mat,lib_size=lib_size)
  Smat<-matrix(0,ncol(mat),ncol(mat))
  sample_matrix<-lapply(1:nbootstrap,function(x) sample(1:nrow(mat),size=nsample))
  for(i in sample_matrix){
    Smat<-Smat+epiConv.matrix(mat=mat[i,],inf_replace=inf_replace,bin=bin)
  }
  Smat<-Smat/nbootstrap

  cell_sample<-500
  temp<-sample(1:nrow(Smat),cell_sample*2)
  retain1<-sort(temp[1:cell_sample])
  retain2<-sort(temp[(cell_sample+1):(cell_sample*2)])
  stablization_pars<-stablization.pars(Smat=Smat[retain1,retain2],
                                       lib_size=lib_size[c(retain1,retain2)])
  Smat<-similarity.stablization(Smat,
                                lib_size=lib_size,
                                stablization_pars=stablization_pars,
                                bin=bin)
  return(Smat)
}
