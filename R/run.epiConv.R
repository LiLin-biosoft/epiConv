#' Calculate similarities between cells
#'
#' Calculate similarities between cells
#'
#' @param mat Matrix object constains the peak by cell matrix.
#' @param lib_size library size of single cells.
#' @param nbootstrap number of bootstraps performed.
#' @param nsample number of peaks sampled in each bootstrap.
#' @param bin Matrix are divided by bin to avoid very long vectors, It does not affect the results.
#' @param inf a small value used to replace -Inf.
#' @param ... Arguments passed to big.matrix when creating the similarity matrix. You can specify a backup file through "backingfile" and "descriptorfile". If not specified, the function will create a in-memory backup file.
#'
#' @examples
#' \dontrun{
#' Smat<-run.epiConv(mat=mat,
#'                   lib_size=colSums(mat),
#'                   nbootstrap=15,
#'                   nsample=floor(nrow(mat)*0.2),
#'                   bin=5000,
#'                   inf=(-8),
#'                   backingfile="backupfile",
#'                   descriptorfile="backupfile.descriptor")
#'}
#'
#' @return Return a similarity matrix.


run.epiConv<-function(mat,lib_size,nbootstrap,nsample,bin=5000,inf=(-8),...){

  sample_matrix<-lapply(1:nbootstrap,function(x) sample(1:nrow(mat),size=nsample))
  mat<-tfidf.norm(mat,lib_size=lib_size)

  ####calculate pars#########
  cell_sample<-500
  temp<-sample(1:ncol(mat),cell_sample*2)
  retain1<-sort(temp[1:cell_sample])
  retain2<-sort(temp[(cell_sample+1):(cell_sample*2)])
  Smat_small<-epiConv.matrix(mat1=mat[,retain1],
                             mat2=mat[,retain2],
                             sample_matrix=sample_matrix,
                             inf=inf)
  adjust_pars<-stablization.pars(Smat=Smat_small,
                                 lib_size=lib_size[c(retain1,retain2)])
  rm(Smat_small)
  gc()
  ####calculate pars#########

  Smat<-big.matrix(nrow=ncol(mat),ncol=ncol(mat),init=0,...)
  for(i in 1:ceiling(ncol(mat)/bin)-1){
    row_index<-(i*bin+1):min((i*bin+bin),ncol(mat))
    for(j in (i+1):ceiling(ncol(mat)/bin)-1){
      col_index<-(j*bin+1):min((j*bin+bin),ncol(mat))
      aa<-paste0(min(row_index),"-",max(row_index))
      bb<-paste0(min(col_index),"-",max(col_index))
      cat("Calculating similarities between cells",aa,"vs",bb,fill=T)
      temp<-epiConv.matrix(mat1=mat[,row_index],
                           mat2=mat[,col_index],
                           sample_matrix=sample_matrix,
                           inf=inf,
                           lib_size1=lib_size[row_index],
                           lib_size2=lib_size[col_index],
                           adjust_pars=adjust_pars)
      Smat[row_index,col_index]<-temp
      if(i!=j){
        Smat[col_index,row_index]<-t(temp)
      }
      rm(temp)
      gc()
    }
  }
  return(Smat)
}
