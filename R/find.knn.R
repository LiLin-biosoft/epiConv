#' Find nearest neighbors for each column
#'
#' Find nearest neighbors for each column
#'
#' @param Smat the similarity matrix.
#' @param knn number of nearest neighbors.
#' @param binarize if binarize=T, the function will return 1 for a neighbor or 0 otherwise. if binarize=F, the function will return the similarity value for a neighbor or 0 otherwise

#'
#' @examples
#' \dontrun{
#' knn<-find.knn(Smat=Smat,
#'               knn=20,
#'               binarize=T)
#'}
#'
#' @return Return a knn matrix.


find.knn<-function(Smat,knn,binarize=T){
  output<-new("dgCMatrix")
  output@Dim<-as.integer(dim(Smat))
  output@p<-as.integer(0)
  nelement<-as.integer(0)
  for(k in 1:ncol(Smat)){
    cutoff<-sort(Smat[,k],decreasing=T)[knn]
    index<-sort(which(Smat[,k]>=cutoff))
    output@i<-c(output@i,as.integer(index-1))
    nelement<-nelement+length(index)
    output@p<-c(output@p,nelement)
    if(binarize){
      output@x<-c(output@x,rep(1,length(index)))
    }else{
      output@x<-c(output@x,Smat[index,k])
    }
  }
  return(output)
}
