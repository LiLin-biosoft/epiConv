#' Calculate frequency of peaks
#'
#' Calculate frequency of peaks
#'
#' @param mat matrix.
#'
#' @examples
#' \dontrun{
#' freq<-freq.estimate(mat)
#'}
#'
#' @return Return a vector with frequency of peaks.

freq.estimate<-function(mat){
  mat@x[mat@x>1]<-1
  return(Matrix::rowSums(mat)/ncol(mat))
}
