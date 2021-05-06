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
  freq<-sapply(1:nrow(mat)-1,function(x){
    sum(mat@i==x)/ncol(mat)
  })
  return(freq)
}
