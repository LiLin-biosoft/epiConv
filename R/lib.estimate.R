#' Calculate number of accessible peaks of cells
#'
#' Calculate number of accessible peaks of cells
#'
#' @param mat matrix.
#'
#' @examples
#' \dontrun{
#' lib_size<-lib.estimate(mat)
#'}
#'
#' @return Return a vector with number of accessible of cells

lib.estimate<-function(mat){
  lib<-mat@p[1:ncol(mat)+1]-mat@p[1:ncol(mat)]
  return(lib)
}
