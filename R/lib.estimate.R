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
  mat@x[mat@x>1]<-1
  return(Matrix::colSums(mat))}
