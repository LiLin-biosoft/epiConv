epiConv.matrix<-function(mat,bin=50000,inf_replace=-10){
  ncell<-ncol(mat)
  output<-log10(symm.crossprod(mat,mat,bin=bin))
  output[is.infinite(output)]<-inf_replace
  return(output)
}
