#' Find nearest neighbors for each column
#'
#' Find nearest neighbors for each column
#'
#' @param Smat the similarity matrix.
#' @param batch the batch information of cells.
#' @param reference the name of reference batch. Cells are sequentially aligned according the order of this vector.
#' @param neigs the number of Eigen vectors to calculate.
#' @param features features used to refine the knn matrix. It should be NULL if no feature is provided or a named list contains the features of cells from some batches.
#' @param knn_mat the pre-computed knn_mat. It should be NULL if no knn_mat is provided or a named list contains the knn_mat of some batches.
#' @param knn_target the number of nearest neighbors for cells from  query batch to find in reference batch.
#' @param knn_reference the number of nearest neighbors for cells from reference batch to find in query batch.
#' @param threshold the Z-score threshold to filter false neighbors.

#'
#' @examples
#' \dontrun{
#' res_joint<-epiConv.joint(Smat=Smat,
#'                          batch=batch,
#'                          reference=c("ref1","ref2"),
#'                          neigs=30,
#'                          features=NULL,
#'                          knn_mat=NULL,
#'                          knn_target=50,
#'                          knn_reference=10,
#'                          threshold=2)
#'
#'}
#'
#' @return Return a list contains Eigen vectors before correction, after correction and knn matrix across batches.


epiConv.joint<-function(Smat,batch,reference,neigs=30,features=list(),knn_mat=list(),
                        knn_target=50,knn_reference=10,threshold=2){
  knn_transfer_correction<-10
  anchor_threshold<-5
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
  cat("decomposition finished",fill=T)
  guide_features<-list()
  for(x in setdiff(levels(batch),reference[1])){
    if(is.null(features[[x]])){
      index<-which(batch==x)
      guide_features[[x]]<-dim.reduce(Smat[index,index],neigs=neigs)$vectors
    }else{
      index<-which(batch==x)
      guide_features[[x]]<-cbind(dim.reduce(Smat[index,index],neigs=neigs)$vectors,
                                 features[[x]])
    }

  }
  #  guide_features<-sapply(setdiff(levels(batch),reference[1]),function(x){
  #    index<-which(batch==x)
  #    eigs$vectors[index,1:neigs]
  #  })
  #  names(guide_features)<-setdiff(levels(batch),reference[1])
  #  if(!is.null(features)){
  #    for(i in names(features)){
  #      guide_features[[i]]<-cbind(guide_features[[i]],
  #      				features[[i]])
  #    }
  #  }
  knn_update<-eigs.knn(Smat=Smat,
                       features=guide_features,
                       batch=batch,
                       reference=reference,
                       knn_target=knn_target,
                       knn_reference=knn_reference,
                       threshold=threshold)
  for(i in names(knn_mat)){
    knn_update[[i]][[i]]<-knn_mat[[i]]
  }
  for(i in setdiff(levels(batch),c(names(knn_mat),reference[1]))){
    index<-which(batch==i)
    cat("calculating knn matrix:",i,fill=T)
    knn_update[[i]][[i]]<-find.knn(Smat[index,index],
                                   knn=floor(length(index)*0.01))
  }

  cat("knn calculation finished",fill=T)
  eigs_corrected<-eigs.scale(eigs=eigs$vectors,
                             knn_mat=knn_update,
                             batch=batch,
                             reference=reference,
                             knn_transfer_correction=knn_transfer_correction)
  cat("eigs correction finished",fill=T)

  if(is.big.matrix(Smat)){
    bin<-10000
    for(i in 1:ceiling(ncol(Smat)/bin)-1){
      row_index<-(i*bin+1):min((i*bin+bin),ncol(Smat))
      for(j in (i+1):ceiling(ncol(Smat)/bin)-1){
        col_index<-(j*bin+1):min((j*bin+bin),ncol(Smat))
        aa<-paste0(min(row_index),"-",max(row_index))
        bb<-paste0(min(col_index),"-",max(col_index))
        cat("Calculating corrected similarities between cells",aa,"vs",bb,fill=T)
        temp<-Smat[row_index,col_index]-
          tcrossprod(eigs$vectors[row_index,],
                     t(t(eigs$vectors[col_index,])*eigs$values))+
          tcrossprod(eigs_corrected[row_index,],
                     t(t(eigs_corrected[col_index,])*eigs$values))
        Smat[row_index,col_index]<-temp
        if(i!=j){
          Smat[col_index,row_index]<-t(temp)
        }
        rm(temp)
        gc()
      }
    }
  }else{
    Smat<-Smat-
      tcrossprod(eigs$vectors,
                 t(t(eigs$vectors)*eigs$values))+
      tcrossprod(eigs_corrected,
                 t(t(eigs_corrected)*eigs$values))
  }
  cat("similarity matrix correction finished",fill=T)
  return(list(eigs=eigs,knn_update=knn_update,eigs_corrected=eigs_corrected))
}
