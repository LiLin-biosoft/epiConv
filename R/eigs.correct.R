#' correct Eigen vectors.
#'
#' correct Eigen vectors.
#'
#' @param eigs the Eigen vectors.
#' @param knn_mat the knn matrix, output from eigs.knn.
#' @param batch the batch information.
#' @param reference datasets used as reference.
#' @param knn_transfer_correction the number of anchors used by non-anchor cells to learn their corrections.
#' @param anchor_threshold the minimum number of neighbors of an anchor cell.
#'
#' @examples
#' \dontrun{
#'eigs_corrected<-eigs.correct(eigs=eigs$vectors,
#'                             knn_mat=knn_update,
#'                             batch=batch,
#'                             reference="Resting",
#'                             knn_transfer_correction=10)
#'}
#'
#' @return Return corrected Eigen vectors.


eigs.correct<-function(eigs,knn_mat,batch,reference,knn_transfer_correction=10,anchor_threshold=5){
  align_order<-c(reference,setdiff(levels(batch),reference))
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    index_self<-which(batch==tag1)
    index_other<-which(batch%in%intersect(reference,align_order[1:(i-1)]))
    index_other<-sapply(intersect(reference,align_order[1:(i-1)]),function(x){
      which(batch==x)
    })
    index_other<-unlist(index_other)

    anchor_mat<-NULL
    for(tag2 in intersect(reference,align_order[1:(i-1)])){
      anchor_mat<-rbind(anchor_mat,knn_mat[[tag1]][[tag2]])
    }

    correction<-sapply(1:length(index_self),function(x){
      apply(eigs,2,function(y){
        wt<-anchor_mat[,x]
        if(sum(wt!=0)<anchor_threshold){
          return(NA)
        }else{
          weighted.mean(x=y[index_other],w=wt,na.rm=T)-y[index_self[x]]
        }
      })
    })
    correction<-Matrix::t(correction)
    snn_mat<-knn_mat[[tag1]][[tag1]]
    correction<-apply(correction[!is.na(correction[,1]),],2,function(x,neighbor_mat){
      apply(neighbor_mat,2,function(y){
        cutoff<-sort(y,decreasing=T)[knn_transfer_correction]
        y[y<cutoff]<-0
        return(weighted.mean(x=x,w=y,na.rm=T))
      })
    },neighbor_mat=snn_mat[!is.na(correction[,1]),])

    eigs[index_self,]<-eigs[index_self,]+correction
  }
  return(eigs)
}
