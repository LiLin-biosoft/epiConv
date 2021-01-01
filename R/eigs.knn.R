#' learn neighbors across multiple datasets
#'
#' learn neighbors across multiple datasets
#'
#' @param Smat the similarity matrix.
#' @param features the guiding features.
#' @param batch the batch information.
#' @param reference which dataset in batch is used as reference. Reference should contain all cell types. If there are more than one references (A, B and C), we can specify reference=c("A","B","C") to align B to A, and align C to A and B.
#' @param knn_target number of neighbors for cells in target datasets (non-reference datasets) to learn in each reference.
#' @param knn_reference number of neighbors for cells in references to learn in each target dataset.
#' @param threshold the Z-score threshold used to filter false neighbors.
#'
#' @examples
#' \dontrun{
#'knn_update<-eigs.knn(Smat=Smat,
#'                     features=guide_features,
#'                     batch=batch,
#'                     reference="Resting",
#'                     knn_target=50,
#'                     knn_reference=10,
#'                     threshold=2)
#'}
#'
#' @return Return a list containing knn matrix. It is the input of eigs.correct.

eigs.knn<-function(Smat,features,batch,reference,knn_target=50,knn_reference=20,threshold=2){
  align_order<-c(reference,setdiff(levels(batch),reference))
  output<-list()
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    output[[tag1]]<-list()
    index1<-which(batch==tag1)
    temp<-cal.snn(Smat[index1,index1],
                  knn=max(floor(length(index1)*0.01),20))
    temp<-Matrix(data=temp,sparse=T,doDiag=F)
    output[[tag1]][[tag1]]<-temp

    if(!is.null(features[[tag1]])){
      for(tag2 in setdiff(reference,align_order[i:length(align_order)])){
        index2<-which(batch==tag2)
        temp<-refine.knn(target=find.knn(Smat[index2,index1],knn=knn_target),
                         reference=find.knn(Smat[index1,index2],knn=knn_reference),
                         features=features[[tag1]],
                         threshold=threshold)
        temp<-Matrix(data=temp,sparse=T,doDiag=F)
        output[[tag1]][[tag2]]<-temp

      }
    }else{
      for(tag2 in setdiff(reference,align_order[i:length(align_order)])){
        index2<-which(batch==tag2)
        temp<-find.knn(Smat[index2,index1],knn=knn_target)
        temp<-Matrix(data=temp,sparse=T,doDiag=F)
        output[[tag1]][[tag2]]<-temp
      }
    }
  }
  return(output)
}
