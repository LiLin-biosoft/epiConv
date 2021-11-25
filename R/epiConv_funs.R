#' Calculate number of accessible peaks of cells
#'
#' Calculate number of accessible peaks of cells
#'
#' @param mat a dgCMatrix object.
#'
#' @examples
#' \dontrun{
#' lib_size<-lib.estimate(mat)
#'}
#'
#' @return Return a vector of library sizes

lib.estimate<-function(mat){
  lib<-mat@p[1:ncol(mat)+1]-mat@p[1:ncol(mat)]
  return(lib)
}


num2int<-function(x,nbin=32,bound=quantile(x,c(0.01,0.99))){
  win<-(bound[2]-bound[1])/(nbin-2)*(0:(nbin-2))+bound[1]
  return(continuous.match(x,win))
}

stablization.pars<-function(Smat,lib_size){
  lib_mat<-matrix(0,nrow(Smat),ncol(Smat))
  lib_mat<-lib_mat+log10(lib_size[1:nrow(Smat)])
  lib_mat<-Matrix::t(Matrix::t(lib_mat)+log10(lib_size[(nrow(Smat)+1):(nrow(Smat)*2)]))

  lm_res<-lm(as.vector(Smat)~as.vector(lib_mat),weight=10^as.vector(lib_mat))
  Smat_residual<-matrix(lm_res$residuals,nrow(Smat),nrow(Smat))
  lib_discrete<-num2int(lib_mat,nbin=100,bound=quantile(lib_mat,c(0.01,0.99)))
  lib_window<-quantile(lib_mat,0.01)+(quantile(lib_mat,0.99)-quantile(lib_mat,0.01))/98*(0:98)
  sd_value<-tapply(Smat_residual,list(factor(lib_discrete)),function(x){
    sqrt(sum(x^2)/length(x))
  })
  return(list(coef=lm_res$coefficients,lib_window=lib_window,sd_value=sd_value))
}

tfidf.norm<-function(mat,lib_size){
  mat@x[mat@x>0]<-1
  freq<-Matrix::rowSums(mat)/ncol(mat)
  wt<-log10(1+1/freq)
  wt[freq==0]<-0
  temp<-sapply(1:length(lib_size),function(i){
    return(mat@x[(mat@p[i]+1):mat@p[i+1]]/lib_size[i])
  })
  mat@x<-unlist(temp)
  return(mat*wt)
}

continuous.match<-function(x,win){
  nbin<-length(win)+1
  bin<-win[2]-win[1]
  x<-floor((x-win[1])/bin)+2
  x[x>nbin]<-nbin
  x[x<1]<-1
  return(x)
}
lib2sd<-function(lib,adjust_pars){
  win<-adjust_pars$lib_window
  nbin<-length(win)+1
  bin<-win[2]-win[1]
  lib<-floor((lib-win[1])/bin)+2
  lib[lib>nbin]<-nbin
  lib[lib<1]<-1
  lib<-adjust_pars$sd_value[lib]
  return(lib)
}

epiConv.matrix<-function(mat1,mat2,sample_matrix,inf,lib_size1,lib_size2,adjust_pars=NULL){
  Smat<-matrix(0,ncol(mat1),ncol(mat2))
  nbootstrap<-length(sample_matrix)
  for(k in sample_matrix){
    temp<-log10(Matrix::crossprod(mat1[k,],mat2[k,]))
    temp@x[is.infinite(temp@x)]<-inf
    Smat<-Smat+as.matrix(temp)/nbootstrap
    rm(temp)
  }
  if(is.null(adjust_pars)){
    return(Smat)
  }else{
    ncell1<-length(lib_size1)
    ncell2<-length(lib_size2)
    Smat<-Smat-log10(lib_size1)*adjust_pars$coef[2]-adjust_pars$coef[1]
    Smat<-sweep(x=Smat,MARGIN=2,STATS=log10(lib_size2)*adjust_pars$coef[2],FUN="-")
    sd_mat<-sweep(matrix(log10(lib_size1),ncell1,ncell2),2,log10(lib_size2),FUN="+")

    sd_mat<-apply(sd_mat,2,lib2sd,adjust_pars=adjust_pars)
    Smat<-Smat/sd_mat
    rm(sd_mat)
    return(Smat)
  }
}

#' Calculate similarities between cells
#'
#' Calculate similarities between cells
#'
#' @param mat Matrix object constains the peak by cell matrix.
#' @param row_sample indices of sampled cells.
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
#'                   row_sample=1:ncol(mat),
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


run.epiConv<-function(mat,row_sample=NULL,lib_size,nbootstrap,nsample,inf=(-8),bin=1000,...){

  if(is.null(row_sample)){
    Smat_index<-1:ncol(mat)
  }else{
    Smat_index<-row_sample
  }
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

  Smat<-big.matrix(nrow=length(Smat_index),ncol=ncol(mat),init=0,type="double",...)
  descriptor<-describe(Smat)

  row_index_list<-lapply(1:ceiling(length(Smat_index)/bin)-1,function(j){
    return((j*bin+1):min((j*bin+bin),length(Smat_index)))
  })
  col_index_list<-lapply(row_index_list,function(j){
    Smat_index[j]
  })
  non_row<-setdiff(1:ncol(mat),Smat_index)
  if(length(non_row)!=0){
    temp<-lapply(1:ceiling(length(non_row)/bin)-1,function(j){
      return(non_row[(j*bin+1):min((j*bin+bin),length(non_row))])
    })
    col_index_list<-c(col_index_list,temp)
  }
  task_index<-lapply(1:length(col_index_list),function(x) rbind(1:min(x,length(row_index_list)),x))
  task_index<-matrix(unlist(task_index),2,)

  nbin_row<-length(row_index_list)

  mat_split<-list()
  for(i in 1:length(col_index_list)){
    mat_split[[i]]<-mat[,col_index_list[[i]]]
  }
  rm(mat)
  gc()

  cat("Total jobs: ",ncol(task_index),sep="",fill=T)

  for(k in 1:ncol(task_index)){
    x<-task_index[,k]
    row_index<-row_index_list[[x[1]]]
    col_index<-col_index_list[[x[2]]]

    cat("Begin job ",k,sep="",fill=T)
    temp<-epiConv.matrix(mat1=mat_split[[x[1]]],
                         mat2=mat_split[[x[2]]],
                         sample_matrix=sample_matrix,
                         inf=inf,
                         lib_size1=lib_size[Smat_index][row_index],
                         lib_size2=lib_size[col_index],
                         adjust_pars=adjust_pars)
    Smat[row_index,col_index]<-temp
    if((x[2]<=nbin_row)&(x[1]!=x[2])){
      Smat[match(col_index,Smat_index),Smat_index[row_index]]<-t(temp)
    }
    rm(temp)
  }

  return(Smat)
}

#' Calculate similarities between cells in parallel
#'
#' Calculate similarities between cells in parallel
#'
#' @param mat Matrix object constains the peak by cell matrix.
#' @param row_sample indices of sampled cells.
#' @param lib_size library size of single cells.
#' @param nbootstrap number of bootstraps performed.
#' @param nsample number of peaks sampled in each bootstrap.
#' @param bin Matrix are divided by bin to avoid very long vectors, It does not affect the results.
#' @param inf a small value used to replace -Inf.
#' @param ncore number of threads.
#' @param ... Arguments passed to big.matrix when creating the similarity matrix. You can specify a backup file through "backingfile" and "descriptorfile". If not specified, the function will create a in-memory backup file.
#'
#' @examples
#' \dontrun{
#' Smat<-run.epiConv.parallel(mat=mat,
#'                            row_sample=1:ncol(mat),
#'                            lib_size=colSums(mat),
#'                            nbootstrap=15,
#'                            nsample=floor(nrow(mat)*0.2),
#'                            bin=5000,
#'                            inf=(-8),
#'                            ncore=5,
#'                            backingfile="backupfile",
#'                            descriptorfile="backupfile.descriptor")
#'}
#'
#' @return Return a similarity matrix.


run.epiConv.parallel<-function(mat,row_sample=NULL,lib_size,nbootstrap,nsample,inf=(-8),bin=1000,ncore=5,...){

  if(is.null(row_sample)){
    Smat_index<-1:ncol(mat)
  }else{
    Smat_index<-row_sample
  }
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

  Smat<-big.matrix(nrow=length(Smat_index),ncol=ncol(mat),init=0,type="double",...)
  descriptor<-describe(Smat)

  row_index_list<-lapply(1:ceiling(length(Smat_index)/bin)-1,function(j){
    return((j*bin+1):min((j*bin+bin),length(Smat_index)))
  })
  col_index_list<-lapply(row_index_list,function(j){
    Smat_index[j]
  })
  non_row<-setdiff(1:ncol(mat),Smat_index)
  if(length(non_row)!=0){
    temp<-lapply(1:ceiling(length(non_row)/bin)-1,function(j){
      return(non_row[(j*bin+1):min((j*bin+bin),length(non_row))])
    })
    col_index_list<-c(col_index_list,temp)
  }
  task_index<-lapply(1:length(col_index_list),function(x) rbind(1:min(x,length(row_index_list)),x))
  task_index<-matrix(unlist(task_index),2,)

  nbin_row<-length(row_index_list)

  mat_split<-list()
  for(i in 1:length(col_index_list)){
    mat_split[[i]]<-mat[,col_index_list[[i]]]
  }
  rm(mat)
  gc()

  cat(ncol(task_index)," jobs running on ",ncore," threads",sep="",fill=T)

  cl<-makeCluster(ncore,type="FORK")
  parApply(cl=cl,task_index,2,function(x){
    row_index<-row_index_list[[x[1]]]
    col_index<-col_index_list[[x[2]]]
    temp<-epiConv.matrix(mat1=mat_split[[x[1]]],
                         mat2=mat_split[[x[2]]],
                         sample_matrix=sample_matrix,
                         inf=inf,
                         lib_size1=lib_size[Smat_index][row_index],
                         lib_size2=lib_size[col_index],
                         adjust_pars=adjust_pars)
    Smat<-attach.big.matrix(descriptor)
    Smat[row_index,col_index]<-temp
    if((x[2]<=nbin_row)&(x[1]!=x[2])){
      Smat[match(col_index,Smat_index),Smat_index[row_index]]<-t(temp)
    }
    rm(Smat)
    gc()
    return(1)
  })

  stopCluster(cl)
  return(Smat)
}

index2MM<-function(index){
  index<-apply(index,2,sort)
  ncell<-ncol(index)
  knn<-nrow(index)
  output<-new("dgCMatrix")
  output@Dim<-as.integer(c(ncell,ncell))
  output@x<-rep(1,knn*ncell)
  output@p<-as.integer((0:ncell)*knn)
  output@i<-as.integer(index-1)
  return(output)
}

indexvalue2MM<-function(index,dim1,dim2){
  knn<-nrow(index)/2
  output<-new("dgCMatrix")
  output@Dim<-c(as.integer(dim1),as.integer(dim2))
  output@p<-as.integer((0:dim2)*knn)
  output@i<-as.integer(index[1:knn,]-1)
  output@x<-as.vector(index[1:knn+knn,])
  return(output)
}

find.knn<-function(Smat,knn,noise=F){
  if(is.big.matrix(Smat)){
    index<-biganalytics::apply(Smat,2,function(x){
      temp<-sort(order(x,decreasing=T)[1:knn])
      return(c(temp-1,1+x[temp]/1000))
    })
  }else{
    index<-apply(Smat,2,function(x){
      temp<-sort(order(x,decreasing=T)[1:knn])
      return(c(temp-1,1+x[temp]/1000))
    })
  }
  output<-new("dgCMatrix")
  output@Dim<-as.integer(dim(Smat))
  output@p<-as.integer((0:ncol(Smat))*knn)
  output@i<-as.integer(index[1:knn,])
  if(noise){
    output@x<-as.vector(index[1:knn+knn,])
  }else{
    output@x<-rep(1,ncol(Smat)*knn)
  }
  return(output)
}

trunc.crossprod<-function(mat1,mat2,knn,ncore=1){

  if(ncore==1){
    output_index<-sapply(X=1:ncol(mat2),FUN=function(k){
      temp<-mat2@i[(mat2@p[k]+1):mat2@p[k+1]]+1
      temp<-Matrix::colSums(mat1[temp,])
      nonzero<-which(temp!=0)
      if(length(nonzero)>knn){
        aa<-sort(nonzero[order(temp[nonzero],decreasing=T)[1:knn]])
      }else{
        aa<-sort(c(nonzero,setdiff(1:length(temp),nonzero))[1:knn])
      }

      return(c(aa,temp[aa]))
    })
  }else{
    cl<-makeCluster(ncore,type="FORK")
    output_index<-parSapply(cl=cl,X=1:ncol(mat2),FUN=function(k){
      temp<-mat2@i[(mat2@p[k]+1):mat2@p[k+1]]+1
      temp<-Matrix::colSums(mat1[temp,])
      nonzero<-which(temp!=0)
      if(length(nonzero)>knn){
        aa<-sort(nonzero[order(temp[nonzero],decreasing=T)[1:knn]])
      }else{
        aa<-sort(c(nonzero,setdiff(1:length(temp),nonzero))[1:knn])
      }
      return(c(aa,temp[aa]))
    })
    stopCluster(cl)
  }

  output<-new("dgCMatrix")
  output@Dim<-c(as.integer(ncol(mat1)),as.integer(ncol(mat2)))
  output@p<-as.integer((0:ncol(mat2))*knn)
  output@i<-as.integer(output_index[1:knn,]-1)
  output@x<-as.vector(output_index[(knn+1):(2*knn),])
  return(output)
}

dim.reduce<-function(Smat,Smat_index,neigs=30){
  Smat_index<-sort(Smat_index)
  eigs1<-PRIMME::eigs_sym(Smat[,Smat_index],NEig=neigs,which="LA")
  eigs2<-PRIMME::eigs_sym(Smat[,Smat_index],NEig=neigs,which="SA")
  values<-c(eigs1$values,eigs2$values)
  vectors<-cbind(eigs1$vectors,eigs2$vectors)
  odr<-order(abs(values),decreasing=T)[1:neigs]
  eigs<-sweep(vectors[,odr],MARGIN=2,STAT=values[odr],FUN="*")
  non_index<-setdiff(1:ncol(Smat),Smat_index)
  if(length(non_index)==0){
    return(list(vectors=vectors[,odr],values=values[odr]))
  }

  left_mat<-crossprod(eigs)
  right_mat<-crossprod(eigs,Smat[,non_index])
  temp<-eigen(left_mat)
  root<-temp$vectors%*%diag(1/temp$values)%*%t(temp$vectors)%*%right_mat
  return(list(vectors=rbind(vectors[,odr],t(root))[match(1:ncol(Smat),c(Smat_index,non_index)),],
              values=values[odr]))

}

#' Find nearest neighbors between batches
#'
#' Find nearest neighbors between batches
#'
#' @param Smat the similarity matrix.
#' @param row_sample indices of sampled cells for rows.
#' @param batch the batch information of cells.
#' @param reference the name of reference batch. Cells are sequentially aligned according the order of this vector.
#' @param neigs the number of Eigen vectors to calculate.
#' @param features features used to refine the knn matrix. It should be NULL if no feature is provided or a named list, where each element is a cell by feature matrix and the name is equal to the corresponding batch.
#' @param knn_target the number of nearest neighbors for cells from  query batch to find in reference batch.
#' @param knn_reference the number of nearest neighbors for cells from reference batch to find in query batch.
#' @param threshold the Z-score threshold to filter false neighbors.

#'
#' @examples
#' \dontrun{
#' expr<-matrix(runif(10000),100,100)
#' PCs<-prcomp(expr)$x[,1:20]
#' res_anchor<-epiConv.anchor(Smat=Smat,
#'                            row_sample=1:ncol(Smat),
#'                            batch=batch,
#'                            reference=c("Batch1","Batch2"),
#'                            neigs=30,
#'                            features=list(Batch1=PCs),
#'                            knn_target=50,
#'                            knn_reference=20,
#'                            threshold=2)
#'
#'}
#'
#' @return Return a list contains knn matrix between batches and associated data.

epiConv.anchor<-function(Smat,row_sample,batch,reference,neigs=30,features=list(),
                         knn_target=50,knn_reference=20,threshold=2){

  Smat_index<-row_sample
  guide_features<-list()
  knn_mat<-list()
  for(x in setdiff(levels(batch),reference[1])){
    index<-which(batch==x)
    aa<-which(Smat_index%in%index)
    bb<-match(Smat_index[aa],index)
    if(is.null(features[[x]])){
      guide_features[[x]]<-dim.reduce(Smat[aa,index],Smat_index=bb,neigs=neigs)$vectors
    }else{
      guide_features[[x]]<-cbind(dim.reduce(Smat[aa,index],Smat_index=bb,neigs=neigs)$vectors,
                                 features[[x]])
      knn_mat[[x]]<-index2MM(t(FNN::knn.index(features[[x]],k=max(floor(nrow(features[[x]])*0.01),20))))
    }
  }

  rm(features)
  cat("decomposition finished",fill=T)
  gc()

  ref_knn<-list()
  for(tag in reference){
    cat("calculate snn of ",tag,sep="",fill=T)
    index<-which(batch==tag)
    index<-index[index%in%Smat_index]
    row_index<-match(index,Smat_index)
    temp_mat<-find.knn(Smat[row_index,index],knn=max(floor(length(row_index)/100),2),noise=T)
    if(tag!=reference[1]&(is.null(knn_mat[[tag]]))){
      knn_mat[[tag]]<-temp_mat
    }

    output_index<-sapply(X=1:ncol(temp_mat),FUN=function(k){
      temp<-temp_mat@i[(temp_mat@p[k]+1):temp_mat@p[k+1]]+1
      temp<-Matrix::colSums(temp_mat[temp,])
      nonzero<-which(temp!=0)
      if(length(nonzero)>knn_reference){
        aa<-sort(nonzero[order(temp[nonzero],decreasing=T)[1:knn_reference]])
      }else{
        aa<-sort(c(nonzero,setdiff(1:length(temp),nonzero))[1:knn_reference])
      }

      return(c(aa,temp[aa]))
    })
    ref_knn[[tag]]<-indexvalue2MM(output_index,dim1=ncol(temp_mat),dim2=ncol(temp_mat))
  }

  align_order<-c(reference,setdiff(levels(batch),reference))
  knn_update<-list()
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    knn_update[[tag1]]<-list()
    index1<-which(batch==tag1)

    for(tag2 in setdiff(reference,align_order[i:length(align_order)])){
      index2<-which(batch==tag2)
      index2<-index2[index2%in%Smat_index]
      row_index2<-match(index2,Smat_index)

      cat("finding anchors:",tag1,"to",tag2,fill=T)

      #################
      ntarget<-length(index1)
      nreference<-length(row_index2)
      nfeatures<-ncol(guide_features[[tag1]])
      unsmoothed<-find.knn(t(Smat[row_index2,index1]),knn=knn_reference,noise=T)
      targetMM<-Matrix::t(find.knn(Smat[row_index2,index1],knn=knn_target))
      gc()
      smoothed_index<-sapply(1:ncol(ref_knn[[tag2]]),function(j){
        temp<-Matrix::rowSums(unsmoothed[,ref_knn[[tag2]]@i[(ref_knn[[tag2]]@p[j]+1):ref_knn[[tag2]]@p[j+1]]+1])
        order(temp,decreasing=T)[1:knn_reference]
      })
      whether_retain<-sapply(which(targetMM@p[1:nreference+1]!=targetMM@p[1:nreference]),function(k){
        temp<-guide_features[[tag1]][smoothed_index[,k],]
        mu<-colMeans(temp)
        sdev<-sqrt(colVars(temp))
        linked_target<-targetMM@i[(targetMM@p[k]+1):targetMM@p[k+1]]+1
        output<-rep(0,length(linked_target))
        if(length(linked_target)>1){
          output[colSums(abs((t(guide_features[[tag1]][linked_target,])-mu)/sdev)<threshold)==nfeatures]<-1
        }else{
          output<-as.numeric(sum(abs((guide_features[[tag1]][linked_target,]-mu)/sdev)<2)==nfeatures)
        }
        return(output)
      })
      targetMM@x<-unlist(whether_retain)
      knn_update[[tag1]][[tag2]]<-Matrix::t(drop0(targetMM))
      rm(targetMM)
      gc()

      #################
    }
  }

  for(i in names(knn_mat)){
    knn_update[[i]][[i]]<-knn_mat[[i]]
  }

  for(i in setdiff(levels(batch),c(names(knn_mat),reference[1]))){
    index<-which(batch==i)
    row_index<-which(Smat_index%in%index)
    cat("calculating knn matrix:",i,fill=T)
    knn_update[[i]][[i]]<-find.knn(Smat[row_index,index],
                                   knn=max(floor(length(row_index)*0.01),20),noise=T)
  }
  rm(knn_mat)
  gc()
  cat("anchor finding finished",fill=T)

  return(list(knn_update=knn_update,
              ref_knn=ref_knn,
              guide_features=guide_features))
}


#' Find nearest neighbors between batches in parallel
#'
#' Find nearest neighbors between batches in parallel
#'
#' @param Smat the similarity matrix.
#' @param row_sample indices of sampled cells for rows.
#' @param batch the batch information of cells.
#' @param reference the name of reference batch. Cells are sequentially aligned according the order of this vector.
#' @param neigs the number of Eigen vectors to calculate.
#' @param features features used to refine the knn matrix. It should be NULL if no feature is provided or a named list, where each element is a cell by feature matrix and the name is equal to the corresponding batch.
#' @param knn_target the number of nearest neighbors for cells from  query batch to find in reference batch.
#' @param knn_reference the number of nearest neighbors for cells from reference batch to find in query batch.
#' @param threshold the Z-score threshold to filter false neighbors.
#' @param ncore number of threads.
#'
#' @examples
#' \dontrun{
#' expr<-matrix(runif(10000),100,100)
#' PCs<-prcomp(expr)$x[,1:20]
#' res_anchor<-epiConv.anchor.parallel(Smat=Smat,
#'                                     row_sample=1:ncol(Smat),
#'                                     batch=batch,
#'                                     reference=c("Batch1","Batch2"),
#'                                     neigs=30,
#'                                     features=list(Batch1=PCs),
#'                                     knn_target=50,
#'                                     knn_reference=20,
#'                                     threshold=2,
#'                                     ncore=5)
#'
#'}
#'
#' @return Return a list contains knn matrix between batches and associated data.


epiConv.anchor.parallel<-function(Smat,row_sample,batch,reference,neigs=30,features=list(),
                                  knn_target=50,knn_reference=20,threshold=2,ncore=1){

  Smat_index<-row_sample
  guide_features<-list()
  knn_mat<-list()
  for(x in setdiff(levels(batch),reference[1])){
    index<-which(batch==x)
    aa<-which(Smat_index%in%index)
    bb<-match(Smat_index[aa],index)
    if(is.null(features[[x]])){
      guide_features[[x]]<-dim.reduce(Smat[aa,index],Smat_index=bb,neigs=neigs)$vectors
    }else{
      guide_features[[x]]<-cbind(dim.reduce(Smat[aa,index],Smat_index=bb,neigs=neigs)$vectors,
                                 features[[x]])
      knn_mat[[x]]<-index2MM(t(FNN::knn.index(features[[x]],k=max(floor(nrow(features[[x]])*0.01),20))))
    }
  }

  rm(features)
  cat("decomposition finished",fill=T)
  gc()

  cl<-makeCluster(ncore,type="FORK")
  ref_knn<-list()
  for(tag in reference){
    cat("calculate snn of ",tag,sep="",fill=T)
    index<-which(batch==tag)
    index<-index[index%in%Smat_index]
    row_index<-match(index,Smat_index)
    temp_mat<-find.knn(Smat[row_index,index],knn=max(floor(length(row_index)/100),2),noise=T)
    if(tag!=reference[1]&(is.null(knn_mat[[tag]]))){
      knn_mat[[tag]]<-temp_mat
    }

    output_index<-parSapply(cl=cl,X=1:ncol(temp_mat),FUN=function(k){
      temp<-temp_mat@i[(temp_mat@p[k]+1):temp_mat@p[k+1]]+1
      temp<-Matrix::colSums(temp_mat[temp,])
      nonzero<-which(temp!=0)
      if(length(nonzero)>knn_reference){
        aa<-sort(nonzero[order(temp[nonzero],decreasing=T)[1:knn_reference]])
      }else{
        aa<-sort(c(nonzero,setdiff(1:length(temp),nonzero))[1:knn_reference])
      }

      return(c(aa,temp[aa]))
    })
    ref_knn[[tag]]<-indexvalue2MM(output_index,dim1=ncol(temp_mat),dim2=ncol(temp_mat))
  }

  align_order<-c(reference,setdiff(levels(batch),reference))
  knn_update<-list()
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    knn_update[[tag1]]<-list()
    index1<-which(batch==tag1)

    for(tag2 in setdiff(reference,align_order[i:length(align_order)])){
      index2<-which(batch==tag2)
      index2<-index2[index2%in%Smat_index]
      row_index2<-match(index2,Smat_index)

      cat("finding anchors:",tag1,"to",tag2,fill=T)

      #################
      ntarget<-length(index1)
      nreference<-length(row_index2)
      nfeatures<-ncol(guide_features[[tag1]])
      unsmoothed<-find.knn(t(Smat[row_index2,index1]),knn=knn_reference,noise=T)
      targetMM<-Matrix::t(find.knn(Smat[row_index2,index1],knn=knn_target))
      gc()
      smoothed_index<-parSapply(cl=cl,1:ncol(ref_knn[[tag2]]),function(j){
        temp<-Matrix::rowSums(unsmoothed[,ref_knn[[tag2]]@i[(ref_knn[[tag2]]@p[j]+1):ref_knn[[tag2]]@p[j+1]]+1])
        order(temp,decreasing=T)[1:knn_reference]
      })
      whether_retain<-parSapply(cl=cl,which(targetMM@p[1:nreference+1]!=targetMM@p[1:nreference]),function(k){
        temp<-guide_features[[tag1]][smoothed_index[,k],]
        mu<-colMeans(temp)
        sdev<-sqrt(colVars(temp))
        linked_target<-targetMM@i[(targetMM@p[k]+1):targetMM@p[k+1]]+1
        output<-rep(0,length(linked_target))
        if(length(linked_target)>1){
          output[colSums(abs((t(guide_features[[tag1]][linked_target,])-mu)/sdev)<threshold)==nfeatures]<-1
        }else{
          output<-as.numeric(sum(abs((guide_features[[tag1]][linked_target,]-mu)/sdev)<2)==nfeatures)
        }
        return(output)
      })
      targetMM@x<-unlist(whether_retain)
      knn_update[[tag1]][[tag2]]<-Matrix::t(drop0(targetMM))
      rm(targetMM)
      gc()

      #################
    }
  }
  stopCluster(cl)

  for(i in names(knn_mat)){
    knn_update[[i]][[i]]<-knn_mat[[i]]
  }

  for(i in setdiff(levels(batch),c(names(knn_mat),reference[1]))){
    index<-which(batch==i)
    row_index<-which(Smat_index%in%index)
    cat("calculating knn matrix:",i,fill=T)
    knn_update[[i]][[i]]<-find.knn(Smat[row_index,index],
                                   knn=max(floor(length(row_index)*0.01),20),noise=T)
  }
  rm(knn_mat)
  gc()
  cat("anchor finding finished",fill=T)

  return(list(knn_update=knn_update,
              ref_knn=ref_knn,
              guide_features=guide_features))
}

#' Correct batch effects on similarity matrix
#'
#' Correct batch effects on similarity matrix
#'
#' @param Smat the similarity matrix.
#' @param row_sample indices of sampled cells for rows.
#' @param batch the batch information of cells.
#' @param reference the name of reference batch. Cells are sequentially aligned according the order of this vector.
#' @param neigs the number of Eigen vectors to calculate.
#' @param knn_update the knn matrix between batches.
#'
#' @examples
#' \dontrun{
#' res_correct<-epiConv.correct(Smat=Smat,
#'                              row_sample=1:ncol(Smat),
#'                              batch=batch,
#'                              reference=c("Batch1","Batch2"),
#'                              neigs=30,
#'                              knn_update=res_anchor$knn_update)
#'
#'}
#'
#' @return Return a list contains Eigen vectors before and after correction.

epiConv.correct<-function(Smat,row_sample,batch,reference,neigs=30,knn_update){

  Smat_index<-row_sample
  eigs<-dim.reduce(Smat=Smat,Smat_index=Smat_index,neigs=neigs)
  eigs$corrected<-eigs$vectors
  cat("decomposition finished",fill=T)
  gc()

  anchor_threshold<-5
  knn_transfer_correction<-10
  align_order<-c(reference,setdiff(levels(batch),reference))
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    index_self<-which(batch==tag1)
    nself<-length(index_self)
    index_other<-sapply(intersect(reference,align_order[1:(i-1)]),function(x){
      which(batch==x)
    })
    index_other<-unlist(index_other)
    index_other<-index_other[index_other%in%Smat_index]

    anchor_mat<-NULL
    for(tag2 in intersect(reference,align_order[1:(i-1)])){
      anchor_mat<-rbind(anchor_mat,knn_update[[tag1]][[tag2]])
    }
    valid_anchor<-which(Matrix::colSums(anchor_mat)>=anchor_threshold)
    anchor_mat<-anchor_mat[,valid_anchor]
    cat(tag1,": ",length(valid_anchor)," anchors in ",length(index_self)," cells",sep="",fill=T)

    output_index<-sapply(X=1:ncol(knn_update[[tag1]][[tag1]]),FUN=function(k){
      temp<-knn_update[[tag1]][[tag1]]@i[(knn_update[[tag1]][[tag1]]@p[k]+1):knn_update[[tag1]][[tag1]]@p[k+1]]+1
      temp<-Matrix::colSums(knn_update[[tag1]][[tag1]][temp,])
      nonzero<-which(temp!=0)
      odr<-c(nonzero[order(temp[nonzero],decreasing=T)],setdiff(1:length(temp),nonzero))
      aa<-sort(odr[1:knn_transfer_correction])
      bb<-odr[odr%in%valid_anchor][1:knn_transfer_correction]
      return(c(aa,temp[aa],match(bb,valid_anchor),temp[bb]))
    })
    snn_mat1<-indexvalue2MM(output_index[1:(2*knn_transfer_correction)+2*knn_transfer_correction,],dim1=length(valid_anchor),dim2=nself)
    snn_mat2<-indexvalue2MM(output_index[1:(2*knn_transfer_correction),],dim1=nself,dim2=nself)
    rm(output_index)

    cat(tag1,": snn between anchors and non-anchors finished",sep="",fill=T)

    correction_valid<-apply(anchor_mat,2,function(x){
      colMeans(eigs$corrected[index_other[x==1],])
    })
    rm(anchor_mat)
    correction_valid<-t(correction_valid)-eigs$corrected[index_self[valid_anchor],]

    correction<-sapply(1:ncol(snn_mat1),function(k){
      index<-snn_mat1@i[(snn_mat1@p[k]+1):(snn_mat1@p[k+1])]+1
      wt<-snn_mat1@x[(snn_mat1@p[k]+1):(snn_mat1@p[k+1])]
      if(sum(wt)==0){
        return(rep(0,ncol(correction_valid)))
      }else{
        return((colSums(correction_valid[index,]*wt))/sum(wt))
      }
    })
    correction<-t(correction)
    rm(correction_valid,snn_mat1)

    correction<-sapply(1:ncol(snn_mat2),function(k){
      index<-snn_mat2@i[(snn_mat2@p[k]+1):(snn_mat2@p[k+1])]+1
      return((colMeans(correction[index,])))
    })
    correction<-t(correction)
    eigs$corrected[index_self,]<-eigs$corrected[index_self,]+correction
    rm(correction,snn_mat2)
    gc()
  }

  cat("eigs correction finished",fill=T)

  bin<-10000
  for(j in 1:ceiling(ncol(Smat)/bin)-1){
    col_index<-(j*bin+1):min((j*bin+bin),ncol(Smat))
    bb<-paste0(min(col_index),"-",max(col_index))
    cat("Calculating corrected similarities between cells",bb,fill=T)
    temp<-Smat[,col_index]-
      tcrossprod(eigs$vectors[Smat_index,],
                 t(t(eigs$vectors[col_index,])*eigs$values))+
      tcrossprod(eigs$corrected[Smat_index,],
                 t(t(eigs$corrected[col_index,])*eigs$values))
    Smat[,col_index]<-temp
    rm(temp)
  }

  cat("similarity matrix correction finished",fill=T)
  return(eigs)
}

#' Correct batch effects on similarity matrix in parallel
#'
#' Correct batch effects on similarity matrix in parallel
#'
#' @param Smat the similarity matrix.
#' @param row_sample indices of sampled cells for rows.
#' @param batch the batch information of cells.
#' @param reference the name of reference batch. Cells are sequentially aligned according the order of this vector.
#' @param neigs the number of Eigen vectors to calculate.
#' @param knn_update the knn matrix between batches.
#' @param ncore number of threads.
#'
#' @examples
#' \dontrun{
#' res_correct<-epiConv.correct.parallel(Smat=Smat,
#'                                       row_sample=1:ncol(Smat),
#'                                       batch=batch,
#'                                       reference=c("Batch1","Batch2"),
#'                                       neigs=30,
#'                                       knn_update=res_anchor$knn_update,
#'                                       ncore=5)
#'
#'}
#'
#' @return Return a list contains Eigen vectors before and after correction.


epiConv.correct.parallel<-function(Smat,row_sample,batch,reference,neigs=30,knn_update,ncore=1){

  Smat_index<-row_sample
  eigs<-dim.reduce(Smat=Smat,Smat_index=Smat_index,neigs=neigs)
  eigs$corrected<-eigs$vectors
  cat("decomposition finished",fill=T)
  gc()

  anchor_threshold<-5
  knn_transfer_correction<-10
  align_order<-c(reference,setdiff(levels(batch),reference))
  cl<-makeCluster(ncore,type="FORK")
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    index_self<-which(batch==tag1)
    nself<-length(index_self)
    index_other<-sapply(intersect(reference,align_order[1:(i-1)]),function(x){
      which(batch==x)
    })
    index_other<-unlist(index_other)
    index_other<-index_other[index_other%in%Smat_index]

    anchor_mat<-NULL
    for(tag2 in intersect(reference,align_order[1:(i-1)])){
      anchor_mat<-rbind(anchor_mat,knn_update[[tag1]][[tag2]])
    }
    valid_anchor<-which(Matrix::colSums(anchor_mat)>=anchor_threshold)
    anchor_mat<-anchor_mat[,valid_anchor]
    cat(tag1,": ",length(valid_anchor)," anchors in ",length(index_self)," cells",sep="",fill=T)


    output_index<-parSapply(cl=cl,X=1:ncol(knn_update[[tag1]][[tag1]]),FUN=function(k){
      temp<-knn_update[[tag1]][[tag1]]@i[(knn_update[[tag1]][[tag1]]@p[k]+1):knn_update[[tag1]][[tag1]]@p[k+1]]+1
      temp<-Matrix::colSums(knn_update[[tag1]][[tag1]][temp,])
      nonzero<-which(temp!=0)
      odr<-c(nonzero[order(temp[nonzero],decreasing=T)],setdiff(1:length(temp),nonzero))
      aa<-sort(odr[1:knn_transfer_correction])
      bb<-odr[odr%in%valid_anchor][1:knn_transfer_correction]
      return(c(aa,temp[aa],match(bb,valid_anchor),temp[bb]))
    })
    snn_mat1<-indexvalue2MM(output_index[1:(2*knn_transfer_correction)+2*knn_transfer_correction,],dim1=length(valid_anchor),dim2=nself)
    snn_mat2<-indexvalue2MM(output_index[1:(2*knn_transfer_correction),],dim1=nself,dim2=nself)
    rm(output_index)

    cat(tag1,": snn between anchors and non-anchors finished",sep="",fill=T)

    correction_valid<-apply(anchor_mat,2,function(x){
      colMeans(eigs$corrected[index_other[x==1],])
    })
    rm(anchor_mat)
    correction_valid<-t(correction_valid)-eigs$corrected[index_self[valid_anchor],]

    correction<-sapply(1:ncol(snn_mat1),function(k){
      index<-snn_mat1@i[(snn_mat1@p[k]+1):(snn_mat1@p[k+1])]+1
      wt<-snn_mat1@x[(snn_mat1@p[k]+1):(snn_mat1@p[k+1])]
      if(sum(wt)==0){
        return(rep(0,ncol(correction_valid)))
      }else{
        return((colSums(correction_valid[index,]*wt))/sum(wt))
      }
    })
    correction<-t(correction)
    rm(correction_valid,snn_mat1)

    correction<-sapply(1:ncol(snn_mat2),function(k){
      index<-snn_mat2@i[(snn_mat2@p[k]+1):(snn_mat2@p[k+1])]+1
      return((colMeans(correction[index,])))
    })
    correction<-t(correction)
    eigs$corrected[index_self,]<-eigs$corrected[index_self,]+correction
    rm(correction,snn_mat2)
    gc()
  }
  stopCluster(cl)

  cat("eigs correction finished",fill=T)

  bin<-10000
  for(j in 1:ceiling(ncol(Smat)/bin)-1){
    col_index<-(j*bin+1):min((j*bin+bin),ncol(Smat))
    bb<-paste0(min(col_index),"-",max(col_index))
    cat("Calculating corrected similarities between cells",bb,fill=T)
    temp<-Smat[,col_index]-
      tcrossprod(eigs$vectors[Smat_index,],
                 t(t(eigs$vectors[col_index,])*eigs$values))+
      tcrossprod(eigs$corrected[Smat_index,],
                 t(t(eigs$corrected[col_index,])*eigs$values))
    Smat[,col_index]<-temp
    rm(temp)
  }

  cat("similarity matrix correction finished",fill=T)
  return(eigs)
}


#' calculate shared-nearest-neighbor(SNN) matrix
#'
#' calculate shared-nearest-neighbor(SNN) matrix
#'
#' @param Smat the similarity matrix.
#' @param knn number of nearest neighbors kept.
#' @param ncore number of threads.
#'
#' @examples
#' \dontrun{
#' snn_mat<-Smat2snn(Smat=Smat,
#'                   knn=20,
#'                   ncore=5)
#'
#'}
#'
#' @return a sparse matrix


Smat2snn<-function(Smat,knn,ncore=1){

  mat<-find.knn(Smat,knn=max(floor(nrow(Smat)/100),20),noise=T)
  if(ncore==1){
    output_index<-sapply(X=1:ncol(mat),FUN=function(k){
      temp<-mat@i[(mat@p[k]+1):mat@p[k+1]]+1
      temp<-Matrix::colSums(mat[temp,])
      nonzero<-which(temp!=0)
      if(length(nonzero)>knn){
        aa<-sort(nonzero[order(temp[nonzero],decreasing=T)[1:knn]])
      }else{
        aa<-sort(c(nonzero,setdiff(1:length(temp),nonzero))[1:knn])
      }

      return(c(aa,temp[aa]))
    })
  }else{
    cl<-makeCluster(ncore,type="FORK")
    output_index<-parSapply(cl=cl,X=1:ncol(mat),FUN=function(k){
      temp<-mat@i[(mat@p[k]+1):mat@p[k+1]]+1
      temp<-Matrix::colSums(mat[temp,])
      nonzero<-which(temp!=0)
      if(length(nonzero)>knn){
        aa<-sort(nonzero[order(temp[nonzero],decreasing=T)[1:knn]])
      }else{
        aa<-sort(c(nonzero,setdiff(1:length(temp),nonzero))[1:knn])
      }
      return(c(aa,temp[aa]))
    })
    stopCluster(cl)
  }

  output<-new("dgCMatrix")
  output@Dim<-c(as.integer(ncol(mat)),as.integer(ncol(mat)))
  output@p<-as.integer((0:ncol(mat))*knn)
  output@i<-as.integer(output_index[1:knn,]-1)
  output@x<-as.vector(output_index[(knn+1):(2*knn),])
  return(output)
}

graph.norm<-function(x){
  knn<-length(x)
  target_fun<-function(x,similarity,knn){
    sum(exp(x*similarity))-log2(knn)
  }
  root<-uniroot(target_fun,interval=c(0,1e10),similarity=x-max(x),knn=knn)$root
  return(exp(root*(x-max(x))))
}

#' Perform louvain clustering on SNN matrix
#'
#' Perform louvain clustering on SNN matrix
#'
#' @param snn the SNN matrix.
#' @param resolution the resoltions for louvain clustering.
#'
#' @examples
#' \dontrun{
#' clust<-epiConv.louvain(snn=snn_mat,
#'                        resolution=c(0.8,0.6,0.4,0.2))
#' head(clust)
#'
#'}
#'
#' @return clustering results

epiConv.louvain<-function(snn,resolution,...){
  temp<-sapply(1:ncol(snn),function(i){
    graph.norm(snn@x[(snn@p[i]+1):snn@p[i+1]])
  })
  snn@x<-as.vector(temp)
  graph<-snn+Matrix::t(snn)-snn*Matrix::t(snn)
  rownames(graph)<-as.character(1:ncol(snn))
  colnames(graph)<-as.character(1:ncol(snn))

  graph<-Seurat::as.Graph(graph)

  clust<-matrix("",ncol(snn),length(resolution))
  for(i in 1:length(resolution)){
    temp<-Seurat::FindClusters(graph,resolution=resolution[i])[,1]
    clust[,i]<-sprintf(paste0("%0",ceiling(log10(nlevels(temp))),"d"),temp)
  }
  colnames(clust)<-c(paste0("resolution_",resolution))
  rownames(clust)<-rownames(snn)
  return(clust)
}
