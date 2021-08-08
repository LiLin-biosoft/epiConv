find.knn2<-function(mat,knn){
  output<-new("dgCMatrix")
  output@Dim<-as.integer(dim(mat))
  output@p<-as.integer(0)
  nelement<-as.integer(0)
  for(k in 1:ncol(mat)){
    temp<-mat@i[(mat@p[k]+1):mat@p[k+1]]+1
    temp<-Matrix::colSums(mat[temp,])
    cutoff<-sort(temp,decreasing=T)[knn]
    index<-which(temp>=cutoff)
    output@i<-c(output@i,as.integer(index-1))
    nelement<-nelement+length(index)
    output@p<-c(output@p,nelement)
    output@x<-c(output@x,rep(1,length(index)))
  }
  return(output)
}

dim.reduce<-function(Smat,neigs){
  eigs1<-PRIMME::eigs_sym(Smat,NEig=neigs,which="LA")
  eigs2<-PRIMME::eigs_sym(Smat,NEig=neigs,which="SA")
  values<-c(eigs1$values,eigs2$values)
  vectors<-cbind(eigs1$vectors,eigs2$vectors)
  odr<-order(abs(values),decreasing=T)[1:neigs]
  return(list(vectors=vectors[,odr],values=values[odr]))
}

tfidf.norm<-function(mat,lib_size){
  mat@x[mat@x>0]<-1
  freq<-Matrix::rowSums(mat)/ncol(mat)
  wt<-log10(1+1/freq)
  wt[freq==0]<-0
  mat<-Matrix::t(Matrix::t(mat)/lib_size)*wt
  return(mat)
}

num2int<-function(x,nbin=32,bound=quantile(x,c(0.05,0.95))){
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

#################################
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
    Smat<-Smat-t(matrix(log10(lib_size2),ncell2,ncell1))*adjust_pars$coef[2]
    sd_mat<-matrix(log10(lib_size1),ncell1,ncell2)+
      t(matrix(log10(lib_size2),ncell2,ncell1))
    sd_mat<-apply(sd_mat,2,lib2sd,adjust_pars=adjust_pars)
    Smat<-Smat/sd_mat
    rm(sd_mat)
    return(Smat)
  }
}


refine.knn<-function(Smat,ref_knn,target_index,reference_index,knn_target,knn_reference,features,threshold=1.5){

  unsmoothed<-find.knn(Smat[target_index,reference_index],knn=knn_reference)
  smoothed<-Matrix::crossprod(Matrix::t(unsmoothed),ref_knn)

  ######calculate target matrix#########
  nfeature<-ncol(features)
  reference_range<-matrix(numeric(),2*nfeature,length(reference_index))
  #  target<-deepcopy(Smat,cols=target_index,rows=reference_index)
  target_MM<-new("dgCMatrix")
  #  target_MM@Dim<-as.integer(dim(target))
  target_MM@Dim<-c(length(reference_index),length(target_index))
  target_MM@p<-as.integer(0)
  nelement<-as.integer(0)
  for(i in 1:length(target_index)){
    edge<-sort(order(Smat[reference_index,target_index[i]],decreasing=T)[1:knn_target])
    nodata<-edge[is.na(reference_range[1,edge])]
    for(j in nodata){
      neighbor_index<-order(smoothed[,j]+Smat[target_index,reference_index[j]]/100,decreasing=T)[1:knn_reference]
      temp<-features[neighbor_index,]
      mu<-colMeans(temp)
      sdev<-apply(temp,2,sd)
      reference_range[,j]<-c(mu-sdev*threshold,mu+sdev*threshold)
    }

    temp<-reference_range[,edge]
    temp<-colSums(features[i,]>temp[1:nfeature,]&
                    features[i,]<temp[(nfeature+1):(nfeature*2),])
    edge<-edge[temp==nfeature]
    target_MM@i<-c(target_MM@i,as.integer(edge-1))
    nelement<-nelement+length(edge)
    target_MM@p<-c(target_MM@p,nelement)
    target_MM@x<-c(target_MM@x,rep(1,length(edge)))
  }
  return(target_MM)
  ######calculate target matrix#########
}

eigs.knn<-function(Smat,features,batch,reference,knn_target=50,knn_reference=20,threshold=2){
  align_order<-c(reference,setdiff(levels(batch),reference))
  ref_knn<-list()
  output<-list()
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    output[[tag1]]<-list()
    index1<-which(batch==tag1)

    for(tag2 in setdiff(reference,align_order[i:length(align_order)])){
      index2<-which(batch==tag2)
      if(is.null(ref_knn[[tag2]])){
        temp<-find.knn(Smat[index2,index2],knn=floor(length(index2)/100))
        ref_knn[[tag2]]<-find.knn2(temp,knn=knn_reference)
      }
      cat("calculating knn matrix:",tag1,"to",tag2,fill=T)
      output[[tag1]][[tag2]]<-refine.knn(Smat=Smat,
                                         ref_knn=ref_knn[[tag2]],
                                         target_index=index1,
                                         reference_index=index2,
                                         knn_target=knn_target,
                                         knn_reference=knn_reference,
                                         features=features[[tag1]],
                                         threshold=threshold)
    }
  }
  return(output)
}

eigs.scale<-function(eigs,knn_mat,batch,reference,knn_transfer_correction=10,anchor_threshold=5){
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
    {
      temp<-which(!is.na(eigs[index_other,1]))
      anchor_mat<-anchor_mat[temp,]
      index_other<-index_other[temp]

      valid_anchor<-which(Matrix::colSums(anchor_mat)>=anchor_threshold)
      correction<-apply(anchor_mat[,valid_anchor],2,function(x){
        colMeans(eigs[index_other[x==1],])
      })
      correction<-t(correction)-eigs[index_self[valid_anchor],]
    }

    snn_mat<-Matrix::crossprod(knn_mat[[tag1]][[tag1]][,valid_anchor],
                               knn_mat[[tag1]][[tag1]])
    correction<-apply(snn_mat,2,function(x){
      cutoff<-sort(x,decreasing=T)[knn_transfer_correction]
      x[x<cutoff]<-0
      if(sum(x!=0)==0)
        x<-rep(1,length(x))
      (colSums(correction*x))/sum(x)
    })
    correction<-t(correction)

    eigs[index_self,]<-eigs[index_self,]+correction
  }
  return(eigs)
}

graph.norm<-function(x,knn){
  idx<-order(x,decreasing=T)[1:knn]
  y<-x[idx]
  target_fun<-function(x,similarity,knn){
    sum(exp(x*similarity))-log2(knn)
  }
  root<-uniroot(target_fun,interval=c(0,1e10),similarity=y-max(y),knn=knn)$root
  x<-rep(0,length(x))
  x[idx]<-exp(root*(y-max(y)))
  return(x)
}
