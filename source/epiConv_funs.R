suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Seurat))
suppressMessages(library(bigmemory))
suppressMessages(library(biganalytics))
suppressMessages(library(PRIMME))
suppressMessages(library(ggplot2))
{
  setClass("epiConv",
           slots=list(mat="list",similarity="list",embedding="list",
                      ncell="numeric",meta.features="data.frame"))
  setMethod("[[", "epiConv", function(x,i){x@similarity[[i]]})
  setMethod("$", "epiConv", function(x,name){x@meta.features[[name]]})
  create.epiconv<-function(mat=list(),similarity=list(),embedding=list(),meta.features){
    output<-new("epiConv",
                mat=mat,
                similarity=similarity,
                embedding=embedding,
                ncell=nrow(meta.features),
                meta.features=meta.features)
    return(output)
  }
  add.similarity<-function(obj,x,name){
    if(is.vector(x)){
      temp<-matrix(0,obj@ncell,obj@ncell)
      temp[upper.tri(temp)]<-x
      temp<-temp+t(temp)
      obj@similarity[[name]]<-temp
    }else{
      obj@similarity[[name]]<-x
    }
    return(obj)
  }
  add.embedding<-function(obj,x,name){
    colnames(x)<-paste(name,"_",1:ncol(x),sep="")
    obj@embedding[[name]]<-x
    return(obj)
  }
  add.feature<-function(obj,x,name){
    obj@meta.features[[name]]<-x
    return(obj)
  }
}

big.crossprod<-function(mat1,mat2,bin=50000){
  if(ncol(mat1)<bin&ncol(mat2)<bin){
    output<-Matrix::crossprod(mat1,mat2)
  }else{
    output<-big.matrix(nrow=ncol(mat1),ncol=ncol(mat2),init=0)
    for(i in 1:ceiling(ncol(mat1)/bin)-1){
      row_index<-(i*bin+1):min((i*bin+bin),ncol(mat1))
      for(j in 1:ceiling(ncol(mat2)/bin)-1){
        col_index<-(j*bin+1):min((j*bin+bin),ncol(mat2))
        output[row_index,col_index]<-as.matrix(Matrix::crossprod(mat1[,row_index],mat2[,col_index]))
      }
    }
  }
  output<-as.matrix(output)
  rownames(output)<-colnames(mat1)
  colnames(output)<-colnames(mat2)
  return(output)
}

symm.crossprod<-function(mat,bin=50000,ncore=1){
  if(ncol(mat)<=bin){
    output<-Matrix::crossprod(mat)
  }else{
    output<-big.matrix(nrow=ncol(mat),ncol=ncol(mat),init=0)
    for(i in 1:ceiling(ncol(mat)/bin)-1){
      row_index<-(i*bin+1):min((i*bin+bin),ncol(mat))
      for(j in (i+1):ceiling(ncol(mat)/bin)-1){
        col_index<-(j*bin+1):min((j*bin+bin),ncol(mat))
        temp<-as.matrix(Matrix::crossprod(mat[,row_index],mat[,col_index]))
        output[row_index,col_index]<-temp
        if(i!=j){
          output[col_index,row_index]<-t(temp)
        }
      }
    }
  }
  output<-as.matrix(output)
  rownames(output)<-colnames(mat)
  colnames(output)<-colnames(mat)
  return(output)
}


big.tcrossprod<-function(mat1,mat2,bin=50000){
  if(nrow(mat1)<bin&nrow(mat2)<bin){
    output<-Matrix::tcrossprod(mat1,mat2)
  }else{
    output<-big.matrix(nrow=nrow(mat1),ncol=nrow(mat2),init=0)
    for(i in 1:ceiling(nrow(mat1)/bin)-1){
      row_index<-(i*bin+1):min((i*bin+bin),nrow(mat1))
      for(j in 1:ceiling(nrow(mat2)/bin)-1){
        col_index<-(j*bin+1):min((j*bin+bin),nrow(mat2))
        output[row_index,col_index]<-as.matrix(Matrix::tcrossprod(mat1[row_index,],mat2[col_index,]))
      }
    }
  }
  output<-as.matrix(output)
  rownames(output)<-rownames(mat1)
  colnames(output)<-rownames(mat2)
  return(output)
}

tfidf.norm<-function(mat,lib_size){
  mat@x[mat@x>0]<-1
  freq<-Matrix::rowSums(mat)/ncol(mat)
  wt<-log10(1+1/freq)
  mat<-t(t(mat)/lib_size)*wt
  return(mat)
}

inf.estimate<-function(mat,sample_size=0.2,nsim=15){
  sample_size<-floor(nrow(mat)*sample_size)
  sample_matrix<-lapply(1:nsim,function(x) sample(1:nrow(mat),size=sample_size))
  conv<-sapply(sample_matrix,function(x){
    temp<-as.matrix(Matrix::crossprod(mat[x,]))
    return(log10(min(temp[temp!=0])))
  })
  return(min(conv))
}

epiConv.matrix<-function(mat,bin=50000,inf_replace=-10){
  ncell<-ncol(mat)
  output<-log10(symm.crossprod(mat,mat,bin=bin))
  output[is.infinite(output)]<-inf_replace
  return(output)
}

dim.reduce<-function(Smat,neigs=50){
  eigs1<-eigs_sym(Smat,NEig=neigs,which="LA")
  eigs2<-eigs_sym(Smat,NEig=neigs,which="SA")
  values<-c(eigs1$values,eigs2$values)
  vectors<-cbind(eigs1$vectors,eigs2$vectors)
  odr<-order(abs(values),decreasing=T)[1:neigs]
  return(list(vectors=vectors[,odr],values=values[odr]))
}

sim.blur<-function(Smat,weight_scale=1,neighbor_frac=0.25,knn=20,bin=10000){
  ncell<-nrow(Smat)
  Wmat<-10^Smat
  Wmat<-biganalytics::apply(Wmat,2,function(x,weight_scale){
    x[order(x)[1:(length(x)-knn)]]<-0
    x<-x*weight_scale/sum(x*weight_scale)*neighbor_frac
    return(x)
  },weight_scale=weight_scale)
  diag(Wmat)<-1-neighbor_frac
  Smat_blurred<-big.crossprod(Smat,Wmat,bin=bin)
  Smat_blurred<-(Smat_blurred+t(Smat_blurred))/2
  return(Smat_blurred)
}

lib.estimate<-function(mat){
  mat@x[mat@x>1]<-1
  return(Matrix::colSums(mat))}
freq.estimate<-function(mat){
  mat@x[mat@x>1]<-1
  return(Matrix::rowSums(mat)/ncol(mat))
}

cal.snn<-function(Smat,knn=50,bin=5000){
  diag(Smat)<-apply(Smat,2,max)
  ncell<-nrow(Smat)
  neighbor_mat<-Matrix(apply(Smat,2,function(x,knn){
    output<-rep(0,length(x))
    output[order(x,decreasing=T)[1:knn]]<-1
    return(output)
  },knn=knn))
  output<-symm.crossprod(mat=neighbor_mat,bin=bin)
  return(output)
}

continuous.match<-function(x,win){
  output<-rep(1,length(x))
  for(i in 1:length(win)){
    output[x>=win[i]]<-output[x>=win[i]]+1
  }
  return(output)
}

num2int<-function(x,nbin=32,bound=quantile(x,c(0.05,0.95))){
  win<-(bound[2]-bound[1])/(nbin-2)*(0:(nbin-2))+bound[1]
  return(continuous.match(x,win))
}

stablization.pars<-function(Smat,lib_size){
  lib_mat<-matrix(0,nrow(Smat),ncol(Smat))
  lib_mat<-lib_mat+log10(lib_size[1:nrow(Smat)])
  lib_mat<-t(t(lib_mat)+log10(lib_size[(nrow(Smat)+1):(nrow(Smat)*2)]))
  
  lm_res<-lm(as.vector(Smat)~as.vector(lib_mat),weight=10^as.vector(lib_mat))
  Smat_residual<-matrix(lm_res$residuals,nrow(Smat),nrow(Smat))
  lib_discrete<-num2int(lib_mat,nbin=100,bound=quantile(lib_mat,c(0.01,0.99)))
  lib_window<-quantile(lib_mat,0.01)+(quantile(lib_mat,0.99)-quantile(lib_mat,0.01))/98*(0:98)
  sd_value<-tapply(Smat_residual,list(factor(lib_discrete)),function(x){
    sqrt(sum(x^2)/length(x))
  })
  return(list(coef=lm_res$coefficients,lib_window=lib_window,sd_value=sd_value))
}

similarity.stablization.small<-function(Smat,lib_size1,lib_size2,stablization_pars){
  lib_mat<-matrix(0,nrow(Smat),ncol(Smat))
  lib_mat<-lib_mat+log10(lib_size1)
  lib_mat<-t(t(lib_mat)+log10(lib_size2))
  lib_mat_discrete<-continuous.match(lib_mat,stablization_pars$lib_window)
  output<-(as.vector(Smat)-as.vector(lib_mat)*stablization_pars$coef[2]-stablization_pars$coef[1])
  output<-output/stablization_pars$sd_value[lib_mat_discrete]
  output<-matrix(output,nrow(Smat),ncol(Smat))
  return(output)
}


similarity.stablization<-function(Smat,lib_size,stablization_pars,bin=5000){
  if(nrow(Smat)<=bin){
    output<-similarity.stablization.small(Smat,
                                          lib_size1=lib_size,
                                          lib_size2=lib_size,
                                          stablization_pars=stablization_pars)
  }else{
    output<-big.matrix(nrow=nrow(Smat),ncol=ncol(Smat),init=0)
    for(i in 1:ceiling(nrow(Smat)/bin)-1){
      row_index<-(i*bin+1):min((i*bin+bin),nrow(Smat))
      for(j in (i+1):ceiling(nrow(Smat)/bin)-1){
        col_index<-(j*bin+1):min((j*bin+bin),nrow(Smat))
        Smat_small<-Smat[row_index,col_index]
        temp<-similarity.stablization.small(Smat[row_index,col_index],
                                            lib_size1=lib_size[row_index],
                                            lib_size2=lib_size[col_index],
                                            stablization_pars=stablization_pars)
        if(i==j){
          output[row_index,col_index]<-temp
        }else{
          output[row_index,col_index]<-temp
          output[col_index,row_index]<-t(temp)
        }
      }
    }
    output<-as.matrix(output)
  }
  return(output)
}

find.knn<-function(Smat,knn,binarize=T){
  apply(Smat,2,function(x,knn){
    output<-rep(0,length(x))
    cutoff<-sort(x,decreasing=T)[knn]
    if(binarize){
      output[x>=cutoff]<-1
    }else{
      output[x>=cutoff]<-x[x>=cutoff]
    }
    return(output)
  },knn=knn)
}

reference.range<-function(features,reference,threshold=1.5){
  output<-apply(reference,2,function(x){
    mu<-mean(features[x!=0])
    sdev<-sd(features[x!=0])
    return(c(mu-threshold*sdev,mu+threshold*sdev))
    # return(quantile(features[x!=0],probs=c(threshold,1-threshold),na.rm=T))
  })
  return(c(output[1,],output[2,]))
} ##accept one feature and knn matrix, return c(lower bound vector, upper bound vector)

neighbor_vector.update<-function(neighbor_vector,features,ref_range){
  index<-which(neighbor_vector!=0)
  ref_range<-ref_range[c(index,index+length(neighbor_vector)),]
  temp<-sapply(1:length(features),function(i){
    rng<-matrix(ref_range[,i],,2)
    if(is.na(features[i])){
      return(rep(1,nrow(rng)))
    }else{
      return(as.numeric(features[i]>rng[,1]&features[i]<rng[,2]))
    }
  })
  neighbor_vector[index[rowSums(temp)!=ncol(temp)]]<-0
  return(neighbor_vector)
} ##accept one knn vector, multiple features of one object and range matrix, returns updated knn vector

refine.knn<-function(target,reference,features,threshold=1.5){
  reference_range<-apply(features,2,reference.range,reference=reference,threshold=threshold)
  target_update<-sapply(1:ncol(target),function(x){
    neighbor_vector.update(target[,x],features=features[x,],ref_range=reference_range)
  })
  return(target_update)
}

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

eigs.snn<-function(Smat,features,batch,reference,knn_target=50,knn_reference=20,threshold=1.5){
  self_knn<-list()
  align_order<-c(reference,setdiff(levels(batch),reference))
  output<-list()
  for(k in levels(batch)){
    index<-which(batch==k)
    self_knn[[k]]<-find.knn(Smat[index,index],
                            knn=max(floor(length(index)*0.01),20))
  }
  for(i in 2:length(align_order)){
    tag1<-align_order[i]
    output[[tag1]]<-list()
    index1<-which(batch==tag1)
    temp<-symm.crossprod(mat=self_knn[[tag1]])
    temp<-Matrix(data=temp,sparse=T,doDiag=F)
    output[[tag1]][[tag1]]<-temp
    
    if(!is.null(features[[tag1]])){
      for(tag2 in setdiff(reference,align_order[i:length(align_order)])){
        index2<-which(batch==tag2)
        temp1<-find.knn(Smat[index2,index1],knn=max(20,floor(length(index2)*0.01)))
        temp1<-big.crossprod(self_knn[[tag2]],temp1)
        temp1<-find.knn(temp1,knn=knn_target,binarize=F)
        
        temp2<-find.knn(Smat[index1,index2],knn=max(20,floor(length(index1)*0.01)))
        temp2<-big.crossprod(self_knn[[tag1]],temp2)
        temp2<-find.knn(temp2,knn=knn_reference,binarize=F)
        
        temp<-refine.knn(target=temp1,
                         reference=temp2,
                         features=features[[tag1]],
                         threshold=threshold)
        temp<-Matrix(data=temp,sparse=T,doDiag=F)
        output[[tag1]][[tag2]]<-temp
        
      }
    }else{
      for(tag2 in setdiff(reference,align_order[i:length(align_order)])){
        index2<-which(batch==tag2)
        
        temp1<-find.knn(Smat[index2,index1],knn=max(20,floor(length(index2)*0.01)))
        temp1<-big.crossprod(self_knn[[tag2]],temp1)
        temp1<-find.knn(temp1,knn=knn_target,binarize=F)
        
        temp<-Matrix(data=temp1,sparse=T,doDiag=F)
        output[[tag1]][[tag2]]<-temp
        
      }
    }
  }
  return(output)
}

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
    correction<-t(correction)
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

run.epiConv<-function(mat,lib_size,nbootstrap=15,nsample=floor(nrow(mat)*0.2),bin=10000,inf_replace=(-8)){
  mat<-tfidf.norm(mat,lib_size=lib_size)
  Smat<-matrix(0,ncol(mat),ncol(mat))
  sample_matrix<-lapply(1:nbootstrap,function(x) sample(1:nrow(mat),size=nsample))
  for(i in sample_matrix){
    Smat<-Smat+epiConv.matrix(mat=mat[i,],inf_replace=inf_replace,bin=bin)
  }
  Smat<-Smat/nbootstrap
  
  cell_sample<-500
  temp<-sample(1:nrow(Smat),cell_sample*2)
  retain1<-sort(temp[1:cell_sample])
  retain2<-sort(temp[(cell_sample+1):(cell_sample*2)])
  stablization_pars<-stablization.pars(Smat=Smat[retain1,retain2],
                                       lib_size=lib_size[c(retain1,retain2)])
  Smat<-similarity.stablization(Smat,
                                lib_size=lib_size,
                                stablization_pars=stablization_pars,
                                bin=bin)
  return(Smat)
}

cal.feature<-function(Smat,pcs,neigs=50){
  eigs<-dim.reduce(Smat-median(Smat),neigs=neigs)
  ratio<-sqrt(sum(apply(pcs,2,function(x) sum(x^2)))/sum(abs(eigs$values))*abs(eigs$values))
  return(cbind(pcs,t(t(eigs$vectors)*ratio)))
}

graph.norm<-function(x,knn){
  cutoff<-sort(x,decreasing=T)[knn]
  x[x<cutoff]<-0
  y<-x[x!=0]
  target_fun<-function(x,similarity,knn){
    sum(exp(x*similarity))-log2(knn)
  }
  root<-uniroot(target_fun,interval=c(0,1e10),similarity=y-max(y),knn=knn)$root
  x[x!=0]<-exp(root*(y-max(y)))
  return(x)
}

run.umap_louvain<-function(Smat,knn=20,umap_settings=NULL,resolution=NULL){
  knn<-max(min(floor(nrow(Smat)*0.01),200),20)
  snn_mat<-cal.snn(Smat,knn=knn)+Smat/1e4
  if(is.null(umap_settings)){
    umap_settings<-umap::umap.defaults
  }
  umap_settings$input<-"dist"
  umap_settings$n_neighbors<-knn
  umap_res<-umap::umap(max(snn_mat)-snn_mat,config=umap_settings)$layout
  
  if(is.null(resolution)){
    return(umap_res)
  }
  
  diag(snn_mat)<-0
  temp<-apply(snn_mat,2,graph.norm,knn=knn)
  temp<-temp+t(temp)-temp*t(temp)
  rownames(temp)<-as.character(1:nrow(Smat))
  colnames(temp)<-as.character(1:nrow(Smat))
  
  graph<-as.Graph(Matrix(temp,sparse=T))
  clust<-matrix("",nrow(Smat),length(resolution))
  for(i in 1:length(resolution)){
    temp<-FindClusters(graph,resolution=resolution[i])[,1]
    clust[,i]<-sprintf(paste0("%0",ceiling(log10(nlevels(temp))),"d"),temp)
  }
  output<-data.frame(umap_res,clust)
  colnames(output)<-c("umap_1","umap_2",paste0(resolution,resolution))
  return(output)
}