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
  output<-bigmemory::as.matrix(output)
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
          output[col_index,row_index]<-Matrix::t(temp)
        }
      }
    }
  }
  output<-bigmemory::as.matrix(output)
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
  output<-bigmemory::as.matrix(output)
  rownames(output)<-rownames(mat1)
  colnames(output)<-rownames(mat2)
  return(output)
}

tfidf.norm<-function(mat,lib_size){
  mat@x[mat@x>0]<-1
  freq<-Matrix::rowSums(mat)/ncol(mat)
  wt<-log10(1+1/freq)
  wt[freq==0]<-0
  mat<-Matrix::t(Matrix::t(mat)/lib_size)*wt
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

similarity.stablization.small<-function(Smat,lib_size1,lib_size2,stablization_pars){
  lib_mat<-matrix(0,nrow(Smat),ncol(Smat))
  lib_mat<-lib_mat+log10(lib_size1)
  lib_mat<-Matrix::t(Matrix::t(lib_mat)+log10(lib_size2))
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
          output[col_index,row_index]<-Matrix::t(temp)
        }
      }
    }
    output<-bigmemory::as.matrix(output)
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
