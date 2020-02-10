suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(bigmemory))
suppressMessages(library(biganalytics))
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
  add.mat<-function(obj,x,name){
    obj@mat[[name]]<-x
    return(obj)
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

  embedding.plot<-function(obj,name,ident,pch=16,pt.size=2,cols=NULL,plot.na=T){
    if(length(ident)==1){
      colby<-obj@meta.features[[ident]]
    }else{
      colby=ident
    }
    if(class(colby)=="character"){
      colby<-factor(colby)
    }
    p<-ggplot(data=as.data.frame(obj@embedding[[name]]))+labs(col=NULL)+
      geom_point(show.legend=T,pch=pch,size=pt.size,
                 mapping=aes_string(x=paste(name,"_1",sep=""),y=paste(name,"_2",sep=""),
                                    col=colby))

  if(class(colby)=="factor"){
    if(is.null(cols))
      cols<-c("red","chocolate","darkorange","gold","darkgreen","green3","darkolivegreen1",
              "blue","lightblue","cyan","lightpink","magenta","purple")
  cols_scale<-colorRampPalette(cols)(length(levels(colby)))
  p<-p+scale_color_manual(values=cols_scale,na.value="gray",na.translate=plot.na)

  }else if(class(colby)=="numeric"){
    if(is.null(cols))
      cols<-c("darkblue","yellow")
    p<-p+scale_colour_gradient(low=cols[1],high=cols[2],na.value="gray")
  }
    return(p)
  }

}

tfidf.norm<-function(mat,lib_size){
  mat@x[mat@x>1]<-1
  freq<-Matrix::rowSums(mat)/ncol(mat)
  wt<-log10(1+1/freq)
  mat<-t(t(mat)/lib_size)*wt
  return(mat)
}

inf.estimate<-function(mat,sample_size=0.125,nsim=30){
  sample_size<-floor(nrow(mat)*sample_size)
  sample_matrix<-lapply(1:nsim,function(x) sample(1:nrow(mat),size=sample_size))
  conv<-sapply(sample_matrix,function(x){
    temp<-as.matrix(t(mat[x,])%*%mat[x,])
    return(log10(min(temp[temp!=0])))
  })
  return(min(conv))
}

epiConv.matrix<-function(mat,bin=10000,inf_replace=-10){
  ncell<-ncol(mat)
  if(ncell<=bin){
    output<-log10(as.matrix(t(mat)%*%mat))
    output[is.infinite(output)]<-inf_replace
  }else{
    block_index<-lapply(1:ceiling(ncell/bin),function(x) {
      temp<-1:ncell
      temp[temp>=((x-1)*bin+1)&temp<=(x*bin)]
    })
    aa<-matrix(,ceiling(ncell/bin),ceiling(ncell/bin))
    block_multiple<-matrix(c(row(aa)[upper.tri(aa,diag=T)],col(aa)[upper.tri(aa,diag=T)]),,2)
    output<-big.matrix(nrow=ncell,ncol=ncell,init=0)
    for(j in 1:nrow(block_multiple)){
      index1<-block_index[[block_multiple[j,1]]]
      index2<-block_index[[block_multiple[j,2]]]
      temp<-log10(as.matrix(t(mat[,index1])%*%mat[,index2]))
      temp[is.infinite(temp)]<-inf_replace
      if(block_multiple[j,1]!=block_multiple[j,2]){
      	output[index1,index2]<-temp
	output[index2,index1]<-t(temp)
      }else{
        output[index1,index2]<-temp
      }
    }
    output<-as.matrix(output)
  }
  return(output)
}

sim.blur<-function(Smat,weight_scale=1,neighbor_frac=0.25,knn=20,bin=10000){
  ncell<-nrow(Smat)
  Wmat<-10^Smat
  diag(Wmat)<-0
  Wmat<-biganalytics::apply(Wmat,2,function(x,weight_scale){
    x[order(x)[1:(length(x)-knn)]]<-0
    x<-x*weight_scale/sum(x*weight_scale)*neighbor_frac
    return(x)
  },weight_scale=weight_scale)
  diag(Wmat)<-1-neighbor_frac
  diag(Smat)<-0
  diag(Smat)<-biganalytics::apply(Smat,2,function(x) quantile(x,0.99))
  if(ncell<=bin){
    Smat_blurred<-Smat%*%Wmat
    Smat_blurred<-(Smat_blurred+t(Smat_blurred))/2
    diag(Smat_blurred)<-0
    return(as.matrix(Smat_blurred))
  }else{
    nbin<-ceiling(ncell/bin)
    block_index<-lapply(1:nbin,function(x) {
      temp<-1:ncell
      temp[temp>=((x-1)*bin+1)&temp<=(x*bin)]
    })
    Smat_blurred<-big.matrix(nrow=ncell,ncol=ncell,init=0)
    for(j in 1:nbin){
      for(k in 1:nbin){
        index1<-block_index[[j]]
        index2<-block_index[[k]]
        temp<-as.matrix(Smat[index1,]%*%Wmat[,index2])/2
        Smat_blurred[index1,index2]<-Smat_blurred[index1,index2]+temp
        Smat_blurred[index2,index1]<-Smat_blurred[index2,index1]+t(temp)
      }
    }
    diag(Smat_blurred)<-0
    return(as.matrix(Smat_blurred))
  }
}

lib.estimate<-function(mat){
  mat@x[mat@x>1]<-1
  return(Matrix::colSums(mat))}
freq.estimate<-function(mat){
  mat@x[mat@x>1]<-1
  return(Matrix::rowSums(mat)/ncol(mat))
}

sd.cal<-function(total,positive,sample){
  var<-sample*(positive/total)*((total-positive)/total)*(total-sample)/(total-1)
  return(sqrt(var))
}
z.cal<-function(x1,x2){
  output<-(sum(x1*x2)-sum(x1)*sum(x2)/length(x1))/sd.cal(total=length(x1),positive=sum(x1),sample=sum(x2))
  return(output)
}
p.cal<-function(x1,x2){
  output<-phyper(q=sum(x1*x2),m=sum(x1),n=length(x1)-sum(x1),k=sum(x2),lower.tail=F)
  return(output)
}

zscore.bulk<-function(mat,clust){
  mat@x[mat@x>0]<-1
  clust<-factor(clust)
  output<-apply(mat,1,function(x,clust){
    zscore<-sapply(levels(clust),function(y,x,clust){
      z.cal(x,as.numeric(clust==y))
    },x=x,clust=clust)
    pvalue<-sapply(levels(clust),function(y,x,clust){
      p.cal(x,as.numeric(clust==y))
    },x=x,clust=clust)
    return(as.vector(rbind(zscore,pvalue)))
  },clust=clust)
  output<-t(output)
  rownames(output)<-rownames(mat)
  colnames(output)<-paste0(sapply(levels(clust),rep,times=2),"_",c("z","p"))
  return(output)
}

zscore.single<-function(mat,Smat,qt=0.05,lib_size=NULL){
  diag(Smat)<-apply(Smat,2,max)
  Smat<-apply(Smat,2,function(x){
    output<-rep(0,length(x))
    output[x>quantile(x,1-qt)]<-1
    return(output)
  })
  mat@x[mat@x>0]<-1
  freq<-Matrix::rowSums(mat)
  zmat<-mat%*%Smat
  colnames(zmat)<-colnames(mat)
  if(is.null(lib_size)){
    n_neighbour<-sum(Smat[,1])
    peak_sd<-sd.cal(total=ncol(mat),positive=n_neighbour,sample=freq)
    zmat<-(zmat-n_neighbour/ncol(mat)*freq)/peak_sd
    return(as.matrix(zmat))
  }  
  neighbour_size<-Matrix::colSums(Smat*lib_size)/mean(lib_size)
  mean_mat<-sapply(neighbour_size,function(x) x/ncol(mat)*freq)
  sd_mat<-sapply(neighbour_size,function(y) sd.cal(total=ncol(mat),positive=y,sample=freq))
  zmat<-(zmat-mean_mat)/sd_mat
  return(as.matrix(zmat))
}
