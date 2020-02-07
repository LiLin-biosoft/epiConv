suppressMessages(library(Matrix))
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

}

