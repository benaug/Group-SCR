dztpois=function(k,lambda){
  return(lambda^k/((exp(lambda)-1)*factorial(k)))
}

init.group=function(data,M=NA,inits=NA){
  buff<- data$buff
  X<-as.matrix(data$X)
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  lambda.P=inits$lambda.P
  lam0.P=inits$lam0.P
  sigma.P=inits$sigma.P
  lambda.I=inits$lambda.I
  y.I=data$y.I.obs
  ni=nrow(y.I)
  J=dim(y.I)[2]
  K=dim(y.I)[3]
  groupID=data$groupID.obs
  
  #initialize group activity centers
  s=cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2]))
  
  #make sure groups with captured individuals turned on
  z=rep(0,M)
  z[unique(groupID)]=1
  
  #improve initial s for captured groups
  capgroups=unique(groupID)
  y.I2D=apply(y.I,c(1,2),sum)
  for(i in 1:length(capgroups)){
    these.guys=which(groupID==capgroups[i])
    y.tmp=y.I2D[these.guys,]
    if(length(these.guys)>1){
      y.tmp=colSums(y.tmp)
    }
    trapcaps=which(y.tmp>0)
    if(length(trapcaps)>1){
      meanloc=colMeans(X[trapcaps,])
    }else{
      meanloc=X[trapcaps,]
    }
    s[capgroups[i],]=meanloc
  }
  
  
  D<- e2dist(s, X)
  lamd=lam0.P*exp(-D*D/(2*sigma.P*sigma.P))
  y.P=array(0,dim=c(M,J,K))
  for(i in 1:M){
    for(j in 1:J){
      for(k in 1:K){
        y.P[i,j,k]=rpois(1,lamd[i,j])
      }
    }
  }
  #make sure y.P>0 if sample assigned to it
  for(i in 1:ni){
    idx=which(y.I[i,,]>0,arr.ind=TRUE)
    if(K>1){
      for(j in 1:nrow(idx)){
        if(y.P[groupID[i],idx[j,1],idx[j,2]]==0){
          y.P[groupID[i],idx[j,1],idx[j,2]]=1
        }
        
      }
    }else{
      for(j in 1:length(idx)){
        if(y.P[groupID[i],idx[j],1]==0){
          y.P[groupID[i],idx[j],1]=1
        }
      }
    }
  }
  
  #group site use likelihood
  D<- e2dist(s, X)
  lamd=lamd.cand=lam0.P*exp(-D*D/(2*sigma.P*sigma.P))
  ll.yp=dpois(y.P,lamd,log=TRUE)
  ll.yp.sum=sum(ll.yp)
  if(!is.finite(sum(ll.yp.sum)))stop("Initial site visitation likelihood not finite. Try raising lam0.P and/or sigma.P")
  
  #individual likelihood given group site use
  lambda.tmp=y.P[groupID,,]*lambda.I #reorder by individual groupID
  ll.yi=dpois(y.I,lambda.tmp,log=TRUE)
  ll.yi.cand=ll.yi
  ll.yi.sum=sum(ll.yi)
  if(!is.finite(sum(ll.yi.sum)))stop("Initial individual detection likelihood not finite. Shouldn't happen.")
  
  counts=counts.detected=ll.group=rep(0,M)
  for(i in 1:M){
    if(z[i]==1){
      counts[i]=sum(groupID==i)
      counts.detected[i]=counts[i]
      ll.group[i]=log(dztpois(counts[i],lambda.P))
    }else{
      counts[i]=rzapois(1,lambda.P)
      ll.group[i]=log(dztpois(counts[i],lambda.P))
    }
  }
  ll.group.sum=sum(ll.group)
  if(!is.finite(sum(ll.group.sum)))stop("Initial group size likelihood not finite. Shouldn't happen.")
  return(list(z=z,s=s,y.I=y.I,y.P=y.P,counts=counts,counts.detected=counts.detected,groupID=groupID))
}