e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


sim.GroupSCR.Bernoulli <- function(Np=10,lambda.P=6,p0.P=0.2,sigma.P=0.50,p.I=0.2,K=10,X=X,buff=3){
  library(VGAM)
  # # simulate a population of group activity centers
  s <- cbind(runif(Np, min(X[,1])-buff,max(X[,1])+buff), runif(Np,min(X[,2])-buff,max(X[,2])+buff))
  D <- e2dist(s,X)
  pd.P <- p0.P*exp(-D*D/(2*sigma.P*sigma.P))
  J <- nrow(X)
  # Simulate latent group site use history
  y.P <- array(0,dim=c(Np,J,K))
  for(i in 1:Np){
    for(j in 1:J){
      for(k in 1:K){
        y.P[i,j,k] <- rbinom(1,1,pd.P[i,j])
      }
    }
  }
  #Group membership stuff
  NperGroup <- rzapois(Np,lambda.P)
  Ni <- sum(NperGroup)
  groupID <- rep(1:Np,times=NperGroup)
  
  #Simulate observed individual history
  y.I <- array(0,dim=c(Ni,J,K))
  for(i in 1:Ni){
    for(j in 1:J){
      for(k in 1:K){
        y.I[i,j,k] <- rbinom(1,1,p.I*y.P[groupID[i],j,k])
      }
    }
  }
  si <- cbind(rep(s[,1],times=NperGroup),rep(s[,2],times=NperGroup))
  keep <- which(rowSums(y.I)>0)
  y.I.obs <- y.I[keep,,]
  if(K==1){
    y.I.obs <- array(y.I.obs,dim=c(dim(y.I.obs),K))
  }
  si.obs <- si[keep,]
  groupID.obs <- groupID[keep]
  xlim <- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim <- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  out <- list(y.P=y.P,s=s,X=X, K=K,y.I=y.I,si=si,y.I.obs=y.I.obs,si.obs=si.obs,NperGroup=NperGroup,
            groupID=groupID,groupID.obs=groupID.obs,buff=buff,Ni=Ni,
            xlim=xlim,ylim=ylim)
  return(out)
}