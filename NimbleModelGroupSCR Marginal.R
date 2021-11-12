NimModel <- nimbleCode({
  #group size prior
  lambda.P~dunif(0,100)
  #group site visitation priors
  lam0.P~dunif(0,10)
  sigma.P~dunif(0,100)
  #individual detection prior
  lambda.I~dunif(0,10)
  #group data augmentation prior
  psi~dunif(0,1)
  #group likelihoods (except for s priors)
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #group size
    counts[i] ~ T(dpois(lambda.P),1,Inf) #zero-truncated Poisson
    counts.zero[i] <- counts[i]-counts.detected[i] #number of undetected group members
    #site visitation expected values
    d2[i,1:J] <- Getd2(s=s[i,1:2],X=X[1:J,1:2],J=J,z=z[i])
    lamd.P[i,1:J] <- GetVisitRate(d2=d2[i,1:J],lam0.P=lam0.P,sigma.P=sigma.P,J=J,z=z[i])
    #Marginal observation likelihood for each group
    #Trickery here to make nimble do what we want
    #y.I.unobs are all 0 capture histories for undetected group members
    #y.I captured ind data fed in as an argument. Subsetting for group i inside function
    y.I.unobs[i,1:J,1:K]~dGroupVisitDetect(y.I[1:ni,1:J,1:K],this.group=i,groupID=groupID[1:ni],
                                          lambda.I=lambda.I,lamd.P=lamd.P[i,1:J],z=z[i],
                                          counts.zero=counts.zero[i],J=J,K=K,maxvisit=maxvisit)
  }
  N.group <- sum(z[1:M])
  N.ind <- sum(counts[1:M]*z[1:M])
})# end model

Getd2 <- nimbleFunction(
  run = function(s = double(1), X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0){
      return(rep(0,J))
    }else{
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      return(d2)
    }
  }
)
GetVisitRate <- nimbleFunction(
  run = function(d2=double(1),lam0.P=double(0), sigma.P=double(0), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0){
      return(rep(0,J))
    }else{
      lamd.P <- lam0.P*exp(-d2/(2*sigma.P^2))
    }
    return(lamd.P)
  }
)

dGroupVisitDetect <- nimbleFunction(
  run = function(x = double(2), y.I = double(3), this.group = double(0), groupID = double(1),
                 lambda.I = double(0), lamd.P = double(1),
                 z = double(0), counts.zero = double(0),J = double(0), K = double(0),
                 maxvisit = double(0), log = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }else{
      if(counts.zero<0){#keeps counts >= counts.detected
        return(-Inf)
      }else{
        group.idx=which(groupID==this.group)
        n.group=length(group.idx)
        ll.y=matrix(0,nrow=J,ncol=K)
        for(j in 1:J){
          for(k in 1:K){
            logProb.detect=rep(0,maxvisit+1)
            logProb.nodetect=rep(0,maxvisit+1)
            logProb.visit=rep(0,maxvisit+1)
            for(v in 0:maxvisit){
              #group site visit likelihood
              logProb.visit[v+1]=dpois(v,lamd.P[j],log=log)
              #undetected ind likelihood
              if(counts.zero>0){
                logProb.nodetect[v+1]=counts.zero*dpois(x[j,k],v*lambda.I,log=log)
              }
              #detected ind likelihoods
              if(n.group>0){#were any group members detected?
                for(i in 1:n.group){#add group member detection likelihoods
                  logProb.detect[v+1]=logProb.detect[v+1]+dpois(y.I[group.idx[i],j,k],v*lambda.I,log=log)
                }
              }
            }
            #total likelihood for site j, occasion k over all v
            logProb.total=logProb.visit+logProb.detect+logProb.nodetect
            maxlp=max(logProb.total)
            #marginal (over v) likelihood for site j, occasion k
            ll.y[j,k]=maxlp+log(sum(exp(logProb.total-maxlp)))
          }
        }
      }
      return(sum(ll.y))
    }
  }
)

#dummy RNG to make nimble happy, not used
rGroupVisitDetect <- nimbleFunction(
  run = function(n=integer(0), y.I = double(3), this.group = double(0), groupID = double(1),
                 lambda.I = double(0), lamd.P = double(1),
                 z = double(0), counts.zero = double(0),J = double(0), K = double(0),
                 maxvisit = double(0)) {
    returnType(double(2))
    return(matrix(0,nrow=J,ncol=K))
  }
)