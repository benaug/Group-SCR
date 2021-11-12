#could do some more vectorization, skip some z[i]==0 calculations...

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
    #site visitation
    d2[i,1:J] <- Getd2(s=s[i,1:2],X=X[1:J,1:2],J=J)
    lamd.P[i,1:J] <- GetVisitRate(d2=d2[i,1:J],lam0.P=lam0.P,sigma.P=sigma.P,J=J)
    for(j in 1:J){
      for(k in 1:K){
        y.P[i,j,k]~dpois(lamd.P[i,j])
      }
    }
  }
  #Individual likelihoods
  #observed individuals
  for(l in 1:ni){#loop over individuals
    for(j in 1:J){
      for(k in 1:K){
        y.I[l,j,k]~dpois(lambda.I*y.P[groupID[l],j,k])
      }
    }
  }
  #unobserved individuals
  for(i in 1:M){#loop over groups
    #these data are all 0 capture histories
    y.I.unobs[i,1:J,1:K]~dIndNoDetect(y.P=y.P[i,1:J,1:K],z=z[i],counts.zero=counts.zero[i],
                                     lambda.I=lambda.I)
  }
  N.group <- sum(z[1:M])
  N.ind <- sum(counts[1:M]*z[1:M])
})# end model

Getd2 <- nimbleFunction(
  run = function(s = double(1), X=double(2), J=double(0)){ 
    returnType(double(1))
    d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
    return(d2)
  }
)

GetVisitRate <- nimbleFunction(
  run = function(d2=double(1),lam0.P=double(0), sigma.P=double(0), J=double(0)){ 
    returnType(double(1))
    lamd.P <- lam0.P*exp(-d2/(2*sigma.P^2))
    return(lamd.P)
  }
)

dIndNoDetect <- nimbleFunction(
  run = function(x = double(2), y.P = double(2), z = double(0), counts.zero = double(0),
                 lambda.I = double(0), log = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }else{
      if(counts.zero<0){#keeps counts >= counts.detected
        return(-Inf)
      }else if(counts.zero==0){ #0 uncaptured individuals in this group
        return(0)
      }else{ #counts.zero uncaptured individuals in this group
        logProb <- counts.zero*sum(dpois(x, lambda = lambda.I*y.P, log = TRUE))
        return(logProb)
      }
    }
  }
)

#dummy RNG to make nimble happy, not used
rIndNoDetect <- nimbleFunction(
  run = function(n=integer(0), y.P = double(2), z = double(0), counts.zero = double(0),
                 lambda.I = double(0)) {
    returnType(double(2))
    J=nimDim(y.P)[1]
    K=nimDim(y.P)[2]
    return(matrix(0,nrow=J,ncol=K))
  }
)