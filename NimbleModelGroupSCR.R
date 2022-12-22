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
  #group likelihoods
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #group size
    counts[i] ~ dztpois(lambda.P) #zero-truncated Poisson
    counts.zero[i] <- counts[i]-counts.detected[i] #number of undetected group members
    #site visitation - don't skip z[i]=0 here so y.P in place when turning z's on
    lamd.P[i,1:J] <- GetVisitRate(s=s[i,1:2],X=X[1:J,1:2],lam0.P=lam0.P,sigma.P=sigma.P,J=J)
    for(j in 1:J){
      for(k in 1:K){
        y.P[i,j,k] ~ dpois(lamd.P[i,j])
      }
    }
  }
  #Individual likelihoods
  #observed individuals
  for(l in 1:ni){#loop over individuals
    for(j in 1:J){
      for(k in 1:K){
        y.I[l,j,k] ~ dpois(lambda.I*y.P[groupID[l],j,k])
      }
    }
  }
  #unobserved individuals
  for(i in 1:M){#loop over groups
    #these data are all 0 capture histories, can skip z[i]=0 here
    y.I.unobs[i,1:J,1:K] ~ dIndNoDetect(y.P=y.P[i,1:J,1:K],counts.zero=counts.zero[i],
                                     lambda.I=lambda.I, z=z[i])
  }
  N.group <- sum(z[1:M])
  N.ind <- sum(counts[1:M]*z[1:M])
})# end model


GetVisitRate <- nimbleFunction(
  run = function(s = double(1), X=double(2), lam0.P=double(0), sigma.P=double(0), J=double(0)){ 
    returnType(double(1))
    d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
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


dztpois <- nimbleFunction(
  run = function(x=double(0), lambda.P = double(0), log = integer(0)) {
    returnType(double(0))  
    if(x==0){
      prob = 0
    }else{
      prob = dpois(x, lambda.P)/(1 - exp(-lambda.P))
    }
    if(log){
      return(log(prob))
    }else{
      return(prob)
    }
  }
)

#From Peter Dalgaard via Simon Bonnor
#https://simon.bonners.ca/bonner-lab/wpblog/?p=102
rztpois <- nimbleFunction(
  run = function(n=integer(0), lambda.P = double(0)) {
    returnType(double(0))
    if(n>1){
      print("rztpois only handles n=1")
    }else{
      tol=1e-10
      if(lambda.P < tol){
        x <- 1 #if lambda practically 0, all obs will be 1
      }else{
        tmp <- runif(n, dpois(0, lambda.P), 1)
        x <- qpois(tmp[1], lambda.P)
      }
      return(x)
    }
  }
)



#Required custom update for number of individuals/group
countSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    i <- control$i
    counts.detected <- control$counts.detected
    count.ups <- control$count.ups
  },
  run = function() {
    if(model$z[i]==0){ #if z is off, propose from prior
      counts.cand=rztpois(1,model$lambda.P[1])
      model$counts[i] <<- counts.cand
      model$calculate(calcNodes)
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    }else{ #if z is on, use Metropolis-Hastings
      #we'll propose to add or subtract 1 individual with equal probability.
      #if adding, we can always choose a individual to turn on
      #if subtracting, we need to consider that if we choose a detected individual
      #it cannot be turned off. This could be incorporated into the proposal probs
      #but I'm just rejecting the update if you select a detected individual to subtract.
      #You must account for this or it won't work.
      for(up in 1:count.ups){ #how many updates per iteration?
        #propose to add/subtract 1
        updown=rbinom(1,1,0.5) #p=0.5 is symmetric
        reject=FALSE #we reject if 1) proposed counts <1 (bc zero-truncation) or 
        #2) you select a detected individual
        if(updown==0){#subtract
          counts.cand=model$counts[i]-1
          if(model$counts[i]<1){
            reject=TRUE
          }else{
            tmp=rcat(1,rep(1/model$counts[i],model$counts[i])) #select one of the individuals in this group
            if(tmp<=counts.detected){ #is it one of the detected individuals?
              reject=TRUE #if so, we reject
            }
          }
        }else{#add
          counts.cand=model$counts[i]+1
        }
        if(!reject){
          model.lp.initial <- model$getLogProb(calcNodes)
          model$counts[i] <<- counts.cand
          model.lp.proposed <- model$calculate(calcNodes)#proposed logProb
          log_MH_ratio <- (model.lp.proposed) - (model.lp.initial)
          accept <- decide(log_MH_ratio)
          if(accept) {
            copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
          } else {
            copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
          }
        }
      }
    }
  },
  methods = list( reset = function () {} )
)


