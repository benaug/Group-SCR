NimModel <- nimbleCode({
  #group size prior
  lambda.P ~ dunif(0,100)
  #group site visitation priors
  p0.P ~ dunif(0,1)
  sigma.P ~ dunif(0,100)
  #individual detection prior
  p.I ~ dunif(0,1)
  #group data augmentation prior
  psi ~ dunif(0,1)
 
  #group likelihoods
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #group size
    counts[i] ~ dztpois(lambda.P) #zero-truncated Poisson
    counts.zero[i] <- counts[i]-counts.detected[i] #number of undetected group members
    #site visitation - can skip z[i]=0 when marginalizing over y.P
    pd.P[i,1:J] <- GetVisitProb(s=s[i,1:2],X=X[1:J,1:2],p0.P=p0.P,sigma.P=sigma.P,J=J,z=z[i])
    #Marginal observation likelihood for each group
    #Trickery here to make nimble do what we want
    #y.I.unobs are all 0 capture histories for undetected group members
    #y.I captured ind data fed in as an argument. Subsetting for group i inside function
    y.I.unobs[i,1:J,1:K] ~ dGroupVisitDetect(y.I[1:ni,1:J,1:K],this.group=i,groupID=groupID[1:ni],
                                             p.I=p.I,pd.P=pd.P[i,1:J],z=z[i],
                                             counts.zero=counts.zero[i],J=J,K=K)
  }
  N.group <- sum(z[1:M])
  N.ind <- sum(counts[1:M]*z[1:M])
})# end model

GetVisitProb <- nimbleFunction(
  run = function(s = double(1), X=double(2),p0.P=double(0), sigma.P=double(0), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0){
      return(rep(0,J))
    }else{
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      pd.P <- p0.P*exp(-d2/(2*sigma.P^2))
    }
    return(pd.P)
  }
)

#this version much faster because for each group and trap, we only need to compute the likelihood
#for occasions with 0 detections once, because it is the same for all 0 detection occasions.
dGroupVisitDetect <- nimbleFunction(
  run = function(x = double(2), y.I = double(3), this.group = double(0), groupID = double(1),
                 p.I = double(0), pd.P = double(1),
                 z = double(0), counts.zero = double(0),J = double(0), K = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }else{
      if(counts.zero<0){#keeps counts >= counts.detected
        return(-Inf)
      }else{
        group.idx <- which(groupID==this.group)
        n.group <- length(group.idx)
        lp.y <- matrix(0,nrow=J,ncol=K)
        #undetected ind likelihood - only need to calculate once, does not depend on j,k
        logProb.nodetect <- rep(0,2) #logProb is 0 if no site visit, only fill in 2nd element for visited
        if(counts.zero>0){
          logProb.nodetect[2] <- counts.zero*dbinom(0,size=1,prob=p.I,log=log)
        }
        for(j in 1:J){
          fill.zeros <- FALSE#have we encountered a zero detection k for this trap, yet?
          for(k in 1:K){
            #count detections for this group at this j-k
            if(n.group>0){
              detects <- 0
              for(i in 1:n.group){
                detects <- detects + y.I[group.idx[i],j,k]
              }
            }
            if(!fill.zeros|detects>0){ #compute if we have detections, or if there are 0 detections and we haven't computed that, yet
              logProb.detect <- rep(0,2)
              logProb.visit <- rep(0,2)
              for(v in 0:1){
                #group site visit likelihood
                logProb.visit[v+1] <- dbinom(v,prob=pd.P[j],size=1,log=log)
                #detected ind likelihoods
                if(n.group>0){#were any group members detected?
                  for(i in 1:n.group){#add group member detection likelihoods
                    logProb.detect[v+1] <- logProb.detect[v+1]+dbinom(y.I[group.idx[i],j,k],size=1,prob=v*p.I,log=log)
                  }
                }
              }
              #total likelihood for site j, occasion k over all v
              logProb.total <- logProb.visit+logProb.detect+logProb.nodetect
              maxlp <- max(logProb.total)
              #marginal (over v) likelihood for site j, occasion k
              lp.y[j,k] <- maxlp + log(sum(exp(logProb.total-maxlp)))
              if(detects==0){
                fill.zeros <- TRUE
                lp.zero <- lp.y[j,k]
              }
            }else{#0 detection occasion likelihood already computed, fill in for the remainder
              lp.y[j,k] <- lp.zero
            }
          }
        }
      }
      return(sum(lp.y))
    }
  }
)

#this one computes every i x j x k, much slower
# dGroupVisitDetect <- nimbleFunction(
#   run = function(x = double(2), y.I = double(3), this.group = double(0), groupID = double(1),
#                  p.I = double(0), pd.P = double(1),
#                  z = double(0), counts.zero = double(0),J = double(0), K = double(0),
#                  log = integer(0)) {
#     returnType(double(0))
#     if(z==0){
#       return(0)
#     }else{
#       if(counts.zero<0){#keeps counts >= counts.detected
#         return(-Inf)
#       }else{
#         group.idx <- which(groupID==this.group)
#         n.group <- length(group.idx)
#         lp.y <- matrix(0,nrow=J,ncol=K)
#         #undetected ind likelihood - only need to calculate once, does not depend on j,k
#         logProb.nodetect <- rep(0,2) #logProb is 0 if no site visit, only fill in 2nd element for visited
#         if(counts.zero>0){
#           logProb.nodetect[2] <- counts.zero*dbinom(0,size=1,prob=p.I,log=log)
#         }
#         for(j in 1:J){
#           for(k in 1:K){
#             logProb.detect <- rep(0,2)
#             # logProb.nodetect <- rep(0,2)
#             logProb.visit <- rep(0,2)
#             for(v in 0:1){
#               #group site visit likelihood
#               logProb.visit[v+1] <- dbinom(v,prob=pd.P[j],size=1,log=log)
#               #detected ind likelihoods
#               if(n.group>0){#were any group members detected?
#                 for(i in 1:n.group){#add group member detection likelihoods
#                   logProb.detect[v+1] <- logProb.detect[v+1]+dbinom(y.I[group.idx[i],j,k],size=1,prob=v*p.I,log=log)
#                 }
#               }
#             }
#             #total likelihood for site j, occasion k over all v
#             logProb.total <- logProb.visit+logProb.detect+logProb.nodetect
#             maxlp <- max(logProb.total)
#             #marginal (over v) likelihood for site j, occasion k
#             lp.y[j,k] <- maxlp + log(sum(exp(logProb.total-maxlp)))
#           }
#         }
#       }
#       return(sum(lp.y))
#     }
#   }
# )

#dummy RNG to make nimble happy, not used
rGroupVisitDetect <- nimbleFunction(
  run = function(n=integer(0),  y.I = double(3), this.group = double(0), groupID = double(1),
                 p.I = double(0), pd.P = double(1),
                 z = double(0), counts.zero = double(0),J = double(0), K = double(0)) {
    returnType(double(2))
    return(matrix(0,J,K))
  }
)


dztpois <- nimbleFunction(
  run = function(x=double(0), lambda.P = double(0), log = integer(0)) {
    returnType(double(0))  
    if(x==0){
      prob <- 0
    }else{
      prob <- dpois(x, lambda.P)/(1 - exp(-lambda.P))
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
      tol <- 1e-10
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
      counts.cand <- rztpois(1,model$lambda.P[1])
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
        updown <- rbinom(1,1,0.5) #p=0.5 is symmetric
        reject <- FALSE #we reject if 1) proposed counts <1 (bc zero-truncation) or 
        #2) you select a detected individual
        if(updown==0){#subtract
          counts.cand <- model$counts[i]-1
          if(model$counts[i]<1){
            reject <- TRUE
          }else{
            tmp <- rcat(1,rep(1/model$counts[i],model$counts[i])) #select one of the individuals in this group
            if(tmp<=counts.detected){ #is it one of the detected individuals?
              reject <- TRUE #if so, we reject
            }
          }
        }else{#add
          counts.cand <- model$counts[i]+1
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