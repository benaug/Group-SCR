#This sampler marginalizes out site visitation. 
#Mixes better than full CDL sampler for site visitation and detection parameters
#This version probably most efficient.

source("sim.GroupSCR.R")
source("init.group.R")
library(nimble)
source("NimbleModelGroupSCR Marginal.R")
Np=15 #number of groups
lambda.P=5 #group size parameter (zero-truncated Poisson lambda)
lam0.P=2 #expected # group site visits at activity center
sigma.P=0.75 #group site visit spatial scale parameter
lambda.I=0.25 #individual detection rate|group site visitation
K=2 #number of occasions. Sampler needs to be modified if only 1 occasion, reducing data dimensions from 3D to 2D.
X=expand.grid(1:10,1:10) #trapping array
buff=2.5 #state space buffer

data=sim.GroupSCR(Np=Np,lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I,K=K,X=X,buff=buff)

par(mfrow=c(1,1),ask=FALSE)
xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
plot(NA,pch=16,xlim=xlim,ylim=ylim,xlab="X",ylab="Y")
visiti=which(rowSums(data$y.P)>0)
y.P2D=apply(data$y.P,c(1,2),sum)
for(i in visiti){
  trapcap=which(y.P2D[i,]>0)
  for(j in trapcap){
    lines(x=c(data$s[i,1],X[j,1]),y=c(data$s[i,2],X[j,2]),lty=3)
  }
}
capi=which(rowSums(data$y.I)>0)
y.I2D=apply(data$y.I,c(1,2),sum)
for(i in capi){
  trapcap=which(y.I2D[i,]>0)
  for(j in trapcap){
    ID=data$groupID[i]
    lines(x=c(data$s[ID,1],X[j,1]),y=c(data$s[ID,2],X[j,2]),lty=1,col="#999999",lwd=2)
  }
}
points(data$s,cex=2,col="#56B4E9",pch=16)
points(data$s,pch=as.character(data$NperGroup),cex=0.75)
points(X,pch=4)


####Fit model
M=50 #group data augmentation level
#to integrate out the site visits, we need to supply a max to sum over.
#this should be large enough that the probability of site visits (per occasion) higher than
#maxvist is negligible. Be careful!
maxvisit=12

#If you know lam0.P, a good value would be
qpois(0.99,lam0.P) #99% of site visits are this number or smaller given lam0.P
#But we need to account for the fact that the lam0.P posterior spans above the point estimate. 
#On the other hand, lam0.P is expected # site visits *at the activity center*, and most ACs
#are not exactly on top of traps. Maybe set maxvisits to 2x expected lam0.P? 
qpois(0.99,2*lam0.P)

#initialize some data objects/latent variables
inits=list(lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I)
nimbuild=init.group(data=data,M=M,inits=inits)

#inits for nimble.
Niminits <- list(z=nimbuild$z,s=nimbuild$s,counts=nimbuild$counts,
                 #Starting at true parameter values to speed convergence. Don't do this in practice.
                 lambda.P=lambda.P,lam0.P=lam0.P,sigma.P=sigma.P,lambda.I=lambda.I)

#constants for Nimble
J=nrow(data$X)
ni=nrow(data$y.I.obs) #number of captured individuals
#supplying y.I data as constants to trick nimble. Must be in constants not data!
constants<-list(y.I=data$y.I.obs,M=M,J=J,K=K,xlim=data$xlim,ylim=data$ylim,ni=nrow(data$y.I.obs),ni=ni,
                counts.detected=nimbuild$counts.detected,groupID=nimbuild$groupID,maxvisit=maxvisit)

#Supply data to Nimble.
z.data=nimbuild$z
z.data[z.data==0]=NA #unobserved group z's are unknown, must be estimated
#y.I.unobs are 0 capture histories for each group to account for uncaptured individuals
Nimdata<-list(y.I.unobs=array(0,dim=c(M,J,K)),counts=rep(NA,M),X=data$X,z=z.data)

# set parameters to monitor
parameters<-c('psi',"lambda.P","lam0.P","sigma.P","lambda.I","N.group","N.ind","psi")


start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1, useConjugacy = TRUE)
# conf$removeSampler(paste("y.P[1:",M,",1:",J,",1:",K,"]", sep=""))
conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
for(i in 1:M){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE),silent = TRUE)
}

#AF slice really improves mixing, but slower than independent RW samplers. Probably a good trade-off.
#But does improved mixing for these params propagate to N.ind and N.group?
#RW_block may be efficient, but takes longer to converge than independent RW samplers.
conf$removeSampler(c("lambda.I","lam0.P","sigma.P"))
conf$addSampler(target = c("lambda.I","lam0.P","sigma.P"),
                type = 'AF_slice',
                control = list(adaptive=TRUE),
                silent = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model
start.time2<-Sys.time()
Cmcmc$run(500,reset=FALSE) #short run for demonstration. Can run again to continue sampling where it stopped.
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples=as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))
