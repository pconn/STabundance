#simulate a small test dataset for debugging
source("./STabundance/R/sim_data_generic.R")
source("./STabundance/R/sim_data_no_error.R")
source("./STabundance/R/mcmc_CPIF.R")
source("./STabundance/R/mcmc_STPC.R")
source("./STabundance/R/mcmc_OPRS.R")
source("./STabundance/R/mcmc_AST.R")

set.seed(12345)
S=400
delta=0
t.steps=20
n.transects=1
line.width=0.05
sim.type="RS2closed"

Sim.data=sim_data_generic(sim.type=sim.type,S=S,t.steps=t.steps,n.transects=n.transects,line.width=line.width,delta=delta)
#Sim.data=sim_data_no_error(sim.type=sim.type,S,t.steps,n.transects,line.width,delta=delta)
#fname="sim_data_test"
#save(Sim.data,file=paste("./sim_generic_data/",fname,sep=''))



#now do some testing

Data=Sim.data$Data

#set knot locations, etc. for process convolution models
x.length=20
y.length=20
Which.knot=NA
cur.pl=1
for(i in 1:x.length){
  for(j in 1:y.length){
    if(i%%4==1 & j%%4==1)Which.knot=c(Which.knot,cur.pl)
    cur.pl=cur.pl+1
  }
}
require(rgeos)
Knot.centers=gCentroid(Data$Grid[[1]][Which.knot[-1],],byid=TRUE)
n.knots=length(Knot.centers)
#calculate kernel densities at grid cell centroids 
Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
Distances=gDistance(Knot.centers,Cell.centroids,byid=TRUE)
K=matrix(dnorm(Distances,0,4),S,n.knots)  #knot sd=4 
K=K/rowSums(K)        
Data$K=K

Distances=gDistance(gCentroid(Data$Grid[[1]],byid=TRUE),gCentroid(Data$Grid[[1]],byid=TRUE),byid=TRUE)
Distances=Matrix(Distances)
#the following neighborhood uses all cells that intersect with a 75km radius circle surrounding a given grid cells centroid
Distances[which(Distances>3.7)]=NA  #replace all distances >3.7 with NA (not part of my kernel)
Distances[which(Distances<0.01)]=9
Distances[which(Distances>.9 & Distances<1.1)]=8
Distances[which(Distances>1.4 & Distances<1.42)]=7
Distances[which(Distances>1.9 & Distances<2.1)]=6
Distances[which(Distances>2.1 & Distances<2.3)]=5
Distances[which(Distances>2.8 & Distances<2.9)]=4
Distances[which(Distances>2.9 & Distances<3.1)]=3
Distances[which(Distances>3.1 & Distances<3.2)]=2
Distances[which(Distances>3.6 & Distances<3.7)]=1
Data$Which.distances=which(is.na(Distances)==FALSE)
Data$Dist.entries=Distances[Data$Which.distances]  



Control=list(iter=10010,burnin=10,thin=1,srr.tol=0.5,predict=TRUE,MH.Z=rep(0.08,t.steps),MH.mu=rep(0.1,S),MH.N=0.05,MH.omega=rep(0.05,t.steps),MH.beta=c(.05,0.1,0.2),MH.tau.d=0.5,adapt=TRUE,fix.tau.epsilon=FALSE,fix.tau.varepsilon=TRUE)
model=~0+matern+matern2
Predictors=c("matern","matern2")
Prior.pars=NULL
Area.adjust=NULL

MCMC=mcmc_OPRS(Predictors=Predictors,Data=Data,Control=Control)

#if(Est.mods[iest]=="RS2")MCMC=NULL   



Control=list(iter=10100,burnin=4100,thin=1,predict=FALSE,MH.N=0.01,MH.omega=rep(0.001,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_CPIF(model=model,Data=Data,Control=Control)
Control=list(iter=40100,burnin=100,thin=20,predict=TRUE,MH.N=MCMC$Control$MH.N,MH.omega=MCMC$Control$MH.omega,adapt=FALSE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_CPIF(model=model,Data=Data,Control=Control)

#if(Est.mods[iest]=="SST")
MCMC=mcmc_AST(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
#if(Est.mods[iest]=="STPC")MCMC=mcmc_STPC(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
Control=list(iter=10100,burnin=100,thin=1,srr.tol=0.5,predict=TRUE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_STPC(model=~~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
Control=list(iter=100100,burnin=100,thin=100,srr.tol=0.5,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_STPC(model=~1,Prior.pars=NULL,Data=Data,Control=Control)


#Obs.data
Obs.data=Data$Count.data[which(Data$Count.data[,"Count"]>0),]
Cov=rep(0,nrow(Obs.data))
for(i in 1:nrow(Obs.data))Cov[i]=Data$Grid[[Obs.data[i,"Time"]]]@data[Obs.data[i,"Cell"],"matern"]

#calculate posterior for N
MCMC$MCMC$N.est=apply(MCMC$MCMC$Pred,c(2,3),'sum')
N.true=apply(Sim.data$N,2,'sum')
N.est=apply(MCMC$MCMC$N.est,1,'mean')
plot(N.est)
lines(N.true)

plot_N_map(1,as.matrix(Data$Grid[[1]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
plot_N_map(1,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
plot_N_map(15,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
plot_N_map(20,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
plot_N_map(20,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")

for(cur.yr in 1:20){
  cur.acf=acf(MCMC$MCMC$N.est[cur.yr,],lag.max=2000)
  m = Control$iter/(1+2*sum(cur.acf$acf)) #effective sample size 
  n.cv01=(1+2*sum(cur.acf$acf))*var(MCMC$MCMC$N.est[cur.yr,])*10000/mean(MCMC$MCMC$N.est[cur.yr,])^2  
  print(n.cv01)
}

#for CPIF
cur.acf=acf(MCMC$MCMC$N,lag.max=1500)
m = Control$iter/(1+2*sum(cur.acf$acf)) #effective sample size 
n.cv01=(1+2*sum(cur.acf$acf))*var(MCMC$MCMC$N)*10000/mean(MCMC$MCMC$N)^2  
print(n.cv01)
