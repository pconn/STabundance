# run_generic_sims.R
# script to run generic spatio-temporal count data simulations
source("./STabundance/R/sim_data_generic.R")
source("./STabundance/R/mcmc_STPC.R")
source("./STabundance/R/mcmc_CPIF.R")
source("./STabundance/R/mcmc_SST.R")

set.seed(12345)
n.sims=100 #number of simulations at each design point
Delta=c(-0.02,-0.02+0.04*c(1:(n.sims-1))/(n.sims-1)) #equally spaced from 2% decrease to 2% increase in habitat
S=400
t.steps=20
N.transects=c(1,5)
line.width=0.05
Model.list=c("RS2closed","CPIF","RS2open")
GENERATE=FALSE
if(GENERATE==TRUE){
  for(igen in 1:3){  #loop over generating model to generate data sets
    for(itrans in 1:2){ #loop over number of transects in each cell
      for(isim in 1:n.sims){
        Sim.data=sim_data_generic(sim.type=Model.list[igen],S=S,t.steps=t.steps,n.transects=N.transects[itrans],line.width=line.width,delta=Delta[isim])
        fname=paste("simdata_gen",Model.list[igen],"_trans",N.transects[itrans],"_sim",isim,sep='')
        save(Sim.data,file=paste("./sim_generic_data/",fname,sep=''))
      }
    }  
  }
}

Distances=gDistance(gCentroid(Data$Grid[[1]],byid=TRUE),gCentroid(Data$Grid[[1]],byid=TRUE),byid=TRUE)
Distances=Matrix(Distances)
#the following neighborhood uses all cells that intersect with a 75km radius circle surrounding a given grid cells centroid
Distances[which(Distances<0.01)]=9
Distances[which(Distances>.5 & Distances<1.1)]=8
Distances[which(Distances>1.4 & Distances<1.42)]=7
Distances[which(Distances>1.9 & Distances<2.1)]=6
Distances[which(Distances>2.1 & Distances<2.3)]=5
Distances[which(Distances>2.8 & Distances<2.9)]=4
Distances[which(Distances>2.9 & Distances<3.1)]=3
Distances[which(Distances>3.1 & Distances<3.2)]=2
Distances[which(Distances>3.6 & Distances<3.7)]=1
Distances[which(Distances>3.7)]=NA  #replace all distances >3.7 with NA (not part of my kernel)


#call estimation routines
Est.mods=c("RS1","RS2","CPIF","SST","STPC")
#for(igen in 1:3){  #loop over generating model to generate data sets
#  for(itrans in 1:2){ #loop over number of transects in each cell
#    for(isim in 1:n.sims){  
      igen=1
      itrans=2
      isim=52
      fname=paste("simdata_gen",Model.list[igen],"_trans",N.transects[itrans],"_sim",isim,sep='')
      load(paste("./sim_generic_data/",fname,sep=''))
      Data=Sim.data$Data
      Data$Which.distances=which(is.na(Distances)==FALSE)
      Data$Dist.entries=Distances[Data$Which.distances]  

      #set knot locations, etc. for process convolution models
      x.length=30
      y.length=30
      Which.knot=NA
      cur.pl=1
      for(i in 1:x.length){
        for(j in 1:y.length){
          if(i%%5==3 & j%%5==3)Which.knot=c(Which.knot,cur.pl)
          cur.pl=cur.pl+1
        }
      }
      require(rgeos)
      Knot.centers=gCentroid(Data$Grid[[1]][Which.knot[-1],],byid=TRUE)
      n.knots=length(Knot.centers)
      #calculate kernel densities at grid cell centroids 
      Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
      Distances=gDistance(Knot.centers,Cell.centroids,byid=TRUE)
      K=matrix(dnorm(Distances,0,5),S,n.knots)  #knot sd=5 
      K=K/rowSums(K)        
      
      Data$K=K
      for(iest in 1:5){ #loop over estimation model
        
        if(Est.mods[iest]=="RS1")MCMC=NULL
        if(Est.mods[iest]=="RS2")MCMC=NULL  
        
        Control=list(iter=5000,burnin=3000,thin=1000,predict=FALSE,MH.N=0.05,MH.omega=rep(0.04,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
        if(Est.mods[iest]=="CPIF")MCMC=mcmc_CPIF(model=~0+matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
        Control=list(iter=50100,burnin=100,thin=20,predict=TRUE,MH.N=MCMC$Control$MH.N,MH.omega=MCMC$Control$MH.omega,adapt=FALSE,fix.tau.epsilon=FALSE)        
        if(Est.mods[iest]=="CPIF")MCMC=mcmc_CPIF(model=~0+matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
        Control=list(iter=5000,burnin=100,thin=1000,srr.tol=0.5,predict=TRUE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
        if(Est.mods[iest]=="SST")MCMC=mcmc_AST(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
        Control=list(iter=25100,burnin=100,thin=10,srr.tol=0.5,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
        if(Est.mods[iest]=="SST")MCMC=mcmc_AST(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
        Control=list(iter=5000,burnin=100,thin=1000,srr.tol=0.5,predict=TRUE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
        if(Est.mods[iest]=="STPC")MCMC=mcmc_STPC(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
        Control=list(iter=100100,burnin=100,thin=100,srr.tol=0.5,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
        if(Est.mods[iest]=="STPC")MCMC=mcmc_STPC(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
        
        
        #Obs.data
        Obs.data=Data$Count.data[which(Data$Count.data[,"Count"]>0),]
        Cov=rep(0,nrow(Obs.data))
        for(i in 1:nrow(Obs.data))Cov[i]=Data$Grid[[Obs.data[i,"Time"]]]@data[Obs.data[i,"Cell"],"matern"]
        
        #calculate posterior for N
        MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
        N.true=apply(Sim.data$N,2,'sum')
        N.est=apply(MCMC$MCMC$N,1,'mean')
        plot(N.est)
        lines(N.true)
        
        plot(apply(MCMC$MCMC$Pred,3,'sum')/30)

        plot_N_map(1,as.matrix(Data$Grid[[1]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
        plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
        plot_N_map(1,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
        plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
        plot_N_map(15,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
        plot_N_map(30,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
        plot_N_map(30,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
        
      }
#    }  
#  }
#}

  
  