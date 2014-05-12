# run_generic_sims.R
# script to run generic spatio-temporal count data simulations
require(STabundance)

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




#call estimation routines
Est.mods=c("CPIF","AST","STPC","OPRS")
Sim.ind.20= c(1:19)*5
for(igen in 1:3){  #loop over generating model to generate data sets
  for(itrans in 1:1){ #loop over number of transects in each cell
    for(isim in 1:(n.sims/2)){  

      if(igen==3){
        fname=paste("simdata_gen",Model.list[igen-1],"_trans",N.transects[itrans],"_sim",isim,sep='')
        load(paste("./sim_generic_data/",fname,sep=''))
        Which.distances=Sim.data$Data$Which.distances  
      }
      
      fname=paste("simdata_gen",Model.list[igen],"_trans",N.transects[itrans],"_sim",isim,sep='')
      load(paste("./sim_generic_data/",fname,sep=''))
      Data=Sim.data$Data
      if(igen==3)Data$Which.distances=Which.distances
      
      n.knots=length(Data$Knot.locations)
      #calculate kernel densities at grid cell centroids 
      Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
      Distances=gDistance(Data$Knot.locations,Cell.centroids,byid=TRUE)
      K=matrix(dnorm(Distances,0,5),S,n.knots)  #knot sd=5 
      K=K/rowSums(K)          
      Data$K=K
      
      #Data$Count.data=Data$Count.data[-which(Data$Count.data[,"Time"]%in%c(3,4,10)),]
        
      for(iest in 1:1){ #loop over estimation model
        
        run.flag=0
        if(Est.mods[iest]=="CPIF"){
          if(itrans==1 & isim%in%Sim.ind.10){
            #Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.N=0.05,MH.omega=rep(0.04,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
            Control=list(iter=5000,burnin=10,thin=10,predict=TRUE,MH.N=0.05,MH.omega=rep(0.04,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_CPIF(model=~0+matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)            
            Control=list(iter=110000,burnin=10000,thin=50,predict=TRUE,MH.N=MCMC$Control$MH.N,MH.omega=MCMC$Control$MH.omega,adapt=FALSE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_CPIF(model=~0+matern+matern2,Prior.pars=NULL,Data=Data,Control=Control) 
            run.flag=1
          } 
          if(itrans==2){
            Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.N=0.05,MH.omega=rep(0.04,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_CPIF(model=~0+matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)            
            Control=list(iter=20000,burnin=10000,thin=5,predict=TRUE,MH.N=MCMC$Control$MH.N,MH.omega=MCMC$Control$MH.omega,adapt=FALSE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_CPIF(model=~0+matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)  
            run.flag=1
          }
        } 
        
        if(Est.mods[iest]=="AST"){
          if(itrans==1){
            Control=list(iter=5000,burnin=10,thin=10,predict=TRUE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_AST(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
            Control=list(iter=30000,burnin=10000,thin=10,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_AST(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)  
            run.flag=1
          }
          if(itrans==2){
            Control=list(iter=5000,burnin=100,thin=1000,predict=FALSE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_AST(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
            Control=list(iter=20000,burnin=10000,thin=5,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_AST(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control) 
            run.flag=1
          }
        }        
        
        
        if(Est.mods[iest]=="STPC"){
          if(itrans==1 & isim%in%Sim.ind.10){
            Control=list(iter=5000,burnin=10,thin=10,srr.tol=0.5,predict=TRUE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_STPC(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
            Control=list(iter=110000,burnin=10000,thin=50,srr.tol=0.5,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_STPC(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)   
            run.flag=1
          }
          if(itrans==2){
            Control=list(iter=5000,burnin=100,thin=1000,srr.tol=0.5,predict=FALSE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_STPC(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)
            Control=list(iter=20000,burnin=10000,thin=50,srr.tol=0.5,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
            MCMC=mcmc_STPC(model=~matern+matern2,Prior.pars=NULL,Data=Data,Control=Control)  
            run.flag=1
          }          
        }

        if(Est.mods[iest]=="OPRS" & itrans==2 & isim%in%Sim.ind.10){
          Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.Z=rep(0.08,t.steps),MH.mu=rep(0.1,S),MH.N=0.05,MH.omega=rep(0.05,t.steps),MH.beta=c(.05,0.1,0.2),MH.tau.d=0.5,adapt=TRUE,fix.tau.epsilon=FALSE,fix.tau.varepsilon=TRUE)
          model=~0+matern+matern2
          Predictors=c("matern","matern2")
          MCMC=mcmc_OPRS(Predictors=Predictors,Data=Data,Control=Control)
          Control=list(iter=20000,burnin=10000,thin=5,predict=TRUE,MH.Z=MCMC$Control$MH.Z,MH.mu=rep(0.1,S),MH.N=0.05,MH.omega=rep(0.05,t.steps),MH.beta=MCMC$Control$MH.beta,MH.tau.d=MCMC$Control$MH.tau.d,adapt=TRUE,fix.tau.epsilon=FALSE,fix.tau.varepsilon=TRUE)
          model=~0+matern+matern2
          Predictors=c("matern","matern2")
          MCMC=mcmc_OPRS(Predictors=Predictors,Data=Data,Control=Control) 
          run.flag=1
        }  
        
        out.file=paste("d:/ST_out/ST_out_gen",igen,"_est_",Est.mods[iest],"_trans_",itrans,"_sim",isim,".Rdata",sep="")
        if(run.flag==1)save(MCMC,file=out.file)
 
        #save results of MCMC run as .Rdata
        
        
        #Obs.data
        #Obs.data=Data$Count.data[which(Data$Count.data[,"Count"]>0),]
        #Cov=rep(0,nrow(Obs.data))
        #for(i in 1:nrow(Obs.data))Cov[i]=Data$Grid[[Obs.data[i,"Time"]]]@data[Obs.data[i,"Cell"],"matern"]
        
        #calculate posterior for N
        MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
        N.true=apply(Sim.data$N,2,'sum')
        N.est=apply(MCMC$MCMC$N,1,'mean')
        plot(N.est)
        lines(N.true)
        
        #plot(apply(MCMC$MCMC$Pred,3,'sum')/30)

        #plot_N_map(1,as.matrix(Sim.data$Data$Grid[[5]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
        #plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
        #plot_N_map(1,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
        #plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
        #plot_N_map(15,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
        #plot_N_map(20,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
        #plot_N_map(30,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
        
        
        
      }
    }  
  }
}

  
  