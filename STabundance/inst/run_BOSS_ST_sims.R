# run_BOSS_ST_sims.R
# script to run generic spatio-temporal count data simulations
require(STabundance)

data(BOSSst_2012data.Rdata)  #boss grid, ice data
data(Effort2012_BOSSst.Rdata)  #load effort data indicating grid cells and times surveyed (data produced with format_effort.R)
data(Knot_cell_distances.Rdata) #load object giving K matrix, Q for knots

Model.list=c("RS2closed","CPIF")
set.seed(12345)
n.sims=10 #number of simulations at each design point

#limit gridded covariate data to April 10 - May 8
t.steps=29
Old.Grid=Data$Grid
Data$Grid=vector('list',t.steps)
for(it in 1:t.steps)Data$Grid[[it]]=Old.Grid[[it+3]]
rm(Old.Grid)

S=nrow(Data$Grid[[1]])

GENERATE=FALSE
if(GENERATE==TRUE){
  for(igen in 1:2){  #loop over generating model to generate data sets
      for(isim in 1:n.sims){
        Sim.data=sim_data_Bering(sim.type=Model.list[igen])
        fname=paste("simBering_gen",Model.list[igen],"_sim",isim,sep='')
        save(Sim.data,file=paste("./sim_BOSS_data/",fname,sep=''))
      }
  }  
}




#call estimation routines
Est.mods=c("CPIF","AST","STPC")
for(igen in 1:2){  #loop over generating model to generate data sets
  for(isim in 1:n.sims){  
    
    fname=paste("simBering_gen",Model.list[igen],"_sim",isim,sep='')
    load(paste("./sim_BOSS_data/",fname,sep=''))
    Data=Sim.data$Data
    if(sum(Effort$Area.hab==0)>0)Effort$Area.hab[which(Effort$Area.hab==0)]=0.0001
     
    n.knots=length(Data$Knot.locations)
    #calculate kernel densities at grid cell centroids 
    #Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
    #Distances=gDistance(Data$Knot.locations,Cell.centroids,byid=TRUE)
     
    Data$K=K.data$K
    for(iest in 1:3){ #loop over estimation model
      
      run.flag=0
      if(Est.mods[iest]=="CPIF" & (isim%%10)==1){
          #Control=list(iter=8000,burnin=10,thin=10,predict=TRUE,MH.N=0.06,MH.omega=rep(0.07,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)                
          Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.N=0.06,MH.omega=rep(0.07,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
          MCMC=mcmc_CPIF(model=~0+ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab)            
          Control=list(iter=60000,burnin=10000,thin=25,predict=TRUE,MH.N=MCMC$Control$MH.N,MH.omega=MCMC$Control$MH.omega,adapt=FALSE,fix.tau.epsilon=FALSE)        
          MCMC=mcmc_CPIF(model=~0+ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
          run.flag=1
      } 
      
      if(Est.mods[iest]=="AST"){
          #Control=list(iter=20000,burnin=10,thin=10,predict=TRUE,MH.mu=rep(1.5,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
          Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
          MCMC=mcmc_AST(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
          Control=list(iter=30000,burnin=10000,thin=10,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
          MCMC=mcmc_AST(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
          run.flag=1
      }        
      
      
      if(Est.mods[iest]=="STPC" & (isim%%10)==1){
          #Control=list(iter=10000,burnin=10,thin=10,predict=TRUE,MH.mu=rep(1.5,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
          Control=list(iter=5000,burnin=4000,thin=1000,srr.tol=0.5,predict=FALSE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
          MCMC=mcmc_STPC(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
          Control=list(iter=160000,burnin=10000,thin=75,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
          MCMC=mcmc_STPC(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
          run.flag=1
      }
      
      
      out.file=paste("d:/ST_out/ST_BOSS_out_gen",igen,"_est_",Est.mods[iest],"_sim",isim,".Rdata",sep="")
      if(run.flag==1)save(MCMC,file=out.file)
      
      #save results of MCMC run as .Rdata
      
      
      #Obs.data
      #Obs.data=Data$Count.data[which(Data$Count.data[,"Count"]>0),]
      #Cov=rep(0,nrow(Obs.data))
      #for(i in 1:nrow(Obs.data))Cov[i]=Data$Grid[[Obs.data[i,"Time"]]]@data[Obs.data[i,"Cell"],"matern"]
      
      #calculate posterior for N
      #MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
      #N.true=apply(Sim.data$N,2,'sum')
      #N.est=apply(MCMC$MCMC$N,1,'mean')
      #plot(N.est)
      #lines(N.true)
      
      #crap=acf(MCMC$MCMC$N[1,200:999],lag.max=140)
      #(1+2*sum(crap$acf)*var(MCMC$MCMC$N[1,200:999])*10000)/(mean(MCMC$MCMC$N[1,200:999]))^2
      
      #plot(apply(MCMC$MCMC$Pred,3,'sum')/30)
      
      #plot_N_map(1,as.matrix(Sim.data$Data$Grid[[5]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
      #plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
      #plot_N_map(1,apply(MCMC$MCMC$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")
      #plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
      #plot_N_map(15,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
      #plot_N_map(20,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
      #plot_N_map(30,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")
      
    }  
  }
}

  
  