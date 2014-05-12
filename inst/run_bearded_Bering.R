# run_BOSS_ST_sims.R
# script to run generic spatio-temporal count data simulations
require(rgeos)
require(sp)
require(Matrix)
require(hierarchicalDS)
require(ggplot2)
require(plyr)
require(grid)

#load spatial-temporal covariate & grid data
source("./STabundance/R/sim_data_Bering.R")
source("./STabundance/R/sim_funcs.R")
source("./STabundance/R/util_funcs.R")
source("./STabundance/R/mcmc_STPC.R")
source("./STabundance/R/mcmc_CPIF.R")
source("./STabundance/R/mcmc_AST.R")

load("AlaskaBeringData2012_17April2014.Rdat")  #boss grid, ice data
load("Effort2012_BOSSst_bearded_29Apr2014.Rdata")  #load effort data indicating grid cells and times surveyed (data produced with format_effort.R)
load("Knot_cell_distances.Rdata") #load object giving K matrix, Q for knots


set.seed(123454)
t.steps=29
Old.Grid=Data$Grid
Data$Grid=vector('list',t.steps)
for(it in 1:t.steps){
  Data$Grid[[it]]=Old.Grid[[it+3]]
  Data$Grid[[it]]@data$ice2=Data$Grid[[it]]@data$ice_conc^2
  Data$Grid[[it]]@data$sqrt_edge=Data$Grid[[it]]@data$dist_edge^(0.5)
}
Data$Count.data=Effort$Count.data
rm(Old.Grid)

S=nrow(Data$Grid[[1]])

 
Effort$Area.hab[which(Effort$Area.hab==0)]=0.00001
n.knots=length(Data$Knot.locations)
Data$K=K.data$K
#Control=list(iter=8000,burnin=10,thin=10,predict=TRUE,MH.N=0.06,MH.omega=rep(0.07,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)                
Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.N=0.06,MH.omega=rep(0.07,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_CPIF(model=~0+ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab)            
Control=list(iter=80000,burnin=10000,thin=150,predict=TRUE,MH.N=MCMC$Control$MH.N,MH.omega=MCMC$Control$MH.omega,adapt=FALSE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_CPIF(model=~0+ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
out.file="./ST_out/spotted_Bering_CPIF3.Rdata"
save(MCMC,file=out.file)

Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.N=0.06,MH.omega=rep(0.07,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_CPIF(model=~0+ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab)            
Control=list(iter=310000,burnin=10000,thin=100,predict=TRUE,MH.N=MCMC$Control$MH.N,MH.omega=MCMC$Control$MH.omega,adapt=FALSE,fix.tau.epsilon=FALSE)        
MCMC=mcmc_CPIF(model=~0+ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
out.file="./ST_out/spotted_Bering_CPIF4.Rdata"
save(MCMC,file=out.file)

#Control=list(iter=8000,burnin=10,thin=10,predict=TRUE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)                
#Control=list(iter=5000,burnin=4000,thin=1000,predict=FALSE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
#MCMC=mcmc_AST(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
#Control=list(iter=30000,burnin=10000,thin=10,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
#MCMC=mcmc_AST(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 

#Control=list(iter=8000,burnin=10,thin=10,predict=TRUE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)                
#Control=list(iter=5000,burnin=4000,thin=1000,srr.tol=0.5,predict=FALSE,MH.mu=rep(0.2,nrow(Data$Count.data)),MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=TRUE,fix.tau.epsilon=FALSE)        
#MCMC=mcmc_STPC(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 
#Control=list(iter=410000,burnin=10000,thin=200,predict=TRUE,MH.mu=MCMC$Control$MH.mu,MH.N=0.05,MH.omega=rep(0.05,t.steps),adapt=FALSE,fix.tau.epsilon=FALSE)        
#MCMC=mcmc_STPC(model=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf,Prior.pars=NULL,Data=Data,Control=Control,Area.adjust=Effort$Area.hab) 

#out.file="./ST_out/spotted_Bering_STPC1.Rdata"
#save(MCMC,file=out.file)

#calculate posterior for N
MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
N.est=apply(MCMC$MCMC$N,1,'mean')
plot(N.est)

#crap=acf(MCMC$MCMC$N[18,200:799],lag.max=45)
#(1+2*sum(crap$acf)*var(MCMC$MCMC$N[18,200:799])*10000)/(mean(MCMC$MCMC$N[18,200:799]))^2

#plot(apply(MCMC$MCMC$Pred,3,'sum')/30)

#plot_N_map(1,as.matrix(Sim.data$Data$Grid[[5]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
#plot_N_map(1,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(29,apply(MCMC$MCMC$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(15,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(17,apply(MCMC$MCMC$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")
#plot_N_map(20,Sim.data$N,Grid=Data$Grid,leg.title="True Abundance")
#plot_N_map(29,apply(MCMC$MCMC$Pred,c(1,2),'mean'),Grid=Data$Grid,leg.title="Abundance")

