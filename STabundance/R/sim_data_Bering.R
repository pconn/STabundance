### sim_data_Bering.R
### function to simulate data for spatio-temporal abundance analysis from BOSS ice data
### in Bering Sea

sim_data_Bering<-function(sim.type){
  require(sp)
  require(rgeos)
  require(Matrix)
  require(hierarchicalDS)
  require(ggplot2)
  require(plyr)
  require(grid)
  source('c:/users/paul.conn/git/STabundance/STabundance/R/util_funcs.R')
  source('c:/users/paul.conn/git/STabundance/STabundance/R/sim_funcs.R')
  
  #source('c:/users/paul.conn/git/STabundance/STabundance/R/spat_funcs.R')
  
  #load spatial-temporal covariate & grid data
  load("BOSSst_2012data.Rdata")  #boss grid, ice data
  #limit to April 10 - May 8
  t.steps=29 
  Old.Grid=Data$Grid
  Data$Grid=vector('list',t.steps)
  for(it in 1:t.steps)Data$Grid[[it]]=Old.Grid[[it+3]]
  rm(Old.Grid)
  
  load("Effort2012_BOSSst.Rdata")  #load effort data indicating grid cells and times surveyed (data produced with format_effort.R)
  load("Knot_cell_distances.Rdata")
  
  Data$Effort=data.frame(matrix(0,nrow(Effort$Mapping),1))
  Data$Effort$Cell=Effort$Mapping[,1]
  Data$Effort$Time=Effort$Mapping[,2]
  Data$Effort$AreaSurveyed=Effort$Area.trans
  Data$Effort=Data$Effort[,2:4]
  Area.adjust=1-Data$Grid[[1]]@data[,"land_cover"]
  
  S=nrow(Data$Grid[[1]])
  
  #initialize abundance for a single species over the "Grid"
  set.seed(12345)

  for(it in 1:t.steps){
    Data$Grid[[it]]@data$sqrt_edge=sqrt(Data$Grid[[it]]@data[,"dist_edge"])
    Data$Grid[[it]]@data$ice2=(Data$Grid[[it]]@data[,"ice_conc"])^2
  }
  formula=~ice_conc+ice2+dist_edge+sqrt_edge+dist_shelf

  Sim.pars=list(Hab.init=c(0.7,20,-13,-.7,-.7,-2.5),Hab.evol=c(-1.5,20,-13,-.7,-.7,-2.5),tau.eta=15,kern.sd=3,tau.epsilon=100)
  if(sim.type=="RS2closed")Lambda=sim_RS2(S=S,Data=Data,Sim.pars=Sim.pars,hab.formula=formula,Area.adjust=Area.adjust)
    
  Sim.pars=list(Hab=c(0.7,20,-13,-.7,-.7,-2.5),lambda=280000,tau.eta=15,rho.ar1=0.5,sig.ar1=0.1,tau.epsilon=20)
  if(sim.type=="CPIF")Lambda=sim_CPIF(S=S,Data=Data,Sim.pars=Sim.pars,hab.formula=formula,Q.knot=K.data$Q.knot,K.cpif=K.data$K,Area.adjust=Area.adjust)
    
  Sim.data=sim_Ncounts(S=S,Data=Data,Lambda=Lambda)
  Data$Count.data=Sim.data$Count.data
  Sim.data$Data=Data
  
  #plot(Data$Grid[[1]])
  #plot(Data$Grid[[1]][Data$Effort[which(Data$Effort[,"Time"]==37),"Cell"],],add=TRUE,col="blue")
  
  #function for plotting abundance map
  #plot_N_map(1,matrix(diag(Iw),ncol=1),Grid=Data$Grid)
  #plot_N_map(1,matrix(Sim.data$N[,30],ncol=1),Grid=Data$Grid)
  return(Sim.data)
}




#final data format:
#Dat - same as previous
#  col 1- transect ID (now, separate for each cell/day combo)
#      2- photo obtained (0/1) ? 
#' 		 3- Observation type (integer - the max integer being 'unknown' if applicable) [NOTE: modeled as factor, but need to be input as integers to account for unknown species observations]
#' 		 x- Group size

