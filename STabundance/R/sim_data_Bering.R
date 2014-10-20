#' function to simulate Bering Sea spotted seal data for spatio-temporal abundance analysis 
#' @param sim.type A character string specifying the type of model used for space-time abundance dynamics 
#'        ("RS2closed" - resource selection on absolute abundance intensity on a closed population; "CPIF" closed population ideal free;
#'        ("RS2open" - open population with restricted dis))
#' @return A list object composed of at least two objects, "Data$Grid" is a list vector holding covariate data by day,
#'         and "Count.data" which holds simulated transect counts 
#' @export
#' @import Matrix
#' @keywords spatio-temporal, simulation, spotted seals
#' @author Paul B. Conn
sim_data_Bering<-function(sim.type){
  #source('c:/users/paul.conn/git/STabundance/STabundance/R/util_funcs.R')
  #source('c:/users/paul.conn/git/STabundance/STabundance/R/sim_funcs.R')
  
  #load spatial-temporal covariate & grid data
  data(AlaskaBeringData2012_17April2014)  #boss grid, ice data
  #limit to April 10 - May 8
  t.steps=29 
  Old.Grid=Data$Grid
  Data$Grid=vector('list',t.steps)
  for(it in 1:t.steps)Data$Grid[[it]]=Old.Grid[[it+3]]
  rm(Old.Grid)
  
  data(Effort2012_BOSSst)  #load effort data indicating grid cells and times surveyed (data produced with format_effort.R)
  data(Knot_cell_distances)
  
  Data$Effort=data.frame(matrix(0,nrow(Effort$Mapping),1))
  Data$Effort$Cell=Effort$Mapping[,1]
  Data$Effort$Time=Effort$Mapping[,2]
  Data$Effort$AreaSurveyed=Effort$Area.trans
  Data$Effort=Data$Effort[,2:4]
  Area.adjust=1-Data$Grid[[1]]@data[,"land_cover"]
  
  S=nrow(Data$Grid[[1]])

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

