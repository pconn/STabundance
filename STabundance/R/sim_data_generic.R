# function to simulate spatio-temporal count data over generic simulated landscapes
#' @param sim.type A character string specifying the type of model used for space-time abundance dynamics 
#'        ("RS2closed" - resource selection on absolute abundance intensity on a closed population; "CPIF" closed population ideal free;
#'        ("RS2open" - open population with restricted dis))
#' @param S Number of cells (must be a square number)
#' @param t.steps Number of time steps
#' @param n.transects Number of transects to simulate at each time step
#' @param line.width Proportional of a cell's diameter that is covered by a transect when it surveys a cell
#' @param delta Expected proportion increase/decrease for habitat covariate at each time step
#' @param burnin Any additional #'s of values from beginning of chain to discard before calculating PPL statistic (default is 0)
#' @return A matrix with posterior variance (P), sums of squares (G) for the posterior mean and median predictions (compared to Observations), and total posterior loss (D)
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_data_generic<-function(sim.type="RS2closed",S,t.steps,n.transects,line.width,delta=1){
  require(sp)
  require(rgeos)
  require(Matrix)
  require(hierarchicalDS)
  require(ggplot2)
  require(plyr)
  require(grid)
  require(spatstat)
  require(animation)
  source('./STabundance/R/util_funcs.R')
  source('./STabundance/R/sim_funcs.R')
  
  #source('c:/users/paul.conn/git/BOSSst/BOSSst/R/spat_funcs.R')
  if(sqrt(S)%%1>0)cat("error in sim_data_generic; S must be a square number \n")
  if(!(sim.type %in% c("RS2closed","RS2open","CPIF")))cat("error in sim_data_generic; sim.type not recognized")
  
  
  sim.type="CPIF"
  x.len=sqrt(S)
  x.len2=x.len+10
  S2=x.len2^2
  tau.eta.hab=20
  
  #parameters for matern habiatat covariate
  kappa=12
  r=.25
  mu=1000
  Dat.matern=rMatClust(kappa,r,mu)
  X=round((x.len2)*Dat.matern$x+0.5)
  Y=round((x.len2)*Dat.matern$y+0.5)
  Grid=matrix(0,x.len2,x.len2)
  for(i in 1:length(X))Grid[X[i],Y[i]]=Grid[X[i],Y[i]]+1
  Grid=Grid/max(as.vector(Grid))
  Grid.topo=GridTopology(c(0,0),c(1,1),c(x.len2,x.len2))
  Grid.SpG=SpatialGrid(Grid.topo)
  Grid.SpP=as(as(Grid.SpG,"SpatialPixels"),"SpatialPolygons")
  Matern.df=as.vector(Grid)
  Matern.df=as.data.frame(cbind(Matern.df,Matern.df^2))
  colnames(Matern.df)=c("matern","matern2")
  Grid.SpPDF=SpatialPolygonsDataFrame(Grid.SpP,data=Matern.df,match.ID=FALSE)
  laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                         "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(Grid.SpPDF)=CRS(laea_180_proj)
  Data=list(Meta="Simulated space-time dataset")  #define generic data structure for space-time estimation
  Data$Grid=vector("list",t.steps)
  Data$Effort=NULL
  for(it in 1:t.steps)Data$Grid[[it]]=Grid.SpPDF

  #Determine which cells should be included, removed when going from 40 by 40 grid to 30 by 30 grid
  Which.remove=c(1:(x.len2*5),(S2-x.len2*5+1):S2) #ends
  for(icell in 1:S2){
    if(icell%%x.len2<=5)Which.remove=c(Which.remove,icell)
    if(icell%%x.len2>(x.len2-5))Which.remove=c(Which.remove,icell)
  }
  Which.remove=sort(unique(Which.remove))
  Which.include=rep(1,S2)
  Which.include[Which.remove]=0
  Which.include=which(Which.include==1)
 
  #set up some parameters for covariate s-t process
  tau.cov=0.2
  rho.ar1=0.8
  #sig.ar1=(1-rho.ar1)/tau.cov #to make initial variance approx equal to subsequent variance
  sig.ar1=0.6 #to make initial variance approx equal to subsequent variance
  beta0=log(0.3) #-0.5/tau.cov  #lognormal bias correction
  beta1=log(1+delta)  
  
  #Define knot centers for evolution of habitat covariate using 8 by 8 grid over the larger 40 by 40 grid
  Which.knot=NA
  cur.pl=1
  for(i in 1:x.len2){
    for(j in 1:x.len2){
      if(i%%5==3 & j%%5==3)Which.knot=c(Which.knot,cur.pl)
      cur.pl=cur.pl+1
    }
  }
  Knot.centers=gCentroid(Data$Grid[[1]][Which.knot[-1],],byid=TRUE)
  n.knots=length(Knot.centers)
  #calculate kernel densities at grid cell centroids 
  Cell.centroids=gCentroid(Data$Grid[[1]],byid=TRUE)
  Distances=gDistance(Knot.centers,Cell.centroids,byid=TRUE)
  K=matrix(dnorm(Distances,0,5),S2,n.knots)  #knot sd=5 
  K=K/rowSums(K)
  #initial log kernel weights are set via an ICAR(10) prcoess
  Knot.Adj=rect_adj(sqrt(n.knots),sqrt(n.knots))
  Q.knot=-Knot.Adj
  diag(Q.knot)=apply(Knot.Adj,2,'sum')
  Q.knot=Matrix(Q.knot)
  Alpha=matrix(0,n.knots,t.steps)
  Alpha[,1]=rrw(tau.cov*Q.knot)
  #future log kernel weights are set via a random walk process
  for(it in 2:t.steps)Alpha[,it]=rho.ar1*Alpha[,it-1]+rnorm(n.knots,0,sig.ar1)
  #Covariate values are obtained through process convolution
  for(it in 1:t.steps){
    Data$Grid[[it]]@data[,1]=exp(beta0+beta1*(it-1)+K%*%Alpha[,it])
    Data$Grid[[it]]@data[,2]=(Data$Grid[[it]]@data[,1])^2
  }
  #plot_N_map(1,as.matrix(Data$Grid[[1]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
  #plot_N_map(1,as.matrix(Data$Grid[[2]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
  #plot_N_map(1,as.matrix(Data$Grid[[5]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
  #plot_N_map(1,as.matrix(Data$Grid[[10]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
  #plot_N_map(1,as.matrix(Data$Grid[[15]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")  
  #plot_N_map(1,as.matrix(Data$Grid[[20]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
  #plot_N_map(1,as.matrix(Data$Grid[[25]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")  
  #plot_N_map(1,as.matrix(Data$Grid[[30]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
  
  
  
  #Go ahead and switch to 30 by 30 for RS2closed, CPIF
  Cur.S=S2
  if(sim.type!="RS2open"){
    for(it in 1:t.steps)Data$Grid[[it]]=Data$Grid[[it]][Which.include,]
    Cur.S=S
  }
  Data$Adj=rect_adj(sqrt(Cur.S),sqrt(Cur.S))
    
  
  # Compute Data$Which.distances and Data$Dist.entries for RS models
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
  Data$Which.distances=which(is.na(Distances)==FALSE)
  Data$Dist.entries=Distances[Data$Which.distances]  
  
  hab.formula=~matern+matern2
  
  #now evolve spatial process depending on simulation type
  if(sim.type %in% c("RS2closed","RS2open")){
    Sim.pars=list(Hab.init=c(-1.5,20,-13),Hab.evol=c(-1.5,20,-13),tau.eta=10,kern.sd=2,tau.epsilon=20)
    Lambda=sim_RS2(S=Cur.S,Data=Data,Sim.pars=Sim.pars,hab.formula=hab.formula)
  }
  K.cpif=K[Which.include,]
  if(sim.type=="CPIF"){
    Sim.pars=list(Hab=c(-1.5,20,-13),lambda=70000,tau.eta=20,rho.ar1=0.5,sig.ar1=0.1,tau.epsilon=20)
    Lambda=sim_CPIF(S=Cur.S,Data=Data,Sim.pars=Sim.pars,hab.formula=hab.formula,Q.knot=Q.knot,K.cpif=K.cpif)
  }
  
  
  if(sim.type=="RS2open"){
    Lambda=Lambda[Which.include,]
    for(it in 1:t.steps)Data$Grid[[it]]=Data$Grid[[it]][Which.include,]
  }
  Effort=sim_effort(S=S,Data=Data,n.transects=n.transects,line.width=line.width)
  Data$Effort=Effort
  Sim.data=sim_Ncounts(S=S,Data=Data,Lambda=Lambda)
  Data$Count.data=Sim.data$Count.data
  Sim.data$Data=Data
      

  #function for plotting abundance map
  i.plot=FALSE
  if(i.plot){
    plot_N_map(1,Sim.data$N,Grid=Data$Grid)
    plot_N_map(2,Sim.data$N,Grid=Data$Grid)
    plot_N_map(3,Sim.data$N,Grid=Data$Grid)
    plot_N_map(4,Sim.data$N,Grid=Data$Grid)
    plot_N_map(5,Sim.data$N,Grid=Data$Grid)
    plot_N_map(6,Sim.data$N,Grid=Data$Grid)
    plot_N_map(7,Sim.data$N,Grid=Data$Grid)
    plot_N_map(8,Sim.data$N,Grid=Data$Grid)
    plot_N_map(9,Sim.data$N,Grid=Data$Grid)
    plot_N_map(10,Sim.data$N,Grid=Data$Grid)
    plot_N_map(11,Sim.data$N,Grid=Data$Grid)
    plot_N_map(12,Sim.data$N,Grid=Data$Grid)
    plot_N_map(13,Sim.data$N,Grid=Data$Grid)
    plot_N_map(14,Sim.data$N,Grid=Data$Grid)
    plot_N_map(15,Sim.data$N,Grid=Data$Grid)
    plot_N_map(16,Sim.data$N,Grid=Data$Grid)
    plot_N_map(17,Sim.data$N,Grid=Data$Grid)
    plot_N_map(18,Sim.data$N,Grid=Data$Grid)
    plot_N_map(19,Sim.data$N,Grid=Data$Grid)
    plot_N_map(20,Sim.data$N,Grid=Data$Grid)
    plot_N_map(21,Sim.data$N,Grid=Data$Grid)
    plot_N_map(22,Sim.data$N,Grid=Data$Grid)
    plot_N_map(23,Sim.data$N,Grid=Data$Grid)
    plot_N_map(24,Sim.data$N,Grid=Data$Grid)
    plot_N_map(25,Sim.data$N,Grid=Data$Grid)
    plot_N_map(26,Sim.data$N,Grid=Data$Grid)
    plot_N_map(27,Sim.data$N,Grid=Data$Grid)
    plot_N_map(28,Sim.data$N,Grid=Data$Grid)
    plot_N_map(29,Sim.data$N,Grid=Data$Grid)
    plot_N_map(30,Sim.data$N,Grid=Data$Grid)
  }
  return(Sim.data)
}
