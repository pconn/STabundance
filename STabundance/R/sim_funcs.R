# function to simulate spatio-temporal latent abundance process given an RS2closed formulation
#' @param S Number of cells (must be a square number)
#' @param Data A list including at least the following
#'    Grid - A vectored list where each element holds a time-specific SpatialPolygonsDataFrame full of habitat covariates for the area being modeled   
#'    Adj - Adjacency matrix for use in ICAR modeling
#'    Effort - If provided, this data frame provides records of the amount of area sampled each
#'       sample unit and time step via columns "Cell", "Time", and "AreaSurveyed".  If NULL, transects
#'       are simulated using the "n.transects" option below
#'    Which.distances - A vector providing the entries of an S by S distance matrix that are nonzero (i.e. those connected in resource selection kernel)
#'    Dist.entries - A vector with specific distance category values (particular to 3 cell radius model in paper)
#' @param Sim.pars A list holding parameters that describe evolution of the spatio-temporal process.  Included are 
#'         Hab.init - spatial regression parameters for initial time step
#'         Hab.evol - parameters guiding habitat preference associated with the RS process
#'         kern.sd - SD of the redistribution kernel
#'         tau.eta - precision for spatial random effects on initial abundance
#'         tau.epsilon - precision for random noise
#' @param hab.formula A formula object holding the regression model formulation
#' @return Lambda  A matrix holding spatio-temporal lambda values
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_RS2<-function(S,Data,Sim.pars,hab.formula){
  require(rgeos)
  require(Matrix)
  DEBUG=FALSE
  if(DEBUG){
    tau.epsilon=100000
    tau.eta=100000
  }
  
  t.steps=length(Data$Grid)
  X=model.matrix(hab.formula,Data$Grid[[1]]@data)  
  Lambda=matrix(0,S,t.steps)
  N=Lambda
  Q=-Data$Adj
  diag(Q)=apply(Data$Adj,2,'sum')
  Q=Matrix(Sim.pars$tau.eta*Q)
  Eta=rrw(Q)
  Lambda[,1]=exp(X%*%Sim.pars$Hab.init+Eta+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon)))
  #plot_N_map(1,as.matrix(Lambda[,1],ncol=1),Grid=Data$Grid,highlight=c(1,2),cell.width=1,leg.title="Covariate")
 
  
  Phi=Matrix(0,S,S)
  Distances=c(sqrt(2^2+3^2),sqrt(1+3^2),3,sqrt(2^2+2^2),sqrt(1+2^2),2,sqrt(2),1,0)
  Dist.pdf=dnorm(Distances,0,Sim.pars$kern.sd)
  One=rep(1,S)
  srr.cov=FALSE #if true, add spatially correlated noise to resource selection process
  if(srr.cov)Q.cov=Matrix(10*Q) #add some spatially correlated noise to resource selection process
  
  for(it in 2:t.steps){
    X=Data$Grid[[it]]@data  
    X=model.matrix(hab.formula,data=X)
    if(srr.cov==TRUE)Iw=Diagonal(x=as.vector(exp(X%*%Sim.pars$Hab.evol+rrw(Q.cov))))
    else Iw=Diagonal(x=as.vector(exp(X%*%Sim.pars$Hab.evol)))
    Phi[Data$Which.distances]=Dist.pdf[Data$Dist.entries]
    M=Diagonal(x=as.vector(1/(Phi%*%Iw%*%One)))%*%Phi%*%Iw
    Lambda[,it]=as.vector(t(M)%*%Lambda[,it-1])
  }
  Lambda
}



# function to simulate latent abundance with a closed pop ideal free superpopulation model
#' @param S Number of cells (must be a square number)
#' @param Data A list including at least the following
#'    Grid - A vectored list where each element holds a time-specific SpatialPolygonsDataFrame full of habitat covariates for the area being modeled   
#'    Adj - Adjacency matrix for use in ICAR modeling
#'    Effort - If provided, this data frame provides records of the amount of area sampled each
#'       sample unit and time step via columns "Cell", "Time", and "AreaSurveyed".  If NULL, transects
#'       are simulated using the "n.transects" option below
#'    Which.distances - A vector providing the entries of an S by S distance matrix that are nonzero (i.e. those connected in resource selection kernel)
#'    Dist.entries - A vector with specific distance category values (particular to 3 cell radius model in paper)
#' @param Sim.pars A list holding parameters that describe evolution of the spatio-temporal process.  Included are 
#'         Hab - spatial regression parameters
#'         lambda - total expected population size
#'         tau.eta - precision for spatial random effects on initial abundance
#'         rho.ar1 - correlation of AR1 process for knot weights
#'         sig.ar1 - standard deviation for AR1 process
#'         tau.epsilon - precision of random error
#' @param hab.formula A formula object holding the regression model formulation
#' @param Q.knot A structure matrix for reduced rank rsr model for spatio-temporal random effects at the first time step
#' @param K.cpif A matrix holding the N by k weights associated with process convolution
#' @return Lambda  A matrix holding spatio-temporal lambda values
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_CPIF<-function(S,Data,Sim.pars,hab.formula,Q.knot,K.cpif){
  t.steps=length(Data$Grid)
  n.knots=ncol(K.cpif)
  X=model.matrix(hab.formula,Data$Grid[[1]]@data)  
  Lambda=matrix(0,S,t.steps)
  Alpha=matrix(0,n.knots,t.steps)
  Cell.probs=Lambda
  N=Lambda
  Q=Matrix(Sim.pars$tau.eta*Q.knot)
  Alpha[,1]=rrw(Q)
  Cell.probs[,1]=exp(X%*%Sim.pars$Hab+K.cpif%*%Alpha[,1]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon)))
  Lambda[,1]=rmultinom(1,Sim.pars$lambda,Cell.probs[,1])
  #plot_N_map(1,as.matrix(Lambda[,1],ncol=1),Grid=Data$Grid,highlight=c(1,2),cell.width=1,leg.title="Covariate")
  for(it in 2:t.steps){
    X=model.matrix(hab.formula,Data$Grid[[it]]@data)  
    Alpha[,it]=Sim.pars$rho.ar1*Alpha[,it-1]+rnorm(n.knots,0,Sim.pars$sig.ar1)  
    Cell.probs[,it]=exp(X%*%Sim.pars$Hab+K.cpif%*%Alpha[,it]+rnorm(S,0,1/sqrt(Sim.pars$tau.epsilon)))
    Lambda[,it]=rmultinom(1,Sim.pars$lambda,Cell.probs[,it])
  } 
  Lambda
}
#plot_N_map(1,as.matrix(Data$Grid[[30]]@data[,1],ncol=1),Grid=Data$Grid,leg.title="Covariate")
#plot_N_map(1,as.matrix(Lambda[,30],ncol=1),Grid=Data$Grid,highlight=c(1,2),cell.width=1,leg.title="Covariate")


# function to simulate effort data 
#' @param S Number of cells (must be a square number)
#' @param Data A list including at least the following
#'    Grid - A vectored list where each element holds a time-specific SpatialPolygonsDataFrame full of habitat covariates for the area being modeled   
#' @param n.transects Number of transects to simulate at each time step (set to NULL if effort information is provided instead)
#' @param line.width Proportional diameter covered by a transect passing through a cell
#' @return Effort A data.frame holding the following columns:  "Cell", "Time", and "AreaSurveyed"
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_effort<-function(S,Data,n.transects=NULL,line.width){
  t.steps=length(Data$Grid)
  Coords.y=unique(coordinates(Data$Grid[[1]])[,1])
  Effort=data.frame(Cell=rep(NA,n.transects*t.steps*sqrt(S)),Time=rep(1:t.steps,each=n.transects*sqrt(S)),AreaSurveyed=rep(line.width,t.steps*n.transects*sqrt(S)))
  cur.pl=1
  for(it in 1:t.steps){
    #first transect
    Cur.y=sample(Coords.y,n.transects)
    Sampled=which(coordinates(Data$Grid[[1]])[,1]%in%Cur.y)
    Effort[cur.pl:(cur.pl+length(Sampled)-1),1]=Sampled
    cur.pl=cur.pl+length(Sampled)
  }
  Effort
}
  
# function to simulate transect counts given an effort dataset and spatio-temporal lambda values
#' @param S Number of cells (must be a square number)
#' @param Data A list including at least the following
#'    Grid - A vectored list where each element holds a time-specific SpatialPolygonsDataFrame full of habitat covariates for the area being modeled   
#'    Adj - Adjacency matrix for use in ICAR modeling
#'    Effort - Data frame providing records of the amount of area sampled each
#'       sample unit and time step via columns "Cell", "Time", and "AreaSurveyed".  
#' @param Lambda  A matrix holding spatio-temporal lambda values (sample units on rows, time on columns) 
#' @return A listed holding two objects (1) "Count.data" A simulated count dataset, and (2) Simulated true abundance
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn  
sim_Ncounts<-function(S,Data,Lambda){
  t.steps=length(Data$Grid)
  #Simulate transect counts 
  Not.sampled=matrix(1,S,t.steps)
  Count=rep(NA,nrow(Data$Effort))
  for(iobs in 1:nrow(Data$Effort)){
    Count[iobs]=rpois(1,Data$Effort[iobs,"AreaSurveyed"]*Lambda[Data$Effort[iobs,"Cell"],Data$Effort[iobs,"Time"]])
    Not.sampled[Data$Effort[iobs,"Cell"],Data$Effort[iobs,"Time"]]=Not.sampled[Data$Effort[iobs,"Cell"],Data$Effort[iobs,"Time"]]-Data$Effort[iobs,"AreaSurveyed"]
  }
  #Simulate remaining abundance
  N=matrix(rpois(S*t.steps,Lambda*Not.sampled),S,t.steps)
  for(iobs in 1:nrow(Data$Effort))N[Data$Effort[iobs,"Cell"],Data$Effort[iobs,"Time"]]=N[Data$Effort[iobs,"Cell"],Data$Effort[iobs,"Time"]]+Count[iobs]
  
  Sim.out<-list(Count.data=cbind(Data$Effort,Count),N=N)
  Sim.out
} 
  