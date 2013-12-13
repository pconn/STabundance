### sim_photo_data.R
### function to simulate data for spatio-temporal abundance analysis from BOSS ice data


require(sp)
require(rgeos)
require(Matrix)
require(hierarchicalDS)
require(ggplot2)
require(plyr)
require(grid)
source('c:/users/paul.conn/git/BOSSst/BOSSst/R/util_funcs.R')
#source('c:/users/paul.conn/git/BOSSst/BOSSst/R/spat_funcs.R')

#load spatial-temporal covariate & grid data
load("BOSSst_2012data.Rdata")  #boss grid, ice data
load("Effort2012_BOSSst.Rdata")  #load effort data indicating grid cells and times surveyed (data produced with format_effort.R)

t.steps=length(Data$Grid) # number of time steps
S=nrow(Data$Grid[[1]])

#initialize abundance for a single species over the "Grid"
set.seed(12345)
X=Data$Grid[[1]]@data[,c(5,7,4)]  #initial design matrix has ice concentration, distance from edge, dist from shelf
X=cbind(X,(X[,1])^2,sqrt(X[,2])) #also include 
colnames(X)=c("ice","edge","shelf","ice2","sqrt_edge")
formula=~ice+ice2+edge+sqrt_edge+shelf
X=model.matrix(formula,data=X)

N=matrix(0,S,t.steps)
Hab.pars=c(-1.5,20,-13,-.7,-.7,-2.5)
Q=-Data$Adj
diag(Q)=apply(Data$Adj,2,'sum')
tau=15
Q=Matrix(tau*Q)
Eta=rrw(Q)
N[,1]=exp(X%*%Hab.pars+Eta+rnorm(S,0,0.2))


#plot initial state distribution
plot_N_map(1,N,Grid=Data$Grid)


Phi=Matrix(0,S,S)
sigma.rk=2
Distances=c(sqrt(2^2+3^2),sqrt(1+3^2),3,sqrt(2^2+2^2),sqrt(1+2^2),2,sqrt(2),1,0)
Dist.pdf=dnorm(Distances,0,sigma.rk)
One=rep(1,S)
srr.cov=FALSE #if true, add spatially correlated noise to resource selection process
Q.cov=Matrix(10*Q) #add some spatially correlated noise to resource selection process
#simulate movement assuming an overdispersed multinomial distribution
for(it in 2:t.steps){
  X=Data$Grid[[it]]@data[,c(5,7,4)]  #initial design matrix has ice concentration, distance from edge, dist from shelf
  if(length(which(X[,1]<0))>0)X[which(X[,1]<0),1]=0
  if(length(which(X[,1]>1))>0)X[which(X[,1]>1),1]=1
  X=cbind(X,(X[,1])^2,sqrt(X[,2])) #also include quadratic for ice, sqrt of dist from ice edge
  colnames(X)=c("ice","edge","shelf","ice2","sqrt_edge")
  X=model.matrix(formula,data=X)
  if(srr.cov==TRUE)Iw=Diagonal(x=as.vector(exp(X%*%Hab.pars+rrw(Q.cov))))
  else Iw=Diagonal(x=as.vector(exp(X%*%Hab.pars)))
  Phi[Data$Which.distances]=Dist.pdf[Data$Dist.entries]
  M=Diagonal(x=as.vector(1/(Phi%*%Iw%*%One)))%*%Phi%*%Iw
  N[,it]=t(M)%*%N[,it-1])
}

#function for plotting abundance map
plot_N_map(1,matrix(diag(Iw),ncol=1),Grid=Data$Grid,highlight=Which.1035)

Which.1035=which(M[1035,]>0)
  
#produce plots of abundance over time
plot_N_map(cur.t=1,N=N,Grid=Data$Grid,highlight=1035)
plot_N_map(cur.t=9,N=N,Grid=Data$Grid)
plot_N_map(cur.t=18,N=N,Grid=Data$Grid)
plot_N_map(cur.t=24,N=N,Grid=Data$Grid)
plot_N_map(cur.t=25,N=N,Grid=Data$Grid)
plot_N_map(cur.t=26,N=N,Grid=Data$Grid)
plot_N_map(cur.t=27,N=N,Grid=Data$Grid)
plot_N_map(cur.t=28,N=N,Grid=Data$Grid)
plot_N_map(cur.t=29,N=N,Grid=Data$Grid)
plot_N_map(cur.t=30,N=N,Grid=Data$Grid)
plot_N_map(cur.t=31,N=N,Grid=Data$Grid)
plot_N_map(cur.t=32,N=N,Grid=Data$Grid)
plot_N_map(cur.t=33,N=N,Grid=Data$Grid)

plot_N_map(cur.t=36,N=N,Grid=Data$Grid)


plot_N_map(1,matrix(Data$Grid[[24]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)
plot_N_map(1,matrix(Data$Grid[[25]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)
plot_N_map(1,matrix(Data$Grid[[26]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)
plot_N_map(1,matrix(Data$Grid[[27]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)
plot_N_map(1,matrix(Data$Grid[[28]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)
plot_N_map(1,matrix(Data$Grid[[29]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)
plot_N_map(1,matrix(Data$Grid[[30]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)
plot_N_map(1,matrix(Data$Grid[[36]]@data[,"ice_conc"],ncol=1),Grid=Data$Grid,highlight=1100)






#final data format:
#Dat - same as previous
#  col 1- transect ID (now, separate for each cell/day combo)
#      2- photo obtained (0/1) ? 
#' 		 3- Observation type (integer - the max integer being 'unknown' if applicable) [NOTE: modeled as factor, but need to be input as integers to account for unknown species observations]
#' 		 x- Group size


