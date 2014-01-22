#' generate initial values for if not already specified by user
#' @param DM.hab.pois   a list of design matrices for the Poisson habitat model (elements are named sp1,sp2, etc.)
#' @param DM.hab.bern   If a hurdle model, a list of design matrices for the Bernoulli habitat model (elements are named sp1,sp2, etc.) (NULL if not hurdle)
#' @param N.hab.pois.par  vector giving number of parameters in the Poisson habitat model for each species
#' @param N.hab.bern.par  vector giving number of parameters in the Bernoulli habitat model for each species (NULL if not hurdle)
#' @param G.transect a matrix of the number of groups of animals in area covered by each transect; each row gives a separate species		
#' @param Area.trans	a vector giving the proportion of a strata covered by each transect
#' @param Area.hab	a vector of the relative areas of each strata
#' @param Mapping	a vector mapping each transect to it's associated strata
#' @param spat.ind  is spatial independence assumed? (TRUE/FALSE)
#' @param grp.mean  a vector giving the pois1 parameter for group size (one entry for each species)
#' @return a list of initial parameter values
#' @export
#' @keywords initial values, mcmc
#' @author Paul B. Conn
generate_inits_BOSS<-function(DM.hab.pois,DM.hab.bern,N.hab.pois.par,N.hab.bern.par,G.transect,Area.trans,Area.hab,Mapping,spat.ind,grp.mean){		
  i.hurdle=1-is.null(DM.hab.bern)
  n.species=nrow(G.transect)
  n.cells=length(Area.hab)
  hab.pois=matrix(0,n.species,max(N.hab.pois.par))
  hab.bern=NULL
  tau.eta.bern=NULL
  Eta.bern=NULL
  if(i.hurdle==1){
    hab.bern=matrix(0,n.species,max(N.hab.bern.par))
    tau.eta.bern=runif(n.species,0.5,2)
    Eta.bern=matrix(rnorm(n.species*n.cells),n.species,n.cells)
  }
  Nu=matrix(0,n.species,n.cells)
  for(isp in 1:n.species){
    Nu[isp,]=log(max(G.transect[isp,])/mean(Area.trans)*exp(rnorm(length(Area.hab),0,0.1)))
  }
  Par=list(hab.pois=hab.pois,hab.bern=hab.bern,
           Nu=Nu,Eta.pois=matrix(rnorm(n.species*n.cells),n.species,n.cells),Eta.bern=Eta.bern,
           tau.eta.pois=runif(n.species,0.5,2),tau.eta.bern=tau.eta.bern,tau.nu=runif(n.species,0.5,2))
  Par$hab.pois[,1]=log(apply(G.transect,1,'mean')/(mean(Area.trans)*mean(Area.hab))*exp(rnorm(n.species,0,1)))
  Par$G=round(exp(Par$Nu)*Area.hab*exp(rnorm(length(Par$Nu))))
  I.error=(G.transect>Par$G[Mapping])
  if(sum(I.error)>0){
    Which.error=which(G.transect>Par$G[Mapping])
    for(i in 1:sum(I.error))Par$G[Mapping[Which.error[i]]]=Par$G[Mapping[Which.error[i]]]+G.transect[Which.error[i]]
  }
  for(isp in 1:n.species)Par$N[isp,]=Par$G[isp,]+rpois(n.cells,grp.mean[isp]*Par$G[isp,])
  if(spat.ind==1){
    Par$Eta.bern=0*Par$Eta.bern
    Par$Eta.pois=0*Par$Eta.pois
  }
  Par
}



#' function to convert MCMC list vector (used in estimation) into an mcmc object (cf. coda package) 
#' @param MCMC list vector providing MCMC samples for each parameter type 
#' @param N.hab.pois.par see help for mcmc_ds.R
#' @param N.hab.bern.par see help for mcmc_ds.R
#' @param Cov.par.n see help for mcmc_ds.R
#' @param Hab.pois.names see help for mcmc_ds.R
#' @param Hab.bern.names see help for mcmc_ds.R
#' @param Cov.names see help for mcmc_ds.R
#' @param fix.tau.nu see help for mcmc_ds.R
#' @param spat.ind see help for mcmc_ds.R
#' @export
#' @keywords MCMC, coda
#' @author Paul B. Conn
convert.BOSS.to.mcmc<-function(MCMC,N.hab.pois.par,N.hab.bern.par,Cov.par.n,Hab.pois.names,Hab.bern.names,Det.names,Cov.names,fix.tau.nu=FALSE,spat.ind=TRUE,point.ind=TRUE){
  require(coda)
  i.ZIP=!is.na(N.hab.bern.par)[1]
  n.species=nrow(MCMC$Hab.pois)
  n.iter=length(MCMC$Hab.pois[1,,1])
  n.col=n.species*2+sum(N.hab.pois.par)+(1-spat.ind)*n.species+(1-fix.tau.nu)*n.species+sum(Cov.par.n)*n.species
  if(i.ZIP)n.col=n.col+sum(N.hab.bern.par)+(1-spat.ind)*n.species #for ZIP model
  n.cells=dim(MCMC$G)[3]
  Mat=matrix(0,n.iter,n.col)
  Mat[,1:n.species]=t(MCMC$N.tot)
  counter=n.species
  col.names=paste("Abund.sp",c(1:n.species),sep='')
  for(isp in 1:n.species){
    Mat[,counter+isp]=rowSums(as.matrix(MCMC$G[isp,,],nrow=n.iter,ncol=n.cells)) #total abundance of groups
    col.names=c(col.names,paste("Groups.sp",isp,sep=''))
  }
  counter=counter+n.species
  for(isp in 1:n.species){  #habitat parameters
    Mat[,(counter+1):(counter+N.hab.pois.par[isp])]=MCMC$Hab.pois[isp,,1:N.hab.pois.par[isp]]
    col.names=c(col.names,paste("Hab.pois.sp",isp,Hab.pois.names[[isp]],sep=''))
    counter=counter+sum(N.hab.pois.par[isp])
  }
  if(i.ZIP){
    for(isp in 1:n.species){  #habitat parameters
      Mat[,(counter+1):(counter+N.hab.bern.par[isp])]=MCMC$Hab.bern[isp,,1:N.hab.bern.par[isp]]
      col.names=c(col.names,paste("Hab.bern.sp",isp,Hab.bern.names[[isp]],sep=''))
      counter=counter+sum(N.hab.bern.par[isp])
    }   
  }
  if(spat.ind==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta.pois)
    col.names=c(col.names,paste("tau.eta.pois.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(spat.ind==FALSE & i.ZIP){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.eta.bern)
    col.names=c(col.names,paste("tau.eta.bern.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(fix.tau.nu==FALSE){
    Mat[,(counter+1):(counter+n.species)]=t(MCMC$tau.nu)
    col.names=c(col.names,paste("tau.nu.sp",c(1:n.species),sep=''))
    counter=counter+n.species
  }
  if(is.null(Cov.par.n)==FALSE){
    max.par=max(Cov.par.n)
    for(isp in 1:n.species){
      for(ipar in 1:length(Cov.par.n)){
        Mat[,(counter+1):(counter+Cov.par.n[ipar])]=MCMC$Cov.par[isp,,((ipar-1)*max.par+1):((ipar-1)*max.par+Cov.par.n[ipar])]
        counter=counter+Cov.par.n[ipar]
        col.names=c(col.names,paste("Cov.sp",isp,".",Cov.names[[ipar]],sep=''))
      }
    }
  }
  colnames(Mat)=col.names
  Mat=mcmc(Mat)
  Mat
}


#' function to calculate posterior predictive loss given the output object from ST analysis
#' @param Out Output object from running spatio-temporal estimation routine
#' @param burnin Any additional #'s of values from beginning of chain to discard before calculating PPL statistic (default is 0)
#' @return A matrix with posterior variance (P), sums of squares (G) for the posterior mean and median predictions (compared to Observations), and total posterior loss (D)
#' @export
#' @keywords Posterior predictive loss
#' @author Paul B. Conn
post_loss_boss<-function(Out,burnin=0){
  dims.Pred=dim(Out$Pred.det)
  median.Pred=array(0,dim=dims.Pred[2:3])
  mean.Pred=median.Pred
  var.Pred=mean.Pred
  for(itrans in 1:dims.Pred[2]){
    for(iobs in 1:dims.Pred[3]){
        median.Pred[itrans,iobs]=median(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,iobs])
        mean.Pred[itrans,iobs]=mean(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,iobs])
        var.Pred[itrans,iobs]=var(Out$Pred.det[(burnin+1):dims.Pred[1],itrans,iobs])
    }  
  }
  sum.sq.mean=sum((Out$Obs.det-mean.Pred)^2)
  sum.sq.median=sum((Out$Obs.det-median.Pred)^2)
  Loss=matrix(0,2,3)
  colnames(Loss)=c("P","G","D")
  rownames(Loss)=c("mean","median")
  Loss[,1]=sum(var.Pred)
  Loss[1,2]=sum.sq.mean
  Loss[2,2]=sum.sq.median
  Loss[,3]=rowSums(Loss[1:2,1:2])
  Loss
}

#' Produce an RW2 structure matrix for a line (e.g. for time series)
#' @param x length of vector
#' @return RW2 precision matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn
linear_adj_RW2 <- function(x){
  require(Matrix)
  Q=Matrix(0,x,x)
  for(i in 1:(x-2)){
    Q[i,i+2]=1
    Q[i+2,i]=1
  }
  for(i in 1:(x-1)){
    Q[i,i+1]=-4
    Q[i+1,i]=-4
  }
  Q[1,2]=-2
  Q[2,1]=-2
  Q[x,x-1]=-2
  Q[x-1,x]=-2
  diag(Q)=-rowSums(Q)
  Q
}  

#' Produce an RW1 adjacency matrix for a rectangular grid for use with areal spatial models (queens move)
#' @param x number of cells on horizontal side of grid
#' @param y number of cells on vertical side of grid
#' @param byrow If TRUE, cell indices are filled along rows (default is FALSE)
#' @return adjacency matrix
#' @export 
#' @keywords adjacency
#' @author Paul Conn \email{paul.conn@@noaa.gov}
rect_adj <- function(x,y,byrow=FALSE){
  Ind=matrix(c(1:(x*y)),y,x,byrow=byrow)
  if(byrow==TRUE)Ind=t(Ind)
  n.row=nrow(Ind)
  n.col=ncol(Ind)
  Adj=matrix(0,x*y,x*y)
  for(i in 1:n.row){
    for(j in 1:n.col){
      if(i==1 & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i==1 & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i==1 & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
      }
      if(i>1 & i<n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row+1]=1
      }
      if(i>1 & i<n.row & j>1 & j<n.col){
        cur.nums=c(Ind[i,j]-n.row-1,Ind[i,j]-n.row,Ind[i,j]-n.row+1,Ind[i,j]-1,Ind[i,j]+1,Ind[i,j]+n.row-1,Ind[i,j]+n.row,Ind[i,j]+n.row+1)
        Adj[Ind[i,j],cur.nums]=1
      }
      if(i>1 & i<n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]+1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row+1]=1
        
      }
      if(i==n.row & j==1){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
      }
      if(i==n.row & j>1 & j<n.col){
        Adj[Ind[i,j],Ind[i,j]+n.row]=1
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]+n.row-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
      if(i==n.row & j==n.col){
        Adj[Ind[i,j],Ind[i,j]-1]=1
        Adj[Ind[i,j],Ind[i,j]-n.row]=1
        Adj[Ind[i,j],Ind[i,j]-n.row-1]=1
      }
    }
  }
  if(byrow==TRUE)Adj=t(Adj)
  return(Adj)
}

#' Produce weights for bivariate normal kernel for upper left triangular quadrant
#' To save computing time, some building blocks are passed into
#' the function
#' @param Tmp.vec A length n vector for holding integrals of bivariate normal
#' @param XY 2xn matrix giving x and y distances
#' @param Sigma vector of length 2 giving sd of bivariate normal in the x and y directions 
#' @return a filled vector that holds unnormalized redistribution kernel values
#' @export 
#' @keywords bivariate normal, kernel weight
#' @author Paul Conn \email{paul.conn@@noaa.gov}
d_biv_normal<-function(Tmp.vec,XY,Sigma){
  return(dnorm(XY[1,],0,Sigma[1])*dnorm(XY[2,],0,Sigma[2]))
}
  
#' Function to plot abundance map for BOSS survey grid
#' @param cur.t  Time step to plot map for
#' @param N A matrix holding abundance point estimates; different rows correspond to different sampling units, columns correspond to time step
#' @param Grid A list object, each element of which holds a spatial polygons data frame for each time step
#' @param highlight If provided, this vector specifies which cells to separately highlight during plotting
#' @param cell.width if highlight is provided, cell.width must provide the width of a composite grid cell
#' @param leg.title Title for legend of plot (if different from the default "Abundance")
#' @return An abundance map
#' @export 
#' @keywords abundance map, plot
#' @author Paul Conn \email{paul.conn@@noaa.gov}
plot_N_map<-function(cur.t,N,Grid,highlight=NULL,cell.width,leg.title="Abundance"){
  require(rgeos)
  require(ggplot2)
  Tmp=Grid[[1]]
  if(is.null(highlight)==FALSE){
    midpoints=data.frame(gCentroid(Tmp[highlight,],byid=TRUE))
    colnames(midpoints)=c("Easting","Northing")
  }
  Abundance=N[,cur.t]
  Cur.df=cbind(data.frame(gCentroid(Tmp,byid=TRUE)),Abundance)
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme+scale_fill_continuous(name=leg.title)
  if(is.null(highlight)==FALSE){
    #p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
    p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-1/2,xmax=Easting+1/2,ymin=Northing-1/2,ymax=Northing+1/2))
  }
  p1
}

#' SIMULATE AN ICAR PROCESS 
#' @param Q Precision matrix for the ICAR process
#' @return Spatial random effects
#' @export 
#' @keywords ICAR, simulation
#' @author Devin Johnson
rrw <- function(Q){
  v <- eigen(Q, TRUE)
  val.inv <- sqrt(ifelse(v$values>sqrt(.Machine$double.eps), 1/v$values, 0))
  P <- v$vectors
  sim <- P%*%diag(val.inv)%*%rnorm(dim(Q)[1], 0, 1)
  X <- rep(1,length(sim))
  if(sum(val.inv==0)==2) X <- cbind(X, 1:length(sim))
  sim <- sim-X%*%solve(crossprod(X), crossprod(X,sim))
  return(sim)
}

#' stack list of matrices into one big matrix
#' @param L list vector of matrices
#' @return X one matrix consisting of stacked elements of L
#' @export 
#' @keywords logit, expit
#' @author Paul Conn \email{paul.conn@@noaa.gov}
stack_list<-function(L){
  X=L[[1]]
  len.L=length(L)
  if(len.L>1){
    for(ilen in 2:len.L)X=rbind(X,L[[ilen]])  #could be better optimized
  }
  X
}


#' inverse logit transformation
#' @param x quantity to be transformed
#' @return real valued quantities
#' @export 
#' @keywords logit, expit
#' @author Paul Conn \email{paul.conn@@noaa.gov}
expit<-function(x)1/(1+exp(-x))

#' logit transformation
#' @param x quantity to be transformed
#' @return logit value(s)
#' @export 
#' @keywords logit
#' @author Paul Conn \email{paul.conn@@noaa.gov}
logit<-function(x)log(x/(1-x))
 