#' function to perform Bayesian analysis of count data using a resource selection model operating on the log of abundance intensity
#' @param Predictors A vector specifying the names of predictors for habitat selection (such predictor names also need to be present in each Data$Grid object) 
#'        If Predictors is missing, an intercept model is assumed
#' @param Data   A list holding the following objects:
#'        Adj - Adjacency matrix of ones and zeros describing spatial connectivity (for ICAR modeling)
#'        Grid - SpatialPolygongsDataFrame holding covariate values for all cells one wants to make predictions on (including sampled cells)
#'        Count.data - A data frame with the following columns: (1) "Cell" - entries this column provides the numerical index (in 1:S) for the cell surveyed
#'              (2) "Time" An integer value providing the time (e.g. day) when sampling occurred, (3) "AreaSurveyed" - provides
#'              the proportion of the target cell that was surveyed, (4) "Count" Number of target species counted
#' @param Area.adjust   A vector allowing for differences in suitable habitat for each cell.  Can be used for different grid cell sizes or different
#'        proportions of suitable habitat (e.g., 1.0 = 100% of habitat is suitable for the focal species)
#' @param Control A list giving MCMC controls and additional model options
#'  "iter": number of MCMC iterations;
#'  "burnin": number of MCMC burnin iterations;
#'	"thin": if specified, how many iterations to skip between recorded posterior samples;
#'  "srr.tol": Threshold eigenvalue level for SRR; only eigenvectors with higher eigenvalues than srr.tol are included in SRR formulation (default is 0.5)
#'  "predict" If TRUE, calculate posterior predictions across the entire grid
#'  "MH.mu" A vector providing continuous Uniform half-range values for Metropolis-Hastings proposals
#'  "adapt" If true, adapts tuning parameters for nu updates; e.g., one could run a chain with adapt=TRUE, and use adapted MH tuning parameters in a longer analysis (default FALSE)
#'  "fix.tau.epsilon" If TRUE, fixes tau.epsilon to 100
#'  "fix.tau.varepsilon" If TRUE, fixes tau.varepsilon to 100
#' @param Prior.pars  An optional list giving prior parameters; these include
#'  "beta.tau" precision for Gaussian prior on regression parameters (default 0.1)
#'  "a.eps" alpha parameter for tau_epsilon~Gamma(alpha,beta) (default = 1.0)
#'  "b.eps" beta parameter for tau_epsilon~Gamma(alpha,beta) (default = 0.01)
#'   If provided, all values for Prior.pars need to be specified
#' @return returns a list with the following objecs: 
#' 	MCMC: A list object containing posterior samples;
#'  Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis-Hastings;
#'  Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#' @export
#' @import Matrix
#' @keywords abundance, mcmc, spatio-temporal model, spatial prediction
#' @author Paul B. Conn \email{paul.conn@@noaa.gov} 
#' @examples print("Later!")


mcmc_OPRS<-function(Predictors=NULL,Data,Prior.pars=NULL,Control,Area.adjust=NULL){
  require(Matrix)
  require(mvtnorm)
  
  DEBUG=FALSE
  S=nrow(Data$Grid[[1]])
  t.steps=length(Data$Grid)
  n.obs=nrow(Data$Count.data)
  n.ST=S*t.steps
  if(is.null(Prior.pars)){
    Prior.pars=list(beta.tau=0.01,
                    a.eps=1,
                    b.eps=0.01,
                    a.vareps=5,
                    b.vareps=0.05,
                    a.eta=1,
                    b.eta=0.01,
                    a.gamma=1,
                    b.gamma=0.01,
                    a.tau.d=1,
                    b.tau.d=0.01)
  }
  beta.sd=1/sqrt(Prior.pars$beta.tau)
  if(is.null(Area.adjust))Area.adjust=rep(1,S)
  #summary statistics to quantify number of counts made in each time step
  Nt.obs=rep(0,t.steps)
  Time.indices=matrix(0,t.steps,2) #beginning and ending entries of "Count" for each time step
  cur.pl=1
  for(it in 1:t.steps){
    Nt.obs[it]=sum(Data$Count.data[,"Time"]==it)
    Time.indices[it,1]=cur.pl
    Time.indices[it,2]=(cur.pl+Nt.obs[it]-1)
    cur.pl=cur.pl+Nt.obs[it]
  }

  Offset=Data$Count.data[,"AreaSurveyed"]  
  Count=Data$Count.data[,"Count"]
  Offset=Offset*Area.adjust[Data$Count.data[,"Cell"]]
  Log.offset=log(Offset)
  Log.area.adjust=log(Area.adjust)
  
  #Construct model formulas from 'Predictors' (need intercept for initial distribution but not for resource selection)
  model=vector("list",2)
  if(is.null(Predictors)){
    model[[1]]=~1
    model[[2]]=NULL
  }
  else{
    model[[1]]=paste('~',Predictors[1],sep='')
    model[[2]]=paste('~0+',Predictors[1],sep='')
    if(length(Predictors)>1){
      for(i in 2:length(Predictors)){
        model[[1]]=paste(model[[1]],'+',Predictors[i],sep='')
        model[[2]]=paste(model[[2]],'+',Predictors[i],sep='')
      }
    }
    model[[1]]=formula(model[[1]])
    model[[2]]=formula(model[[2]])
  }
  
  #set up design matrices
  X.pred.list=vector("list",t.steps)
  X.pred.list[[1]]=model.matrix(model[[1]],data=Data$Grid[[1]]@data)
  for(it in 2:t.steps){
    X.pred.list[[it]]=model.matrix(model[[2]],data=Data$Grid[[it]]@data)
    X.pred.list[[it]]=cbind(rep(0,S),X.pred.list[[it]])  #for missing intercept
  }
  X.pred=stack_list(X.pred.list)
  n.beta=ncol(X.pred)
  Which.obs=(Data$Count.data[,"Time"]-1)*S+(Data$Count.data[,"Cell"])
  Which.no.obs=c(1:nrow(X.pred))
  Which.no.obs=Which.no.obs[-Which.obs]
  X.obs=matrix(X.pred[Which.obs,],n.obs,n.beta)
  n.no=length(Which.no.obs)
  #X.no.obs=matrix(X.pred[Which.no.obs,],n.obs,n.beta)
  #XpXinv=solve(crossprod(X.obs))
  #XpXinvXp=XpXinv%*%t(X.obs)
    
  
  ###### Declare MCMC structures, set initial values
  mcmc.length=(Control$iter-Control$burnin)/Control$thin  
  sigma=runif(1,0.1,1.0)
  Sigma.mc=rep(0,mcmc.length)
  Beta=rep(0,n.beta)
  #set Beta to reasonable starting values
  X.beta=data.frame(X.obs)
  colnames(X.beta)=colnames(X.pred.list[[1]])
  X.beta$Count1=log(Count+1)
  beta.formula=as.formula(paste("Count1 ~",model[[1]])[2])
  Beta = lm(beta.formula,data=X.beta)$coefficients
  Beta.mc=matrix(0,n.beta,mcmc.length)
  
  
  tau.d=runif(1,1,3)
  if(DEBUG)tau.d=1/4
  tau.varepsilon=100
  tau.epsilon=100
  tau.eta=runif(1,10,100)
  tau.gamma=runif(1,10,100)
  Tau.varepsilon.mc=rep(0,mcmc.length)
  Tau.epsilon.mc=rep(0,mcmc.length)
  Tau.d.mc=Tau.epsilon.mc
  Tau.eta.mc=Tau.epsilon.mc
  Tau.gamma.mc=Tau.epsilon.mc
  Pred.mc=array(0,dim=c(S,t.steps,mcmc.length))
  Accept=list(Z=rep(0,t.steps),Beta=rep(0,n.beta),tau.d=0)
  Accept.old=Accept
  
  if(DEBUG){
    Beta=c(3,10,-10)
    tau.epsilon=15
    tau.varepsilon=100
  }
  
  #provide initial estimate of Mu  
  Mu=rep(0,n.ST)
  Mu=exp(X.pred%*%Beta)
  Mu[Which.obs]=Data$Count.data[,"Count"]/Offset
  Mu[which(Mu==0)]=1
  #Mu[which(is.na(Mu))]=median(Mu[which(Mu>0)])
  Mu=log(Mu)
  Mu.list=vector("list",t.steps)
  Z.list=Mu.list
  Z=Mu+0.1*rnorm(n.ST)
  Z_func<-function(z,count,mu,tau,offset)return((offset*exp(z)+z*tau-count-tau*mu)^2)
  #initialize Z in sampled cells so gradient of Posterior w.r.t. Z is zero (better starting place for Langevin updates)
  for(iobs in 1:n.obs){
    Z[Which.obs[iobs]]=optimize(f=Z_func,interval=c(Mu[Which.obs[iobs]]-3,Mu[Which.obs[iobs]]+3),offset=Offset[iobs],tau=tau.varepsilon,count=Count[iobs],mu=Mu[Which.obs[iobs]])$minimum
  }  
  counter=1  
  for(it in 1:t.steps){
    Mu.list[[it]]=Mu[counter:(counter+S-1)]
    Z.list[[it]]=Z[counter:(counter+S-1)]
    counter=counter+S
  }
  
  #setup RW2 time series model
  Gamma<-c(0,rnorm(t.steps-1,0,1/sqrt(tau.gamma)))
  
  #set up process conv model for spatial effects
  #initial log kernel weights are set via an ICAR(10) prcoess
  n.knots=ncol(Data$K)
  Alpha=rnorm(n.knots,0,sqrt(1/tau.eta)) #initial kernel weights/random effects
  Eta=Data$K%*%Alpha
  K.t=t(Data$K)
  cross.K<-crossprod(Data$K,Data$K)
  I.knot=diag(n.knots)
  if(DEBUG)Eta=Eta*0
  
  #set up resource selection, forward-backward machinery
  Distances=c(sqrt(2^2+3^2),sqrt(1+3^2),3,sqrt(2^2+2^2),sqrt(1+2^2),2,sqrt(2),1,0)
  Dist.pdf=dnorm(Distances,0,sqrt(1/tau.d))
  One=rep(1,S)
  M=MM=Iw=Theta=Which.list=Which.list2=Count.list=Log.offset.list=H=H.t=R=vector("list",t.steps)
  I=Diagonal(S+1)
  Q=Matrix(0,S+1,S+1)
  diag(Q)=c(rep(tau.epsilon^{-1},S),1/tau.gamma)
  Phi=Matrix(0,S,S)
  Phi[Data$Which.distances]=Dist.pdf[Data$Dist.entries]  
  counter=1
  for(it in 1:t.steps){  #Note some of the first time step values (e.g. M[[1]]) won't be used
    Which.list[[it]]=Data$Count.data[which(Data$Count.data[,"Time"]==it),"Cell"]  #which list entries are observed at each time step
    Iw[[it]]=Diagonal(x=as.vector(exp(X.pred.list[[it]]%*%Beta)))
    M[[it]]=Diagonal(x=as.vector(1/(Phi%*%Iw[[it]]%*%One)))%*%Phi%*%Iw[[it]]
    Theta[[it]]=c(Mu.list[[it]],Gamma[it])
    MM[[it]]=Matrix(0,S+1,S+1) #\mathcal{M} in ms
    MM[[it]][1:S,1:S]=M[[it]]  
    MM[[it]][,S+1]=1
    MM[[it]][S+1,S+1]=sigma
    Count.list[[it]]=Count[which(Data$Count.data[,"Time"]==it)]
    Log.offset.list[[it]]=Log.offset[which(Data$Count.data[,"Time"]==it)]
    if(Nt.obs[it]>0){
      Which.list2[[it]]=c(counter:(counter+Nt.obs[it]-1))
      counter=counter+Nt.obs[it]
      H[[it]]=matrix(0,Nt.obs[it],S+1)
      for(iobs in 1:Nt.obs[it])H[[it]][iobs,Which.list[[it]][iobs]]=1
      H.t[[it]]=t(H[[it]])
      R[[it]]=Matrix(0,Nt.obs[it],Nt.obs[it])
      diag(R[[it]])=1/tau.varepsilon
    }
  }  
  P=vector("list",t.steps)
  P[[1]]=Q
  W.plus.z=matrix(0,n.obs,t.steps)
  W.plus.y=matrix(0,S+1,t.steps)
  Y.plus=W.plus.y
  Z.plus=W.plus.z
  Nu.list=Z.list
  Finv=Z.list
  G=vector("list",t.steps)
  L=G
  A=L
  R.disturb=L
  R.disturb[[t.steps]]=Matrix(0,S+1,1)
  Hat=L
  Iw.prop=Iw
  Phi.prop=Phi
  M.prop=M
  
  if(DEBUG){
    for(it in 1:t.steps){
      Z.list[[it]]=log(Sim.data$N[,it]+0.01)
      Mu.list[[it]]=Z.list[[it]]
    }
    Z=unlist(Z)
    Mu=unlist(Mu.list)
    Eta=0
  }

  
  for(iiter in 1:Control$iter){
    #if(iiter%%1000==0)
    cat(paste('iteration ',iiter,' of ',Control$iter,'\n'))
    
    #update tau.varepsilon
    if(Control$fix.tau.varepsilon==FALSE){
      if(iiter>100){  #start out tau.varepsilon low until Z, Mu get compatible
        Z.ssq=sum((Z[Which.obs]-Mu[Which.obs])^2)
        tau.varepsilon<-rgamma(1,0.5*n.obs + Prior.pars$a.vareps, Z.ssq*0.5 + Prior.pars$b.vareps)
        for(it in 1:t.steps)diag(R[[it]])=1/tau.varepsilon     
      }
    }
    
    #update Z (sampled cells) 
    if(DEBUG==FALSE){
    for(it in 1:t.steps){
      if(Nt.obs[it]>0){
        Dold=Control$MH.Z[it]^2*0.5*d_logP_Z(Z=Z.list[[it]][Which.list[[it]]],Log.offset=Log.offset.list[[it]],Counts=Count.list[[it]],Mu=Mu.list[[it]][Which.list[[it]]],tau=tau.varepsilon) 
        Prop=Z.list[[it]][Which.list[[it]]]+Dold+rnorm(Nt.obs[it],0,Control$MH.Z[it])
        post.new=sum(Count.list[[it]]*Prop-exp(Log.offset.list[[it]]+Prop)-0.5*tau.varepsilon*(Prop-Mu.list[[it]][Which.list[[it]]])^2)
        post.old=sum(Count.list[[it]]*Z.list[[it]][Which.list[[it]]]-exp(Log.offset.list[[it]]+Z.list[[it]][Which.list[[it]]])-0.5*tau.varepsilon*(Z.list[[it]][Which.list[[it]]]-Mu.list[[it]][Which.list[[it]]])^2)
        Temp=Prop-Z.list[[it]][Which.list[[it]]]-Dold
        jump.old.to.new=-0.5/Control$MH.Z[it]^2*Temp%*%Temp
        Dstar=Control$MH.Z[it]^2*0.5*d_logP_Z(Z=Prop,Log.offset=Log.offset.list[[it]],Counts=Count.list[[it]],Mu=Mu.list[[it]][Which.list[[it]]],tau=tau.varepsilon) 
        Temp=Z.list[[it]][Which.list[[it]]]-Prop-Dstar
        jump.new.to.old=-0.5/Control$MH.Z[it]^2*Temp%*%Temp
        if(runif(1)<exp(post.new-post.old+jump.new.to.old-jump.old.to.new)){
          Z.list[[it]][Which.list[[it]]]=Prop
          Accept$Z[it]=Accept$Z[it]+1
        }  
      }
    }   
    Z=unlist(Z.list)
    Z[Which.no.obs]=rnorm(n.no,Mu[Which.no.obs],tau.varepsilon^(-0.5))
    
  

    
    #update mu using disturbance filter/smoother algorithm of Durbin and Koopman 2002
    W.plus.y[1:S,]=sqrt(1/tau.epsilon)*rnorm(S*t.steps)
    W.plus.y[S+1,]=sqrt(1/tau.gamma)*rnorm(t.steps)
    Y.plus[1:S,1]=X.pred.list[[1]]%*%Beta+Eta+sqrt(1/tau.epsilon)*rnorm(S) #holding Y-plus
    Y.plus[S+1,1]=rnorm(1,0,sqrt(1/tau.gamma))
    for(it in 2:t.steps)Y.plus[,it]=as.numeric(MM[[it]]%*%Y.plus[,it-1]+W.plus.y[,it])  
    Z.star=Z[Which.obs]-(Y.plus[1:S,][Which.obs]+sqrt(1/tau.varepsilon)*rnorm(n.obs)) #Z* = Z - Z^+
    #now use Kalman filter and disturbance smoother to obtain \hat{Z}* = E(Y|Z*) 
    #Init.pred=c(X.pred.list[[1]]%*%Beta+Eta,0)
    Init.pred=rep(0,S+1) #since running the KF on z-z.plus, think this needs to be zero
    Nu.list[[1]]=Z.star[Which.list2[[1]]]-H[[1]]%*%Init.pred
    P[[1]]=Q
    Finv[[1]]=solve(H[[1]]%*%P[[1]]%*%H.t[[1]]+R[[1]])
    G[[1]]=MM[[2]]%*%(P[[1]]%*%(H.t[[1]]%*%Finv[[1]]))
    L[[1]]=MM[[2]]-G[[1]]%*%H[[1]]
    A[[1]]=Init.pred
    for(it in 2:(t.steps-1)){  #forward
      P[[it]]=tcrossprod(MM[[it]]%*%P[[it-1]],L[[it-1]])+Q
      Finv[[it]]=solve(H[[it]]%*%P[[it]]%*%H.t[[it]]+R[[it]])
      G[[it]]=MM[[it+1]]%*%(P[[it]]%*%(H.t[[it]]%*%Finv[[it]]))
      A[[it]]=MM[[it]]%*%A[[it-1]]+G[[it-1]]%*%Nu.list[[it-1]]
      Nu.list[[it]]=Z.star[Which.list2[[it]]]-H[[it]]%*%A[[it]]
      L[[it]]=MM[[it+1]]-G[[it]]%*%H[[it]]    
    }
    P[[t.steps]]=tcrossprod(MM[[t.steps]]%*%P[[t.steps-1]],L[[t.steps-1]])+Q
    Finv[[t.steps]]=solve(H[[t.steps]]%*%P[[t.steps]]%*%H.t[[t.steps]]+R[[t.steps]])
    G[[t.steps]]=MM[[t.steps]]%*%(P[[t.steps]]%*%(H.t[[t.steps]]%*%Finv[[t.steps]]))
    L[[t.steps]]=MM[[t.steps]]-G[[t.steps]]%*%H[[t.steps]]  #use MM[[t.steps]] here since we don't really care about t.steps+1
    A[[t.steps]]=MM[[t.steps]]%*%A[[t.steps-1]]+G[[t.steps-1]]%*%Nu.list[[t.steps-1]]    
    Nu.list[[t.steps]]=Z.star[Which.list2[[t.steps]]]-H[[t.steps]]%*%A[[t.steps]]
    for(it in t.steps:2)R.disturb[[it-1]]=H.t[[it]]%*%(Finv[[it]]%*%Nu.list[[it]])+crossprod(L[[it]],R.disturb[[it]])
    R.disturb.0=H.t[[1]]%*%(Finv[[1]]%*%Nu.list[[1]])+crossprod(L[[1]],R.disturb[[1]])
    Hat[[1]]=A[[1]]+P[[1]]%*%R.disturb.0
    for(it in 1:(t.steps-1)){
      Hat[[it+1]]=MM[[it+1]]%*%Hat[[it]]+Q%*%R.disturb[[it]]
      Theta[[it]]=Hat[[it]]+Y.plus[,it]
      Mu.list[[it]]=Theta[[it]][1:S]
      Gamma[it]=Theta[[it]][S+1]
    }
    Gamma[1]=0
    Theta[[t.steps]]=Hat[[t.steps]]+Y.plus[,t.steps]
    Mu.list[[t.steps]]=Theta[[t.steps]][1:S]
    Gamma[[t.steps]]=Theta[[t.steps]][S+1]
    Mu=unlist(Mu.list)
    }

    #if(DEBUG){
    #  icell=1
    #  cur.mu=rep(0,20)
    #  cur.disturb=cur.mu
    #  cur.a=cur.mu
    #   for(it in 1:t.steps){
    #    cur.mu[it]=Mu.list[[it]][icell]
    #    cur.disturb[it]=R.disturb[[it]][icell]
    #    cur.a[it]=A[[it]][icell]
    #  }
    # plot(log(Sim.data$N[icell,]),type="l",col="red",ylim=c(0,7))
     # crap=Data$Count.data[which(Data$Count.data[,"Cell"]==icell),]
    #  points(crap$Time,log(crap$Count*5))
    #  lines(cur.mu)
    #  lines(Y.plus[icell,],lty=2)
    #  lines(cur.a+Y.plus[icell,],col='blue',lty=2)
      #lines(cur.disturb,col="red")
    #}
    
    
    #update beta
    Prop=Beta
    Old.resid.ssq=crossprod(Mu.list[[1]]-X.pred.list[[1]]%*%Beta-Eta)
    for(it in 2:t.steps)Old.resid.ssq=Old.resid.ssq+crossprod(Mu.list[[it]]-M[[it]]%*%Mu.list[[it-1]]-Gamma[it])
    old.logL=as.numeric(-0.5*tau.epsilon*Old.resid.ssq) 
    for(ipar in 1:n.beta){  
      Prop[ipar]=Beta[ipar]+rnorm(1,0,Control$MH.beta[ipar])
      New.resid.ssq=crossprod(Mu.list[[1]]-X.pred.list[[1]]%*%Prop-Eta)
      for(it in 2:t.steps){
        diag(Iw.prop[[it]])=as.vector(exp(X.pred.list[[it]]%*%Prop))
        M.prop[[it]]=Diagonal(x=as.vector(1/(Phi%*%Iw.prop[[it]]%*%One)))%*%Phi%*%Iw.prop[[it]]
        New.resid.ssq=New.resid.ssq+crossprod(Mu.list[[it]]-M.prop[[it]]%*%Mu.list[[it-1]]-Gamma[it])        
      }
      new.logL=as.numeric(-0.5*tau.epsilon*New.resid.ssq)
      if(runif(1)<exp(new.logL+dnorm(Prop[ipar],0,beta.sd,log=TRUE)-old.logL-dnorm(Beta[ipar],0,beta.sd,log=TRUE))){
        Beta[ipar]=Prop[ipar]
        Accept$Beta[ipar]=Accept$Beta[ipar]+1
        old.logL=new.logL
        M=M.prop
        Iw=Iw.prop
        Old.resid.ssq=New.resid.ssq
      }
      else Prop[ipar]=Beta[ipar]
    }
        
    #update precision for system state errors (epsilon)
    if(Control$fix.tau.epsilon==FALSE){
      tau.epsilon <- rgamma(1,0.5*n.ST + Prior.pars$a.eps, as.numeric(Old.resid.ssq)*0.5 + Prior.pars$b.eps)
    }
    
    #update tau_d (precision for redistribution kernel)
    prop=tau.d+rnorm(1,0,Control$MH.tau.d)
    if(prop>0){
      Old.resid.ssq=Old.resid.ssq-crossprod(Mu.list[[1]]-X.pred.list[[1]]%*%Beta-Eta)
      old.logL=as.numeric(0.5*tau.epsilon*Old.resid.ssq)
      New.resid.ssq=0
      Dist.pdf=dnorm(Distances,0,sqrt(1/prop))
      Phi.prop[Data$Which.distances]=Dist.pdf[Data$Dist.entries]
      for(it in 2:t.steps){
        M.prop[[it]]=Diagonal(x=as.vector(1/(Phi.prop%*%Iw[[it]]%*%One)))%*%Phi.prop%*%Iw[[it]]
        New.resid.ssq=New.resid.ssq+crossprod(Mu.list[[it]]-M.prop[[it]]%*%Mu.list[[it-1]]-Gamma[it])        
      }
      if(runif(1)<exp(new.logL+(Prior.pars$a.tau.d-1)*log(prop)-Prior.pars$b.tau.d*prop-(old.logL+(Prior.pars$a.tau.d-1)*log(tau.d)-Prior.pars$b.tau.d*tau.d))){  #accounting for kernel of gamma prior
        tau.d=prop
        M=M.prop
        Accept$tau.d=Accept$tau.d+1
      }
    }
    
    #update kernel weights/REs for spatial model
    Dat.minus.Exp=Mu.list[[1]]-X.pred.list[[1]]%*%Beta
    V.eta.inv <- cross.K*tau.epsilon+tau.eta*I.knot
    M.eta <- solve(V.eta.inv,tau.epsilon*K.t%*%Dat.minus.Exp)
    Alpha <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(n.knots,0,1))
    Eta=Data$K%*%Alpha
    if(DEBUG)Eta=Eta*0
    
    #update tau.eta
    tau.eta <- rgamma(1, n.knots*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Alpha,Alpha)*0.5) + Prior.pars$b.eta)    
    
    #update AR1 temporal autocorrelation parameter (sigma)
    tau.sigma <- crossprod(Gamma[2:(t.steps-1)])
    mu.sigma <- 1/tau.sigma * crossprod(Gamma[2:(t.steps-1)],Gamma[3:t.steps])
    temp=rnorm(1,mu.sigma,sqrt(1/(tau.sigma*tau.gamma)))
    if(temp<1 & temp>0)sigma=temp
    
    #update tau.gamma
    Delta=c(Gamma[2],Gamma[3:t.steps]-sigma*Gamma[2:(t.steps-1)])
    tau.gamma <- rgamma(1,0.5*(t.steps-1)+Prior.pars$a.gamma, as.numeric(crossprod(Delta)*0.5) + Prior.pars$b.gamma)
    
    #reset some objects needed for KF
    diag(Q)=c(rep(tau.epsilon^{-1},S),1/tau.gamma)
    for(it in 1:t.steps){
      MM[[it]][1:S,1:S]=M[[it]]
      MM[[it]][S+1,S+1]=sigma
      diag(R[[it]])=1/tau.varepsilon
    }
    
    
       
    #adapt proposal distributions if Control$adapt=TRUE
    if(Control$adapt==TRUE & iiter%%100==0){
      Diff=Accept$Z-Accept.old$Z
      for(iZ in 1:t.steps){
        if(Diff[iZ]<30)Control$MH.Z[iZ]=Control$MH.Z[iZ]*.95
        if(Diff[iZ]>40)Control$MH.Z[iZ]=Control$MH.Z[iZ]*1.0526
      }
      Diff=Accept$Beta-Accept.old$Beta
      for(ibeta in 1:n.beta){
        if(Diff[ibeta]<30)Control$MH.beta[ibeta]=Control$MH.beta[ibeta]*.95
        if(Diff[ibeta]>40)Control$MH.beta[ibeta]=Control$MH.beta[ibeta]*1.0526
      }
      diff=Accept$tau.d-Accept.old$tau.d
      if(diff<30)Control$MH.tau.d=Control$MH.tau.d*.95
      if(diff>40)Control$MH.tau.d=Control$MH.tau.d*1.0526
      Accept.old=Accept
    }
    
    #store MCMC values when appropriate (including predictions)
    if(iiter>Control$burnin & iiter%%Control$thin==0){
      Beta.mc[,(iiter-Control$burnin)/Control$thin]=Beta
      Tau.varepsilon.mc[(iiter-Control$burnin)/Control$thin]=tau.varepsilon
      Tau.epsilon.mc[(iiter-Control$burnin)/Control$thin]=tau.epsilon
      Tau.d.mc[(iiter-Control$burnin)/Control$thin]=tau.d
      Tau.eta.mc[(iiter-Control$burnin)/Control$thin]=tau.eta
      Tau.gamma.mc[(iiter-Control$burnin)/Control$thin]=tau.gamma
      Sigma.mc[(iiter-Control$burnin)/Control$thin]=sigma
      if(Control$predict==TRUE){ #make predictions
        #simulate mu
        #Mu[Which.no.obs]=rnorm(n.no,X.no.obs%*%Beta+Eta[Sites.no.obs]+Gamma[Times.no.obs],sqrt(1/tau.epsilon))
        #posterior predictions of abundance across landscape
        if(iiter==50){
          crap=1
        }
        if(iiter==950){
          crap=1
        }
          Pred.mc[,,(iiter-Control$burnin)/Control$thin]=rpois(S*t.steps,exp(Log.area.adjust+Z))    
        #}
      }
    }
  } 
  Out=list(MCMC=list(tau.eta=Tau.eta.mc,tau.gamma=Tau.gamma.mc,tau.d=Tau.d.mc,Beta=Beta.mc,tau.epsilon=Tau.epsilon.mc,tau.varepsilon=Tau.varepsilon.mc,sigma=Sigma.mc,Pred=Pred.mc),Accept=Accept,Control=Control)
  Out
}


