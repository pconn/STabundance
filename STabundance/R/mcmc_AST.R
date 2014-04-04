#' function to perform Bayesian analysis of count data using an additive spatio-temporal log Gaussian Cox process model
#' @param model A formula object specifying the linear predictor for abundance intensity
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
#'  "predict" If TRUE (default), calculate posterior predictions across the entire grid
#'  "MH.mu" A vector providing continuous Uniform half-range values for Metropolis-Hastings proposals
#'  "adapt" If true, adapts tuning parameters for nu updates; e.g., one could run a chain with adapt=TRUE, and use adapted MH tuning parameters in a longer analysis (default FALSE)
#'  "fix.tau.epsilon" If TRUE, fixes tau.epsilon to 100
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


mcmc_AST<-function(model,Data,Prior.pars=NULL,Control,Area.adjust=NULL){
  require(Matrix)
  require(mvtnorm)
  
  S=nrow(Data$Grid[[1]])
  t.steps=length(Data$Grid)
  n.obs=nrow(Data$Count.data)
  if(is.null(Prior.pars)){
    Prior.pars=list(beta.tau=0.01,
                    a.eps=1,
                    b.eps=0.01,
                    a.eta=1,
                    b.eta=0.01,
                    a.gamma=1,
                    b.gamma=0.01)
  }
  if(is.null(Area.adjust))Area.adjust=rep(1,S)
  
  Data$Count.data=Data$Count.data[order(Data$Count.data$Time,Data$Count.data$Cell),]
  row.names(Data$Count.data)=c(1:nrow(Data$Count.data))

  Offset=Data$Count.data[,"AreaSurveyed"]  
  Count=Data$Count.data[,"Count"]
  Offset=Offset*Area.adjust[Data$Count.data[,"Cell"]]
  Log.offset=log(Offset)
  Log.area.adjust=log(Area.adjust)
  
  #set up design matrices
  X.pred.list=vector("list",t.steps)
  for(it in 1:t.steps){
    X.pred.list[[it]]=model.matrix(model,data=Data$Grid[[it]]@data)
  }
  All.times=rep(c(1:t.steps),each=S)
  All.sites=rep(c(1:S),t.steps)
  X.pred=stack_list(X.pred.list)
  n.beta=ncol(X.pred)
  Which.obs=(Data$Count.data[,"Time"]-1)*S+(Data$Count.data[,"Cell"])
  Which.no.obs=c(1:nrow(X.pred))
  Which.no.obs=Which.no.obs[-Which.obs]
  Times.no.obs=All.times[Which.no.obs]
  Sites.no.obs=All.sites[Which.no.obs]
  X.obs=matrix(X.pred[Which.obs,],n.obs,n.beta)
  n.no=length(Which.no.obs)
  X.no.obs=matrix(X.pred[Which.no.obs,],n.no,n.beta)
  XpXinv=solve(crossprod(X.obs))
  XpXinvXp=XpXinv%*%t(X.obs)
    
  #Which.not.sampled=c(1:S)
  #Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]
  #n.no=S-n
  #X.no=as.matrix(X.pred[Which.not.sampled,],ncol=ncol(X.obs))
  
  ###### Declare MCMC structures, set initial values
  mcmc.length=(Control$iter-Control$burnin)/Control$thin
  #provide initial estimate of Mu
  n.ST=S*t.steps
  Mu=rep(0,n.ST)
  Mu[Which.obs]=Data$Count.data[,"Count"]/Offset
  Mu[which(Mu==0)]=min(Mu[which(Mu>0)])
  Mu[which(is.na(Mu))]=median(Mu[which(Mu>0)])
  Mu=log(Mu)
  
  Beta=rep(0,n.beta)
  Beta[1]=mean(Mu) #set intercept to reasonable starting value
  Beta.mc=matrix(0,n.beta,mcmc.length)
  tau.epsilon=runif(1,10,100)
  tau.eta=runif(1,10,100)
  tau.gamma=runif(1,10,100)
  Tau.epsilon.mc=rep(0,mcmc.length)
  Tau.eta.mc=Tau.epsilon.mc
  Tau.gamma.mc=Tau.epsilon.mc
  Pred.mc=array(0,dim=c(S,t.steps,mcmc.length))
  Accept=rep(0,nrow(Data$Count.data))
  Accept.old=Accept
  
  #setup RW2 time series model
  QT=linear_adj_RW2(t.steps)  #precision matrix
  QT.t=t(QT)
  Tmp.dat=Data$Count.data
  Tmp.dat[,"Time"]=factor(Tmp.dat[,"Time"],levels=c(1:t.steps))
  XT=model.matrix(~0+Time,data=Tmp.dat)
  XT.t=t(XT)
  cross.XT=crossprod(XT,XT)
  X.trend=matrix(1,t.steps,2)
  X.trend[,2]=c(1:t.steps)
  A=Matrix(solve(crossprod(X.trend,X.trend),t(X.trend)))  #for resolving identifiability - see e.g. eqn 2.30 of Rue and Held
  A.t=t(A)
  n.gamma.ts=t.steps-2 #2 less effective re's due to RW2 model  
  Gamma<-rrw(tau.gamma*QT)
  I.T=diag(t.steps)
  
  #setup process conv model for spatial effects
  #initial log kernel weights are set via an ICAR(10) prcoess
  n.knots=ncol(Data$K)
  Alpha=rnorm(n.knots,0,sqrt(1/tau.eta)) #initial kernel weights/random effects
  Eta=Data$K%*%Alpha
  K.obs=Data$K[Data$Count.data[,"Cell"],]
  K.obs.t=t(K.obs)
  cross.K<-crossprod(K.obs,K.obs)
  I.knot=diag(n.knots)
  
  for(iiter in 1:Control$iter){
    if(iiter%%1000==0)cat(paste('iteration ',iiter,' of ',Control$iter,'\n'))
    #update mu (sampled cells)
    Mu.pred=X.obs%*%Beta+Eta[Data$Count.data[,"Cell"]]+Gamma[Data$Count.data[,"Time"]]
    sd=sqrt(1/tau.epsilon)
    full.cond.old=dnorm(Mu[Which.obs],Mu.pred,sd,log=1)+dpois(Data$Count.data[,"Count"],exp(Log.offset+Mu[Which.obs]),log=1)
    Prop=Mu[Which.obs]+rnorm(n.obs,0,Control$MH.mu)
    full.cond.new=dnorm(Prop,Mu.pred,sd,log=1)+dpois(Data$Count.data[,"Count"],exp(Log.offset+Prop),log=1)
    I.accept=(runif(n.obs)<exp(full.cond.new-full.cond.old))
    Mu[Which.obs[I.accept==1]]=Prop[I.accept==1]
    Accept=Accept+I.accept
    
    #update beta
    Beta=t(rmvnorm(1,XpXinvXp%*%(Mu[Which.obs]-Eta[Data$Count.data[,"Cell"]]-Gamma[Data$Count.data[,"Time"]]),XpXinv/(tau.epsilon+Prior.pars$beta.tau)))
    
    #update precision for exchangeable errors
    if(Control$fix.tau.epsilon==FALSE){
      Mu.pred=X.obs%*%Beta+Eta[Data$Count.data[,"Cell"]]+Gamma[Data$Count.data[,"Time"]]
      Diff=Mu[Which.obs]-Mu.pred
      tau.epsilon <- rgamma(1,0.5*n.obs + Prior.pars$a.eps, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.eps)
    }
    
    #update kernel weights/REs for spatial model
    Dat.minus.Exp=Mu[Which.obs]-X.obs%*%Beta-Gamma[Data$Count.data[,"Time"]]
    V.eta.inv <- cross.K*tau.epsilon+tau.eta*I.knot
    M.eta <- solve(V.eta.inv,tau.epsilon*K.obs.t%*%Dat.minus.Exp)
    Alpha <- M.eta + solve(chol(as.matrix(V.eta.inv)), rnorm(n.knots,0,1))
    Eta=Data$K%*%Alpha
    #update tau.eta
    tau.eta <- rgamma(1, n.knots*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Alpha,Alpha)*0.5) + Prior.pars$b.eta)    
    
    #update REs, precision for time series RW2 ICAR model
    Dat.minus.Exp=Mu[Which.obs]-X.obs%*%Beta-Eta[Data$Count.data[,"Cell"]]
    V.gamma.inv <- tau.epsilon*cross.XT + tau.gamma*QT
    M.eta <- solve(V.gamma.inv,tau.epsilon*XT.t%*%(Dat.minus.Exp))
    Gamma <- M.eta + solve(chol(V.gamma.inv),rnorm(t.steps,0,1))
    #center using eq 2.30 of Rue and Held
    #Gamma <- Gamma - V.gamma.inv %*% A.t %*% solve(A %*% V.gamma.inv %*% A.t,A %*% Gamma)
    Gamma=Gamma-mean(Gamma)  #just center first moment so no confounding with fixed effect intercept  
    #update tau.gamma
    tau.gamma <- rgamma(1,t.steps*0.5 + Prior.pars$a.gamma, as.numeric(crossprod(Gamma,QT %*% Gamma)*0.5) + Prior.pars$b.gamma)
    
    #adapt proposal distributions if Control$adapt=TRUE
    if(Control$adapt==TRUE & iiter%%100==0){
      Diff=Accept-Accept.old
      for(imu in 1:n.obs){
        if(Diff[imu]<30)Control$MH.mu[imu]=Control$MH.mu[imu]*.95
        if(Diff[imu]>40)Control$MH.mu[imu]=Control$MH.mu[imu]*1.0526
      }
      Accept.old=Accept
    }
    
    #store MCMC values when appropriate (including predictions)
    if(iiter>Control$burnin & iiter%%Control$thin==0){
      Beta.mc[,(iiter-Control$burnin)/Control$thin]=Beta
      Tau.epsilon.mc[(iiter-Control$burnin)/Control$thin]=tau.epsilon
      Tau.eta.mc[(iiter-Control$burnin)/Control$thin]=tau.eta
      Tau.gamma.mc[(iiter-Control$burnin)/Control$thin]=tau.gamma
      if(Control$predict==TRUE){ #make predictions
        #simulate mu
        Mu[Which.no.obs]=rnorm(n.no,X.no.obs%*%Beta+Eta[Sites.no.obs]+Gamma[Times.no.obs],sqrt(1/tau.epsilon))
        #posterior predictions of abundance across landscape
        Pred.mc[,,(iiter-Control$burnin)/Control$thin]=rpois(S*t.steps,exp(Log.area.adjust+Mu))    
      }
    }
  } 
  Out=list(MCMC=list(tau.eta=Tau.eta.mc,tau.gamma=Tau.gamma.mc,Beta=Beta.mc,tau.epsilon=Tau.epsilon.mc,Pred=Pred.mc),Accept=Accept,Control=Control)
  Out
}

