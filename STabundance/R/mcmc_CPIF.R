#' function to perform Bayesian analysis of count data using a separable spatio-temporal log Gaussian Cox process model
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
#'  "MH.lambda" A vector providing continuous Uniform half-range values for Metropolis-Hastings proposals
#'  "adapt" If true, adapts tuning parameters for nu updates; e.g., one could run a chain with adapt=TRUE, and use adapted MH tuning parameters in a longer analysis (default FALSE)
#'  "fix.tau.epsilon" If TRUE, fixes tau.epsilon to 100
#' @param Prior.pars  An optional list giving prior parameters; these include
#'  "beta.tau" precision for Gaussian prior on regression parameters (default 0.1)
#'  "a.eps" alpha parameter for tau_epsilon~Gamma(alpha,beta) (default = 1.0)
#'  "b.eps" beta parameter for tau_epsilon~Gamma(alpha,beta) (default = 0.01)
#'  "a.eta" alpha param for gamma prior precision of spatio-temporal model
#'  "b.eta" beta param for gamma prior precision of spatio-temporal process
#'  "scale.prior"  If TRUE (default), uses a prior for lambda of 1/lambda
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


mcmc_CPIF<-function(model,Data,Prior.pars=NULL,Control,Area.adjust=NULL){
  require(Matrix)
  require(mvtnorm)
  source("c:/users/paul.conn/git/STabundance/STabundance/R/util_funcs.R")
  
  S=nrow(Data$Grid[[1]])
  t.steps=length(Data$Grid)
  n.obs=nrow(Data$Count.data)
  if(is.null(Prior.pars)){
    Prior.pars=list(beta.tau=0.01,
                    a.eps=1,
                    b.eps=0.01,
                    a.eta=1,
                    b.eta=0.01)
  }
  if(is.null(Area.adjust))Area.adjust=rep(1,S)
  
  Offset=Data$Count.data[,"AreaSurveyed"]  
  Count=Data$Count.data[,"Count"]
  Offset=Offset*Area.adjust[Data$Count.data[,"Cell"]]
  Log.offset=log(Offset)
  Area.adjust=rep(Area.adjust,t.steps)
  Log.area.adjust=log(Area.adjust)
  
  #set up design matrices
  X.pred.list=vector("list",t.steps)
  for(it in 1:t.steps){
    X.pred.list[[it]]=model.matrix(model,data=Data$Grid[[it]]@data)
  }
  All.times=rep(c(1:t.steps),each=S)
  All.sites=rep(c(1:S),t.steps)
  X.pred=stack_list(X.pred.list)
  Which.obs=(Data$Count.data[,"Time"]-1)*S+(Data$Count.data[,"Cell"])
  Which.no.obs=c(1:nrow(X.pred))
  Which.no.obs=Which.no.obs[-Which.obs]
  Times.no.obs=All.times[Which.no.obs]
  Sites.no.obs=All.sites[Which.no.obs]
  X.obs=X.pred[Which.obs,]
  X.no.obs=X.pred[Which.no.obs,]
  n.no=length(Which.no.obs)
  XpXinv=solve(crossprod(X.obs))
  XpXinvXp=XpXinv%*%t(X.obs)
  
  n.beta=ncol(X.obs)
  
  ###### Declare MCMC structures, set initial values
  mcmc.length=(Control$iter-Control$burnin)/Control$thin
  #provide initial estimate of Omega
  n.ST=S*t.steps
  Omega=rep(0,n.ST)
  Omega[Which.obs]=(Data$Count.data[,"Count"]/Offset)/mean((Data$Count.data[,"Count"]/Offset))
  Omega[which(is.na(Omega))]=median(Omega[which(Omega>0)])
  Omega[which(Omega==0)]=min(Omega[which(Omega>0)])
  Omega=log(Omega)
  One=matrix(1,t.steps,1)
  Y=Matrix(0,n.ST,t.steps)
  for(it in 1:t.steps)Y[((it-1)*S+1):(it*S),it]=1
  Omega.exp=exp(Omega)
  D=Diagonal(x=Omega.exp)
  Pi=t(One)%*%solve((t(Y)%*%D%*%Y),t(Y))*Omega.exp
  
  Beta=rep(0,n.beta)
  Beta=c(7.6,-5)
  Beta.mc=matrix(0,n.beta,mcmc.length)
  Lambda.mc=rep(0,mcmc.length)
  log.lambda=log(S*sum(Data$Count.data[,"Count"])/sum(Offset))
  Lambda=exp(log.lambda)*Pi
  #tau.epsilon=runif(1,10,100)
  tau.epsilon=100
  tau.eta=runif(1,10,100)
  Tau.epsilon.mc=rep(0,mcmc.length)
  Tau.eta.mc=Tau.epsilon.mc
  Pred.mc=array(0,dim=c(S,t.steps,mcmc.length))
  Accept=rep(0,ncol(Data$Count.data))
  Accept.old=Accept
  accept.lambda=0
  
  #setup process conv/RW2 space-time model
  Q1=linear_adj_RW2(t.steps) #precision matrix
  Q=Q1
  n.knots=ncol(Data$K)
  K=Data$K
  for(i in 2:n.knots)Q=bdiag(Q,Q1)  #precision matrix a block diagonal matrix
  for(i in 2:t.steps)K=bdiag(K,Data$K) 
  Q.t=t(Q)
  Alpha=rrw(tau.eta*Q) #initial values for space-time random effects
  Eta=K%*%Alpha
  K.obs=K[Which.obs,]
  K.obs.t=t(K.obs)
  cross.K<-crossprod(K.obs,K.obs)
  
  for(iiter in 1:Control$iter){
    #if(iiter%%1000==0)
      cat(paste('iteration ',iiter,' of ',Control$iter,'\n'))

    #update lambda
    log.lambda.prop=log.lambda+runif(1,-Control$MH.lambda,Control$MH.lambda)
    Lambda.prop=exp(log.lambda.prop)*Pi[Which.obs]
    old.lik=sum(dpois(Data$Count.data[,"Count"],Offset*Lambda[Which.obs],log=TRUE))
    new.lik=sum(dpois(Data$Count.data[,"Count"],Offset*Lambda.prop,log=TRUE))
    if(runif(1)<exp(new.lik-old.lik)){
      Lambda[Which.obs]=Lambda.prop
      accept.lambda=accept.lambda+1
      log.lambda=log.lambda.prop
    }
    #log.lambda=log(70000)
    
    #update omega (sampled cells)
    Omega.pred=X.obs%*%Beta+Eta[Which.obs]+Log.area.adjust[Which.obs]
    sd=sqrt(1/tau.epsilon)
    full.cond.old=dnorm(Omega[Which.obs],Omega.pred,sd,log=1)+dpois(Data$Count.data[,"Count"],Offset*Lambda[Which.obs],log=1)
    Prop=Omega[Which.obs]+runif(n.obs,-Control$MH.omega,Control$MH.omega)
    Omega.exp[Which.obs]=exp(Prop)
    diag(D)=Omega.exp
    Pi.prop=t(One)%*%solve((t(Y)%*%D%*%Y),t(Y))*Omega.exp
    full.cond.new=dnorm(Prop,Omega.pred,sd,log=1)+dpois(Data$Count.data[,"Count"],Offset*exp(log.lambda)*Pi.prop[Which.obs],log=1)
    I.accept=(runif(n.obs)<exp(full.cond.new-full.cond.old))
    Omega[Which.obs[I.accept==1]]=Prop[I.accept==1]
    Omega.exp=exp(Omega)
    diag(D)=Omega.exp
    Pi=t(One)%*%solve((t(Y)%*%D%*%Y),t(Y))*Omega.exp
    Lambda[Which.obs]=exp(log.lambda)*Pi[Which.obs]
    Accept=Accept+I.accept
    
    #update beta
    if(iiter>2000)Beta=t(rmvnorm(1,XpXinvXp%*%(Omega[Which.obs]-Eta[Which.obs]),XpXinv/(tau.epsilon+Prior.pars$beta.tau)))
    
    #update precision for exchangeable errors
    #if(Control$fix.tau.epsilon==FALSE & iiter>3000){
    #  Omega.pred=X.obs%*%Beta+Eta[Which.obs]+Log.area.adjust[Which.obs]
    #  Diff=Omega[Which.obs]-Omega.pred
    #  tau.epsilon <- rgamma(1,n.obs/2 + Prior.pars$a.eps, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.eps)
    #}
    
    #update kernel weights/REs for spatial model
    if(iiter>3000){
      Dat.minus.Exp=Omega[Which.obs]-X.obs%*%Beta
      V.eta.inv <- cross.K*tau.epsilon+tau.eta*Q
      M.eta <- solve(V.eta.inv,tau.epsilon*K.obs.t%*%Dat.minus.Exp)
      Alpha <- M.eta + solve(Cholesky(V.eta.inv), rnorm(length(M.eta),0,1))
      Alpha=Alpha-mean(Alpha)  #so intercept of fixed effects identifiable... note there's still implicit linear trend parameters for each alpha t-series
      Eta=K%*%Alpha
      #Eta=0*Eta
      #update tau.eta
      tau.eta <- rgamma(1, length(M.eta)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Alpha,Q %*% Alpha)*0.5) + Prior.pars$b.eta)    
      #tau.eta=100
    }
    
    #adapt proposal distributions if Control$adapt=TRUE
    if(Control$adapt==TRUE & iiter%%100==0){
      if((accept.lambda-accept.lambda.old)<30)Control$MH.labmda=Control$MH.lambda*0.95
      if((accept.lambda-accept.lambda.old)>40)Control$MH.labmda=Control$MH.lambda*1.0526
      accept.lambda.old=accept.lambda
      Diff=Accept-Accept.old
      for(iomega in 1:n.obs){
        if(Diff[iomega]<30)Control$MH.omega[iomega]=Control$MH.omega[iomega]*.95
        if(Diff[iomega]>40)Control$MH.omega[iomega]=Control$MH.omega[iomega]*1.0526
      }
      Accept.old=Accept
    }
    
    #store MCMC values when appropriate (including predictions)
    if(iiter>Control$burnin & iiter%%Control$thin==0){
      Beta.mc[,(iiter-Control$burnin)/Control$thin]=Beta
      Lambda.mc[(iiter-Control$burnin)/Control$thin]=exp(log.lambda)
      Tau.epsilon.mc[(iiter-Control$burnin)/Control$thin]=tau.epsilon
      Tau.eta.mc[(iiter-Control$burnin)/Control$thin]=tau.eta
      if(Control$predict==TRUE){ #make predictions
        #simulate omega
        Omega[Which.no.obs]=rnorm(n.no,X.no.obs%*%Beta+Eta[Which.no.obs]+Log.area.adjust[Which.no.obs],sqrt(1/tau.epsilon))
        #posterior predictions of abundance across landscape
        Omega.exp=exp(Omega)
        diag(D)=Omega.exp
        Pi=t(One)%*%solve((t(Y)%*%D%*%Y),t(Y))*Omega.exp
        Pred.mc[,,(iiter-Control$burnin)/Control$thin]=rpois(n.ST,as.numeric(exp(log.lambda)*Pi))    
      }
    }
  } 
  Out=list(MCMC=list(tau.eta=Tau.eta.mc,Beta=Beta.mc,Lambda=Lambda.mc,tau.epsilon=Tau.epsilon.mc,Pred=Pred.mc),Accept=list(Omega=Accept,lambda=accept.lambda),Control=Control)
  Out
}

