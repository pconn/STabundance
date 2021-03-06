#' function to perform Bayesian analysis of count data using spatio-temporal interactive log Gaussian Cox process model using dynamic process convolution
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
#'  "a.eta" alpha param for gamma prior precision of spatio-temporal model
#'  "b.eta" beta param for gamma prior precision of spatio-temporal process
#'  "beta0.tau.rw2" prior precision for knot intercepts in rw2 model
#'  "beta1.tau.rw2" prior precision for knot slopes in rw2 model
#'   If provided, all values for Prior.pars need to be specified
#' @return returns a list with the following objecs: 
#' 	MCMC: A list object containing posterior samples;
#'  Accept: A list object indicating the number of proposals that were accepted for parameters updated via Metropolis-Hastings;
#'  Control: A list object giving MCMC tuning parameters (which are updated if the 'adapt' alorithm is used)
#' @export
#' @import Matrix
#' @keywords abundance, mcmc, spatio-temporal model, spatial prediction
#' @author Paul B. Conn \email{paul.conn@@noaa.gov} 
mcmc_STPC<-function(model,Data,Prior.pars=NULL,Control,Area.adjust=NULL){
  
  S=nrow(Data$Grid[[1]])
  t.steps=length(Data$Grid)
  n.obs=nrow(Data$Count.data)
  if(is.null(Prior.pars)){
    Prior.pars=list(beta.tau=0.01,
                    a.eps=1,
                    b.eps=0.01,
                    a.eta=1,
                    b.eta=0.01,
                    beta0.tau.rw2=1,
                    beta1.tau.rw2=10
    )
  }
  if(is.null(Area.adjust))Area.adjust=rep(1,S)
  if(sum(Area.adjust<=0)>0)cat("ERROR: sample unit areas must be >0")

  #sort count data by time and cell so that some of the list to vector code will work right
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
  X.pred=stack_list(X.pred.list)
  n.beta=ncol(X.pred)
  Which.obs=(Data$Count.data[,"Time"]-1)*S+(Data$Count.data[,"Cell"])
  Which.no.obs=c(1:nrow(X.pred))
  Which.no.obs=Which.no.obs[-Which.obs]
  X.obs=matrix(X.pred[Which.obs,],n.obs,n.beta)
  n.no=length(Which.no.obs)
  X.no.obs=matrix(X.pred[Which.no.obs,],n.no,n.beta)
  XpXinv=solve(crossprod(X.obs))
  XpXinvXp=XpXinv%*%t(X.obs)
    
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
  Tau.epsilon.mc=rep(0,mcmc.length)
  Tau.eta.mc=Tau.epsilon.mc
  Pred.mc=array(0,dim=c(S,t.steps,mcmc.length))
  Accept=rep(0,nrow(Data$Count.data))
  Accept.old=Accept
  
  #setup process conv/RW2 space-time model
  Q1=linear_adj_RW2(t.steps) #precision matrix
  Q=Q1
  n.knots=ncol(Data$K)
  K=matrix(0,S*t.steps,n.knots*t.steps)
  for(i in 2:n.knots)Q=bdiag(Q,Q1)  #precision matrix a block diagonal matrix
  #for(i in 2:t.steps)K=bdiag(K,Data$K) 
  #rearrage K to grab right elements of alpha
  cur.row=0
  for(it in 1:t.steps){
    for(iknot in 1:n.knots){
      K[(cur.row+1):(cur.row+S),(iknot-1)*t.steps+it]=Data$K[,iknot]
    }
    cur.row=cur.row+S
  }
  K=Matrix(K)
  Q.t=t(Q)
  Alpha=rrw(tau.eta*Q) #initial values for space-time random effects
  Eta=K%*%Alpha
  K.t=t(K)
  cross.K=crossprod(K,K)
  #K.obs=K[Which.obs,]
  #cross.K.obs=crossprod(K.obs,K.obs)
  #K.obs.t=t(K.obs)
  Diag=diag(nrow(cross.K))*0.001
  X.rw2=matrix(0,n.knots*t.steps,2*n.knots)
  for(iknot in 1:n.knots){
    X.rw2[((iknot-1)*t.steps+1):((iknot-1)*t.steps+t.steps),iknot]=1
    X.rw2[((iknot-1)*t.steps+1):((iknot-1)*t.steps+t.steps),n.knots+iknot]=c(1:t.steps)
  }
  X.rw2.t=t(X.rw2)
  A=solve(crossprod(X.rw2),X.rw2.t) #for conditioning by kriging
  A.t=t(A)
  KX.rw2=K%*%X.rw2
  #KX.rw2.obs=K.obs%*%X.rw2
  KX.rw2.t=t(KX.rw2)
  #KX.rw2.obs.t=t(KX.rw2.obs)
  #KXpX.rw2.obs=crossprod(KX.rw2.obs)
  KXpX.rw2=crossprod(KX.rw2)
  Sigma.inv.rw2=Matrix(diag(c(rep(Prior.pars$beta0.tau.rw2,n.knots),rep(Prior.pars$beta1.tau.rw2,n.knots))))
  #KXpXinv.rw2=solve(crossprod(KX.rw2))
  #KXpXinvXp.rw2=KXpXinv.rw2%*%t(KX.rw2)
  Beta.rw2=rep(0,ncol(X.rw2))
  
  for(iiter in 1:Control$iter){
    if(iiter%%1000==0)cat(paste('iteration ',iiter,' of ',Control$iter,'\n'))
    #update mu (sampled cells)
    Mu.pred=X.obs%*%Beta+Eta[Which.obs]
    sd=sqrt(1/tau.epsilon)
    full.cond.old=dnorm(Mu[Which.obs],Mu.pred,sd,log=1)+dpois(Data$Count.data[,"Count"],exp(Log.offset+Mu[Which.obs]),log=1)
    Prop=Mu[Which.obs]+runif(n.obs,-Control$MH.mu,Control$MH.mu)
    full.cond.new=dnorm(Prop,Mu.pred,sd,log=1)+dpois(Data$Count.data[,"Count"],exp(Log.offset+Prop),log=1)
    I.accept=(runif(n.obs)<exp(full.cond.new-full.cond.old))
    Mu[Which.obs[I.accept==1]]=Prop[I.accept==1]
    Accept=Accept+I.accept
 
    
    Mu[Which.no.obs]=rnorm(n.no,X.no.obs%*%Beta+Eta[Which.no.obs],sqrt(1/tau.epsilon))
        
    #update precision for exchangeable errors - note dependence on Mu.pred (defined in previous mu update)
    if(Control$fix.tau.epsilon==FALSE){
      Mu.pred=X.obs%*%Beta+Eta[Which.obs]
      Diff=Mu[Which.obs]-Mu.pred
      tau.epsilon <- rgamma(1,n.obs/2 + Prior.pars$a.eps, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.eps)
    }

    #update beta
    Beta=t(rmvnorm(1,XpXinvXp%*%(Mu[Which.obs]-Eta[Which.obs]),XpXinv/(tau.epsilon+Prior.pars$beta.tau)))
    
    #update kernel weights/REs for spatial model
    #first, update intercepts and slopes for RW2 models for each knot 
    Dat.minus.Exp=Mu-X.pred%*%Beta-K%*%Alpha
    V.inv.rw2=KXpX.rw2*tau.epsilon+Sigma.inv.rw2
    Beta.rw2=t(rmvnorm(1,solve(V.inv.rw2,KX.rw2.t%*%Dat.minus.Exp*tau.epsilon),as.matrix(solve(V.inv.rw2))))
    #now, update alpha and standardize so has mean and slope of 0 
    Dat.minus.Exp=Mu-X.pred%*%Beta-KX.rw2%*%Beta.rw2
    V.eta.inv <- cross.K*tau.epsilon+tau.eta*Q
    M.eta <- solve(V.eta.inv,tau.epsilon*K.t%*%Dat.minus.Exp)
    Alpha <- M.eta + solve(chol(V.eta.inv), rnorm(length(M.eta),0,1))
    #center using eq 2.30 of Rue and Held
    #Alpha=Alpha-mean(Alpha)
    Alpha=Alpha-V.eta.inv %*% A.t %*% solve(A %*% V.eta.inv %*% A.t,A%*%Alpha)
    
    Eta=K%*%(X.rw2%*%Beta.rw2+Alpha)
    #Eta=K%*%Alpha
    #update tau.eta
    tau.eta <- rgamma(1, (length(M.eta)-2)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Alpha,Q %*% Alpha)*0.5) + Prior.pars$b.eta) 
        
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
      if(Control$predict==TRUE){ #make predictions
        #simulate mu
        Mu[Which.no.obs]=rnorm(n.no,X.no.obs%*%Beta+Eta[Which.no.obs],sqrt(1/tau.epsilon))
        #posterior predictions of abundance across landscape
        Pred.mc[,,(iiter-Control$burnin)/Control$thin]=rpois(S*t.steps,exp(Log.area.adjust+Mu))    
      }
      if(sum(is.na(Pred.mc[,,(iiter-Control$burnin)/Control$thin]))>0){
        temp=1
      }
    }
  } 
  Out=list(MCMC=list(tau.eta=Tau.eta.mc,Beta=Beta.mc,tau.epsilon=Tau.epsilon.mc,Pred=Pred.mc),Accept=Accept,Control=Control)
  Out
}

