#' function to perform Bayesian analysis of count data using a spatio-temporal multinomial model
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
#'  "MH.omega" A vector providing standard deviations for Langevin-Hastings proposals (length = number of time steps)
#'  "MH.N" Standard deviation for total abundance MH updates
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
mcmc_CPIF<-function(model,Data,Prior.pars=NULL,Control,Area.adjust=NULL){
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
                    beta1.tau.rw2=10)
  }
  if(is.null(Area.adjust))Area.adjust=rep(1,S)
  
  if(sum(Area.adjust<=0)>0)cat("ERROR: sample unit areas must be >0")
  
  #sort count data by time and cell so that some of the list to vector code will work right
  Data$Count.data=Data$Count.data[order(Data$Count.data$Time,Data$Count.data$Cell),]
  row.names(Data$Count.data)=c(1:nrow(Data$Count.data))
  
  P=Data$Count.data[,"AreaSurveyed"]  
  Count=Data$Count.data[,"Count"]
  Area.adjust=rep(Area.adjust,t.steps)
  Log.area.adjust=log(Area.adjust)
  Ct=rep(0,t.steps)
  for(it in 1:t.steps){
    Which.Ct.gt0=which(Data$Count.data[,"Time"]==it)
    if(length(Which.Ct.gt0>0))Ct[it]=sum(Count[Which.Ct.gt0])
  }
  log.Mtp1=log(max(Ct)) #lower bound on log(abundance)
  
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

  ###### Declare MCMC structures, set initial values
  mcmc.length=(Control$iter-Control$burnin)/Control$thin
  Beta=optim(rep(0,n.beta),fn=multinom_logL_inits,Counts=Count,P=P,X=X.obs)$par
  Beta.mc=matrix(0,n.beta,mcmc.length)
  N.mc=rep(0,mcmc.length)
  log.N=log(S*sum(Data$Count.data[,"Count"])/sum(P))
  #tau.epsilon=runif(1,10,100)
  tau.epsilon=100
  tau.eta=runif(1,10,100)
  Tau.epsilon.mc=rep(0,mcmc.length)
  Tau.eta.mc=Tau.epsilon.mc
  Pred.mc=array(0,dim=c(S,t.steps,mcmc.length))
  Accept=rep(0,t.steps)
  Accept.old=Accept
  accept.N=accept.N.old=0
  
  
  #provide initial estimate of Omega
  n.ST=S*t.steps
  Omega.pred=X.pred%*%Beta+Log.area.adjust
  Omega=Omega.pred+rnorm(n.ST,0,1/sqrt(tau.epsilon))
  One=matrix(1,t.steps,1)
  Y=Matrix(0,n.ST,t.steps)
  for(it in 1:t.steps)Y[((it-1)*S+1):(it*S),it]=1
  Omega.exp=exp(Omega)
  D=Diagonal(x=Omega.exp)
  Pi=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Omega.exp
  Pi.mat=Diagonal(x=as.numeric(Pi))%*%Y
  Pi.obs=Pi.mat[Which.obs,]*P
  Pi.obs.stnd=Pi.obs%*%Diagonal(x=1/colSums(Pi.obs))
  
  #convert some of these to list format
  Count.list=Omega.list=Pi.list=Dt.old=vector('list',t.steps)
  for(it in 1:t.steps){
    if(Nt.obs[it]>0){
      Count.list[[it]]=Count[which(Data$Count.data[,"Time"]==it)]
      Omega.list[[it]]=Omega[Which.obs[which(Data$Count.data[,"Time"]==it)]]
      Pi.list[[it]]=Pi.obs.stnd[Time.indices[it,1]:Time.indices[it,2],it]
      Dt.old[[it]]=0.5*Control$MH.omega[it]^2*d_logP_omega(Omega.list[[it]],Count.list[[it]],Mu=Omega.pred[Which.obs[Time.indices[it,1]:Time.indices[it,2]]],tau=tau.epsilon)
    }
  }
    
  #setup process conv/RW2 space-time model
  n.knots=ncol(Data$K)
  Q1=linear_adj_RW2(t.steps) #precision matrix
  Q=Q1
  for(i in 2:n.knots)Q=bdiag(Q,Q1)  #precision matrix a block diagonal matrix
  Q.t=t(Q)
  K=matrix(0,S*t.steps,n.knots*t.steps)
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
  Alpha=rrw(tau.eta*Q) #initial values for space-time random effects
  #Eta=K%*%Alpha
  Eta=rep(0,n.ST) #start Eta at 0
  #K.obs=K[Which.obs,]
  #K.obs.t=t(K.obs)
  K.t=t(K)
  cross.K<-crossprod(K,K)
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
  KX.rw2.t=t(KX.rw2)
  KXpX.rw2=crossprod(KX.rw2)
  Sigma.inv.rw2=Matrix(diag(c(rep(Prior.pars$beta0.tau.rw2,n.knots),rep(Prior.pars$beta1.tau.rw2,n.knots))))
  Beta.rw2=rep(0,ncol(X.rw2))
  
  for(iiter in 1:Control$iter){
    cat(paste('iter ',iiter,'\n'))
    if(iiter%%1000==0)cat(paste('iteration ',iiter,' of ',Control$iter,'\n'))

    #update total abundance
    log.N.prop=log.N+rnorm(1,0,Control$MH.N)
    if(log.N.prop>log.Mtp1){
      Cur.p=colSums(Pi.obs)
      N.old=exp(log.N)
      N.prop=exp(log.N.prop)
      old.lik=t.steps*lfactorial(N.old)-sum(lfactorial(N.old-Ct)-(N.old-Ct)*log(1-Cur.p))
      new.lik=t.steps*lfactorial(N.prop)-sum(lfactorial(N.prop-Ct)-(N.prop-Ct)*log(1-Cur.p))
      if(runif(1)<exp(new.lik-old.lik)){
        accept.N=accept.N+1
        log.N=log.N.prop
      }
    }
    
    #update omega (sampled cells)
    Omega.pred=X.obs%*%Beta+Eta[Which.obs]+Log.area.adjust[Which.obs]
    sd=sqrt(1/tau.epsilon)
    for(it in 1:t.steps){
      if(Nt.obs[it]>0){
        Prop=Omega.list[[it]]+Dt.old[[it]]+rnorm(Nt.obs[it],0,Control$MH.omega[it])
        Prop.exp=exp(Prop)
        Pi.prop=Prop.exp/sum(Prop.exp)
        post.new=sum(dnorm(Prop,Omega.pred[Time.indices[it,1]:Time.indices[it,2]],sd,log=TRUE))+sum(Count.list[[it]]*log(Pi.prop))
        post.old=sum(dnorm(Omega.list[[it]],Omega.pred[Time.indices[it,1]:Time.indices[it,2]],sd,log=TRUE))+sum(Count.list[[it]]*log(Pi.list[[it]]))
        Temp=Prop-Omega.list[[it]]-Dt.old[[it]]
        jump.old.to.new=-0.5/Control$MH.omega[it]^2*Temp%*%Temp
        Dstar=0.5*Control$MH.omega[it]^2*d_logP_omega(Prop,Count.list[[it]],Mu=Omega.pred[Time.indices[it,1]:Time.indices[it,2]],tau=tau.epsilon)
        Temp=Omega.list[[it]]-Prop-Dstar
        jump.new.to.old=-0.5/Control$MH.omega[it]^2*Temp%*%Temp
        if(runif(1)<exp(post.new-post.old+jump.new.to.old-jump.old.to.new)){
          Omega.list[[it]]=Prop
          Pi.list[[it]]=Pi.prop
          Accept[it]=Accept[it]+1
          Dt.old[[it]]=Dstar
        }  
      }
    }
    Omega[Which.obs]=unlist(Omega.list)
    #Omega[Which.obs]=stack_list_vector(Omega.list) #convert back to vector format
    #simulate omega for unobserved times and places
    Omega[Which.no.obs]=rnorm(n.no,X.no.obs%*%Beta+Eta[Which.no.obs]+Log.area.adjust[Which.no.obs],sd)
    #compute Pi.obs (need for abundance updates)
    Omega.exp=exp(Omega)
    D=Diagonal(x=Omega.exp)
    Pi=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Omega.exp
    Pi.mat=Diagonal(x=as.numeric(Pi))%*%Y
    Pi.obs=Pi.mat[Which.obs,]*P

    #plot_N_map(1,matrix(Omega[1:400,1],400,1),Grid=Data$Grid,leg.title="Abundance")
    
    #update precision for exchangeable errors
    if(Control$fix.tau.epsilon==FALSE){
      Omega.pred=X.obs%*%Beta+Eta[Which.obs]+Log.area.adjust[Which.obs]
      Diff=Omega[Which.obs]-Omega.pred
      tau.epsilon <- rgamma(1,n.obs/2 + Prior.pars$a.eps, as.numeric(crossprod(Diff,Diff))*0.5 + Prior.pars$b.eps)
    }    
    #update beta
    Beta=t(rmvnorm(1,XpXinvXp%*%(Omega[Which.obs]-Eta[Which.obs]-Log.area.adjust[Which.obs]),XpXinv/(tau.epsilon+Prior.pars$beta.tau)))
    
    #update kernel weights/REs for space-time model
    #first update mean and slope for each knot
    Dat.minus.Exp=Omega-X.pred%*%Beta-Log.area.adjust-K%*%Alpha
    V.inv.rw2=KXpX.rw2*tau.epsilon+Sigma.inv.rw2
    Beta.rw2=t(rmvnorm(1,solve(V.inv.rw2,KX.rw2.t%*%Dat.minus.Exp*tau.epsilon),as.matrix(solve(V.inv.rw2))))
    #now update alpha and correct for mean=0 and slope=0 constraints
    Dat.minus.Exp=Omega-X.pred%*%Beta-Log.area.adjust-KX.rw2%*%Beta.rw2
    V.eta.inv <- cross.K*tau.epsilon+tau.eta*Q
    M.eta <- solve(V.eta.inv,tau.epsilon*K.t%*%Dat.minus.Exp)
    Alpha <- M.eta + solve(chol(V.eta.inv), rnorm(length(M.eta),0,1))
    Alpha=Alpha-V.eta.inv %*% A.t %*% solve(A %*% V.eta.inv %*% A.t,A%*%Alpha)    
    Eta=K%*%(X.rw2%*%Beta.rw2+Alpha)
 
    #update tau.eta
    tau.eta <- rgamma(1, length(M.eta)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Alpha,Q %*% Alpha)*0.5) + Prior.pars$b.eta)    
    #tau.eta=100
    
    #adapt proposal distributions if Control$adapt=TRUE
    if(Control$adapt==TRUE & iiter%%100==0){
      if((accept.N-accept.N.old)<30)Control$MH.N=Control$MH.N*0.95
      if((accept.N-accept.N.old)>40)Control$MH.N=Control$MH.N*1.0526
      accept.N.old=accept.N
      Diff=Accept-Accept.old
      for(iomega in 1:t.steps){
        if(Diff[iomega]<30)Control$MH.omega[iomega]=Control$MH.omega[iomega]*.95
        if(Diff[iomega]>40)Control$MH.omega[iomega]=Control$MH.omega[iomega]*1.0526
      }
      Accept.old=Accept
    }
    
    #store MCMC values when appropriate (including predictions)
    if(iiter>Control$burnin & iiter%%Control$thin==0){
      Beta.mc[,(iiter-Control$burnin)/Control$thin]=Beta
      N.mc[(iiter-Control$burnin)/Control$thin]=exp(log.N)
      Tau.epsilon.mc[(iiter-Control$burnin)/Control$thin]=tau.epsilon
      Tau.eta.mc[(iiter-Control$burnin)/Control$thin]=tau.eta
      if(Control$predict==TRUE){ #make predictions
        #simulate omega
        #Omega[Which.no.obs]=rnorm(n.no,X.no.obs%*%Beta+Eta[Which.no.obs]+Log.area.adjust[Which.no.obs],sqrt(1/tau.epsilon))
        #posterior predictions of abundance across landscape
        #Omega.exp=exp(Omega)
        #diag(D)=Omega.exp
        #Pi=t(t(One)%*%solve((t(Y)%*%D%*%Y),t(Y)))*Omega.exp
        #Pi.mat=Diagonal(x=as.numeric(Pi))%*%Y
        for(it in 1:t.steps)Pred.mc[,it,(iiter-Control$burnin)/Control$thin]=rmultinom(1,exp(log.N),Pi.mat[(S*(it-1)+1):(S*it),it])    
      }
    }
  } 
  Out=list(MCMC=list(tau.eta=Tau.eta.mc,Beta=Beta.mc,N=N.mc,tau.epsilon=Tau.epsilon.mc,Pred=Pred.mc),Accept=list(Omega=Accept,N=accept.N),Control=Control)
  Out
}

