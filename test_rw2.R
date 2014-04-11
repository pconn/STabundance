#test RW2 model
require(Matrix)
require(mvtnorm)

#1 generate data

source("./STabundance/R/sim_funcs.R")
source("./STabundance/R/util_funcs.R")
n.knots=4
S=9
t.steps=10
Q1=linear_adj_RW2(t.steps) #precision matrix
Q=Q1

for(i in 2:n.knots)Q=bdiag(Q,Q1)  #precision matrix a block diagonal matrix

Prior.pars=list(beta.tau=0.01,
                a.eps=1,
                b.eps=0.01,
                a.eta=1,
                b.eta=0.01,
                beta.tau.rw2=1
)

tau.eta=5
tau.epsilon=100  #treated as known here

Alpha=rrw(tau.eta*Q)
Alpha.mat.t=t(matrix(0,n.knots,t.steps))
Alpha.mat.t[]=Alpha  #make sure to fill by time step first
Alpha.mat=t(Alpha.mat.t)

K=matrix(0,S,n.knots)  #K matrix in spatial dimension only (3 by 3 grid with knots in the middle of corner cells)
K[1,]=c(0,2,2,sqrt(8))
K[2,]=c(1,sqrt(5),1,sqrt(5))
K[3,]=c(2,sqrt(8),0,2)
K[4,]=c(1,1,sqrt(5),sqrt(5))
K[5,]=sqrt(2)
K[6,]=c(sqrt(5),sqrt(5),1,1)
K[7,]=c(2,0,sqrt(8),2)
K[8,]=c(sqrt(5),1,sqrt(5),1)
K[9,]=c(sqrt(8),2,2,0)

K=dnorm(K,0,2)
K=K/rowSums(K)

Mu=matrix(0,S,t.steps)
for(it in 1:t.steps){
  Mu[,it]=K%*%Alpha.mat[,it]
}
Mu=as.vector(Mu)

#now rearrange K into a big matrix to grab the right elements of alpha
K2=matrix(0,S*t.steps,n.knots*t.steps)
cur.row=0
for(it in 1:t.steps){
  for(iknot in 1:n.knots){
    K2[(cur.row+1):(cur.row+S),(iknot-1)*t.steps+it]=K[,iknot]
  }
  cur.row=cur.row+S
}
K.old=K
K=K2
Mu2=K%*%Alpha  #same thing as Mu... so this method of constructing the big K matrix checks out.

Alpha.true=Alpha
set.seed(12345)
#Which.obs=sort(sample(c(1:(S*t.steps)),S*t.steps/2)) #not missing any time steps entirely (for the moment)
Which.obs=c(1:(S*t.steps))

#now set up some quantities for estimation
K=Matrix(K)
Q.t=t(Q)
Alpha=rrw(tau.eta*Q) #initial values for space-time random effects
Eta=K%*%Alpha
K.obs=K[Which.obs,]
cross.K.obs=crossprod(K.obs,K.obs)
K.obs.t=t(K.obs)
Diag=diag(nrow(cross.K.obs))*0.001
X.rw2=matrix(0,n.knots*t.steps,2*n.knots)
for(iknot in 1:n.knots){
  X.rw2[((iknot-1)*t.steps+1):((iknot-1)*t.steps+t.steps),iknot]=1
  X.rw2[((iknot-1)*t.steps+1):((iknot-1)*t.steps+t.steps),n.knots+iknot]=c(1:t.steps)
}
X.rw2.t=t(X.rw2)
A=solve(crossprod(X.rw2),X.rw2.t) #for conditioning by kriging
A.t=t(A)
#KX.rw2=K%*%X.rw2
KX.rw2.obs=K.obs%*%X.rw2
#KX.rw2.t=t(KX.rw2)
KX.rw2.obs.t=t(KX.rw2.obs)
KXpX.rw2.obs=crossprod(KX.rw2.obs)
#KXpX.rw2=crossprod(KX.rw2)
Sigma.inv.rw2=diag(x=Prior.pars$beta.tau.rw2,nrow=2*n.knots)  
#KXpXinv.rw2=solve(crossprod(KX.rw2))
#KXpXinvXp.rw2=KXpXinv.rw2%*%t(KX.rw2)
Beta.rw2=rep(0,ncol(X.rw2))

#generate data
n.obs=length(Which.obs)
Data=rpois(n.obs,lambda=exp(Mu))

#now, conduct MCMC inference on Alpha, Beta.rw2, and tau.eta
n.iter=1000
tau.eta=100
#Alpha=rep(0,length(Alpha))
Mu=Mu+rnorm(length(Mu),0,sqrt(1/tau.epsilon))

tau.eta.mc=rep(0,n.iter)
Beta.rw2.mc=matrix(0,length(Beta.rw2),n.iter)

#alternate matrix for just making the sums of time series for each knot = 0
#to use, uncomment code A, A.t below
DM.mean=matrix(0,n.knots*t.steps,n.knots)
DM.mean[1:10,1]=1
DM.mean[11:20,2]=1
DM.mean[21:30,3]=1
DM.mean[31:40,4]=1
#A=solve(crossprod(DM.mean),t(DM.mean))
#A.t=t(A)

for(iiter in 1:n.iter){
  #update mu
  
  #update beta.rw2
  Dat.minus.Exp=Mu[Which.obs]-K.obs%*%Alpha
  V.inv.rw2=KXpX.rw2.obs*tau.epsilon+Sigma.inv.rw2
  Beta.rw2=t(rmvnorm(1,solve(V.inv.rw2,KX.rw2.obs.t%*%Dat.minus.Exp*tau.epsilon),as.matrix(solve(V.inv.rw2))))
  Beta.rw2.mc[,iiter]=Beta.rw2
  #Beta.rw2=rep(0,length(Beta.rw2))
  
  #update alpha
  Dat.minus.Exp=Mu[Which.obs]-KX.rw2.obs%*%Beta.rw2 
  V.eta.inv <- cross.K.obs*tau.epsilon+tau.eta*Q #+Diag
  M.eta <- solve(V.eta.inv,tau.epsilon*K.obs.t%*%Dat.minus.Exp)
  Alpha <- M.eta + solve(chol(V.eta.inv), rnorm(length(M.eta),0,1))
  #center using eq 2.30 of Rue and Held
  #Alpha=Alpha-mean(Alpha)
  Alpha=Alpha-V.eta.inv %*% A.t %*% solve(A %*% V.eta.inv %*% A.t,A%*%Alpha)
   
  #update tau.eta
  tau.eta <- rgamma(1, (length(Alpha)-2)*0.5 + Prior.pars$a.eta, as.numeric(crossprod(Alpha,Q %*% Alpha)*0.5) + Prior.pars$b.eta) 
  tau.eta.mc[iiter]=tau.eta
  
}





