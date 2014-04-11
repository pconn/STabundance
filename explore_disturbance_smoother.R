#test out kalman smoother forward backward algorithm
library(mvtnorm)
library(Matrix)
t.steps=100
S=10
Y=matrix(0,S,t.steps)
Y[,1]=2+rnorm(S)
sigma=0.95
for(it in 2:t.steps){
  Y[,it]=sigma*Y[,it-1]+rnorm(S)
}
Z=Y+rnorm(S*t.steps)


M=sigma*diag(S)
H=diag(S)
R=diag(S)
Q=diag(S)
Z.list=vector("list",t.steps)
Nu.list=Z.list
Finv=Z.list
G=Z.list
L=G
A=L
R.disturb=L
R.disturb[[t.steps]]=Matrix(0,S,1)
P=Z.list
Hat=P
H.t=H

Z.true=Z
Y.true=Y
P[[1]]=1*diag(S)

tau.epsilon=100
tau.varepsilon=100000
diag(Q)=1/sqrt(tau.epsilon)
diag(R)=1/sqrt(tau.varepsilon)

W.plus.z=matrix(rnorm(S*t.steps,0,1/sqrt(tau.varepsilon)),S,t.steps)
W.plus.y=matrix(rnorm(S*t.steps,0,1/sqrt(tau.epsilon)),S,t.steps)
Y.plus=matrix(0,S,t.steps)
Z.plus=Y.plus


#update mu using disturbance filter/smoother algorithm of Durbin and Koopman 2002
ptm=proc.time()
Y.plus[,1]=2+rnorm(S) #holding Y-plus
for(it in 2:t.steps)Y.plus[,it]=as.numeric(M%*%Y.plus[,it-1]+W.plus.y[,it])
Z.plus=Y.plus+W.plus.z
Z.star=Z-Z.plus
#now use Kalman filter and disturbance smoother to obtain \hat{Z}* = E(Y|Z*) 
#Init.pred=c(X.pred.list[[1]]%*%Beta+Eta,0)
Init.pred=rep(0,S) #since running the KF on z-z.plus, think this needs to be zero
Nu.list[[1]]=Z.star[1:S]-H%*%Init.pred
P[[1]]=Q
Finv[[1]]=solve(H%*%P[[1]]%*%H.t+R)
G[[1]]=M%*%(P[[1]]%*%(H.t%*%Finv[[1]]))
L[[1]]=M-G[[1]]%*%H
A[[1]]=Init.pred
for(it in 2:t.steps){  #forward
  P[[it]]=tcrossprod(M%*%P[[it-1]],L[[it-1]])+Q
  Finv[[it]]=solve(H%*%P[[it]]%*%H.t+R)
  G[[it]]=M%*%(P[[it]]%*%(H.t%*%Finv[[it]]))
  A[[it]]=M%*%A[[it-1]]+G[[it-1]]%*%Nu.list[[it-1]]
  Nu.list[[it]]=Z.star[((it-1)*S+1):(it*S)]-H%*%A[[it]]
  L[[it]]=M-G[[it]]%*%H
}
for(it in t.steps:2)R.disturb[[it-1]]=H.t%*%(Finv[[it]]%*%Nu.list[[it]])+crossprod(L[[it]],R.disturb[[it]])
R.disturb.0=H.t%*%(Finv[[1]]%*%Nu.list[[1]])+crossprod(L[[1]],R.disturb[[1]])
Hat[[1]]=A[[1]]+P[[1]]%*%R.disturb.0
for(it in 1:(t.steps-1)){
  Hat[[it+1]]=M%*%Hat[[it]]+Q%*%R.disturb[[it]]
  Y[,it]=as.vector(Hat[[it]]+Y.plus[,it])
}
Y[,t.steps]=as.vector(Hat[[t.steps]]+Y.plus[,t.steps])
cat(proc.time()-ptm)


Y.hat=matrix(0,S,t.steps)
for(it in 1:t.steps)
plot(Y.true[3,])
lines(Y[3,],col='red')
lines(Z[3,],col='blue')


plot(Y[1,],col='red',type="l")
lines(Z[1,],col='blue')