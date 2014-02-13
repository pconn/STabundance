#test out kalman smoother forward backward algorithm
library(mvtnorm)
t.steps=100
S=10
Y=matrix(0,S,t.steps)
Y[,1]=rnorm(S)
sigma=0.95
for(it in 2:t.steps){
  Y[,it]=sigma*Y[,it-1]+rnorm(S)
}
Z=Y+rnorm(S*t.steps)


M=sigma*diag(S)
H=diag(S)
R=diag(S)
Q=diag(S)

Z.true=Z
Y.true=Y
Pinit=1000*diag(S)


Ytt=matrix(0,S,t.steps)
Yttmin1=Ytt
Ytt[,1]=Y[,1]
Pttmin1=array(0,dim=c(S,S,t.steps)) #probably eventually want to do this with list of 'Matrix' types 
Ptt=Pttmin1
Ptt[,,1]=Pinit
I=diag(S)
#forward
for(it in 2:t.steps){
  Yttmin1[,it]=M%*%Ytt[,it-1]
  Pttmin1[,,it]=Q+M%*%Ptt[,,it-1]%*%t(M)
  K=Pttmin1[,,it]%*%t(H)%*%solve(t(H)%*%Pttmin1[,,it]%*%H+R)
  Ytt[,it]=Yttmin1[,it]+K%*%(Z[,it]-H%*%Yttmin1[,it])
  Ptt[,,it]=(I-K%*%H)%*%Pttmin1[,,it]
}

Y[,t.steps]=rmvnorm(1,Ytt[,t.steps],Ptt[,,t.steps])  #probably want to use cholesky here
#backward
for(it in (t.steps-1):1){
  J=Ptt[,,it]%*%t(M)%*%solve(Pttmin1[,,it+1])
  Y[,it]=rmvnorm(1,Ytt[,it]+J%*%(Y[,it+1]-Yttmin1[,it+1]),Ptt[,,it]-J%*%Pttmin1[,,it+1]%*%t(J))
}

plot(Y.true[3,])
lines(Y[3,])
lines(Ytt[3,],col='red')
lines(Z[3,],col='blue')
