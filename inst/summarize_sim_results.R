### analyze and visualize 'generic' simulation data



Trans10.sims=c(1:10)*10-5
Trans20.sims= c(1:19)*5

DF=data.frame(matrix(0,1000,9))
colnames(DF)=c("Gen.mod","Transects","Sim.num","Est.mod","MSE","Bias","CV","Coverage","Trend.bias")

Gen.names=c("RS2closed","CPIF","RS2open")
Est.names=c("AST","STPC","CPIF","OPRS")
Trans.names=c("1","5")
ipl=0

st_dev=function(x)sqrt(var(x))

get_coverage=function(Sample,N.true){
  Coverage=rep(0,length(N.true))
  for(it in 1:length(N.true)){
    upper=quantile(Sample[it,],0.95)
    lower=quantile(Sample[it,],0.05)
    if(N.true[it]<upper & N.true[it]>lower)Coverage[it]=1
  }
  mean(Coverage)
}

MSE_N=function(n.bar,N){
  return(sum((N-n.bar)^2))
}

get_avg_abundance_estimator<-function(N){ #N is vector of total abundance by iteration
  Avg.N=optim(par=200000,fn=MSE_N,method="BFGS",N=N)$par
}

for(igen in 1:3){
  for(iest in 1:4){
    for(itrans in 1:2){
      for(isim in 1:100){
        fname=paste("d:/ST_out/ST_out_gen",igen,"_est_",Est.names[iest],"_trans_",itrans,"_sim",isim,".Rdata",sep="")
        fname2=paste("./sim_generic_data/simdata_gen",Gen.names[igen],"_trans",Trans.names[itrans],"_sim",isim,sep='')
        if((Est.names[iest]=="OPRS" & itrans==2 & isim%in%Trans10.sims) |
             Est.names[iest]=="AST" | (Est.names[iest]=="STPC" & isim%in%Trans20.sims) |  
             (Est.names[iest]=="CPIF" & isim%in%Trans20.sims)){
                  
          load(fname)
          ipl=ipl+1
          DF[ipl,"Sim.num"]=isim
          DF[ipl,"Gen.mod"]=Gen.names[igen]
          DF[ipl,"Transects"]=Trans.names[itrans]
          DF[ipl,"Est.mod"]=Est.names[iest]
          
          MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
          N.bar=get_avg_abundance_estimator(MCMC$MCMC$N)
          N.est=apply(MCMC$MCMC$N,1,'mean')
          N.se=apply(MCMC$MCMC$N,1,'st_dev')
          load(fname2)
          N.true=apply(Sim.data$N,2,'sum')
          t.steps=length(N.true)
          trend.est=lm(N.est~c(1:t.steps))$coef
          trend.est=1+trend.est[2]/trend.est[1]*t.steps  #fitted lambda for whole study
          trend.true=lm(N.true~c(1:length(N.true)))$coef
          trend.true=1+trend.true[2]/trend.true[1]*t.steps
          if(igen==1 | igen==2)n.est=N.bar
          if(igen==3)n.est=N.est
          DF[ipl,"Bias"]=mean((n.est-N.true)/N.true)
          DF[ipl,"MSE"]=mean((n.est-N.true)^2)
          DF[ipl,"CV"]=median(N.se/N.est)
          DF[ipl,"Coverage"]=get_coverage(MCMC$MCMC$N,N.true)
          DF[ipl,"Trend.bias"]=(trend.est-trend.true)/trend.true
        }
      }
    }
  }
}

#move points with bias>2.0 into an auxiliary dataset
DF[which(DF[,"Transects"]=="1"),"Transects"]="1 Transect"
DF[which(DF[,"Transects"]=="5"),"Transects"]="5 Transects"
DF[which(DF[,"Gen.mod"]=="CPIF"),"Gen.mod"]="CPUD"
DF[which(DF[,"Gen.mod"]=="RS2closed"),"Gen.mod"]="CPRD"
DF[which(DF[,"Gen.mod"]=="RS2open"),"Gen.mod"]="OPRD"

DF=DF[1:858,]


#DF.trunc.1 = DF[which(DF[,"Bias"]>2.0),]
#DF.trunc.5 = DF[which(DF[,"Bias"]>0.25 & DF[,"Transects"]=="5 Transects"),]

#DF2=DF[-which(DF[,"Bias"]>2.0),]
#DF2=DF2[-which(DF2[,"Bias"]>0.25 & DF2[,"Transects"]=="5 Transects"),]

#DF.trunc.1[,"Bias"]=2.0
#DF.trunc.5[,"Bias"]=0.25
#delete any values from primary 

#remove extreme outlier
#outlier=which(DF[,"Bias"]>10)
#DF2=DF[-outlier,]

#plot proportion relative bias
library(ggplot2)
bias.plot = ggplot(DF,aes(Est.mod,Bias))+geom_boxplot()+facet_grid(Transects~Gen.mod,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")
pdf("bias_generic.pdf")
bias.plot
dev.off()

save(DF,file="sim_data_generic_DF.Rdata")
#Now create tables of results
require(dplyr)

poop = group_by(DF,Est.mod,Transects,Gen.mod)
Sim.generic.table=summarise(poop,Median.bias=median(Bias),Coverage=mean(Coverage),CV=median(CV),MSE=mean(MSE),Trend.bias=median(Trend.bias))

require(xtable)
my.table=xtable(Sim.generic.table)
print(my.table,file="sim_generic_results_table.txt",include.rownames=FALSE)


#analyze bias as a function of delta
crap=c(1:100)*0.04
crap2=c(1:19)*0.04
crap3=c(1:10)*0.04
summary(lm(DF[1:100,"Bias"]~crap))  #AST, CPRD, 1
summary(lm(DF[101:200,"Bias"]~crap)) #AST, CPRD, 5
summary(lm(DF[201:219,"Bias"]~crap2)) #STPC, CPRD, 1
summary(lm(DF[220:238,"Bias"]~crap2)) #STPC, CPRD, 5
summary(lm(DF[239:257,"Bias"]~crap2)) #CPIF, CPRD, 1
summary(lm(DF[258:276,"Bias"]~crap2)) #CPIF, CPRD, 5
summary(lm(DF[277:286,"Bias"]~crap3)) #OPRS, CPRD, 5

summary(lm(DF[287:386,"Bias"]~crap))  #AST, CPUD, 1
summary(lm(DF[387:486,"Bias"]~crap)) #AST, CPUD, 5
summary(lm(DF[487:505,"Bias"]~crap2)) #STPC, CPUD, 1
summary(lm(DF[506:524,"Bias"]~crap2)) #STPC, CPUD, 5
summary(lm(DF[525:543,"Bias"]~crap2)) #CPIF, CPUD, 1
summary(lm(DF[544:562,"Bias"]~crap2)) #CPIF, CPUD, 5
summary(lm(DF[563:572,"Bias"]~crap3)) #OPRS, CPUD, 5

summary(lm(DF[573:672,"Bias"]~crap))  #AST, OPRD, 1
summary(lm(DF[673:772,"Bias"]~crap)) #AST, OPRD, 5
summary(lm(DF[773:791,"Bias"]~crap2)) #STPC, OPRD, 1
summary(lm(DF[792:810,"Bias"]~crap2)) #STPC, OPRD, 5
summary(lm(DF[811:829,"Bias"]~crap2)) #CPIF, OPRD, 1
summary(lm(DF[830:848,"Bias"]~crap2)) #CPIF, OPRD, 5
summary(lm(DF[849:858,"Bias"]~crap3)) #OPRS, OPRD, 5

summary(lm(DF[1:100,"Trend.bias"]~crap))  #AST, CPRD, 1
summary(lm(DF[101:200,"Trend.bias"]~crap)) #AST, CPRD, 5
summary(lm(DF[201:219,"Trend.bias"]~crap2)) #STPC, CPRD, 1
summary(lm(DF[220:238,"Trend.bias"]~crap2)) #STPC, CPRD, 5
summary(lm(DF[239:257,"Trend.bias"]~crap2)) #CPIF, CPRD, 1
summary(lm(DF[258:276,"Trend.bias"]~crap2)) #CPIF, CPRD, 5
summary(lm(DF[277:286,"Trend.bias"]~crap3)) #OPRS, CPRD, 5

summary(lm(DF[287:386,"Trend.bias"]~crap))  #AST, CPUD, 1
summary(lm(DF[387:486,"Trend.bias"]~crap)) #AST, CPUD, 5
summary(lm(DF[487:505,"Trend.bias"]~crap2)) #STPC, CPUD, 1
summary(lm(DF[506:524,"Trend.bias"]~crap2)) #STPC, CPUD, 5
summary(lm(DF[525:543,"Trend.bias"]~crap2)) #CPIF, CPUD, 1
summary(lm(DF[544:562,"Trend.bias"]~crap2)) #CPIF, CPUD, 5
summary(lm(DF[563:572,"Trend.bias"]~crap3)) #OPRS, CPUD, 5

summary(lm(DF[573:672,"Trend.bias"]~crap))  #AST, OPRD, 1
summary(lm(DF[673:772,"Trend.bias"]~crap)) #AST, OPRD, 5
summary(lm(DF[773:791,"Trend.bias"]~crap2)) #STPC, OPRD, 1
summary(lm(DF[792:810,"Trend.bias"]~crap2)) #STPC, OPRD, 5
summary(lm(DF[811:829,"Trend.bias"]~crap2)) #CPIF, OPRD, 1
summary(lm(DF[830:848,"Trend.bias"]~crap2)) #CPIF, OPRD, 5
summary(lm(DF[849:858,"Trend.bias"]~crap3)) #OPRS, OPRD, 5
