### analyze and visualize 'generic' simulation data


t.steps=29
Trans10.sims=c(5,15,25,35,45,55,65,75,81,91)

DF=data.frame(matrix(0,240,8))  #this one holds performance for "average" N
colnames(DF)=c("Gen.mod","Sim.num","Est.mod","MSE","Bias","CV","Coverage","Trend.bias")

Results=array(0,dim=c(240,t.steps,4))  #This one holds performance results for individual years

Gen.names=c("RS2closed","CPIF")
Est.names=c("AST","STPC","CPIF")

st_dev=function(x)sqrt(var(x))

get_coverage=function(Sample,N.true){
  Coverage=rep(0,length(N.true))
  if(length(N.true)>1){
    for(it in 1:length(N.true)){
      upper=quantile(Sample[it,],0.95)
      lower=quantile(Sample[it,],0.05)
      if(N.true[it]<upper & N.true[it]>lower)Coverage[it]=1
    }
  }
  else{
    upper=quantile(Sample,0.95)
    lower=quantile(Sample,0.05)
    if(N.true<upper & N.true>lower)Coverage=1
  }
  mean(Coverage)
}

MSE_N=function(n.bar,N){
  return(sum((N-n.bar)^2))
}

get_avg_abundance_estimator<-function(N){ #N is vector of total abundance by iteration
  Avg.N=optim(par=200000,fn=MSE_N,method="BFGS",N=N)$par
}

ipl=0
for(igen in 1:2){
  for(iest in 1:3){
    #for(itrans in 1:2){
      for(isim in 1:100){
        cat(paste('ipl ',ipl,'\n'))
        fname=paste("d:/ST_out/ST_BOSS_out_gen",igen,"_est_",Est.names[iest],"_sim",isim,".Rdata",sep="")
        fname2=paste("./sim_BOSS_data/simBering_gen",Gen.names[igen],"_sim",isim,sep='')
        if(Est.names[iest]=="AST" | (Est.names[iest]=="STPC" & isim%in%Trans10.sims) |  
             (Est.names[iest]=="CPIF" & isim%in%Trans10.sims)){
                  
          load(fname)
          ipl=ipl+1
          DF[ipl,"Sim.num"]=isim
          DF[ipl,"Gen.mod"]=Gen.names[igen]
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
          
          #now for the year specific sims
          Results[ipl,,1]= (N.est -N.true)^2
          Results[ipl,,2] = (N.est-N.true)/N.true
          Results[ipl,,3]=N.se/N.est
          for(iyr in 1:t.steps){
            Results[ipl,iyr,4]=get_coverage(MCMC$MCMC$N[iyr,],N.true[iyr])          
          }
        }
      }
    #}
  }
}

#limit analysis to CPIF generating model - RS has too many issues with fastly 
#dissipating ice in Bristol Bay
Which.CPIF=which(DF[,"Gen.mod"]=="CPIF")
Results=Results[Which.CPIF,,]
DF=DF[which(DF[,"Gen.mod"]=="CPIF"),]
DF[which(DF[,"Gen.mod"]=="CPIF"),"Gen.mod"]="CPUD"
DF.18=DF
DF.18[,4:7]=Results[,18,]  #day surveys went over southeastern Bering




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

DF=rbind(DF,DF.18)
rownames(DF)=c(1:240)
Limit.to.27Apr=c(rep("Posterior abundance (average)",120),rep("Posterior abundance (18Apr)",120))
DF=cbind(DF,Limit.to.27Apr)

#plot proportion relative bias
library(ggplot2)
bias.plot = ggplot(DF,aes(Est.mod,Bias))+geom_boxplot()+facet_wrap(~Limit.to.27Apr,ncol=1,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")
pdf("bias_Bering.pdf")
bias.plot
dev.off()

save(DF,file="sim_data_Bering_DF.Rdata")
#Now create tables of results
require(dplyr)

poop = group_by(DF,Limit.to.27Apr,Est.mod)
Sim.Bering.table=summarise(poop,Median.bias=median(Bias),Coverage=mean(Coverage),CV=median(CV),MSE=mean(MSE))

require(xtable)
my.table=xtable(Sim.Bering.table)
print(my.table,file="sim_Bering_results_table.txt",include.rownames=FALSE)

