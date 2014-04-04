### analyze and visualize 'generic' simulation data



Trans10.sims=c(1:10)*10-5

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

for(igen in 1:3){
  for(iest in 1:4){
    for(itrans in 1:2){
      for(isim in 1:100){
        fname=paste("d:/ST_out/ST_out_gen",igen,"_est_",Est.names[iest],"_trans_",itrans,"_sim",isim,".Rdata",sep="")
        fname2=paste("./sim_generic_data/simdata_gen",Gen.names[igen],"_trans",Trans.names[itrans],"_sim",isim,sep='')
        if((Est.names[iest]=="OPRS" & itrans==2 & isim%in%Trans10.sims) |
             Est.names[iest]=="AST" | (Est.names[iest]=="STPC" & itrans==2) 
             | (Est.names[iest]=="STPC" & isim%in%Trans10.sims) | (Est.names[iest]=="CPIF" & itrans==2)
                | (Est.names[iest]=="CPIF" & isim%in%Trans10.sims)){
                  
          load(fname)
          ipl=ipl+1
          DF[ipl,"Sim.num"]=isim
          DF[ipl,"Gen.mod"]=Gen.names[igen]
          DF[ipl,"Transects"]=Trans.names[itrans]
          DF[ipl,"Est.mod"]=Est.names[iest]
          
          MCMC$MCMC$N=apply(MCMC$MCMC$Pred,c(2,3),'sum')
          N.est=apply(MCMC$MCMC$N,1,'mean')
          N.se=apply(MCMC$MCMC$N,1,'st_dev')
          load(fname2)
          N.true=apply(Sim.data$N,2,'sum')
          t.steps=length(N.true)
          trend.est=lm(N.est~c(1:t.steps))$coef
          trend.est=1+trend.est[2]/trend.est[1]*t.steps  #fitted lambda for whole study
          trend.true=lm(N.true~c(1:length(N.true)))$coef
          trend.true=1+trend.true[2]/trend.true[1]*t.steps
          DF[ipl,"Bias"]=mean((N.est-N.true)/N.true)
          DF[ipl,"MSE"]=mean((N.est-N.true)^2)
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


DF.trunc.1 = DF[which(DF[,"Bias"]>2.0),]
DF.trunc.5 = DF[which(DF[,"Bias"]>0.25 & DF[,"Transects"]=="5 Transects"),]

DF2=DF[-which(DF[,"Bias"]>2.0),]
DF2=DF2[-which(DF2[,"Bias"]>0.25 & DF2[,"Transects"]=="5 Transects"),]

DF.trunc.1[,"Bias"]=2.0
DF.trunc.5[,"Bias"]=0.25
#delete any values from primary 

#remove extreme outlier
#outlier=which(DF[,"Bias"]>10)
#DF2=DF[-outlier,]

#plot proportion relative bias
bias.plot = ggplot(DF2,aes(Est.mod,Bias))+geom_boxplot()+facet_grid(Transects~Gen.mod,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
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

