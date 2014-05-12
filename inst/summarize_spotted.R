### analyze and visualize 'generic' simulation data
library(grid)
library(ggplot2)
library(animation)
library(gridExtra)
t.steps=29

#input effort and grid info for plotting
#load('c:/users/paul.conn/git/STabundance/BOSS_2012Effort_22Apr14.Rdata')  
load("Effort2012_BOSSst_22Apr2014.Rdata")
load("AlaskaBeringData2012_17April2014.Rdat")  #boss grid, ice data

#load MCMC results from chain 1
load('d:/ST_out/spotted_Bering_CPIF3.Rdata')
MCMC.CPIF1=MCMC
load('d:/ST_out/spotted_Bering_CPIF4.Rdata')
MCMC.CPIF2=MCMC

load('d:/ST_out/spotted_Bering_STPC1.Rdata')
MCMC.STPC1=MCMC
load('d:/ST_out/spotted_Bering_STPC2.Rdata')
MCMC.STPC2=MCMC
MCMC.STPC1$MCMC$N=apply(MCMC.STPC1$MCMC$Pred,c(2,3),'sum')
MCMC.STPC2$MCMC$N=apply(MCMC.STPC2$MCMC$Pred,c(2,3),'sum')
MCMC.STPC.N=cbind(MCMC.STPC1$MCMC$N,MCMC.STPC2$MCMC$N)

n.linex.iter=100 #how many MCMC samples to calculate a distribution of the linex 'a' parameter
n.mcmc.iter=dim(MCMC.STPC1$MCMC$Pred)[3]
n.obs=nrow(Effort$Count.data)
Which.iter1=sample(c(1:(n.linex.iter/2)),n.linex.iter/2)
Which.iter2=sample(c(1:(n.linex.iter/2)),n.linex.iter/2)
Cur.pred=matrix(0,n.obs,n.linex.iter)
set.seed(12345)
for(iiter in 1:(n.linex.iter/2)){  #half from each chain
  for(iobs in 1:n.obs){
    Cur.pred[iobs,iiter]=MCMC.STPC1$MCMC$Pred[Effort$Count.data[iobs,"Cell"],Effort$Count.data[iobs,"Time"],Which.iter1[iiter]]*Effort$Count.data[iobs,"AreaSurveyed"]*Effort$Area.hab[Effort$Count.data[iobs,"Cell"]]    
    Cur.pred[iobs,n.linex.iter/2+iiter]=MCMC.STPC2$MCMC$Pred[Effort$Count.data[iobs,"Cell"],Effort$Count.data[iobs,"Time"],Which.iter2[iiter]] *  Effort$Count.data[iobs,"AreaSurveyed"]*Effort$Area.hab[Effort$Count.data[iobs,"Cell"]]         
  }
}
calc_linex_a(Pred.G=Cur.pred,Obs.G=Effort$Count.data[,"Count"],min.a=0.00001,max.a=1.0)$minimum
#plot_N_map(17,apply(MCMC.CPIF1$MCMC$Pred,c(1,2),'median'),Grid=Data$Grid,leg.title="Abundance")

#plot observed, predicted
pdf('Spotted.STPC.plots.pdf')
par(mfrow=c(2,2))
plot(c(MCMC.STPC1$MCMC$tau.epsilon,MCMC.STPC1$MCMC$tau.epsilon),xlab="MCMC iteration",ylab="tau_epsilon",main='')
plot(Effort$Count.data[,"Count"],apply(Cur.pred,1,'mean'),xlab="Observed count",ylab="Predicted count")
crap=c(0:14)
lines(crap,crap,col='red')
plot(MCMC.STPC.N[18,],xlab="MCMC iteration",ylab="Abundance, 27Apr")
N.18=MCMC.STPC.N[18,]
hist(N.18[which(N.18<1000000)],xlab="Abundance",ylab="Posterior density",freq=FALSE,main='')
dev.off()

#plot trend in posterior median
pdf('Spotted_STPC_trend.pdf')
plot(apply(MCMC.STPC.N,1,'median',na.rm=TRUE),type="l",xlab="Day",ylab="Spotted seal abundance (STPC)",lwd=2,cex.lab=1.2,cex.axis=1.2)
dev.off()

#define iterations to use and combine data from both chains
Iter=c(500:1999)
dim.Pred=dim(MCMC.CPIF1$MCMC$Pred)
N.CPIF=c(MCMC.CPIF1$MCMC$N[Iter],MCMC.CPIF2$MCMC$N[Iter])
Pred.CPIF=array(0,dim=c(dim.Pred[1],dim.Pred[2],length(Iter)*2))
Pred.CPIF[,,1:length(Iter)]=MCMC.CPIF1$MCMC$Pred[,,Iter]
Pred.CPIF[,,(length(Iter)+1):(2*length(Iter))]=MCMC.CPIF2$MCMC$Pred[,,Iter]
Pred.mean=apply(MCMC$MCMC$Pred,c(1,2),'mean')

pdf('N.chains.spotted.pdf')
plot(MCMC.CPIF1$MCMC$N[Iter],xlab="MCMC iteration (thinned)",ylab="Spotted seal abundance")
points(MCMC.CPIF2$MCMC$N[Iter],col='red')
dev.off()

###########
# plot ice covariate, surveys, estimates as a function of 4 different time steps
#########
S=nrow(Data$Grid[[1]]@data)
Cells=c(2,11,18,29)
DF=data.frame(Type=c(rep("Covariate",S*4),rep("Count",S*4),rep("Estimate",S*4)),Time=rep(c(rep(Cells[1],S),rep(Cells[2],S),rep(Cells[3],S),rep(Cells[4],S)),3))
DF$Cell=rep(1:S,12)

#fill in covariate values
DF$Response=rep(NA,S*12)
DF$Response[1:S]=Data$Grid[[Cells[1]]]@data[,"ice_conc"]
DF$Response[(S+1):(2*S)]=Data$Grid[[Cells[2]]]@data[,"ice_conc"]
DF$Response[(2*S+1):(3*S)]=Data$Grid[[Cells[3]]]@data[,"ice_conc"]
DF$Response[(3*S+1):(4*S)]=Data$Grid[[Cells[4]]]@data[,"ice_conc"]
# Count - unobserved cells remain NAs
Cur.count.data=Effort$Count.data[which(Effort$Count.data[,"Time"]%in%Cells),]
for(irow in 1:nrow(Cur.count.data)){
  DF$Response[which(DF$Time==Cur.count.data[irow,"Time"] & DF$Type=="Count" & DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Count"]
}
# Abundance
#N.est=apply(MCMC$MCMC$Pred,c(1,2),'mean')
N.est=N.CPIF
DF$Response[(8*S+1):(12*S)]=as.vector(Pred.mean[,c(Cells[1],Cells[2],Cells[3],Cells[4])])

require(rgeos)
require(ggplot2)
Centroids=data.frame(gCentroid(Data$Grid[[1]],byid=TRUE))
DF$Easting=rep(Centroids[,1],12)
DF$Northing=rep(Centroids[,2],12)

tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),legend.title=element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))
#p1=ggplot(DF)+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme #+scale_fill_continuous(name=leg.title)
#p1=p1+facet_grid(Type~Time,scales='free')
#if(is.null(highlight)==FALSE){
#p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
#  p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-1/2,xmax=Easting+1/2,ymin=Northing-1/2,ymax=Northing+1/2))
#}
#p1
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

DF[which(DF[,"Time"]==Cells[1]),"Time"]=paste("t = ",Cells[1],sep='')
DF[which(DF[,"Time"]==Cells[2]),"Time"]=paste("t = ",Cells[2],sep='')
DF[which(DF[,"Time"]==Cells[3]),"Time"]=paste("t = ",Cells[3],sep='')
DF[which(DF[,"Time"]==Cells[4]),"Time"]=paste("t = ",Cells[4],sep='')
DF[,"Time"]=factor(DF[,"Time"],levels=c(paste("t = ",Cells[1],sep=''),paste("t = ",Cells[2],sep=''),paste("t = ",Cells[3],sep=''),paste("t = ",Cells[4],sep='')))
library(gridExtra)
p1 <- ggplot(subset(DF,Type=="Covariate"))+aes(Easting,Northing,fill=Response)+geom_raster()+facet_grid(Type~Time,scales='free')+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),plot.margin=unit(c(1,0.75,-1.25,1),"line"))+scale_fill_gradientn(colours=myPalette(100))
p2 <- ggplot(subset(DF,Type=="Count"))+aes(Easting,Northing,fill=Response)+geom_raster()+facet_grid(Type~Time,scales='free')+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),plot.margin=unit(c(-.5,1,.5,1),"line"))+scale_fill_gradientn(colours=myPalette(100))
p3 <- ggplot(subset(DF,Type=="Estimate"))+aes(Easting,Northing,fill=Response)+geom_raster()+facet_grid(Type~Time,scales='free')+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),plot.margin=unit(c(-1.5,0.5,1,1),"line"))+scale_fill_gradientn(limits=c(0,1250),colours=myPalette(100))

pdf(file="spotted_maps.pdf")
grid.arrange(arrangeGrob(p1,p2,p3,nrow=3))
dev.off()


#create time lapsed video of covariate, transect counts, abundance
Dates=data.frame(matrix(0,29,3))
colnames(Dates)=c("x","y","date")
Dates[,"date"]=c(paste(c(10:30),"April 2012"),paste(c(1:8),"May 2012"))

saveHTML(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid[[it]]@data[,"ice_conc"]
    # Count - unobserved cells remain NAs
    Cur.count.data=Effort$Count.data[which(Effort$Count.data[,"Time"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Time"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Count"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=as.vector(Pred.mean[,it])
    
    Centroids=data.frame(gCentroid(Data$Grid[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=myPalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=myPalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1250),colours=myPalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,img.name="Spotted_Bering_ST_html",outdir=getwd()
)

saveLatex(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid[[it]]@data[,"ice_conc"]
    # Count - unobserved cells remain NAs
    Cur.count.data=Effort$Count.data[which(Effort$Count.data[,"Time"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Time"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Count"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=as.vector(Pred.mean[,it])
    
    Centroids=data.frame(gCentroid(Data$Grid[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=myPalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=myPalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1250),colours=myPalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,img.name="Spotted_Bering_ST_pdf",latex.filename="Spotted_Bering_ST_ani.tex"  
  
)

saveGIF(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid[[it]]@data[,"ice_conc"]
    # Count - unobserved cells remain NAs
    Cur.count.data=Effort$Count.data[which(Effort$Count.data[,"Time"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Time"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Count"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=as.vector(Pred.mean[,it])
    
    Centroids=data.frame(gCentroid(Data$Grid[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=myPalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=myPalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1250),colours=myPalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,movie.name="Spotted_Bering_ST.gif",img.name="Spotted_Bering_ST",clean=TRUE 
  
)

saveVideo(
  for(it in 1:t.steps){
    Tmp.DF=data.frame(Type=c(rep("Covariate",S),rep("Count",S),rep("Estimate",S)),Time=rep(rep(it,S),3))
    Tmp.DF$Cell=rep(1:S,3)
    
    #fill in covariate values
    Tmp.DF$Response=rep(NA,S*3)
    Tmp.DF$Response[1:S]=Data$Grid[[it]]@data[,"ice_conc"]
    # Count - unobserved cells remain NAs
    Cur.count.data=Effort$Count.data[which(Effort$Count.data[,"Time"]==it),]
    for(irow in 1:nrow(Cur.count.data)){
      Tmp.DF$Response[which(Tmp.DF$Time==Cur.count.data[irow,"Time"] & Tmp.DF$Type=="Count" & Tmp.DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Count"]
    }
    # Abundance
    Tmp.DF$Response[(2*S+1):(3*S)]=as.vector(Pred.mean[,it])
    
    Centroids=data.frame(gCentroid(Data$Grid[[1]],byid=TRUE))
    Tmp.DF$Easting=rep(Centroids[,1],3)
    Tmp.DF$Northing=rep(Centroids[,2],3)
    Cur.date=Dates[it,]
    
    p1 <- ggplot(subset(Tmp.DF,Type=="Covariate"))+labs(title="Sea ice concentration")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(2.5,"mm"),axis.title.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1),colours=myPalette(100))
    p2 <- ggplot(subset(Tmp.DF,Type=="Count"))+labs(title="Spotted seal count")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(legend.margin=unit(5.5,"mm"),axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,14),colours=myPalette(100))
    p3 <- ggplot(subset(Tmp.DF,Type=="Estimate"))+labs(title="Apparent abundance estimate")+aes(Easting,Northing,fill=Response)+geom_raster()+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),title=element_text(size=10))+scale_fill_gradientn(limits=c(0,1250),colours=myPalette(100))
    p4 <- ggplot(Cur.date)+geom_text(aes(x=x,y=y,label=date),size=8)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),plot.margin=unit(c(0,3,0,0),"cm"))
    grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
  }
  ,video.name="Spotted_Bering_ST.avi",img.name="Spotted_Bering_ST",ffmpeg="C:/ffmpeg/bin/ffmpeg.exe",clean=TRUE 
  
)

#plot histogram and time series of abundance estimators
#N.18.trunc=N.18[which(N.18<1000000)]
#N=c(N.18.trunc,N.CPIF)
#DM.N=data.frame(matrix(0,length(N),2))
#colnames(DM.N)=c("Estimator","Abundance")
#DM.N[,"Estimator"]=c(rep("STPC-27Apr",length(N.18.trunc)),rep("CPIF",length(N.CPIF)))
#DM.N[,2]=N
#h1 <- ggplot(DM.N)+geom_density(aes(x=Abundance,linetype=factor(Estimator)))  #posteriors for avg N, N_18
#h2 <- ggplot2()  #time series
DM.N=data.frame(matrix(0,length(N.CPIF),1))
colnames(DM.N)="Abundance"
DM.N[,1]=N.CPIF
h1 <- ggplot(DM.N)+geom_density(aes(x=Abundance),size=1)+labs(x="Spotted seal apparent abundance",y="Posterior density")+theme(text=element_text(size=16))
pdf("Spotted_N_hist.pdf")
h1
dev.off()

#plot covariate effects
DF.ice=data.frame(matrix(0,length(Iter)*4,2))
colnames(DF.ice)=c("Type","Value")
DF.ice[,1]=c(rep("Linear",length(Iter)*2),rep("Quadratic",length(Iter)*2))
DF.ice[,2]=c(MCMC.CPIF1$MCMC$Beta[1,Iter],MCMC.CPIF2$MCMC$Beta[1,Iter],MCMC.CPIF1$MCMC$Beta[2,Iter],MCMC.CPIF2$MCMC$Beta[2,Iter])
C1=ggplot(DF.ice)+geom_density(aes(x=Value,linetype=Type),size=1)+labs(x="Sea ice effect",y="Posterior density")+theme(text=element_text(size=16))

DF.cov=data.frame(matrix(0,length(Iter)*6,2))
colnames(DF.cov)=c("Type","Value")
DF.cov[,1]=c(rep("Ice edge",length(Iter)*2),rep("Sqrt(Ice edge)",length(Iter)*2),rep("Shelf break",length(Iter)*2))
DF.cov[,2]=c(MCMC.CPIF1$MCMC$Beta[3,Iter],MCMC.CPIF2$MCMC$Beta[3,Iter],MCMC.CPIF1$MCMC$Beta[4,Iter],MCMC.CPIF2$MCMC$Beta[4,Iter],MCMC.CPIF1$MCMC$Beta[5,Iter],MCMC.CPIF2$MCMC$Beta[5,Iter])
C2=ggplot(DF.cov)+geom_density(aes(x=Value,linetype=Type),size=1)+labs(x="Covariate effect",y="Posterior density")+theme(text=element_text(size=16))

DF.cov=data.frame(matrix(0,length(Iter)*10,2))
colnames(DF.cov)=c("Covariate","Value")
DF.cov[,1]=c(rep("Ice conc",length(Iter)*2),rep("Ice conc^2",length(Iter)*2),rep("Edge",length(Iter)*2),rep("Sqrt(Edge)",length(Iter)*2),rep("Shelf",length(Iter)*2))
DF.cov[,2]=c(MCMC.CPIF1$MCMC$Beta[1,Iter],MCMC.CPIF2$MCMC$Beta[1,Iter],MCMC.CPIF1$MCMC$Beta[2,Iter],MCMC.CPIF2$MCMC$Beta[2,Iter],MCMC.CPIF1$MCMC$Beta[3,Iter],MCMC.CPIF2$MCMC$Beta[3,Iter],MCMC.CPIF1$MCMC$Beta[4,Iter],MCMC.CPIF2$MCMC$Beta[4,Iter],MCMC.CPIF1$MCMC$Beta[5,Iter],MCMC.CPIF2$MCMC$Beta[5,Iter])
C3=ggplot(DF.cov)+geom_density(aes(x=Value,linetype=Covariate),size=1)+labs(x="Covariate effect",y="Posterior density")+theme(text=element_text(size=16))


###now, overall effect of covariates
Mean.covs=c(0,0,0,0,0)
Max.covs=apply(Data$Grid[[1]]@data[,c("ice_conc","dist_edge","dist_shelf")],2,'max')
Min.covs=apply(Data$Grid[[1]]@data[,c("ice_conc","dist_edge","dist_shelf")],2,'min')
for(it in 1:t.steps){
  Tmp.cov=Data$Grid[[it]]@data[,c("ice_conc","ice_conc","dist_edge","dist_edge","dist_shelf")]
  #Tmp.cov[,2]=Tmp.cov[,2]^2
  #Tmp.cov[,4]=Tmp.cov[,4]^0.5
  Mean.covs=Mean.covs+apply(Tmp.cov,2,'mean')
  Tmp.min=apply(Data$Grid[[it]]@data[,c("ice_conc","dist_edge","dist_shelf")],2,'min')
  Tmp.max=apply(Data$Grid[[it]]@data[,c("ice_conc","dist_edge","dist_shelf")],2,'max')
  for(icov in 1:3){
    if(Tmp.min[icov]<Min.covs[icov])Min.covs[icov]=Tmp.min[icov]
    if(Tmp.max[icov]>Max.covs[icov])Max.covs[icov]=Tmp.max[icov]
    
  }
}
Mean.covs=Mean.covs/t.steps
Mean.covs[2]=Mean.covs[1]^2
Mean.covs[4]=Mean.covs[3]^0.5

Ice=c(0:100)/100
Edge=c(0:325)/100
Shelf=c(0:240)/100

Beta.post=cbind(MCMC.CPIF1$MCMC$Beta,MCMC.CPIF2$MCMC$Beta)
Post.means=apply(Beta.post,1,'mean')
DF.eff=data.frame(matrix(0,length(Ice)+length(Edge)+length(Shelf),3))
colnames(DF.eff)=c("Covariate","Value","Effect")
DF.eff[,"Covariate"]=c(rep("Ice conc.",length(Ice)),rep("Edge",length(Edge)),rep("Shelf",length(Shelf)))
DF.eff[,"Value"]=c(Ice,Edge,Shelf)
Ice.pred=exp(Ice*Post.means[1]+Ice^2*Post.means[2]+Mean.covs[3:5]%*%Post.means[3:5])
Edge.pred=exp(Edge*Post.means[3]+sqrt(Edge)*Post.means[4]+Mean.covs[c(1,2,5)]%*%Post.means[c(1,2,5)])
Shelf.pred=exp(Shelf*Post.means[5]+Mean.covs[1:4]%*%Post.means[1:4])
DF.eff[,"Effect"]=c(Ice.pred,Edge.pred,Shelf.pred)

E1=ggplot(DF.eff)+geom_line(aes(x=Value,y=Effect,linetype=Covariate),size=1)+theme(text=element_text(size=16))

C1=C1+ggtitle("A.")+theme(plot.title=element_text(hjust=0))
C2=C2+ggtitle("B.")+theme(plot.title=element_text(hjust=0))
h1=h1+ggtitle("D.")+theme(plot.title=element_text(hjust=0))
E1=E1+ggtitle("C.")+theme(plot.title=element_text(hjust=0))

pdf(file="spotted_estimates.pdf")
grid.arrange(arrangeGrob(C1,C2,E1,h1,nrow=2))
dev.off()

C3=C3+ggtitle("A.")+theme(plot.title=element_text(hjust=0))
E1=E1+ggtitle("B.")+theme(plot.title=element_text(hjust=0))    
pdf(file="spotted_estimates.pdf")
grid.arrange(arrangeGrob(C3,E1,nrow=2))
dev.off()
                          

