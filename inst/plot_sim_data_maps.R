#plot simulated data maps
library(ggplot2)
library(grid)
load('d:/ST_out/ST_out_gen1_est_STPC_trans_2_sim5.Rdata')

load("./sim_generic_data/simdata_genRS2closed_trans5_sim5")
Data=Sim.data$Data

#construct data frame for ggplot2...  facet by variable type (covariate, abundance, counts, estimates)
#and time step (4 choices)
Cells=c(2,6,14,20)
DF=data.frame(Type=c(rep("Covariate",S*4),rep("Abundance",S*4),rep("Count",S*4),rep("Estimate",S*4)),
              Time=rep(c(rep(Cells[1],S),rep(Cells[2],S),rep(Cells[3],S),rep(Cells[4],S)),4))
DF$Cell=rep(1:400,16)

#fill in covariate values
DF$Response=rep(NA,S*16)
DF$Response[1:S]=Data$Grid[[Cells[1]]]@data[,1]
DF$Response[(S+1):(2*S)]=Data$Grid[[Cells[2]]]@data[,1]
DF$Response[(2*S+1):(3*S)]=Data$Grid[[Cells[3]]]@data[,1]
DF$Response[(3*S+1):(4*S)]=Data$Grid[[Cells[4]]]@data[,1]
# true abundance
DF$Response[4*S+1:(4*S)]=as.vector(Sim.data$N[,c(Cells[1],Cells[2],Cells[3],Cells[4])])
# Count - unobserved cells remain NAs
Cur.count.data=Sim.data$Count.data[which(Sim.data$Count.data[,"Time"]%in%Cells),]
for(irow in 1:nrow(Cur.count.data)){
  DF$Response[which(DF$Time==Cur.count.data[irow,"Time"] & DF$Type=="Count" & DF$Cell==Cur.count.data[irow,"Cell"])]=Cur.count.data[irow,"Count"]
}
# Abundance
N.est=apply(MCMC$MCMC$Pred,c(1,2),'mean')
DF$Response[(12*S+1):(16*S)]=as.vector(N.est[,c(Cells[1],Cells[2],Cells[3],Cells[4])])

require(rgeos)
require(ggplot2)
Centroids=data.frame(gCentroid(Data$Grid[[1]],byid=TRUE))
DF$Easting=rep(Centroids[,1],16)
DF$Northing=rep(Centroids[,2],16)

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
p1 <- ggplot(subset(DF,Type=="Covariate"))+aes(Easting,Northing,fill=Response)+geom_raster()+facet_grid(Type~Time,scales='free')+tmp.theme+theme(axis.title.x=element_blank(),strip.background=element_rect(fill="white"),plot.margin=unit(c(1,0.75,-1.25,1),"line"))+scale_fill_gradientn(colours=myPalette(100))
p2 <- ggplot(subset(DF,Type=="Abundance"))+aes(Easting,Northing,fill=Response)+geom_raster()+facet_grid(Type~Time,scales='free')+tmp.theme+theme(axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),plot.margin=unit(c(.25,0.45,-.25,1),"line"))+scale_fill_gradientn(limits=c(0,440),colours=myPalette(100))
p3 <- ggplot(subset(DF,Type=="Count"))+aes(Easting,Northing,fill=Response)+geom_raster()+facet_grid(Type~Time,scales='free')+tmp.theme+theme(axis.title.x=element_blank(),strip.text.x=element_blank(),strip.background=element_rect(fill="white"),plot.margin=unit(c(-.5,1,.5,1),"line"))+scale_fill_gradientn(colours=myPalette(100))
p4 <- ggplot(subset(DF,Type=="Estimate"))+aes(Easting,Northing,fill=Response)+geom_raster()+facet_grid(Type~Time,scales='free')+tmp.theme+theme(strip.text.x=element_blank(),strip.background=element_rect(fill="white"),plot.margin=unit(c(-1.5,0.5,1,1),"line"))+scale_fill_gradientn(limits=c(0,440),colours=myPalette(100))

pdf(file="sim_generic_maps.pdf")
grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=4))
dev.off()

