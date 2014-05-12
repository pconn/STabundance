#plot simulated data maps
#library(ggplot2)
#library(grid)
library(sp)
library(rgeos)

load('Data_for_ST_plot.Rdata') #includes grid, other spatial objects associated w landscape

#reconstitute knots
Coords=coordinates(Data$Grid[[1]])
x.min=min(Coords[,1])-100000
x.max=max(Coords[,1])+100000
y.min=min(Coords[,2])-100000
y.max=max(Coords[,2])+100000

X=x.min+(x.max-x.min)/6*c(0:6)
Y=y.min+(y.max-y.min)/6*c(6:0)
XY=expand.grid(x=X,y=Y)

Knots=SpatialPoints(coords=XY,proj4string=CRS(proj4string(Data$Grid[[1]])))

Distances=gDistance(Knots,Data$Grid[[1]],byid=TRUE)
Distances=apply(Distances,2,'min')
my.buffer=150000
Which.include=which(Distances<my.buffer)
Knots=Knots[Which.include,]

#load on effort tracks
load('2012_Tracks_for_ST_analysis.Rdata')

pdf(file="BOSS_survey_map.pdf")
plot(Knots,col='red',pch=20)
plot(Data$Grid[[1]],add=TRUE)
plot(gBoundary(Shelf_break[1,]),col='brown',add=TRUE,lwd=2)
plot(gBoundary(EEZ_Alaska),col='orange',add=TRUE,lwd=2)
plot(Tracks,add=TRUE,col='blue')
plot(alaska_dcw,add=TRUE,lwd=2)
plot(alaska_dcw,add=TRUE,col='gray')
plot(russia_dcw,add=TRUE,lwd=2)
plot(russia_dcw,add=TRUE,col='gray')
plot(Knots,add=TRUE,col='red',pch=20)
dev.off()

