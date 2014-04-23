# create spatial points object to hold knot locations for process convolution models of BOSS data

library(sp)
library(rgeos)
library(Matrix)

#load BOSS grid 
load('AlaskaBeringData2012_17April2014.Rdat') #read in "Data" holding grid info

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

Knot.cell.distances=gDistance(Knots[Which.include,],Data$Grid[[1]],byid=TRUE)
diff.x=(x.max-x.min)/6
diff.y=(y.max-y.min)/6
sigma=(diff.x+diff.y)/2

Knot.Adj=rect_adj(7,7)
Knot.Adj=Knot.Adj[Which.include,Which.include]
Q.knot=-Knot.Adj
diag(Q.knot)=apply(Knot.Adj,2,'sum')
Q.knot=Matrix(Q.knot)

K=dnorm(Knot.cell.distances,0,sigma)
K=K/apply(K,1,'sum')
K.data=list(K=K,Q.knot=Q.knot)
save(K.data,file="Knot_cell_distances.Rdata")

pdf('BOSS_Grid_wKnots.pdf')
plot(Knots[Which.include,],,col='red',pch=20)
plot(Data$Grid[[1]],add=TRUE)
dev.off()