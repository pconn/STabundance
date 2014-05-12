### format_effort.R
### function to format effort and hotspot/count data for spatio-temporal abundance analysis


require(sp)
require(rgeos)
require(rgdal)
require(Matrix)
require(lubridate)
source('c:/users/paul.conn/git/BOSSst/BOSSst/R/util_funcs.R')

date.start=as.Date("2012-04-10")
date.end=as.Date("2012-05-08")
t.steps=as.numeric(date.end-date.start)

#load spatio-temporal covariate & grid data
#read in Sp lines DF for transects ("on_effort_tracks.sldf") and Sp points Df for FMC points ("fmclogs.spdf")
load('c:/users/paul.conn/git/STabundance/BOSS_2012Effort_22Apr14.Rdata')  
load("AlaskaBeringData2012_17April2014.Rdat")  #boss grid, ice data
#rename IDs for grid cells to be 1:1299 instead of old IDs associated with 'full' Bering grid
# (needed for some of the intersection stuff below)
S=nrow(Data$Grid[[1]])
for(it in 1:t.steps){
  for(ipoly in 1:S){
    Data$Grid[[it]]@polygons[[ipoly]]@ID=as.character(ipoly)
  }
}

laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
#reproject onto BOSS Grid
Tracks=spTransform(on_effort_tracks.sldf, CRS(laea_180_proj))
Points=spTransform(fmclogs.spdf, CRS(laea_180_proj))
Hotspots=spTransform(hotspots.spdf,CRS(laea_180_proj))
Hotspots=Hotspots[which(Hotspots@data[,"species"]=="bd"),]  #limit to bearded seals

Tracks$Day=rep(0,nrow(Tracks))
Flt.ids=unique(Tracks$flightid)
#This is really just to help visualize data
Flt.table=data.frame(id=Flt.ids,min.time=rep(0,length(Flt.ids)))
for(irow in 1:nrow(Flt.table))Flt.table[irow,"min.time"]=min(Points[which(Points$flightid==Flt.ids[irow]),]$dt_utc)
plot(Data$Grid[[1]])
#plot(Tracks[which(Tracks$flightid=="12_OtterFl01"),],add=TRUE,col="blue")

#limit to 2012
myfun<-function(x)strsplit(x,'_')[[1]][1]
Year=sapply(as.character(Flt.table[,"id"]),myfun)
Yr12.ind=which(Year=="12")
Flights=Flt.ids[Yr12.ind]

Year=sapply(as.character(Hotspots@data[,"flightid"]),myfun)
Yr12.ind=which(Year=="12")
Hotspots=Hotspots[Yr12.ind,]

# limit to 4/10-5/8
Flights=Flights[-c(30,20,28,29,32,34,36)]
Tracks=Tracks[which(Tracks$flightid %in% Flights),]
Points=Points[which(Points$flightid %in% Flights),]
Hotspots=Hotspots[which(Hotspots@data$flightid %in% Flights),]
#decrease to 1 record every 4 seconds for computation of average swath diameter, time of survey
I.point=rep(0,nrow(Points))
#I.point[1:(ceiling(length(I.point)/10))*10-9]=1
I.point[1:(ceiling(length(I.point)/4))*4-3]=1
Points=Points[which(I.point==1),]
Points=Points[which(Points$effort=="On"),]

plot(Tracks,add=TRUE,col="blue")

save(Tracks,file="2012_Tracks_for_ST_analysis.Rdata")



#intersect on effort spatial points with grid
Flt.anal=Flights
Study.area=gUnionCascaded(Data$Grid[[1]])
int=gIntersects(Points,Study.area,byid=TRUE)
Points=Points[apply(int,2,'sum')==1,]
int=gIntersects(Points,Data$Grid[[1]],byid=TRUE)
which.cell=function(x)which(x==1)
Cell.id=unlist(apply(int,2,which.cell))
Date=Points$dt_utc

Day=as.numeric(as.Date(Date,tz="PST8PDT")-date.start)+1 #make date.start = day 1

#now get hour in UTC
Date=format(Date,tz="UTC",usetz=TRUE)
Hour=round(hour(Date)+minute(Date)/60) #round to nearest hour

#similar intersection, date, time calculations with hotspots
int=gIntersects(Hotspots,Study.area,byid=TRUE) #remove hotspots off of grid
Hotspots=Hotspots[-which(int==FALSE),]
#Hotspots=Hotspots[apply(int,2,'sum')==1,]
int=gIntersects(Hotspots,Data$Grid[[1]],byid=TRUE)
which.cell=function(x)which(x==1)
Cell.id.hs=unlist(apply(int,2,which.cell))
Date.hs=Hotspots$dt_utc


Day.hs=as.numeric(as.Date(Date.hs,tz="PST8PDT")-date.start)+1 #make date.start = day 1
pdf("Bearded_seal_count.pdf")
hist(Day.hs,breaks=(c(0:29)+0.5),main='',xlab="Day",ylab="Bearded seal count")
dev.off()

#now get hour in UTC
Date.hs=format(Date.hs,tz="UTC",usetz=TRUE)
Hour.hs=round(hour(Date.hs)+minute(Date.hs)/60) #round to nearest hour


#calculate swath widths for each on effort point (center, right, left depending on camera)
myfun<-function(x)strsplit(x,'_')[[1]][2]
Temp=sapply(as.character(Points$flightid),myfun)
myfun<-function(x)strsplit(x,'F')[[1]][1]
Airplane=sapply(Temp,myfun)
I_Aero=(Airplane=="Aero")


fov_horz<-24.9*pi/180
offset_aero<-12.5*pi/180
ft_to_meters=0.3048
diameter=rep(0,length(I_Aero))
diameter[I_Aero==TRUE]=ft_to_meters*Points[["gpsalt"]][I_Aero==TRUE]*(tan(offset_aero+fov_horz/2)-tan(offset_aero-fov_horz/2))
diameter[I_Aero==FALSE]=ft_to_meters*Points[["gpsalt"]][I_Aero==FALSE]*2*tan(fov_horz/2)
Points$diameter=diameter

#now compute area for each unique combination of grid cell and area surveyed


#fix some issues with OtterFl06 having effort off grid
#which.rows=which(flight_segs@data[,"flightid"]=="OtterFl06")
#plot(flight_segs[which.rows[c(1:2,3,4:7)],])
#flight_segs=flight_segs[-which.rows[3],]

#overlay flight tracks over grid cells, calculating cumulative distance flown in each
int<-gIntersects(Tracks,Data$Grid[[1]],byid=TRUE)
vec <- vector(mode="list", length=nrow(Tracks))
for (i in seq(along=vec)) vec[[i]] <- try(gIntersection(Tracks[i,],Data$Grid[[1]][int[,i],], byid=TRUE))
out <- do.call("rbind", vec)
rn <- row.names(out)
nrn <- do.call("rbind", strsplit(rn, " "))
Length.df <- data.frame(Fl=nrn[,1], poly=as.numeric(as.character(nrn[,2])), len=gLength(out,byid=TRUE))
Length.df$Day=rep(NA,nrow(Length.df))
Length.df$Hour=Length.df$Day
Length.df$diameter=Length.df$Day
#calculate area surveyed for each of these flight_segment * grid cell combos
Row.index=rep(0,nrow(Length.df))
Seg_ID = rownames(Tracks@data)
for(i in 1:length(Row.index))Row.index[i]=which(Seg_ID==as.character(Length.df[i,"Fl"]))
#get mean date & time, altitude for each flight/grid cell combination
for(i in 1:nrow(Length.df)){
  Which.points=which(Points@data[,"flightid"]==Tracks@data[Row.index[i],"flightid"] & Cell.id==Length.df[i,"poly"])
  if(length(Which.points)>0){
    #attach day and hour
    Length.df[i,"Day"]=mean(Day[Which.points])
    Alpha=2*pi*(Hour[Which.points]/24)
    x.bar=mean(sin(Alpha))
    y.bar=mean(cos(Alpha))
    r=sqrt(x.bar^2+y.bar^2)
    if(x.bar>0)alpha.bar=acos(y.bar/r)
    else alpha.bar=2*pi-acos(y.bar/r)
    Length.df[i,"Hour"]=round(alpha.bar*0.5/pi*24)
    #attach mean diameter
    Length.df[i,"diameter"]=mean(Points@data[Which.points,"diameter"])
  }
}
Which.missing=which(is.na(Length.df[,"Day"]))  #grid/time combos missing hour, day, altitude
if(length(Which.missing>0)){
  for(i in 1:length(Which.missing)){ #for these few missing values, simply input mean values from flight
    Which.points=which(Points@data[,"flightid"]==Tracks@data[Row.index[i],"flightid"])
    Length.df[Which.missing[i],"Day"]=mean(Day[Which.points])
    Alpha=2*pi*(Hour[Which.points]/24)
    x.bar=mean(sin(Alpha))
    y.bar=mean(cos(Alpha))
    r=sqrt(x.bar^2+y.bar^2)
    if(x.bar>0)alpha.bar=acos(y.bar/r)
    else alpha.bar=2*pi-acos(y.bar/r)
    Length.df[Which.missing[i],"Hour"]=round(alpha.bar*0.5/pi*24)
    #attach mean diameter
    Length.df[Which.missing[i],"diameter"]=mean(Points@data[Which.points,"diameter"])
  }
}



#total area surveyed in each flight * grid cell combo
Area.hab=1-Data$Grid[[1]]@data[,"land_cover"]
Length.df$area=Length.df[,"len"]*0.001*(Length.df$diameter*.001)/(625*Area.hab[Length.df[,"poly"]])  #proportional area surveyed for each cell

#Combine area from cells surveyed by more than one flight segment in same day
# (this shouldn't happen with above code - the following includes some artifacts
# from BOSS analysis with constant abundance over time)
Mapping=matrix(0,nrow(Length.df),3)  #includes cell, day, hour
Area.trans=rep(0,nrow(Length.df))
icounter=0
icounter2=0
I.processed=rep(0,nrow(Length.df))
while(sum(I.processed)<nrow(Length.df)){
  icounter=icounter+1
  if(I.processed[icounter]==0){
    icounter2=icounter2+1
    cur.poly=Length.df[icounter,"poly"]
    cur.day=Length.df[icounter,"Day"]
    cur.area=Length.df[icounter,"area"]
    cur.hour=Length.df[icounter,"Hour"]
    Area.trans[icounter2]=cur.area
    Mapping[icounter2,]=c(cur.poly,cur.day,cur.hour)
    Which.eq=which(Length.df[,"poly"]==cur.poly & Length.df[,"Day"]==cur.day)
    if(length(Which.eq)>1){
      Area.trans[icounter2]=Area.trans[icounter2]+sum(Length.df[Which.eq[2:length(Which.eq)],"area"])
      I.processed[Which.eq[2:length(Which.eq)]]=1
    }
    I.processed[icounter]=1
  }
}
Area.trans=Area.trans[1:icounter2]
Mapping=Mapping[1:icounter2,]

Count.data=data.frame(matrix(0,nrow(Mapping),4))
colnames(Count.data)=c("Cell","Time","AreaSurveyed","Count")
Count.data[,1:2]=Mapping[,1:2]
Count.data[,3]=Area.trans
for(irow in 1:nrow(Count.data)){
  Which.entries=which(Day.hs==Mapping[irow,2] & Cell.id.hs==Mapping[irow,1])
  if(length(Which.entries)>0)Count.data[irow,"Count"]=sum(as.numeric(Hotspots[Which.entries,]@data[,"numseals"]))
}
#cat(paste("Total hotspots: ", length(Cell.id.hs),". Total included in output data: ",sum(Count.data[,"Count"])))

Effort=list(Mapping=Mapping,Area.trans=Area.trans,Area.hab=Area.hab,Count.data=Count.data)
save(Effort,file="Effort2012_BOSSst_bearded_29Apr2014.Rdata")


