require(RPostgreSQL)
require(rgeos)
require(rgdal)

# use .Rprofile to store credentials
# for example:
#
# options(pepgeo_usr = "username")
# options(pepgeo_pwd = "password")
# options(pepgeo_ip = "ip_address")

# get our credentials from our .Rprofile
#pepgeo_usr <- getOption('pepgeo_usr')
#pepgeo_pwd <- getOption('pepgeo_pwd')
#pepgeo_ip <- getOption('pepgeo_ip')

pepgeo_usr='connp'
pepgeo_pwd='$e@l'
pepgeo_ip <- '161.55.120.122'

# connect to the PostgreSQL db, pepgeo
conn = dbConnect(
  dbDriver("PostgreSQL"), dbname="pepgeo", host=pepgeo_ip, port="5432", 
  user=pepgeo_usr, password=pepgeo_pwd
)

# create our query to get all hotspot points
strSQL <- paste("SELECT objectid,hotspotid,sealid,assignid,hotspot_found,",
                "hotspot_type,dt_utc,flightid,cameraid,numseals,",
                "species,species_conf,ST_AsText(geom) as wkt_geometry", 
                "FROM surveys.hotspots")

# pull out the SRID so we can set the proper projection
srid <- dbGetQuery(conn,"SELECT Find_SRID('surveys','hotspots','geom')")[1,]
EPSG = make_EPSG()
proj = EPSG[which(EPSG$code == srid), "prj4"]

# get hotspots data and store in a df
# we set the rownames to objectid to assist in the rbind later
hotspots.df = dbGetQuery(conn,strSQL)
rownames(hotspots.df) <- hotspots.df$objectid

# create a simple pre-defined list to hold our sp objects
sp.list <- vector(mode="list",length=nrow(hotspots.df))

# fill the sp.list with objects created via readWKT in the rgeos package
for (n in 1:length(sp.list)) {
  sp.list[[n]] <- SpatialPoints(readWKT(hotspots.df$wkt_geometry[n],
                                             id=hotspots.df$objectid[n])) 
}

# rbind all the spatial objects into one
hotspots.sp <- do.call("rbind",sp.list)

# create a spdf by merging with the source df
# we'll clean up hotspots.df, first, to remove uneeded columns
hotspots.df <- subset(hotspots.df,select = -c(objectid,assignid,wkt_geometry))
hotspots.spdf <- SpatialPointsDataFrame(hotspots.sp,hotspots.df)
proj4string(hotspots.spdf) <- CRS(proj)

# get on-effort track lines; these on-effort track lines include filter
# parameters to NOT include periods where altitude is between
# 750 and 1100 feet and roll is +/2.5 degrees
strSQL <- paste("SELECT objectid,flightid,ST_AsText(track_geom) as wkt_geometry", 
                "FROM surveys.boss_oneffort_tracks")

srid <- dbGetQuery(conn,"SELECT Find_SRID('surveys','boss_oneffort_tracks','track_geom')")[1,]

EPSG = make_EPSG()
proj = EPSG[which(EPSG$code == srid), "prj4"]

on_effort_tracks.df = dbGetQuery(conn,strSQL)
rownames(on_effort_tracks.df) <- on_effort_tracks.df$objectid

# create a simple pre-defined list to hold our sl objects
sl.list <- vector(mode="list",length=nrow(on_effort_tracks.df))

# fill the sl.list with objects created via readWKT
# slightly different approach than with points to accomodate SpatialLines
for (n in 1:length(sl.list)) {
  sl.list[[n]] <- readWKT(on_effort_tracks.df$wkt_geometry[n])@lines[[1]]
  sl.list[[n]]@ID <- as.character(on_effort_tracks.df$objectid[n])
}

# create a sldf by merging with the source df
# we'll clean up on_effort_tracks.df, first, to remove uneeded columns
on_effort_tracks.df <- subset(on_effort_tracks.df,select=-c(objectid,wkt_geometry))
on_effort_tracks.sl <- SpatialLines(sl.list)
on_effort_tracks.sldf <- SpatialLinesDataFrame(on_effort_tracks.sl,on_effort_tracks.df)

proj4string(on_effort_tracks.sldf) <- CRS(proj)

# create our query to get all fmclog points where altitude is between
# 750 and 1100 feet and roll is +/2.5 degrees
strSQL <- paste("SELECT objectid,flightid,gpsalt,gpsspd,gpshead,",
                "roll,pitch,baroalt,effort,dt_utc,effort,",
                "ST_AsText(geom) as wkt_geometry", 
                "FROM surveys.fmclogs",
                "WHERE effort = 'On' AND gpsalt BETWEEN 750 AND 1100 AND (roll >=357.5 OR roll <=2.5)")

# pull out the SRID so we can set the proper projection
srid <- dbGetQuery(conn,"SELECT Find_SRID('surveys','fmclogs','geom')")[1,]
EPSG = make_EPSG()
proj = EPSG[which(EPSG$code == srid), "prj4"]

# get fmclogs data and store in a df
# we set the rownames to objectid to assist in the rbind later
fmclogs.df = dbGetQuery(conn,strSQL)
rownames(fmclogs.df) <- fmclogs.df$objectid

# create a simple pre-defined list to hold our sp objects
sp.list <- vector(mode="list",length=nrow(fmclogs.df))

# fill the sp.list with objects created via readWKT in the rgeos package
# given the size of this dataset, it will take a while
for (n in 1:length(sp.list)) {
  sp.list[[n]] <- SpatialPoints(readWKT(fmclogs.df$wkt_geometry[n],
                                        id=fmclogs.df$objectid[n])) 
}

# rbind all the spatial objects into one
fmclogs.sp <- do.call("rbind",sp.list)

# create a spdf by merging with the source df
# we'll clean up fmclogs.df, first, to remove uneeded columns
fmclogs.df <- subset(fmclogs.df,select = -c(objectid,wkt_geometry))
fmclogs.spdf <- SpatialPointsDataFrame(fmclogs.sp,fmclogs.df)
proj4string(fmclogs.spdf) <- CRS(proj)
save.image(file="BOSS_2012Effort_22Apr14.Rdata")
