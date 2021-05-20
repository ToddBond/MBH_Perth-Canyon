
#### MBH Design Griffen ----

library( rgdal)
library( sp)
library( raster)
library( MBHdesign)
library( pdist)
library( rgeos)

# clear environment ----
rm(list = ls())


# Set names ----
study <- "Griffen_MBH"

platform <- "BOSS"

design.version <- "v4"


# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
o.dir <- paste(w.dir, "outputs", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')



# Read in the inclusion probs ----
inclProbs <- raster(paste(d.dir, "inclProbs-Griffen_MBH-BOSS-v2.tif", sep='/'))
plot(inclProbs)
# check sun of incl probs --
cellStats(inclProbs, 'sum')


inclProbs <- setValues( inclProbs, values( inclProbs) / sum( values( inclProbs), na.rm=TRUE))
plot(inclProbs)
# check sun of incl probs --
cellStats(inclProbs, 'sum')

#inclProbs <- setValues( inclProbs, values( inclProbs)*120)


# Read in data ----

dat <- readRDS(paste(d.dir,"Data-Griffen_MBH.RDS", sep ='/'))
rast <- readRDS(paste(d.dir, "rasters-Griffen_MBH-fine.RDS", sep='/'))
zones <- readRDS(paste(d.dir, "Zones-Griffen_MBH.RDS", sep='/'))
straw.nums <- readRDS(paste(d.dir, "StrawmanNumbers-Griffen_MBH-BOSS-v4.RDS", sep='/'))
strata <- raster(paste(d.dir, "Bathy_cuts-Griffen_MBH-BOSS-v2.tif", sep='/'))

# increase staw nums to oversample --
#straw.nums <- straw.nums + 3 # for design 3
straw.nums <- straw.nums  # for design 4
straw.nums


newSites <- list(g1=NULL, g2=NULL, g3=NULL, g4=NULL, g5=NULL, g6=NULL, g7=NULL, g8=NULL, g9=NULL, cs=NULL, cz=NULL, os=NULL, oz=NULL)

for( zz in c("g1","g2","g3","g4","g5","g6","g7","g8","g9","cs","cz","os","oz")){
  print( zz)
  #the number of samples to take (specified minus the legacy number)
  #numby <- floor( (straw.nums[zz])/4)  # for clustered cluster - without legacy sites
  numby <- floor( (straw.nums[zz])) # for not clustered sites
  #numby <- floor( (straw.nums[zz] - numRef[zz])/2)
  #numby <- floor( (straw.nums[zz] - numRef[zz])) # with legacy sites 
  #set up spatial domain
  myZone <- zones[[zz]]
  #if( zz == "AMP"){
  # myZone = zones$AMP - zones$IUCN2
  #set.seed( 747)
  #}
  #tmpIP <- mask( rootInclProbs, myZone)
  tmpIP <- mask( inclProbs, myZone)
  tmpIP <- crop( tmpIP, myZone)
  #take the sample of clusters based on root incl probs
  newSites[[zz]] <- quasiSamp( n=numby, potential.sites=coordinates( tmpIP), inclusion.probs=values(tmpIP), nSampsToConsider=5000)
  
  #plotting (maybe remove at a later date?)
  tmpIPFull <- mask( inclProbs, myZone)
  tmpIPFull <- crop( tmpIPFull, myZone)
  plot( tmpIPFull)
  #plot( legacySites, add=TRUE, pch=1, col='red')
  points( newSites[[zz]][,c("x","y")], pch=20, col='black')
}
newSites <- do.call( "rbind", newSites)
head(newSites)
newSites.sp <- SpatialPointsDataFrame( coords=newSites[,c("x","y")], data=newSites, proj4string=CRS(proj4string(inclProbs)))


plot(inclProbs)
plot(newSites.sp, add=T)



## Check they are 400m appart ----
### Make sure the clusters centres are ~ 1 km apart ----

## Get CRS in utm ----
crs1 <- CRS("+init=epsg:32750") # WGS 84 / UTM zone 50S


## transform the points into UTM --
p1u <- spTransform(newSites.sp, crs1)

## calculate if 2 points fall within 1500 m of eachother ----
# https://gis.stackexchange.com/questions/102796/remove-points-within-x-distance

dist1 <- gDistance(p1u, byid =T)
dist1
max(dist1)
min(dist1[dist1 > 0]) # minimum distance other than 0

## p1 ----
p1_matrix <- gWithinDistance(p1u, dist = 100, byid = TRUE)
diag(p1_matrix) <- NA
p1_matrix

# extract the upper triangular part of matrix and use the column sums as a criterion to remove the points:

p1_matrix[lower.tri(p1_matrix, diag=TRUE)] <- NA
p1_matrix

colSums(p1_matrix, na.rm=TRUE) == 0
v1 <- colSums(p1_matrix, na.rm=TRUE) == 0
p1u[v1, ] # 98 features left

remaining.sites <- p1u[v1, ]
remaining.sites <- spTransform(remaining.sites, proj4string(inclProbs))

# plot --
plot(inclProbs)
plot(strata)
plot(zones$gAll, add=T)
plot(zones$cs, add=T)
plot(zones$cz, add=T)
plot(zones$os, add=T)
plot(zones$oz, add=T)
plot(remaining.sites, pch = 20, add=T) # 41
#remaining.sites$zone
head(remaining.sites)
id <- paste(1:120)
remaining.sites$uniqueID <- id
head(remaining.sites)
proj4string(remaining.sites)

rem.sites <- spTransform(remaining.sites, proj4string(zones$All))

rem.sites <- raster::extract(zones$All, remaining.sites, xy = TRUE, sp = TRUE)
head(rem.sites)
class(rem.sites)
str(rem.sites)
rem.sites$name <- as.factor(rem.sites$name)
names(rem.sites)
rem.sites <- rem.sites[,c(6,55)]
head(rem.sites)
length(rem.sites$Infra.1)

remaining.sites$infra <- rem.sites$Infra.1
remaining.sites$name <- rem.sites$name
head(remaining.sites)
# save new sites ----

writeOGR(remaining.sites, o.dir, paste(platform, study, design.version, sep='-'), driver = "ESRI Shapefile", overwrite = TRUE)



###   ###   ###   ###

# check distance from structure ----

# read Griffin structure ---
gr <- readOGR(paste(s.dir, "Griffin.shp", sep='/'))
plot(gr)

gru <- spTransform(gr, crs1)

rem.sites.utm <- spTransform(remaining.sites, crs1)

dist.to.str <- apply(gDistance(rem.sites.utm, gru, byid = TRUE), 2, min)

dfs <- as.data.frame(dist.to.str)
min(dfs) 
max(dfs)
hist(dfs[,1], main = paste(platform, study, design.version, sep=' '), xlab = "distance (m)")
hist(dfs[,1], breaks = c(0, 500, 1000, 10000), freq = T)
