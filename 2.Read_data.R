###   ###   ###   Read data for MBH analysis    ###   ###   ###


# clear environment ----
rm(list = ls())


# libraries ----
library( rgdal)
library( sp)
library( raster)

# Set names ----
study <- "Griffen_MBH"



# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')



# Read polygons of NPZs and adjacent areas ----
g1 <- readOGR(paste(s.dir, "EZ_G1.shp", sep = '/'))
g2 <- readOGR(paste(s.dir, "EZ_G2.shp", sep = '/'))
g3 <- readOGR(paste(s.dir, "EZ_G3.shp", sep = '/'))
g4 <- readOGR(paste(s.dir, "EZ_G4.shp", sep = '/'))
g5 <- readOGR(paste(s.dir, "EZ_G5.shp", sep = '/'))
g6 <- readOGR(paste(s.dir, "EZ_G6.shp", sep = '/'))
g7 <- readOGR(paste(s.dir, "EZ_G7.shp", sep = '/'))
g8 <- readOGR(paste(s.dir, "EZ_G8.shp", sep = '/'))
g9 <- readOGR(paste(s.dir, "EZ_G9.shp", sep = '/'))
cz <- readOGR(paste(s.dir, "CautionaryZone.shp", sep = '/'))
oz <- readOGR(paste(s.dir, "OpenAccessZone.shp", sep = '/'))

# check crs --
proj4string(g1)
proj4string(g2)
proj4string(g3)
proj4string(g4)
proj4string(g5)
proj4string(g6)
proj4string(g7)
proj4string(g8)
proj4string(g9)
proj4string(cz)
proj4string(oz)


# Prepare list of polygons ----
zones <- list()
zones$g1 <- g1
zones$g2 <- g2
zones$g3 <- g3
zones$g4 <- g4
zones$g5 <- g5
zones$g6 <- g6
zones$g7 <- g7
zones$g8 <- g8
zones$g9 <- g9
zones$cz <- cz
zones$oz <- oz

# join both polygons
zones$gAll <- raster::union(zones$g1 <- g1, zones$g2)
zones$gAll <- raster::union(zones$gAll, zones$g3)
zones$gAll <- raster::union(zones$gAll, zones$g4)
zones$gAll <- raster::union(zones$gAll, zones$g5)
zones$gAll <- raster::union(zones$gAll, zones$g6)
zones$gAll <- raster::union(zones$gAll, zones$g7)
zones$gAll <- raster::union(zones$gAll, zones$g8)
zones$gAll <- raster::union(zones$gAll, zones$g9)



#intial look to see area
plot( zones$oz, border='black')
plot( zones$cz, add=TRUE, col='orange')
plot( zones$gAll, add=TRUE, col='green')


# Save zones rds ----

saveRDS(zones, file= paste0(paste(d.dir, paste("Zones" , study, sep='-'), sep='/'), ".RDS"))


## Read raster data ----
ders <- stack(paste(r.dir, "Griffen_sea-terrain.tif", sep ='/'))
plot(ders)
names(ders) <-  c("depth","slope", "tpi", "flowdir", "roughness", "aspect")

# Save rasters rds ----
griffen_rasters <- list()
griffen_rasters$depth <- ders$depth
griffen_rasters$slope <- ders$slope

# Save raster data ----
saveRDS(griffen_rasters, file= paste0(paste(d.dir, paste("rasters" , study, sep='-'), sep='/'), ".RDS"))



## Converting polygons to a common raster ----
b <- ders$depth
plot(b)
b1 <- b*(-1)
plot(b1)

t <- readOGR(paste(s.dir, "test-rectangle.shp", sep='/'))

###       ###       ### this takes a while for fine res data  ###      ###       ###
g1_raster <- rasterize(x=zones$g1, y=b, field=zones$g1@data[,1], bkg.value=NA, fun="first")
plot(g1_raster)
oz_raster <- rasterize(x=zones$oz, y=b1, field=zones$oz@data[,1], bkg.value=NA, fun="first")
plot(oz_raster)
g1_raster <- rasterize(x=zones$g1, y=b1)
plot(g1_raster)
t_raster <- rasterize(x=t, y=b1)
plot(t_raster)
out9_raster <- rasterize(x=zones$out9, y=b9, field=zones$out9@data[,1], bkg.value=NA, fun="first")
plot(out9_raster)


#convert and combine --
tmp1 <- as.data.frame( npz6_raster, xy=TRUE)
tmp2 <- as.data.frame( npz9_raster, xy=TRUE)
tmp3 <- as.data.frame( out6_raster, xy=TRUE)
tmp4 <- as.data.frame( out9_raster, xy=TRUE)
tmp5 <- as.data.frame( b6, xy=TRUE)
tmp6 <- as.data.frame( b9, xy=TRUE)
tmp7 <- as.data.frame( s6, xy=TRUE)
tmp8 <- as.data.frame( s9, xy=TRUE)

# Join data for NPZ6 and adjacent analysis --

# NPZ 6 --
npz6Dat <- cbind( tmp1, tmp3[,3])
npz6Dat <- cbind( npz6Dat, tmp5[,3])
npz6Dat <- cbind( npz6Dat, tmp7[,3])
head(npz6Dat)

# NPZ9 --
npz9Dat <- cbind( tmp2, tmp4[,3])
npz9Dat <- cbind( npz9Dat, tmp6[,3])
npz9Dat <- cbind( npz9Dat, tmp6[,3])
head(npz9Dat)

# Set column names --
df.names <- c("Eastern", "Northing", "npz6", "out6", "depth", "slope")

names(npz6Dat) <- df.names
names(npz9Dat) <- df.names


# Save raster dfs rds ----
saveRDS(npz6Dat, file= paste(d.dir, "npz6Dat_forInNOutNPZ.RDS", sep='/'))
saveRDS(npz9Dat, file=paste(d.dir, "npz9Dat_forInNOutNPZ.RDS", sep ='/'))




