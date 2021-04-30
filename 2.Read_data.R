###   ###   ###   Read data for MBH analysis    ###   ###   ###



# libraries ----
library( rgdal)
library( sp)
library( raster)


# clear environment ----
rm(list = ls())



# Set names ----
study <- "Griffen_MBH"



# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')



# Read polygons ----
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
ders <- stack(paste(r.dir, "Griffen_sea-terrain_fine.tif", sep ='/'))
plot(ders)
names(ders) <-  c("depth","slope", "tpi")
plot(ders)

# Save rasters rds ----
griffen_rasters <- list()
griffen_rasters$depth <- ders$depth
griffen_rasters$slope <- ders$slope

# Save raster data ----
saveRDS(griffen_rasters, file= paste0(paste(d.dir, paste("rasters" , study, "fine-complete", sep='-'), sep='/'), ".RDS"))



## Converting polygons to a common raster ----
b <- ders$depth
plot(b)


###       ###       ### this takes a while for fine res data  ###      ###       ###
#g1_raster <- rasterize(x=g1, y=b, field=g1@data[,1], background=NA, fun="max")
g1_raster <- rasterize(x=g1, y=b, field=1, background=NA, fun="first")
plot(g1_raster)
g2_raster <- rasterize(x=g2, y=b, field=1, background=NA, fun="first")
plot(g2_raster)
g3_raster <- rasterize(x=g3, y=b, field=1, background=NA, fun="first")
plot(g3_raster)
g4_raster <- rasterize(x=g4, y=b, field=1, background=NA, fun="first")
plot(g4_raster)
g5_raster <- rasterize(x=g5, y=b, field=1, background=NA, fun="first")
plot(g5_raster)
g6_raster <- rasterize(x=g6, y=b, field=1, background=NA, fun="first")
plot(g6_raster)
g7_raster <- rasterize(x=g7, y=b, field=1, background=NA, fun="first")
plot(g7_raster)
g8_raster <- rasterize(x=g8, y=b, field=1, background=NA, fun="first")
plot(g8_raster)
g9_raster <- rasterize(x=g9, y=b, field=1, background=NA, fun="first")
plot(g9_raster)
cz_raster <- rasterize(x=cz, y=b, field=1, background=NA, fun="first")
plot(cz_raster)
oz_raster <- rasterize(x=oz, y=b, field=1, bkg.value=NA, fun="first")
plot(oz_raster)



#convert and combine --
tmp1 <- as.data.frame( g1_raster, xy=TRUE)
tmp2 <- as.data.frame( g2_raster, xy=TRUE)
tmp3 <- as.data.frame( g3_raster, xy=TRUE)
tmp4 <- as.data.frame( g4_raster, xy=TRUE)
tmp5 <- as.data.frame( g5_raster, xy=TRUE)
tmp6 <- as.data.frame( g6_raster, xy=TRUE)
tmp7 <- as.data.frame( g7_raster, xy=TRUE)
tmp8 <- as.data.frame( g8_raster, xy=TRUE)
tmp9 <- as.data.frame( g9_raster, xy=TRUE)
tmp10 <- as.data.frame( cz_raster, xy=TRUE)
tmp11 <- as.data.frame( oz_raster, xy=TRUE)

# Join data for NPZ6 and adjacent analysis --

# Griffen data --
GriffenDat <- cbind( tmp1, tmp2[,3])
GriffenDat <- cbind( GriffenDat, tmp3[,3])
GriffenDat <- cbind( GriffenDat, tmp4[,3])
GriffenDat <- cbind( GriffenDat, tmp5[,3])
GriffenDat <- cbind( GriffenDat, tmp6[,3])
GriffenDat <- cbind( GriffenDat, tmp7[,3])
GriffenDat <- cbind( GriffenDat, tmp8[,3])
GriffenDat <- cbind( GriffenDat, tmp9[,3])
GriffenDat <- cbind( GriffenDat, tmp10[,3])
GriffenDat <- cbind( GriffenDat, tmp11[,3])
head(GriffenDat)


# Set column names --
df.names <- c("Eastern", "Northing", "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "cz", "oz")

names(GriffenDat) <- df.names



# Save raster dfs rds ----
saveRDS(GriffenDat, file= paste0(paste(d.dir, paste("Data" , study, sep='-'), sep='/'), ".RDS"))




