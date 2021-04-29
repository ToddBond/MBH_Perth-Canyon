
###   ###   ###    Prepare spatial data for analysis    ###   ###   ###

# Libraries ----
library(raster)
library(rgdal)
library(sp)
library(ggplot2)
#library(pals)
library(RColorBrewer)


# clear environment ----
rm(list = ls())


# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')


## Load bathy and park data ----

# load West coast bathy --
b <- raster(paste(r.dir, "Griffen_Bathy.tif", sep='/'))
plot(b)

# Remove the bathy deeper than 200m ----
b[b < -200] <- NA
plot(b)

depth <- b

slope <- terrain(depth, "slope")
plot(slope)

tpi <- terrain(depth, "TPI")
plot(tpi)

flowdir <- terrain(depth, "flowdir")
plot(flowdir)

roughness <- terrain(depth, "roughness")
plot(roughness)

aspect <- terrain(depth, "aspect")
plot(aspect)

sea.terrain <- stack(depth, slope, tpi, flowdir, roughness, aspect)
plot(sea.terrain)
names(sea.terrain) <- c("depth","slope", "tpi", "flowdir", "roughness", "aspect")

# save derivatives of sea terrain ----
writeRaster(sea.terrain, paste(r.dir, "Griffen_sea-terrain.tif", sep='/'))


## Resample bathy ----

b50 <- raster::disaggregate(b, 5)
plot(b50)

slope50 <- raster::disaggregate(slope, 5)
plot(slope50)

tpi50 <- raster::disaggregate(tpi, 5)
plot(tpi50)

sea.terrain <- stack(b50, slope50, tpi50)
plot(sea.terrain)
names(sea.terrain) <- c("depth","slope", "tpi")

# save derivatives of sea terrain ----
writeRaster(sea.terrain, paste(r.dir, "Griffen_sea-terrain_fine.tif", sep='/'))
