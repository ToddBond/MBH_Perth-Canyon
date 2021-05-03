### check bathy of points ###

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
o.dir <- paste(w.dir, "outputs", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')

# Read structures  ----
g1 <- readOGR(paste(s.dir, "G1.shp", sep = '/'))
g2 <- readOGR(paste(s.dir, "G2.shp", sep = '/'))
g3 <- readOGR(paste(s.dir, "G3.shp", sep = '/'))
g4 <- readOGR(paste(s.dir, "G4.shp", sep = '/'))
g5 <- readOGR(paste(s.dir, "G5.shp", sep = '/'))
g6 <- readOGR(paste(s.dir, "G6.shp", sep = '/'))
g7 <- readOGR(paste(s.dir, "G7.shp", sep = '/'))
g8 <- readOGR(paste(s.dir, "G8.shp", sep = '/'))
g9 <- readOGR(paste(s.dir, "G9.shp", sep = '/'))


# make one shpfile with all strucutures
structures <- raster::union(g1, g2)
structures <- raster::union(structures, g3)
structures <- raster::union(structures, g4)
structures <- raster::union(structures, g5)
structures <- raster::union(structures, g6)
structures <- raster::union(structures, g7)
structures <- raster::union(structures, g8)
structures <- raster::union(structures, g9)

plot(structures)

head(structures)

# Add polygon IDs --
poly.id <- paste("g", 1:9, sep='')
poly.id <- c("g1","g2","g2", "g3","g3", "g4", "g4", "g4", "g4", "g5", "g6", "g7", "g7", "g7", "g7", "g8", "g8", "g9")
structures$polyID <- poly.id

# Read bathymetry ----
sterr <- stack(paste(r.dir, "Griffen_sea-terrain_fine.tif", sep='/'))
b <- sterr$Griffen_sea.terrain_fine.1
plot(b)
names(b) <- "depth"

# Get depth of structures ----
str.depth <- raster::extract(b, structures, sp=TRUE)
head(str.depth)

hist(str.depth$depth, main = "Depth")

# Plot ----
str.df <- as.data.frame(str.depth)
str.df$depth <- str.df$depth*(-1)
min(str.df$depth) # 119
max(str.df$depth) # 142

theme_set(theme_bw())

p<-ggplot(data=str.df, aes(x=polyID, y=depth)) +
  geom_bar(stat="identity") +
  coord_cartesian(ylim=c(110,145))
p
