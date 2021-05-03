

library( rgdal)
library( sp)
library( raster)
library( MBHdesign)
library( pdist)
library( rgeos)

# clear environment ----
rm(list = ls())





# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
o.dir <- paste(w.dir, "outputs", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')


# Set names ----
study <- "Griffen_MBH"

platform <- "BOSS"

design.version <- "v3"


# read design ----

remaining.sites <- readOGR(paste(o.dir, "BOSS-Griffen_MBH-v3.shp", sep='/'))


# read Griffin structure ---
gr <- readOGR(paste(s.dir, "Griffin.shp", sep='/'))
plot(gr)


# transform to utm --
## Get CRS in utm ----
crs1 <- CRS("+init=epsg:32750") # WGS 84 / UTM zone 50S

gru <- spTransform(gr, crs1)

rem.sites.utm <- spTransform(remaining.sites, crs1)

# Calculate distance to structure

dist.to.str <- apply(gDistance(rem.sites.utm, gru, byid = TRUE), 2, min) # each point to structure minumum distance

dfs <- as.data.frame(dist.to.str)
min(dfs) 
max(dfs)
hist(dfs[,1], main = paste(platform, study, design.version, sep=' '), xlab = "distance (m)")
#hist(dfs[,1], breaks = c(0, 500, 1000, 10000), freq = T)
