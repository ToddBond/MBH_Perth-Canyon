library(raster)


rm(list = ls())


# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
p.dir <- paste(w.dir, "plots", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')


b <- raster(paste(r.dir, "Perth_Canyon_Marine_Park_Bathymetry_2018_40m_cog.tif", sep='/'), proj4string = CRS ( "+proj=utm+zone=49 +datum=WGS84" ))
plot(b)

b.deep<-b
b.deep[b.deep > -3500] <- NA
plot(b.deep)

b.mid<-b
b.mid[b.mid < -3500] <- NA
b.mid[b.mid > -1500] <- NA
plot(b.mid)


###########################################################################
#### Read in Data from spatial data (.asc here) and Organise ####
#### Foster et al. NESP Biodiversity Hub Field Manuals ####
###########################################################################
##if you don't have MBHdesign installed, please do so using
# install.packages( "MBHdesign")
#Load required packages
library (MBHdesign) #For spatial sampling
library (fields) #for lots of things, but for plotting in this example
library(maps)
library (sp) #for reading the ascii file of cropped depths for the reserve
#Set a seed for reproducability
set.seed ( 666 )
#Read in depth as an ESRI asc file as requested by the sp package.
#This file contains long, lat and depth
#This path/file only exists on the first author's system
# you will need to change it if running this code
#the projection will need to be changed for each region too
#bth.orig.grid <- read.asciigrid("./ExampleGovIsland/gov_bth.asc", proj4string = CRS("+proj=utm +zone=55 +datum=WGS84"))

#bth.orig.grid <- read.asciigrid ( "gov_bth.asc" , proj4string = CRS ( "+proj=utm+zone=55 +datum=WGS84" ))

#convert to a data.frame for ease
#DepthMat <- as.matrix ( bth.orig.grid)
DepthMat <- as.matrix(b.deep)

# 
# bth.orig.grid <- as.data.frame (
#   cbind ( coordinates ( bth.orig.grid), as.numeric ( DepthMat)))


b.deep.1 <- as.data.frame (
  cbind ( coordinates (b.deep), as.numeric ( DepthMat)))



colnames(b.deep.1) <- c ( "Easting" , "Northing" , "Depth" )
b.deep.2 <- b.deep.1[ order(b.deep.1$Northing,b.deep.1$Easting),]

#Setting up plotting for now and later
uniqueEast <- unique(b.deep.2$Easting)
uniqueNorth <- unique(b.deep.2$Northing)
ELims <- range(na.exclude(b.deep.2)$Easting)
NLims <- range(na.exclude(b.deep.2)$Northing)

#Fix up ordering issue
DepthMat <- DepthMat[, rev(1:ncol(DepthMat))]

#plot it to see what we are dealing with.
fields::image.plot(uniqueEast, uniqueNorth, DepthMat,
             xlab= "Easting" , ylab= "Northing" , main= "Perth Canyon" ,
             legend.lab= "Depth (m)" , asp=1 , ylim= NLims, xlim= ELims,
             col= rev(fields::tim.colors()))


bth.orig.grid <- b.deep.2


#number of samples
n <- 30
#take the sample
samp_spatialOnly <- quasiSamp ( n= n, dimension=2 ,
                                potential.sites = bth.orig.grid[, c ( "Easting" , "Northing" )],
                                inclusion.probs= ! is.na ( bth.orig.grid$Depth))
with ( bth.orig.grid, fields::image.plot ( uniqueEast, uniqueNorth, DepthMat,
                                   xlab= "Easting" , ylab= "Northing" , main= "Spatially Balanced Sample" ,
                                   legend.lab= "Depth (m)" , asp=1 , ylim= NLims, xlim= ELims))



#



sd <- quasiSamp ( n= n,potential.sites = bth.orig.grid[, c ( "Easting" , "Northing" )])








#number of samples
n <- 15
#take the sample
samp_spatialOnly <- quasiSamp(n= n, dimension=2,
                                potential.sites = b.deep[,c("Easting","Northing" )],
                                inclusion.probs= !is.na (b.deep$Depth))


with ( bth.orig.grid, image.plot (uniqueEast, uniqueNorth, DepthMat,
                                   xlab= "Easting" , ylab= "Northing" , main= "Spatially Balanced Sample" ,
                                   legend.lab= "Depth (m)" , asp=1 , ylim= NLims, xlim= ELims,
                                   col= rev ( tim.colors ())))
points ( samp_spatialOnly[, c ( "Easting" , "Northing" )], pch=20 , cex=2 )
write.csv (samp_spatialOnly, file= "spatialOnly.csv" , row.names= FALSE )
