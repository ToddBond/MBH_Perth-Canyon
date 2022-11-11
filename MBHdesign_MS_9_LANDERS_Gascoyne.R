## ----libraryLoad, echo=FALSE, cache=FALSE, eval=TRUE---------------------

##############################################################################
####  load packages used in this document                                 ####
##############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library( MBHdesign) #for spatially-balanced designs
  library( raster)  #for spatial data manipulation
  library( colorRamps)  #for nice colours in maps
  library(tidyverse)
})

w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
r.dir <- paste(w.dir, "rasters", sep='/')

##############################################################################
####  Setup the bathymetry data and the inclusion probs from that.        ####
##############################################################################

## Bathymetry data for Hippolyte Islands, Tasmania, Australia
bathy <- raster(paste(r.dir, "Gascoyne_GEBCO_UTM_survey.tif", sep='/'), col=heat.colors)
#bathy[bathy < 0.5] <- NA
plot(bathy)
bathy.df <- as.data.frame( rasterToPoints( bathy))%>%
  filter(!Gascoyne_GEBCO_UTM_survey == 0, !Gascoyne_GEBCO_UTM_survey > -1000)

## Define 4 depth categories to sample within
depthBins <- c(-5700, -3500,- 2000, -1000)
bathy.df$category <- cut( bathy.df$Gascoyne_GEBCO_UTM_survey, depthBins, include.lowest=TRUE)

## Define inclusion probabilities s so that there is the expectation of the same number of samples witin each category
incProb.tmp <- 1 / table( bathy.df$category)  #equal expected sample size for each category.
incProb.tmp <- incProb.tmp / length( levels( bathy.df$category))  #Sum to 1 over all cells
bathy.df$inclusion.prob <- NA
for( ii in 1:length( incProb.tmp)) 
  bathy.df[bathy.df$category == names( incProb.tmp)[ii],"inclusion.prob"] <- incProb.tmp[ii]
incProb <- rasterFromXYZ( bathy.df[,c("x","y","inclusion.prob")], crs=crs( bathy))

## Plot for visual interpretation
#par( mfrow=c(1,2))
#bathymetry
plot( bathy, col=matlab.like2(n=100))
mtext( "A)", adj=0.02, line=-2)
title( main="Gascoyne Bathymetry (m)", xlab="x", ylab="y")
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
#points( legacySites, pch=17, col='darkmagenta', cex=2)
#inclusion probs / shown only to depth bins (no scale)
plot(log(incProb), col=matlab.like2(n=100), legend=FALSE)
mtext( "B)", adj=0.02, line=-2)
title( main="Depth Bins", xlab="x", ylab="")





## ----PointDesigns1, cache=TRUE, echo=FALSE, tidy=FALSE-------------------

#repeatable code
set.seed( 747)  #a 747 is a big plane

#the number of samples to take
nSamp <- 45


# ## ----showEvenCall, cache=TRUE, echo=TRUE, tidy=FALSE---------------------
# #### A spatially-balanced sample within the study area (not depth-related)
# evenSample <- quasiSamp( n=nSamp, potential.sites=bathy.df[,c("x","y")])
# 
# ## ----showSurveySites, cache=TRUE, echo=FALSE, tidy=FALSE-----------------
# knitr::kable( evenSample[sample( 1:nrow( evenSample),4),], booktabs=TRUE, row.names=FALSE, align='c', caption="Four survey sites selected by the function quasiSamp for the even inclusion probability design. In order, the columns are: 'x' coordinate (e.g. longitude), 'y' coordinate (e.g. latitude), the inclusion probability that the site was selected with, and the row number from the potential.sites input argument. Inclusion probabilities, over all potential sites, will sum to 1 by construction.")


## ----showUnevenCall, cache=TRUE, echo=TRUE, tidy=FALSE-------------------
#### A spatially-balanced sample with shallow locations preferred
unevenSample <- quasiSamp( n=nSamp, potential.sites=bathy.df[,c("x","y")],
                           inclusion.probs=bathy.df[,"inclusion.prob"])

#a 3x1 figure with margins that aren't too large
# #plot even design on background of depth
# plot( bathy, col=matlab.like2(n=100), main="Even Sample", legend=FALSE)
# contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
# points( evenSample[,c("x","y")], pch=20, col='black')
# mtext( "A)", adj=0.02, line=-2)
# 
# 


plot( bathy, col=matlab.like2(n=100), main="Uneven Sample", legend=FALSE)
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( unevenSample[,c("x","y")], pch=20, col='black')
#mtext( "A)", adj=0.02, line=-2)


#export samples as shapefile
coordinates(unevenSample)=~x+y
proj4string(unevenSample)<- CRS( "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs")
#LLcoor<-spTransform(WGScoor,CRS("+proj=longlat"))
raster::shapefile(unevenSample, "Gasc_stratified-single.shp", overwrite = TRUE)

#save contours as shapefile
bathy.contour <- rasterToContour(bathy, levels =depthBins)
proj4string(bathy.contour)<- CRS( "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs")
raster::shapefile(bathy.contour, "Gasc_bathy_contour.shp", overwrite = TRUE)

crs(bathy.contour)

##############################################################################
####  Setup inclusion prob raster (and data.frame) with NAs               ####
####  Aggregates first to avoid excessive computation                     ####
##############################################################################

nSamp <- 15

#transectSamp() function requires a grid with missing values given by NA 
#     (not by omission as is often done with rasters)
fullInclusionProbs <- rasterFromXYZ( bathy.df[,c("x","y","inclusion.prob")])
#think carefully about na.rm. FALSE (TRUE) means that the sample area shrinks (grows)
fullInclusionProbs <- aggregate( fullInclusionProbs, fact=5, fun=mean, na.rm=FALSE) 
#spatial subsample to reduce computation (for this example)
incl.prob <- cbind( coordinates( fullInclusionProbs), 
                    inclusion.prob=values( fullInclusionProbs))
#transectSamp() also requires increasing ordering with x moving quicker than y
incl.prob <- incl.prob[ order( incl.prob[,"y"], incl.prob[,"x"]),]


## ----showTransectCall, cache=TRUE, echo=TRUE, tidy=FALSE-----------------
#The representation of transects and other algorithm controls
control <- list( transect.pattern='line', transect.nPts=3, nRotate=21, 
                 line.length=10000, mc.cores=8)
#a spatially-balanced sample with shallow locations preferred
transSample <- transectSamp( n=nSamp, potential.sites=incl.prob[,c("x","y")], 
                             inclusion.probs=incl.prob[,"inclusion.prob"], control=control)


## ----showTransPoints, cache=TRUE, echo=FALSE, tidy=FALSE, caption="First six points (of 120) on transects selected by transSamp. Stored in the second element of the return object. In order, the columns are: transect number, the coordinates of the transects' midpoints (2 directions), the points on the transect, the user-defined inclusion probability and the (internal) probability of selection to maintain the user-defined probabilities."----
#knitr::kable( transSample[[2]][sample(1:nrow(transSample[[2]]),4),], booktabs=TRUE, row.names=FALSE, align='c', caption="Four points (of 120) on transects selected by transSamp. Stored in the second element of the return object. In order, the columns are: transect number, the coordinates of the transects' midpoints, the compass bearing of the transect, the points on the transect, the user-defined inclusion probability and the (internal) probability of selection to maintain the user-defined probabilities.")


#plot transect sample design (no constraints) on background of depth
plot( bathy, col=matlab.like2(n=100), main="Random Transect Sample", legend=FALSE)
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( transSample$points[,c("x","y")], pch=20, col='black')
mtext( "A)", adj=0.02, line=-2)


#export samples as shapefile
WGScoor<-  transSample$points
coordinates(WGScoor)=~x+y
proj4string(WGScoor)<- CRS( "+proj=utm +zone=49 +south +datum=WGS84 +units=m +no_defs")
#LLcoor<-spTransform(WGScoor,CRS("+proj=longlat"))
raster::shapefile(WGScoor, "Gasc_TranSurveyPoints.shp", overwrite = TRUE)






