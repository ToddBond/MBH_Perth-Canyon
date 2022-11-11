## ----libraryLoad, echo=FALSE, cache=FALSE, eval=TRUE---------------------

##############################################################################
####  load packages used in this document                                 ####
##############################################################################




suppressPackageStartupMessages({
  library( MBHdesign) #for spatially-balanced designs
  library( raster)  #for spatial data manipulation
  library( colorRamps)  #for nice colours in maps
})


## ----Setup, cache=TRUE, echo=FALSE, fig.width=10, fig.height=4.5, fig.cap = 'A) Bathymetry of the Hippolyte Rocks survey area. Purple triangles represent the locations of legacy sites that could be incorporated into a survey design. B) The depth bins within which the inclusion probabilities are constant.'----

##############################################################################
####  Setup the bathymetry data and the inclusion probs from that.        ####
##############################################################################

## Bathymetry data for Hippolyte Islands, Tasmania, Australia
bathy <- raster( "HippolyteDepth.tif", col=heat.colors)
bathy.df <- as.data.frame( rasterToPoints( bathy))

## Legacy sites
legacySites <- readRDS("LegacySites.RDS")

## Define 4 depth categories to sample within
depthBins <- c(-80,-30,-45,-60,-5)
bathy.df$category <- cut( bathy.df$HippolyteDepth, depthBins, include.lowest=TRUE)

## Define inclusion probabilities s so that there is the expectation of the same number of samples witin each category
incProb.tmp <- 1 / table( bathy.df$category)  #equal expected sample size for each category.
incProb.tmp <- incProb.tmp / length( levels( bathy.df$category))  #Sum to 1 over all cells
bathy.df$inclusion.prob <- NA
for( ii in 1:length( incProb.tmp)) 
  bathy.df[bathy.df$category == names( incProb.tmp)[ii],"inclusion.prob"] <- incProb.tmp[ii]
incProb <- rasterFromXYZ( bathy.df[,c("x","y","inclusion.prob")], crs=crs( bathy))

## Plot for visual interpretation
par( mfrow=c(1,2))
#bathymetry
plot( bathy, col=matlab.like2(n=100))
mtext( "A)", adj=0.02, line=-2)
title( main="Hippolyte Rocks Bathymetry (m)", xlab="x", ylab="y")
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( legacySites, pch=17, col='darkmagenta', cex=2)
#inclusion probs / shown only to depth bins (no scale)
plot( log( incProb), col=matlab.like2(n=100), legend=FALSE)
mtext( "B)", adj=0.02, line=-2)
title( main="Depth Bins", xlab="x", ylab="")


## ----PointDesigns1, cache=TRUE, echo=FALSE, tidy=FALSE-------------------

#repeatable code
set.seed( 747)  #a 747 is a big plane

#the number of samples to take
nSamp <- 100


## ----showEvenCall, cache=TRUE, echo=TRUE, tidy=FALSE---------------------
#### A spatially-balanced sample within the study area (not depth-related)
evenSample <- quasiSamp( n=nSamp, potential.sites=bathy.df[,c("x","y")])

## ----showSurveySites, cache=TRUE, echo=FALSE, tidy=FALSE-----------------
knitr::kable( evenSample[sample( 1:nrow( evenSample),4),], booktabs=TRUE, row.names=FALSE, align='c', caption="Four survey sites selected by the function quasiSamp for the even inclusion probability design. In order, the columns are: 'x' coordinate (e.g. longitude), 'y' coordinate (e.g. latitude), the inclusion probability that the site was selected with, and the row number from the potential.sites input argument. Inclusion probabilities, over all potential sites, will sum to 1 by construction.")


## ----showUnevenCall, cache=TRUE, echo=TRUE, tidy=FALSE-------------------
#### A spatially-balanced sample with shallow locations preferred
unevenSample <- quasiSamp( n=nSamp, potential.sites=bathy.df[,c("x","y")],
                                inclusion.probs=bathy.df[,"inclusion.prob"])


## ----showAlteredIPCall, cache=TRUE, echo=TRUE, tidy=FALSE----------------
#adjust the inclusion probabilities for the locations of the legacy sites
bathy.df$altered.inclusion.prob <- alterInclProbs( as.matrix( legacySites), 
                potential.sites=bathy.df[,c("x","y")], 
                inclusion.probs=(nSamp-nrow(legacySites))*bathy.df[,"inclusion.prob"], 
                mc.cores=8)
#a spatially-balanced sample with legacy sites and a preference for shallow locations
unevenLegacySample <- quasiSamp( n=nSamp-nrow(legacySites), 
                potential.sites=bathy.df[,c("x","y")], 
                inclusion.probs=bathy.df[,"altered.inclusion.prob"])


## ----PlotPointDes, cache=TRUE, echo=FALSE, tidy=FALSE, fig.width=5, fig.height=8, fig.cap = 'Example point-based designs. A) a spatially-balanced design with equal inclusion probabilities. B) a spatially-balanced design with unequal inclusion probabilities. C) a spatially-balanced design with unequal inclusion probabilities that incorporates 5 legacy sites.'----

##############################################################################
####  Plot even, unequal and legacy site designs                          ####
##############################################################################

#a 3x1 figure with margins that aren't too large
par( mfrow=c(3,1), mar=c(2.1,2.1,3.1,1.1))

#plot even design on background of depth
plot( bathy, col=matlab.like2(n=100), main="Even Sample", legend=FALSE)
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( evenSample[,c("x","y")], pch=20, col='black')
mtext( "A)", adj=0.02, line=-2)

#plot uneven design on background of depth
plot( bathy, col=matlab.like2(n=100), main="Uneven Sample", legend=FALSE)
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( unevenSample[,c("x","y")], pch=20, col='black')
mtext( "B)", adj=0.02, line=-2)

#plot uneven + legacy site design on background of depth
plot( bathy, col=matlab.like2(n=100), 
                main="Uneven Sample with Legacy Sites", legend=FALSE)
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
#the legacy sites
points( legacySites, pch=17, cex=2, col='darkmagenta')
points( unevenLegacySample[,c("x","y")], pch=20, col='black')
mtext( "C)", adj=0.02, line=-2)
legend("bottomleft", legend=c("Legacy Sites","Sample Locations"), 
                pch=c(17,20), pt.cex=c(2,1), col=c("darkmagenta","black"), bty='n')


## ----SetupTrans, cache=TRUE, echo=FALSE, tidy=FALSE----------------------

##############################################################################
####  Setup inclusion prob raster (and data.frame) with NAs               ####
####  Aggregates first to avoid excessive computation                     ####
##############################################################################

#repeatable code
set.seed( 767)  #Still a big plane

nSamp <- 12

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
control <- list( transect.pattern='line', transect.nPts=10, nRotate=21, 
                                              line.length=150, mc.cores=8)
#a spatially-balanced sample with shallow locations preferred
transSample <- transectSamp( n=nSamp, potential.sites=incl.prob[,c("x","y")], 
                            inclusion.probs=incl.prob[,"inclusion.prob"], control=control)


## ----showTransPoints, cache=TRUE, echo=FALSE, tidy=FALSE, caption="First six points (of 120) on transects selected by transSamp. Stored in the second element of the return object. In order, the columns are: transect number, the coordinates of the transects' midpoints (2 directions), the points on the transect, the user-defined inclusion probability and the (internal) probability of selection to maintain the user-defined probabilities."----
knitr::kable( transSample[[2]][sample(1:nrow(transSample[[2]]),4),], booktabs=TRUE, row.names=FALSE, align='c', caption="Four points (of 120) on transects selected by transSamp. Stored in the second element of the return object. In order, the columns are: transect number, the coordinates of the transects' midpoints, the compass bearing of the transect, the points on the transect, the user-defined inclusion probability and the (internal) probability of selection to maintain the user-defined probabilities.")


## ----DownhillTransect, cache=TRUE, echo=FALSE, tidy=FALSE----------------

##############################################################################
####  Generate a downhill transect sample.                                ####
##############################################################################

#repeatable code
set.seed( 787)  #Another a big plane

nSamp <- 12

#transectSamp() function requires a grid with missing values given by NA 
#   (not by omission as is often done for rasters)
fullInclusionProbs <- rasterFromXYZ( bathy.df[,c("x","y","inclusion.prob")])
fullInclusionProbs <- merge( fullInclusionProbs, bathy)
fullRasters <- brick( bathy, fullInclusionProbs)
#decrease resolution to reduce computation time for this example
#think carefully about na.rm -- FALSE will mean that the survey area shrinks
fullRasters <- aggregate( fullRasters, fact=5, fun=mean, na.rm=FALSE)
df <- cbind( coordinates( fullRasters), values( fullRasters))
colnames( df)[3:4] <- c("Bathy","inclusion.prob")
#transectSamp() also requires increasing ordering with x moving quicker than y
df <- df[ order( df[,"y"], df[,"x"]),]

#find those transects that are down-slope. This is a matrix of nCells by nRotate.
descend.type <- findDescendingTrans( potential.sites=df[,c( "x", "y")], 
                bathy=df[,"Bathy"],
                in.area=!is.na( df[,"Bathy"]),
                descend.cutoff=0, control=control)
descend.constr <- ifelse( descend.type=="descend", TRUE, FALSE)
#a spatially-balanced sample with shallow locations preferred
transSample_downhill <- transectSamp( n=nSamp, 
                potential.sites=df[,c("x","y")], 
                inclusion.probs=df[,"inclusion.prob"], 
                constrainedSet=descend.constr, control=control)


## ----WagonwheelTransect, cache=TRUE, tidy=FALSE, echo=FALSE--------------

##############################################################################
####  Generate a 'wagon-wheel' transect sample.                           ####
##############################################################################

#repeatable code
set.seed( 777)  #Another a big plane

nSamp <- 12

#transectSamp() function requires a grid with missing values given by NA
fullInclusionProbs <- rasterFromXYZ( bathy.df[,c("x","y","inclusion.prob")])
fullInclusionProbs <- merge( fullInclusionProbs, bathy)
fullRasters <- brick( bathy, fullInclusionProbs)
#decrease resolution to reduce computation time for this example
#think carefully about na.rm -- FALSE means that the area shrinks, TRUE it grows.
fullRasters <- aggregate( fullRasters, fact=5, fun=mean, na.rm=FALSE)
df <- cbind( coordinates( fullRasters), values( fullRasters))
colnames( df)[3:4] <- c("Bathy","inclusion.prob")
#transectSamp() also requires increasing ordering with x moving quicker than y
df <- df[ order( df[,"y"], df[,"x"]),]

#number of locations to start transects from
nPts <- 20
#lengthen the transects for this example
control$line.length <- 300

#locations to start transects (e.g. top of a seamount).  Choose the nPts highest locations
highLocations <- df[order( df[,"Bathy"], decreasing = TRUE)[1:nPts],c("x","y")]
#find those transects that originate from 1 of the nPts highest points
wagon.trans <- findTransFromPoint( potential.sites=df[,c( "x", "y")], 
                originPoints=highLocations,
                in.area=!is.na( df[,"Bathy"]), control=control)
wagon.constr <- ifelse( wagon.trans=="startsFromPoint", TRUE, FALSE)
#a spatially-balanced sample with shallow locations preferred
transSample_wagon <- transectSamp( n=nSamp, 
                potential.sites=df[,c("x","y")], 
                inclusion.probs=df[,"inclusion.prob"], 
                constrainedSet=wagon.constr, control=control)


## ----PlotTransDes, cache=TRUE, echo=FALSE, tidy=FALSE, fig.width=5, fig.height=8, fig.cap = "Example transect-based designs. A) A spatially-balanced design with unequal inclusion probabilities. B) Like A) but with transects that are downhill. C) Like A) but only with transects that start in the shallowest 20 grid cells. In all panels, black dots represent way-points on the transect and in C) pink diamonds represent the set of potential start locations."----

##############################################################################
####  Plot all three transect samples                                     ####
##############################################################################

#a 3x1 figure with margins that aren't too large
par( mfrow=c(3,1), mar=c(2.1,2.1,3.1,1.1))

#plot transect sample design (no constraints) on background of depth
plot( bathy, col=matlab.like2(n=100), main="Random Transect Sample", legend=FALSE)
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( transSample$points[,c("x","y")], pch=20, col='black')
mtext( "A)", adj=0.02, line=-2)

#plot transect sample design (downhill transects only) on background of depth
plot( bathy, col=matlab.like2(n=100), 
                main="Random Down-Slope Transect Sample", legend=FALSE)
#plot( fullRasters[[1]])
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( transSample_downhill$points[,c("x","y")], pch=20, col='black')
mtext( "B)", adj=0.02, line=-2)

#plot transect sample design (only transects from certain locations) on background of depth
plot( bathy, col=matlab.like2(n=100), 
      main="Random Transect Sample from Defined Starting Points", legend=FALSE)
#plot( fullRasters[[1]])
contour( bathy, levels=depthBins, add=TRUE, col=grey(0.6))
points( highLocations, pch=23, col="black", bg="pink", cex=1.5)
points( transSample_wagon$points[,c("x","y")], pch=20, col='black')
mtext( "C)", adj=0.02, line=-2)

