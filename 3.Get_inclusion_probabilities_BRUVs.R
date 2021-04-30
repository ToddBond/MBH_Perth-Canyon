### Create inclusion probabilities ####

library( rgdal)
library( sp)
library( raster)


# clear environment ----
rm(list = ls())


# Set names ----
study <- "Griffen_MBH"

platform <- "Bruvs"

design.version <- "v1"


# Directories ----
w.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
#w.dir <- "~/MBH_AbroNPZs"
o.dir <- paste(w.dir, "outputs", sep = '/')
d.dir <- paste(w.dir, "data", sep='/')
s.dir <- paste(w.dir, "shapefiles", sep='/')
r.dir <- paste(w.dir, "rasters", sep='/')



# Read in data ----

dat <- readRDS(paste(d.dir,"Data-Griffen_MBH.RDS", sep ='/'))
rast <- readRDS(paste(d.dir, "rasters-Griffen_MBH-fine.RDS", sep='/'))
zones <- readRDS(paste(d.dir, "Zones-Griffen_MBH.RDS", sep='/'))


#  Straw man for numbers of samples in each region ----
# 72 BRUVs in total

straw.nums <- c(3,3,3,3,3,3,4,3,3,10,14,6,14)  # numbers of drops rest w structure + caut with structure,cau wout str, open w str, open wout str
straw.props <- straw.nums / sum( straw.nums)
straw.props
names( straw.nums) <- names( straw.props) <- c( "g1","g2","g3","g4","g5","g6","g7","g8","g9","cs","cz","os","oz")
saveRDS(straw.nums, file = paste0(paste(d.dir, paste("StrawmanNumbers" , study, platform, design.version, sep='-'), sep='/'), ".RDS"))


# Get bathy cut points ----
# And their numbers of drops
Bathy.quant <- c(0,0.20, 0.6, 1)
Bathy.cuts <- quantile(rast$depth, Bathy.quant)#c( -Inf,0.02,0.04,0.08,0.16,Inf)
Bathy.cuts # -200 -149 -120  -86 
#trying to make it so there is no hand-picking (except for the hand-picked function)
tmp <- cumsum( Bathy.quant)
#Bathy.targetNums <- rep( floor( 18/8), 4)#floor( ( tmp / sum( tmp))[-1] * 200)#rep( 40, 5)#c( 20,20,30,65,65)
Bathy.targetNums <- rep(floor( ( tmp / sum( tmp))[-1] * sum(straw.nums)))
Bathy.targetProps <-  Bathy.targetNums / sum( Bathy.targetNums)
# might need to alter these proportions by hand, so that the middle interval is better represented, as all structures fall within that interval


# Depth limiting and some elementary cleaning ----

#minDepthInAMP <- max( NingalooDat$BATHY, na.rm=TRUE)

#NingalooDat[ is.na( NingalooDat$BATHY) | NingalooDat$BATHY < -195 & NingalooDat$BATHY > minDepthInAMP, c("BATHY","TPI_GF")] <- NA

GBDat[ !is.na( GBDat$BATHY) & GBDat$BATHY < -50 , "BATHY"] <- NA
GBDat[ !is.na( GBDat$BATHY) & GBDat$BATHY > 0 , "BATHY"] <- NA


pdf( "Bathy_distribution_zones.pdf")
{
  
  par( mfrow=c(2,2))
  for( ii in c("MUZ", "SPZ", "HPZ", "NPZ"))
    hist( GBDat[GBDat[,ii]==1,"BATHY"], nclass=50, main=ii)
  
  par( mfrow=c(1,1))
  kount <- 1
  plot( ecdf( GBDat[GBDat[,"MUZ"]==1,"BATHY"]), col=2, main="Bathy in Zones", xlim=c(-0.4,0.4))
  for( ii in c("MUZ", "SPZ", "HPZ", "NPZ")){
    plot( ecdf( GBDat[GBDat[,ii]==1,"BATHY"]), col=kount, add=TRUE)
    kount <- kount + 1
  }
  legend( "bottomright", legend=c("MUZ", "SPZ", "HPZ", "NPZ"), lty=1, 
          col=c(1:4), lwd=c(1,1,1,1), bty='n')
  #so surprisingly, the southern control zone has more flat area than the IUCN2 area
  #there may be some issues with the north control having quite a few more 'bumps' than the southern control and IUCN2 area.  Bias.
  dev.off()
}
#########################
#proportion of potential sites in each zone

GBDat_small <- GBDat[!is.na( GBDat$BATHY),]
tmp <- colSums( GBDat_small[,c("MUZ", "SPZ", "HPZ", "NPZ")], na.rm=TRUE) 
tmp[2] <- tmp[2] - tmp[1] # so similar amount of sites in SPZ and MUZ
props <- tmp / nrow( GBDat_small)
props <- props / sum( props) # 1 UP TO HERE

###################################
####  TPI to get cut points
###################################

catB <- cut( gb_rasters$bathy, breaks=Bathy.cuts, na.rm=TRUE)

plot( zones$MUZ); plot( catB, add=TRUE); plot( zones$MUZ, add=TRUE)

writeRaster( catB, file='Bathy_cuts.tif', overwrite=TRUE)

plot(catB)

##################################
####  Within each zone (incl probs)
####  Weight according to straw.props
##################################



inclProbs <- catB
for( zz in c( "NPZ", "HPZ", "SPZ", "MUZ")){
  print( zz)
  #if( zz == "MUZ")
  #zoneID <- extract( x=catB, y=zones$MUZ, cellnumbers=TRUE)
  #zoneID <- extract( x=catB, y=zones$MUZ-zones$NPZ, cellnumbers=TRUE)
  #else
  zoneID <- extract( x=catB, y=zones[[zz]], cellnumbers=TRUE)
  propsOfTPI <- table( catB@data@values[zoneID[[1]][,"cell"]])
  propsOfTPI <- propsOfTPI / sum( propsOfTPI)
  tmp <- Bathy.targetProps / propsOfTPI #the desired inclusion probs (unstandardised)
  for( ii in 1:length( propsOfTPI)){
    inclProbs[zoneID[[1]][,"cell"]][zoneID[[1]][,"value"]==ii] <- tmp[ii]
  }
  inclProbs[zoneID[[1]][,"cell"]][is.na( inclProbs[zoneID[[1]][,"cell"]])] <- 0
  inclProbs[zoneID[[1]][,"cell"]] <- inclProbs[zoneID[[1]][,"cell"]] / sum( inclProbs[zoneID[[1]][,"cell"]])
}
inclProbs@data@values[inclProbs@data@values %in% c(0,1,2,3,4,5,6,7,8)] <- NA  #cheats way to crop
plot( inclProbs)

#standardising so that the zone totals are correct according to straw.props | straw.nums
MUZzone <- extract( x=catB, y=zones$MUZ, cellnumbers=TRUE)
SPZzone <- extract( x=catB, y=zones$SPZ, cellnumbers=TRUE)
HPZZone <- extract( x=catB, y=zones$HPZ, cellnumbers=TRUE)
NPZZone <- extract( x=catB, y=zones$NPZ, cellnumbers=TRUE)

inclProbs@data@values[MUZzone[[1]][,'cell']] <- inclProbs@data@values[MUZzone[[1]][,'cell']] * straw.props["MUZ"]
inclProbs@data@values[SPZzone[[1]][,'cell']] <- inclProbs@data@values[SPZzone[[1]][,'cell']] * straw.props["SPZ"]
inclProbs@data@values[HPZZone[[1]][,'cell']] <- inclProbs@data@values[HPZZone[[1]][,'cell']] * straw.props["HPZ"]
inclProbs@data@values[NPZZone[[1]][,'cell']] <- inclProbs@data@values[NPZZone[[1]][,'cell']] * straw.props["NPZ"]

plot(inclProbs)

writeRaster( inclProbs, file='inclProbs_design1.tif', overwrite=TRUE)

rm( list=lso()$OTHERS)
rm( AMPzone, catTPI, ii, IUCNzone, kount, minDepthInAMP, NingalooDat_small, NthZone, props, propsOfTPI, SthZone, 
    tmp, zoneID, zz, TPI.cuts, TPI.quant, TPI.targetNums, TPI.targetProps)

