# ------------------------------------------------------------------------------
# 1. LOADING LIBRARIES

# loading some useful libraries
library(dismo)
library(XML)
library(maptools)
library(sp)
library(raster)
library(foreign)
library(rgdal)

# for the proj4 strings see www.spatialreference.org
us.atlas.proj <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs")
wgs1984.proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")


# ------------------------------------------------------------------------------
# 2. USA BOUNDARY SHAPEFILE
us <- readShapePoly("US_boundary/US_boundary2.shp")

par(mfcol=c(1,2))
plot(us)
# projecting the us shapefile
proj4string(us) <- wgs1984.proj
# reprojecting to US National Atlas equal-area projection
us.equal <- spTransform(us, us.atlas.proj)
plot(us.equal)

# ------------------------------------------------------------------------------
# 3. EXPERT RANGE-MAP DATA - Lynx canadensis
lynx.shape <- readShapePoly("LYNX/lynx.shp")
# plot shapefile:
par(mfcol=c(1,2))
plot(us)
plot(lynx.shape, add=T, col="pink")
# projecting the range map
proj4string(lynx.shape) <- wgs1984.proj
# reprojecting to US National Atlas equal-area projection
lynx.shape.equal <- spTransform(lynx.shape, us.atlas.proj)
plot(us.equal)
plot(lynx.shape.equal, add=T, col="pink")


# ------------------------------------------------------------------------------
# 4. POINT DATA - Picoides dorsalis


# get GBIF data with function gbif():

picoides.gbif <- gbif("Castor", "fiber", geo=T)

picoides.gbif <- gbif("Picoides", "dorsalis", geo=T)
# plot occurrences:
par(mfcol=c(1,2))
plot(us)
points(picoides.gbif$lon, picoides.gbif$lat,
       pch=19, cex=0.3, col="blue")
# projecting
# important: long goes first and hence 8:7
picoides <- SpatialPoints(coords=picoides.gbif[,8:7],
                                 proj4string=wgs1984.proj)
# and reprojecting
picoides.equal <- spTransform(picoides, us.atlas.proj)
plot(us.equal)
points(picoides.equal, pch=19, cex=0.3, col="blue")

# ------------------------------------------------------------------------------
# 5. CREATING EMPTY GRID
r <- raster(us.equal)
res(r) <- 50000
r[] <- rnorm(ncell(r))
plot(r)
plot(us.equal, add=T)
r[] <- 0

# ------------------------------------------------------------------------------
# 6. RASTERIZING

plot(us.equal)

# rasterizing the point data on Picoides dorsalis
picoides.raster <- rasterize(picoides.equal, r)
picoides.raster[picoides.raster>=1] <- 1
picoides.raster[picoides.raster==0] <- NA
# limiting the data only to US
picoides.raster <- picoides.raster*us.raster
plot(picoides.raster, add=T)
plot(us.equal, add=T)

# rasterizing the shapefile on Lynx arcticus
lynx.raster <- rasterize(lynx.shape.equal, r, getCover=TRUE)
lynx.raster[lynx.raster>=1] <- 1
lynx.raster[lynx.raster==0] <- NA
plot(lynx.raster, add=T)

# rasterizing US boundary
us.raster <- rasterize(us.equal, r, getCover=TRUE)
# conditioning on the amount of mainland in a grid cell
us.raster[us.raster>=1] <- 1
us.raster[us.raster==0] <- NA
plot(us.raster)

# rasterizing environmental data
temp <- raster("ENVI/bio_1")
temp.proj <- projectRaster(temp, crs=wgs1984.proj)
temp.equal <- projectRaster(from=temp.proj, to=r)
plot(temp.equal)
plot(us.equal, add=T)

# ------------------------------------------------------------------------------
# 7. AGGREGATE - picoides
par(mfrow=c(2,2), mai=c(0.1, 0.1, 0.5, 0.1))

# 50 x 50 resolution
plot(picoides.raster, axes=FALSE,
     legend=FALSE, main="50 km x 50 km")
plot(us.equal, add=T)
#points(picoides.equal, col="red", cex=0.4, pch=19)

# 100 x 100 resolution
pic.coarser <- aggregate(picoides.raster, fact=2, fun=max)
plot(pic.coarser, axes=FALSE,
     legend=FALSE, main="100 km x 100 km")
plot(us.equal, add=T)

# 200 x 200 resolution
pic.coarser2 <- aggregate(pic.coarser, fact=2, fun=max)
plot(pic.coarser2, axes=FALSE,
     legend=FALSE, main="200 km x 200 km")
plot(us.equal, add=T)

# 50 x 50 resolution
pic.coarser4 <- aggregate(pic.coarser2, fact=2, fun=max)
plot(pic.coarser4, , axes=FALSE,
     legend=FALSE, main="400 km x 400 km")
plot(us.equal, add=T)


# ------------------------------------------------------------------------------
# 8. AGGREGATE - lynx
par(mfrow=c(2,2), mai=c(0.1, 0.1, 0.5, 0.1))

# 50 x 50 resolution
lynx.raster <- lynx.raster*us.raster
plot(lynx.raster, axes=FALSE,
     legend=FALSE, main="50 km x 50 km")
plot(us.equal, add=T)
#plot(lynx.shape.equal, add=TRUE)

# 100 x 100 resolution
lynx.coarser <- aggregate(lynx.raster, fact=2, fun=max)
plot(lynx.coarser, axes=FALSE,
     legend=FALSE, main="100 km x 100 km")
plot(us.equal, add=T)

# 200 x 200 resolution
lynx.coarser2 <- aggregate(lynx.coarser, fact=2, fun=max)
plot(lynx.coarser2, axes=FALSE,
     legend=FALSE, main="200 km x 200 km")
plot(us.equal, add=T)

# 400 x 400 resolution
lynx.coarser4 <- aggregate(lynx.coarser2, fact=2, fun=max)
plot(lynx.coarser4, , axes=FALSE,
     legend=FALSE, main="400 km x 400 km")
plot(us.equal, add=T)


# ------------------------------------------------------------------------------
# 9. AGGREGATE - temperature
par(mfrow=c(2,2), mai=c(0.1, 0.1, 0.5, 0.1))

temp.equal <- temp.equal*us.raster
plot(temp.equal, axes=FALSE,
     main="50 km x 50 km")
plot(us.equal, add=T)

temp.coarser <- aggregate(temp.equal, fact=2, fun=mean)
plot(temp.coarser, axes=FALSE,
     main="100 km x 100 km")
plot(us.equal, add=T)

temp.coarser2 <- aggregate(temp.coarser, fact=2, fun=mean)
plot(temp.coarser2, axes=FALSE,
     main="200 km x 200 km")
plot(us.equal, add=T)

temp.coarser4 <- aggregate(temp.coarser2, fact=2, fun=mean)
plot(temp.coarser4,
     axes=FALSE, main="400 km x 400 km")
plot(us.equal, add=T)


# ------------------------------------------------------------------------------
# 10. EXTRACTING CARNIVORANS SPECIES RICHNESS

# reading the complete IUCN mammal dataset
# which is downloadable at:
# http://www.iucnredlist.org/technical-documents/spatial-data#mammals
mammterr <- readShapePoly("MAMMTERR/MAMMTERR.shp")

# reading the list of carnivoran species
carniv.list <- read.table("carnivorans.txt",
               sep="\t", header=T)
carniv.list <- paste(carniv.list[,1], carniv.list[,2], sep=" ")

# subsetting the mammal data only to carnivorans
carniv <- mammterr[mammterr$BINOMIAL %in% carniv.list,]
# projecting and reprojecting the dataset
proj4string(carniv) <- wgs1984.proj
carniv.equal <- spTransform(carniv, us.atlas.proj)

# calculation of species richness
# THIS IS A BETA VERSION - IS PROBABLY WRONG!
carniv.raster <- rasterize(carniv.equal, r,
                           fun=function(x)length(unique(x)))
# eliminating cells that do not overlap US
carniv.raster <- carniv.raster*us.raster

# creating empty raster at coarse scale (200 x 200 km)
r4 <- aggregate(r, fact=4, fun=mean)
# calculation of species richness at the coarse scale
carniv.raster.4 <- rasterize(carniv.equal, r4,
                             fun=function(x)length(unique(x)))
# eliminating cells that do not overlap US
carniv.raster.4 <- m2.raster.4*us.raster.4

# plotting richness
par(mfrow=c(1,2), mai=c(0.1, 0.1, 0.5, 0.1))
# 50 x 50 km
plot(carniv.raster, axes=FALSE,
     main="50 km x 50 km",
     col=heat.colors(20))
plot(us.equal, add=T)
# 200 x 200 km
plot(carniv.raster.4, axes=FALSE,
     main="200 km x 200 km",
     col=heat.colors(20))
plot(us.equal, add=T)
