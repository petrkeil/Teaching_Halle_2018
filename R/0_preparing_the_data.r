library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(dismo)


# define projections
LAMBERT <- CRS("+proj=cea +lon_0=Central Meridian +lat_ts=Standard Parallel +x_0=False Easting +y_0=False Northing")
LAM_EUR <- CRS("+proj=laea +lat_0=55 +lon_0=40 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs ")
WGS84 <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# load the raw data
beaver.shape <- readOGR(dsn = "../data/beaver_shapefile", layer = "species_4007")
beaver.shape@data <- dplyr::select(beaver.shape@data, binomial)
world <- readOGR(dsn = "../data/global_boundaries", layer="TM_WORLD_BORDERS-0.3")

beaver.point <- gbif("Castor", "fiber", geo=T) ### BEWARE: 1min slow
beaver.point <- beaver.point[,c("lon","lat")]
beaver.point <- beaver.point[beaver.point$lon > -20,]
beaver.point <- beaver.point[beaver.point$lon < 60,]
beaver.point <- na.omit(beaver.point)
beaver.point.sp <- SpatialPoints(beaver.point)
crs(beaver.point.sp) <- WGS84


# plot the raw data
plot(beaver.shape, col="green")
plot(world, add=T)
points(beaver.point.sp, add=T, cex=0.1, col="red")

# transform the data to an equal-area Lambert projection
beaver.point.eq <- spTransform(beaver.point.sp, CRSobj=LAM_EUR)
coords.eq <- coordinates(beaver.point.eq)
colnames(coords.eq) <- c("X", "Y")

world.eq <- spTransform(world, CRSobj=LAM_EUR)
beaver.shape.eq <- spTransform(beaver.shape, CRSobj=LAM_EUR)

# plot the transformed data
plot(beaver.shape.eq, col="green")
plot(world.eq, add=T)
points(beaver.point.eq, cex=0.1, col="red")

# export the data
writeOGR(obj=beaver.shape.eq, dsn="../data/beaver_shapefile", layer="beaver_shape_equal_area", driver="ESRI Shapefile")
writeOGR(obj=beaver.shape, dsn="../data/beaver_shapefile", layer="beaver_shape_wgs84", driver="ESRI Shapefile")

writeOGR(obj=world.eq, dsn="../data/global_boundaries", layer="world_equal_area", driver="ESRI Shapefile")


write.csv(coords.eq, file="../data/beaver_points/beaver_points_equal_area.csv", row.names = FALSE)
write.csv(beaver.point, file="../data/beaver_points/beaver_points_wgs84.csv", row.names = FALSE)






