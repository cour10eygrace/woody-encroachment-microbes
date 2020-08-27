library(rgdal)
library(raster)
library(sp)

mat<-raster("wc2.1_30s_bio/wc2.1_30s_bio_1.tif")
map<-raster("wc2.1_30s_bio/wc2.1_30s_bio_12.tif")

coords <- data.frame(x=site.table$long,y=site.table$lat)

points <- SpatialPoints(coords, proj4string = mat@crs)
points2 <- SpatialPoints(coords, proj4string = map@crs)

values <- data.frame(raster::extract(mat, points))
values2 <- data.frame(raster::extract(map, points))

site.clim <- cbind.data.frame(coordinates(points),values)
site.clim<-cbind.data.frame(site.clim, values2)

site.clim<-rename(site.clim, long=x, lat=y, 
MAP=raster..extract.map..points., MAT=raster..extract.mat..points.)



