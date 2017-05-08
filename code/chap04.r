############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 4
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 06
## Date of last changes: 04.08.2015
## Data: Megaliths, Bronze Age Barrows and Medieval villages
## Author of data: gshdl
## Purpose: didactic
## Content: 1. one dimension, 2. two dimensions
## Description: applies some basic techniques of
##              density calculation
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #
#
# pdf export einbauen
#
# 0. Preparation ===============================
#  load data
df_vil <- read.table(
    "data/villages.csv",
    header = TRUE, 
    sep = ",",
    stringsAsFactors = TRUE
)

meg <- read.table(
    "data/meg_dw.csv",
    header = TRUE, 
    sep = ";",
    stringsAsFactors = TRUE
)

library(sp)
file_srtm   <- "data/dw_gk3_50_ag.asc"
sgdf_srtm <- read.asciigrid(file_srtm) 

library(raster)
crs1 <- "+init=epsg:31467"
proj4string(sgdf_srtm) <- CRS(as.character(crs1)) 

# 1. one dimension ===============================
vil_fd <- df_vil[,3]

cb <- c(1200,1250,1300,1350,1400)
count <- 1:4
for (i in 1:4) {higher <- which(df_vil[,3] > cb[i]) 
                lower  <- which(df_vil[,3] < cb[i+1]) 
                hl <- intersect(higher,lower) 
                length(hl) 
                count[i] <- length(hl)}
years <- c("1200-1250","1250-1300","1300-1350", "1350-1400")
data.frame(years,count)

pdf("pictures/c4_emVil1.pdf", height=4, width=6, bg = "white") 
    hist(x = df_vil[,3], breaks = 6, col = "gray", border = "white", xlab = "Time A.D.", main = "Histogram of village foundations in different periods")
dev.off() 

library("KernSmooth")
ks_vil <- bkde(df_vil[,3], kernel="normal",  bandwidth=5, gridsize=201, range.x =  c(1200,1400))
plot(ks_vil, pch=20, col = "gray", xlim = c(1240,1360), xlab = "Time A.D.", ylab = "density of village foundation", panel.first = grid())
lines(ks_vil)

plot(df_vil[,3],df_vil[,1], col="gray", pch=16, xlab = "Time A.D.", ylab = "id", panel.first = grid())
lines(df_vil[,3],df_vil[,1])

interval <- c(df_vil[,3],df_vil[13,3]) - c(df_vil[1,3],df_vil[,3])
plot(interval[2:12],col="gray", pch=16, xlab = "Index", ylab = "Interval", panel.first = grid())
lines(interval[2:12])

ts_vil <- ts(ks_vil$y, start=c(1200), end=c(1400),frequency = 1) 
plot(ts_vil, ylab = "density", panel.first = grid())

acf(ts_vil, lag.max = 100, main = "")


# 2. two dimensions ===============================
library(spatstat)

ppp_meg <- ppp(
    x = meg$x, y = meg$y,
    window = owin(
        xrange = range(meg[,2]),
        yrange = range(meg[,3])
    )
)

count <- ppp_meg$n
dx <- (ppp_meg$window$xrange[2] -  ppp_meg$window$xrange[1]) / 1000
dy <- (ppp_meg$window$yrange[2] -  ppp_meg$window$yrange[1]) / 1000
density_1 <- count / (dx*dy)
density_1

library(spatstat)  
ch <- convexhull.xy(ppp_meg$x, ppp_meg$y) 
fl <- area.owin(ch)   

fl <- area.owin(ch)      
density_2 <- count / fl  
(density_2 * 1000000) / density_1 

rw      <- 3000    
xmin    <- ppp_meg$window$xrange[1] - rw/2
xmax    <- ppp_meg$window$xrange[2] + rw/2
ymin    <- ppp_meg$window$yrange[1] - rw/2
ymax    <- ppp_meg$window$yrange[2] + rw/2
rows    <- round((ymax-ymin)/rw, 0) + 1
columns <- round((xmax-xmin)/rw, 0) + 1
z  <- cbind(1:(columns*rows))
df <- data.frame(z)

gt  <- GridTopology(cellcentre.offset=c(xmin,ymin), cellsize=c(rw,rw), cells.dim=c(columns,rows))
sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(crs1))

for (i in seq(along=coordinates(gt)[,1])){
    x <- coordinates(gt)[i,1] - rw/2
    y <- coordinates(gt)[i,2] - rw/2
    xi <- which(ppp_meg$x>x & ppp_meg$x<x+rw)  
    yi <- which(ppp_meg$y>y & ppp_meg$y<y+rw)  
    pz <- length(intersect(xi,yi))   
    sgdf@data$z[i]<- pz / (rw/1000)^2 
}

library(raster)
plot(raster(sgdf), col = gray.colors(25, start = 0.97, end = 0.4), cex.axis = .9)
points(ppp_meg$x, ppp_meg$y, pch=20)

sgdf_kde <- sgdf
sd       <- 3000
for (i in seq(along=coordinates(gt)[,1])){
    x <- coordinates(gt)[i,1]
    y <- coordinates(gt)[i,2]
    g2 <- 0
    for (j in seq(along=ppp_meg$x)){  
        distance <- sqrt((ppp_meg$x[j] - x)^2 +  (ppp_meg$y[j] - y)^2)
        g1 <- dnorm(distance, mean=0, sd=sd)
        g2 <-g2 + g1}
    sgdf_kde@data$z[i]<- g2}
plot(raster(sgdf_kde), col = gray.colors(25, start = 0.97, end = 0.4), cex.axis = .9)
points(ppp_meg$x, ppp_meg$y, pch=20)

library(spatstat)
rw <- 1000  
sd <- 2000
dens_p  <- density(ppp_meg, sd, edge=TRUE, at="points")   
dens_r5 <- density(ppp_meg, sd, eps=rw, edge=TRUE, at="pixels") 
plot(raster(dens_r5), col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r5, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=20)

sdev    <- 3*mean(nndist(ppp_meg))  
dens_r6 <- density(ppp_meg, sdev, eps=rw, edge=TRUE, at="pixels")
plot(raster(dens_r6), col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r6, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=20)

dens_r <- density(ppp_meg, bw = "nrd", eps=rw, edge=TRUE, at="pixels")
plot(raster(dens_r), col = gray.colors(25, start = 0.97, end = 0.4))
contour(dens_r, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=20)



spdf_meg <- SpatialPointsDataFrame(
    coords = meg[,2:3],
    data = meg,
    proj4string = CRS(crs1)
)

rw <- 1000
fs   <- cbind(x=spdf_meg@coords[,1], y=spdf_meg@coords[,2])
rows <- round((bbox(spdf_meg)[2,2]-bbox(spdf_meg)[2,1])/rw, 0) + 2
cols <- round((bbox(spdf_meg)[1,2]- bbox(spdf_meg)[1,1])/rw, 0) + 2
z    <- cbind(1:(rows*cols))
df   <- data.frame(cbind(1:((round((bbox(spdf_meg)[2,2]-bbox(spdf_meg)[2,1])/rw, 0) + 2)*(round((bbox(spdf_meg)[1,2]- bbox(spdf_meg)[1,1])/rw, 0) + 2))))
gt   <- GridTopology(cellcentre.offset=c(bbox(spdf_meg)[1,1] - rw/2,bbox(spdf_meg)[2,1] - rw/2), cellsize=c(rw,rw), cells.dim=c(cols,rows))
sgdf <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))

library(tripack)
fsv   <- voronoi.mosaic(spdf_meg$x, spdf_meg$y, duplicate = 'remove') 
rad   <- fsv$radius
fsvsp <- SpatialPointsDataFrame(cbind(fsv$x, fsv$y), as.data.frame(rad), proj4string= CRS(as.character(crs1)))  
fspv  <- ppp(fsvsp@coords[,1], fsvsp@coords[,2], window = owin(
    xrange = range(fsvsp@coords[,1]),
    yrange = range(fsvsp@coords[,2])
)) 
fs_vd <- cbind(fspv$x,fspv$y,nncross(fspv,ppp_meg)$dist)            
fs_vd_spdf <- SpatialPointsDataFrame(cbind(fs_vd[,1],fs_vd[,2]), as.data.frame(fs_vd[,3]), proj4string=CRS(as.character(crs1)))

library(gstat)
g  <- gstat(formula=fsvsp@data$rad ~ 1, data=fsvsp,  nmin = 5, maxdist = 10000, nmax = 15)
vt <- variogram(g)
v.fit <- fit.variogram(vt, vgm(1, "Gau", 10000, 1), fit.sills = TRUE, fit.ranges = TRUE,fit.method = 1)
g <- gstat(g, id="var1", model=v.fit )
k <- predict(g, model=v.fit, newdata=sgdf_srtm)             

image(raster(k), col = gray.colors(25, start = 0.4, end = 0.97))
contour(k, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

diff6_5   <- dens_r5
diff6_5$v <- dens_r6$v - dens_r5$v

save.image("ws/ws04.rws")




