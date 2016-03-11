t############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 12
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 06
## Date of last changes: 26.09.2015
## Data: Megaliths, Bronze Age Barrows and Medieval villages
## Author of data: gshdl, srtm
## Purpose: didactic
## Content: 1. xxxxx
## Description: perception
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #

# 0. Preparation ===============================
#  working direktory
wd <- "~/modproj_qaam"
wd <- "/home/fon/daten/analyse/modproj_qaam"
setwd(wd)
load("4ws/ws10.rws")
#load("4ws/ws11.rws")
#load("4ws/ws07.rws")


# 1. intro - coordinates ===============================

 trans.pol <- function(a, b=c(0,0)){
     x   <- a[1]
     y   <- a[2]
     xt  <- b[1]
     yt  <- b[2]
     r <- (((x-xt)^2)+((y-yt)^2))^0.5
     if ((x-xt) >= 0 & (y-yt) >= 0)  phi <- atan((y-yt)/(x-xt))
     if ((x-xt) < 0  & (y-yt) >= 0)  phi <- atan((y-yt)/(x-xt))  + pi
     if ((x-xt) < 0  & (y-yt) < 0)   phi <- atan((y-yt)/(x-xt))  - pi
     if ((x-xt) >= 0 & (y-yt) < 0)   phi <- atan((y-yt)/(x-xt))  + 2 * pi
    return(c(r, phi))
     }
 
 trans.cartes <- function(a, b=c(0,0)){
         r   <- a[1]
         phi <- a[2]
         xt  <- b[1]
         yt  <- b[2]
         x   <- r*cos(phi) + xt
         y   <- r*sin(phi) + yt
         return(c(x, y))
         }
 
a <- c(3559376, 6027178)
b <- c(3564474, 6032765)
 
c <- trans.pol(a,b)
c
trans.cartes(c,b)


# 2.  sensual perception - visibility ===========================


# ## define variables
# file_meg <- "1data/meg.dw2.csv";file_tum <- "1data/tum.dw2.csv";file_vil <- "1data/villages.xls"
# sgdf_srtm <- "2geodata/dw_gk3_50_ag.asc"
# 
# #crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1+x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs"
# #crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# crs3 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1+x_0=3500000 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
# 
# ## load library
# library(sp);library(gdata);library(rgdal);library(proj4);library(raster);library(rasterVis)
# 
# ## load data
# df_meg <- read.table(file_meg, sep=";", header=TRUE);df_tum <- read.table(file_tum, sep=";", header=TRUE)
# spdf_meg <- read.table(file_meg, sep=";", header=TRUE);spdf_tum <- read.table(file_tum, sep=";", header=TRUE)
# coordinates(spdf_meg)= ~x+y;coordinates(spdf_tum)= ~x+y
# 
# df_vil_wgs84 <- read.xls(file_vil, 1) # requires gdata
# spdf_vil_wgs84 <- df_vil_wgs84
# coordinates(spdf_vil_wgs84)= ~x+y
# 
# 
# ### srtm
# sgdf_srtm <- readGDAL(sgdf_srtm) # requires rgdalc
# names(sgdf_srtm@data) <- "srtm" # change the name from "band1" to "srtm"
# sgdf_srtm@proj4string@projargs <- crs3
# 
# proj4string(spdf_meg) <- CRS(as.character(crs3))
# proj4string(spdf_tum) <- CRS(as.character(crs3))
# 
# ## reproject villages.xls --> requires proj4string
# df_vil_coord <- proj4::project(cbind(spdf_vil_wgs84@coords[,1], spdf_vil_wgs84@coords[,2]), crs3) # requires proj4
# df_vil_k <- cbind(x=df_vil_coord[,1]+3500000,y=df_vil_coord[,2]) # +3500000 due to GK zoning
# df_vil <- data.frame(id=df_vil_wgs84[,1], village=as.character(df_vil_wgs84[,2]), AD=df_vil_wgs84[,3])
# spdf_vil <- SpatialPointsDataFrame(df_vil_k, as.data.frame(df_vil), proj4string=CRS(as.character(crs3)))


## Visibility
## ========================================
#install.packages("spgrass7", repos="http://R-Forge.R-project.org")
library(spgrass7)
library(rgdal)

# INITIALIZE GRASS
# set up environment --  no spatial ref now
#loc <-  initGRASS("/usr/local/grass-7.0.0svn/",home=tempdir(),mapset = "PERMANENT",override = TRUE)
loc <-  initGRASS("/usr/lib/grass70",home=tempdir(), mapset = "PERMANENT", override = TRUE)
#loc <-  initGRASS("/opt/grass64", home=tempdir(), mapset = "PERMANENT", override = TRUE)

## check the location and define it according to our data
execGRASS("g.proj", flags = c("p"))
execGRASS("g.proj", flags = c("c"), parameters = list(proj4=crs1))

## load data from R to grass
writeRAST(x = sgdf_srtm, vname = "dem")

## GRASS - show raster in mapset
execGRASS("g.list", type = "rast")
#execGRASS("r.info ", parameters = list(map = "dem"))

## adjust Regions resolution!
execGRASS("g.region", flags = c("p"))
execGRASS("g.region", parameters = list(raster = "dem",res = "50"))

## PERFORM VIEWSHED ANALYSIS
## ====================
# display the possible commands in grass
parseGRASS('r.viewshed')

## viewshed for one point
co.meg <- spdf_meg@coords
execGRASS("r.viewshed", flags = c("overwrite","b"), parameters = list(input = "dem",output = "view.meg",coordinates = co.meg[160,]))
single.viewshed <- readRAST("view.meg")


## viewshed for all points
## ------------------------------
## load basic raster from grass
dem <- readRAST(vname = "dem")
## create coordinate vector and take care that all points are within the study area
co.meg <- cbind(spdf_meg@coords[,1][spdf_meg@coords[,1]>dem@bbox[1,1] & spdf_meg@coords[,1]<dem@bbox[1,2]],spdf_meg@coords[,2][spdf_meg@coords[,2]>dem@bbox[2,1] & spdf_meg@coords[,1]<dem@bbox[2,2]])
co.tum <- cbind(spdf_tum@coords[,1][spdf_tum@coords[,1]>dem@bbox[1,1] & spdf_tum@coords[,1]<dem@bbox[1,2]],spdf_tum@coords[,2][spdf_tum@coords[,2]>dem@bbox[2,1] & spdf_tum@coords[,1]<dem@bbox[2,2]])
co.vil <- cbind(spdf_vil@coords[,1][spdf_vil@coords[,1]>dem@bbox[1,1] & spdf_vil@coords[,1]<dem@bbox[1,2]],spdf_vil@coords[,2][spdf_vil@coords[,2]>dem@bbox[2,1] & spdf_vil@coords[,1]<dem@bbox[2,2]])

cum.view.meg <- raster(dem)
cum.view.tum <- raster(dem)
cum.view.vil <- raster(dem)
cum.view     <- raster(dem)

for (i in seq(1, length(co.meg[,1]))) {
    execGRASS("r.viewshed"
              ,flags = c("overwrite","b")
              ,parameters = list(input = "dem",output = "view.meg",coordinates = co.meg[i,])
    )
    viewshed <- readRAST("view.meg")
    if (i==1) cum.view.meg <- raster(viewshed)
    else  cum.view.meg <- raster(viewshed) + cum.view.meg 
    cat("iteration ", i, " of ", length(co.meg[,1]),"\n")
}

for (i in seq(1, length(co.tum[,1]))) {
    execGRASS("r.viewshed"
              ,flags = c("overwrite","b")
              ,parameters = list(input = "dem",output = "view.tum",coordinates = co.tum[i,])
    )
    viewshed <- readRAST("view.tum")
    if (i==1) cum.view.tum <- raster(viewshed)
    else  cum.view.tum <- raster(viewshed)  +  cum.view.tum 
    cat("iteration ", i, " of ", length(co.tum[,1]),"\n")
}

for (i in seq(1, length(co.vil[,1]))) {
    execGRASS("r.viewshed"
              ,flags = c("overwrite","b")
              ,parameters = list(input = "dem",output = "view.vil",coordinates = co.vil[i,])
    )
    viewshed <- readRAST("view.vil")
    if (i==1) cum.view.vil <- raster(viewshed)
    else cum.view.vil <- raster(viewshed) + cum.view.vil
    cat("iteration ", i, " of ", length(co.vil[,1]),"\n")
}

cum.view <- cum.view.meg+cum.view.tum+cum.view.vil
cum.view.b <- brick(cum.view.meg, cum.view.tum, cum.view.vil, cum.view)


pdf("6pictures/c11_view1.pdf", height=3.8, width=6, bg = "white") 
par( mai = c(0, 0, 0, 0))
plot(cum.view.b, col = grey(20:0/30), main=c("megaliths","tumuli","villages","all"), legend.shrink=1.0)
dev.off() 


# 3.  cognitive perception ===========================
# 3.1.  cognitive perception - fuzzy ===========================
memb.low <- function(x){
    y <- 12.53 * dnorm(x, mean=5, sd=5)
    return(y)}
memb.high <- function(x){
    y <- 12.53 * dnorm(x, mean=15, sd=5)
    return(y)}
seq1 <- seq(0,20,0.1)
mlow <- memb.low(seq1)
mhigh <- memb.high(seq1)
pdf("6pictures/c11_membfunctions.pdf", height=3, width=6, bg = "white") 
par (mai = c(1, 1, 0.2, 0.2))
    matplot(seq1, mlow, type="l", lty=1, ylab="membership degree", xlab="altitude (m)")
    matplot(seq1, mhigh, type="l", lty=2, ylab="membership degree", add=T)
dev.off() 

sgdf_srtm_low <- sgdf_srtm_high <- sgdf_srtm
sgdf_srtm_low@data[,1] <- memb.low(sgdf_srtm@data[,1])
sgdf_srtm_high@data[,1] <- memb.high(sgdf_srtm@data[,1])

pdf("6pictures/c11_altitudeclasses.pdf", height=2.3, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    image(sgdf_srtm_low, col = gray.colors(25, start = 0.90, end = 0.2))  
    image(sgdf_srtm_high, col = gray.colors(25, start = 0.90, end = 0.2))  
par(mfcol=c(1,1))
dev.off() 

# 3.2.  cognitive perception - cogmap ===========================
# prepare data.frame
vp <- c(3565015, 6030963)
cm_denspt <- data.frame(coordinates(cent_meg), 10000000 * over(cent_meg, sgdf_meg_dens), r=0, phi=0, r2=0, phi2=0, x2=0, y2=0)

# transform to polar coordinates
tp2 <- function(a,b) trans.pol(a,vp)
cm_pc <- apply(cm_denspt[,1:2], 1, tp2)
cm_denspt[,4] <- cm_pc[1,]
cm_denspt[,5] <- cm_pc[2,]

# manipulate distances and angles
cm_denspt[,6] <- cm_denspt[,4] / (0.9 * cm_denspt[,3]^0.2)
cm_denspt[which(cm_denspt[,5] < pi & cm_denspt[,5] > 0),7] <- cm_denspt[which(cm_denspt[,5] < pi & cm_denspt[,5] > 0),5] * 0.95
cm_denspt[which(cm_denspt[,5] > pi),7] <- 2 * pi - (2 * pi - cm_denspt[which(cm_denspt[,5] > pi),5]) * 0.95 
cm_denspt[which(cm_denspt[,5] > -pi & cm_denspt[,5] < 0),7] <- cm_denspt[which(cm_denspt[,5] > -pi & cm_denspt[,5] < 0),5] * 0.95

# transform back
tc2 <- function(a,b)  trans.cartes (a,vp)
cm_cc <- apply(cm_denspt[,6:7], 1, tc2)
cm_denspt[,8] <- cm_cc[1,]
cm_denspt[,9] <- cm_cc[2,]

# produce a spatialpointsdataframe
cm_despt2 <- SpatialPointsDataFrame(cbind(cm_cc[1,], cm_cc[2,]), cm_denspt, proj4string = CRS(as.character(crs1)))


# Lösung 1

grid_x <- gstat::idw(x2~1, cm_despt2, newdata=sgdf_srtm,  nmax=12, maxdist=20000, idp=4.0)#; image(grid_x)
grid_y <- gstat::idw(y2~1, cm_despt2, newdata=sgdf_srtm,  nmax=12, maxdist=20000, idp=4.0)#; image(grid_y)
#Werte aus den rastern extrahieren
cm_spdf <- data.frame(x=grid_x@data[,1], y=grid_y@data[,1], z=sgdf_srtm@data[,1])
# na entfernen und in spdf umwandeln
cm_spdf <- cm_spdf[intersect(intersect(which(!is.na(cm_spdf[,1])), which(!is.na(cm_spdf[,2]))), which(!is.na(cm_spdf[,3]))),]
coordinates(cm_spdf)  = ~x+y
proj4string(cm_spdf)  <- CRS(as.character(crs1)) 

library(gstat)
vt3    <- variogram(cm_spdf@data$z ~ 1, cm_spdf)
v.fit3 <- fit.variogram(vt3, vgm(1, "Mat", 5000, 1))
plot(vt3,v.fit3)
grid_cm <- krige(cm_spdf@data$z ~ 1, cm_spdf, newdata=sgdf_srtm, v.fit3, nmin = 4, maxdist = 200, nmax = 15)
grid_cm <- raster(grid_cm)
grid_cm[is.na(grid_cm)] <- 0

pdf("6pictures/c11_cm1.pdf", height=2.0, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.1))  
    for (i in 1:length(cm_despt2@data[,1])) {
        xl <- c(cm_despt2@data[i,1], cm_despt2@data[i,8])
        yl <- c(cm_despt2@data[i,2], cm_despt2@data[i,9])
        lines(xl,yl, lwd=2)
    }
    points(cm_despt2@data[,8], cm_despt2@data[,9], pch=16, col='black', cex=sqrt(1.35 * cm_despt2@data[,3]))
    points(vp[1], vp[2], pch=8, col='black')
    
    image(as(grid_cm, 'SpatialGridDataFrame'), col = gray.colors(25, start = 0.97, end = 0.1))  
    for (i in 1:length(cm_despt2@data[,1])) {
        xl <- c(cm_despt2@data[i,1], cm_despt2@data[i,8])
        yl <- c(cm_despt2@data[i,2], cm_despt2@data[i,9])
        lines(xl,yl, lwd=2)
    }
    points(cm_despt2@data[,8], cm_despt2@data[,9], pch=16, col='black', cex=sqrt(1.35 * cm_despt2@data[,3]))
    points(vp[1], vp[2], pch=8, col='black')
par(mfcol=c(1,1))
dev.off() 



# Lösung 2

#library(deldir)
#cm_del <- deldir(cm_denspt[,1], cm_denspt[,2], plot=F, wl='tr')  
#fs_nb_del <- tri2nb(fs, row.names=ids)    
#library(tripack)
#cm_del <- tri.mesh(cm_denspt[,1], cm_denspt[,2], duplicate = 'remove')   # ermittelt die delauny trinangulation




# srtm in data.frame umwandeln
srtm <- data.frame(x=coordinates(sgdf_srtm)[,1], y=coordinates(sgdf_srtm)[,2], z=sgdf_srtm@data[,1], u=0, v=0)
cmdf <- cbind(x=cm_denspt[,1], y=cm_denspt[,2], u=cm_denspt[,8], v=cm_denspt[,9])
n <- length(cmdf[,1])
p <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)
q <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)
b1 <- (sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,3] - cmdf[,1])) - q  * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))) /  (sum((cmdf[,2] - mean(cmdf[,2]))^2)  - q * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,2] - mean(cmdf[,2]))))            
b3 <- (sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,4] - cmdf[,2])) - p  * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))) /  (sum((cmdf[,1] - mean(cmdf[,1]))^2)  - p * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,1] - mean(cmdf[,1]))))    
b4 <- 1 - p * b3 + sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)   
b0 <- 1 - q * b1 + sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)   
b2 <- (1/n) * sum(cmdf[,3]-cmdf[,1]) +  mean(cmdf[,1]) - b0 * mean(cmdf[,1]) - b1 * mean(cmdf[,2])
b5 <- (1/n) * sum(cmdf[,4]-cmdf[,2]) +  mean(cmdf[,2]) - b4 * mean(cmdf[,2]) - b3 * mean(cmdf[,1])
trans.x <- function(x) u <- b0*x[1] +b1*x[2] + b2
trans.y <- function(x) u <- b3*x[1] +b4*x[2] + b5
srtm[,4] <- apply(srtm[,1:2], 1, trans.x)
srtm[,5] <- apply(srtm[,1:2], 1, trans.y)
coordinates(srtm)  = ~u+v
proj4string(srtm)  <- CRS(as.character(crs1)) 
cm_raster <- rasterize(srtm, raster(sgdf_srtm), field='z', update=TRUE, proj4string = CRS(as.character(crs1)))

# pdf("6pictures/c11_cm2.pdf", height=2.0, width=6, bg = "white") 
# par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
#     image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.1))  
#     for (i in 1:length(cm_despt2@data[,1])) {
#         xl <- c(cm_despt2@data[i,1], cm_despt2@data[i,8])
#         yl <- c(cm_despt2@data[i,2], cm_despt2@data[i,9])
#         lines(xl,yl, lwd=2)
#     }
#     points(cm_despt2@data[,8], cm_despt2@data[,9], pch=16, col='black', cex=sqrt(1.35 * cm_despt2@data[,3]))
#     points(vp[1], vp[2], pch=8, col='black')
#     
#     image(as(cm_raster, 'SpatialGridDataFrame'), col = gray.colors(25, start = 0.97, end = 0.1))  
#     for (i in 1:length(cm_despt2@data[,1])) {
#         xl <- c(cm_despt2@data[i,1], cm_despt2@data[i,8])
#         yl <- c(cm_despt2@data[i,2], cm_despt2@data[i,9])
#         lines(xl,yl, lwd=2)
#     }
#     points(cm_despt2@data[,8], cm_despt2@data[,9], pch=16, col='black', cex=sqrt(1.35 * cm_despt2@data[,3]))
#     points(vp[1], vp[2], pch=8, col='black')
# par(mfcol=c(1,1))
# dev.off() 
# 



# Lösung 3
sgdf <- sgdf_srtm 
sgdf@data[,1]  <- 0
srtm <- data.frame(x=coordinates(sgdf_srtm)[,1], y=coordinates(sgdf_srtm)[,2], z=sgdf_srtm@data[,1], u=0, v=0)
library(spatstat)
cm_pp <- ppp(cm_denspt[,1], cm_denspt[,2], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del <- delaunay(cm_pp)
cm_pp2 <- ppp(cm_denspt[,8], cm_denspt[,9], window=owin(c(bb[1,1],bb[1,2]),c(bb[2,1],bb[2,2])))
cm_del2 <- delaunay(cm_pp2)
trans.x <- function(x) u <- b0*x[1] +b1*x[2] + b2
trans.y <- function(x) u <- b3*x[1] +b4*x[2] + b5
for (i in 1:cm_del$n)   {
    x1 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][1]
    y1 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][1]
    x2 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][2]
    y2 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][2]
    x3 <- tiles(cm_del)[[i]]$bdry[[1]][[1]][3]
    y3 <- tiles(cm_del)[[i]]$bdry[[1]][[2]][3]
    u1 <- cm_despt2@data[which(cm_despt2@data[,1] == x1), 8][1]
    v1 <- cm_despt2@data[which(cm_despt2@data[,2] == y1), 9][1]
    u2 <- cm_despt2@data[which(cm_despt2@data[,1] == x2), 8][1]
    v2 <- cm_despt2@data[which(cm_despt2@data[,2] == y2), 9][1]
    u3 <- cm_despt2@data[which(cm_despt2@data[,1] == x3), 8][1]
    v3 <- cm_despt2@data[which(cm_despt2@data[,2] == y3), 9][1]
    cmdf <- cbind(x=c(x1,x2,x3), y=c(y1,y2,y3), u=c(u1,u2,u3), v=c(v1,v2,v3))
    n <- length(cmdf[,1])
    p <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)
    q <- sum((cmdf[,1] - mean(cmdf[,1]))  * (cmdf[,2] - mean(cmdf[,2])))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)
    b1 <- (sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,3] - cmdf[,1])) - q  * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))) /  (sum((cmdf[,2] - mean(cmdf[,2]))^2)  - q * sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,2] - mean(cmdf[,2]))))            
    b3 <- (sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,4] - cmdf[,2])) - p  * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))) /  (sum((cmdf[,1] - mean(cmdf[,1]))^2)  - p * sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,1] - mean(cmdf[,1]))))    
    b4 <- 1 - p * b3 + sum((cmdf[,2] - mean(cmdf[,2])) * (cmdf[,4] - cmdf[,2]))  /  sum((cmdf[,2] - mean(cmdf[,2]))^2)   
    b0 <- 1 - q * b1 + sum((cmdf[,1] - mean(cmdf[,1])) * (cmdf[,3] - cmdf[,1]))  /  sum((cmdf[,1] - mean(cmdf[,1]))^2)   
    b2 <- (1/n) * sum(cmdf[,3]-cmdf[,1]) +  mean(cmdf[,1]) - b0 * mean(cmdf[,1]) - b1 * mean(cmdf[,2])
    b5 <- (1/n) * sum(cmdf[,4]-cmdf[,2]) +  mean(cmdf[,2]) - b4 * mean(cmdf[,2]) - b3 * mean(cmdf[,1])
    cm_ind <- which(inside.owin(srtm[,1], srtm[,2], tiles(cm_del)[[i]]) == T)
    srtm[cm_ind,4] <- apply(srtm[cm_ind,1:2], 1, trans.x)
    srtm[cm_ind,5] <- apply(srtm[cm_ind,1:2], 1, trans.y)
    }
coordinates(srtm)  = ~u+v
proj4string(srtm)  <- CRS(as.character(crs1)) 
cm_raster2 <- rasterize(srtm, raster(sgdf), field='z', update=TRUE, proj4string = CRS(as.character(crs1)))
cm_raster2[which(getValues(cm_raster2) == 0)] <- NA
cm_raster2 <- focal(cm_raster2, w= matrix(rep(1,25), nrow=5),  na.rm=TRUE, mean, NAonly=TRUE) 
cm_raster2[is.na(cm_raster2)] <- 0

pdf("6pictures/c11_cm3.pdf", height=2.0, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    image(as(cm_raster, 'SpatialGridDataFrame'), col = gray.colors(25, start = 0.97, end = 0.1))  
    for (i in 1:length(cm_despt2@data[,1])) {
        xl <- c(cm_despt2@data[i,1], cm_despt2@data[i,8])
        yl <- c(cm_despt2@data[i,2], cm_despt2@data[i,9])
        lines(xl,yl, lwd=2)
    }
    points(cm_despt2@data[,8], cm_despt2@data[,9], pch=16, col='black', cex=sqrt(1.35 * cm_despt2@data[,3]))
    points(vp[1], vp[2], pch=8, col='black')
    
    image(as(cm_raster2, 'SpatialGridDataFrame'), col = gray.colors(25, start = 0.97, end = 0.1))  
    for (i in 1:length(cm_despt2@data[,1])) {
        xl <- c(cm_despt2@data[i,1], cm_despt2@data[i,8])
        yl <- c(cm_despt2@data[i,2], cm_despt2@data[i,9])
        lines(xl,yl, lwd=2)
    }
    points(cm_despt2@data[,8], cm_despt2@data[,9], pch=16, col='black', cex=sqrt(1.35 * cm_despt2@data[,3]))
    points(vp[1], vp[2], pch=8, col='black')
    plot(cm_del2, add=T, lty=2)
par(mfcol=c(1,1))
dev.off() 


save.image("4ws/ws11.rws")





##########################################################




> a <- matrix(c(1,2,3,4),2,2)
> a
[,1] [,2]
[1,]    1    3
[2,]    2    4
> b <- c(b1,b2)
Fehler: Objekt 'b1' nicht gefunden
> ?solve
> 
    > b <- c(5,6)
> solve(a, b)
[1] -1  2
> x <- solve(a, b)
> a %*% x 
[,1]
[1,]    5
[2,]    6



