############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 7
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 06
## Date of last changes: 09.08.2015
## Data: Megaliths, Bronze Age Barrows and Medieval villages
## Author of data: gshdl
## Purpose: didactic
## Content: 2. 1th order prop., 3. 2nd order prop., 
##          4. 3rd order prob.
## Description: applies some basic techniques of
##             point pattern analysis
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #

# 0. Preparation ===============================
load("ws/ws05.rws")
load("ws/ws06.rws")

library(spatstat)
library(maptools)
ppp_tum <- as.ppp(spdf_tum)
ppp_vil <- as.ppp(spdf_vil)

# 1. Point Processes ===============================
bb   <- bbox(sgdf_srtm) 
ppspec <- list(cif="strauss",par=list(beta=2,gamma=0.2,r=0.7), w=c(bb[1,1],bb[1,2],bb[2,1],bb[2,2]))
ppsim <- rmh(model=ppspec,start=list(n.start=200), control=list(nrep=10,nverb=5))
pdf("pictures/c7_ppsim.pdf", height=4, width=6, bg = "white") 
    plot(ppsim)
dev.off() 

# 2. first order properties ===============================
sdev <- 2*mean(nndist(ppp_meg)+mean(nndist(ppp_tum)))
bb   <- bbox(sgdf_srtm)    
win  <- owin(xrange=c(bb[1,1],bb[1,2]), yrange= c(bb[2,1],bb[2,2]), unitname="m")
meg_dens <- density(ppp_meg, kernel="gaussian", sigma=sdev, dimyx=c(36,56), w=win,  edge=TRUE, at="pixels")
tum_dens <- density(ppp_tum, kernel="gaussian", sigma=sdev, dimyx=c(36,56), w=win,  edge=TRUE, at="pixels")
vil_dens <- density(ppp_vil, kernel="gaussian", sigma=sdev, dimyx=c(36,56), w=win,  edge=TRUE, at="pixels")

samppt   <- spsample(sgdf_srtm, 500,  type="random")

library(sp)
library(raster)
sgdf_meg_dens     <- as.SpatialGridDataFrame.im(meg_dens)
projection(sgdf_meg_dens) <- crs1
meg_dens_samp     <- sp::over(samppt, sgdf_meg_dens)
sgdf_tum_dens     <- as.SpatialGridDataFrame.im(tum_dens)
projection(sgdf_tum_dens) <- crs1
tum_dens_samp     <- sp::over(samppt, sgdf_tum_dens)
sgdf_vil_dens     <- as.SpatialGridDataFrame.im(vil_dens)
projection(sgdf_vil_dens) <- crs1
vil_dens_samp     <- sp::over(samppt, sgdf_vil_dens)
elev_dens_samp     <- sp::over(samppt, sgdf_srtm)

dens_samp         <- meg_dens_samp
names(dens_samp)[names(dens_samp) == 'v'] <- 'meg'
dens_samp    <- cbind(dens_samp, tum_dens_samp)
names(dens_samp)[2] <- 'tum'
dens_samp   <- cbind(dens_samp, vil_dens_samp)
names(dens_samp)[3] <- 'vil'
dens_samp   <- cbind(dens_samp, elev_dens_samp)
names(dens_samp)[4] <- 'elev'

cor.test(dens_samp$meg, dens_samp$tum,  method="p")
cor.test(dens_samp$vil, dens_samp$tum,  method="p")
cor.test(dens_samp$elev, dens_samp$tum,  method="p")


prcomp(na.omit(dens_samp))

meg_elev_samp     <- sp::over(spdf_meg, sgdf_srtm)
tum_elev_samp     <- sp::over(spdf_tum, sgdf_srtm)
vil_elev_samp     <- sp::over(spdf_vil, sgdf_srtm)
ks.test(meg_elev_samp$srtm, tum_elev_samp$srtm)
ks.test(meg_elev_samp$srtm, vil_elev_samp$srtm)
ks.test(tum_elev_samp$srtm, vil_elev_samp$srtm)

# 3. second order properties ===============================
library(spatstat)
meg_env_g <- envelope(Y = ppp_meg, fun = Gest, nrank = 2, nsim = 99)   
pdf("pictures/c7_meg_G.pdf", height=4, width=6, bg = "white") 
    plot(meg_env_g)
dev.off() 

meg_env_f <- envelope(Y = ppp_meg, fun = Fest, nrank = 2, nsim = 99)   
pdf("pictures/c7_meg_F.pdf", height=4, width=6, bg = "white") 
    plot(meg_env_f)
dev.off() 

meg_env_k <- envelope(Y = ppp_meg, fun = Kest, nrank = 2, nsim = 99)   
pdf("pictures/c7_meg_K.pdf", height=4, width=6, bg = "white") 
    plot(meg_env_k)
dev.off() 

meg_env_L <- envelope(Y = ppp_meg, fun = Lest, nrank = 2, nsim = 99)   
pdf("pictures/c7_meg_L.pdf", height=4, width=6, bg = "white") 
    plot(meg_env_L)
dev.off() 

meg_env_J <- envelope(Y = ppp_meg, fun = Jest, nrank = 2, nsim = 99)   
pdf("pictures/c7_meg_J.pdf", height=4, width=6, bg = "white") 
    plot(meg_env_J)
dev.off() 


# mowing window
rw      <- 1000   
xmin    <- sgdf_meg_dens@bbox[1,1]
xmax    <- sgdf_meg_dens@bbox[1,2]
ymin    <- sgdf_meg_dens@bbox[2,1]
ymax    <- sgdf_meg_dens@bbox[2,2]
dx      <- xmax - xmin
dy      <- ymax - ymin
xmin    <- xmin - (dx / 5)
xmax    <- xmax + (dx / 5)
ymin    <- ymin - (dy / 5)
ymax    <- ymax + (dy / 5)
rows  <- round((ymax-ymin)/rw, 0) + 1
columns <- round((xmax-xmin)/rw, 0) + 1
z <- cbind(1:(columns*rows))
df <- data.frame(z)
gt <- GridTopology(c(xmin - rw/2,ymin - rw/2), c(rw,rw), c(columns,rows))
ras <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))

fs_nn <- nndist(ppp_tum)                    
radius <- 10000                        
r <- seq(0, radius, 250)             
win <- owin(xrange=c(xmin=xmin,xmax=xmax), yrange=c(ymin=ymin,ymax=ymax), unitname="m")   
for (i in seq(along=ras@data$z))  {   
    xr      <- coordinates(ras)[i,1]
    yr      <- coordinates(ras)[i,2]
    distances <- sqrt((ppp_tum$x -xr)^2 + (ppp_tum$y -yr)^2)     
    indiz <- which(distances<radius)       
    if (length(indiz) > 2 ) {        
        x <- ppp_tum$x[indiz]
        y <- ppp_tum$y[indiz]
        name <- ppp_tum$x[indiz]      
        fsgf <- SpatialPointsDataFrame(cbind(x, y), as.data.frame(name), proj4string= CRS(as.character(crs1)))
        fspt <- ppp(fsgf@coords[,1], fsgf@coords[,2], window=win)  
        gfs <- Gest(fspt, r=r,correction="km") 
        value <-  mean(gfs$theo-gfs$km)       
        ras@data$z[i]   <- value             
    }
    else {ras@data$z[i] <-0}
}

pdf("pictures/c7_movingWindowG.pdf", height=5, width=6, bg = "white") 
    image(ras, col = gray.colors(20))    
    points(ppp_tum$x, ppp_tum$y, pch=16, cex=0.4)     
dev.off()

writeAsciiGrid(ras,  "pictures/c7_movingWindowG.asc",   attr = 1, na.value = -999999, dec=".")                                 

# # complicated point patterns
# library(rgdal)
# dsn <- paste("data/",sep="")
# pp1 <- readOGR(dsn, layer='pp1') 
# pp2 <- readOGR(dsn, layer='pp2x') 
# pp3 <- readOGR(dsn, layer='pp3x') 
# pp4 <- readOGR(dsn, layer='pp4x') 
#
# library(spatstat)
# bb = bbox(sgdf_srtm)    
# win <- owin(xrange=c(bb[1,1],bb[1,2]), yrange=c( bb[2,1],bb[2,2]), unitname="m")
# ppp1 <- ppp(pp1@coords[,1], pp1@coords[,2], window=win)
# ppp2 <- ppp(pp2@coords[,1], pp2@coords[,2], window=win)
# ppp3 <- ppp(pp3@coords[,1], pp3@coords[,2], window=win)
# ppp4 <- ppp(pp4@coords[,1], pp4@coords[,2], window=win)
# 
# svg("6pictures/c7_ppp1g.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp1, fun = Gest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp1f.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp1, fun = Fest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp1k.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp1, fun = Kest, nrank = 2, nsim = 50) ) 
# dev.off()
# 
# svg("6pictures/c7_ppp2g.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp2, fun = Gest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp2f.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp2, fun = Fest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp2k.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp2, fun = Kest, nrank = 2, nsim = 50) ) 
# dev.off()
# 
# svg("6pictures/c7_ppp3g.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp3, fun = Gest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp3f.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp3, fun = Fest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp3k.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp3, fun = Kest, nrank = 2, nsim = 50) ) 
# dev.off()
# 
# svg("6pictures/c7_ppp4g.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp4, fun = Gest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp4f.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp4, fun = Fest, nrank = 2, nsim = 50) ) 
# dev.off()
# svg("6pictures/c7_ppp4k.svg", height=5, width=6, bg = "white") 
# plot(envelope(Y = ppp4, fun = Kest, nrank = 2, nsim = 50) ) 
# dev.off()
# 
# svg("6pictures/c7_ppp1.svg", height=5, width=6, bg = "white") 
# plot(ppp1, pch=16)
# dev.off()
# 
# svg("6pictures/c7_ppp2.svg", height=5, width=6, bg = "white") 
# plot(ppp2, pch=16)
# dev.off()
# 
# svg("6pictures/c7_ppp3.svg", height=5, width=6, bg = "white") 
# plot(ppp3, pch=16)
# dev.off()
# 
# svg("6pictures/c7_ppp4.svg", height=5, width=6, bg = "white") 
# plot(ppp4, pch=16)
# dev.off()


# 4. third order properties ===============================
meg_env_t <- envelope(ppp_meg, fun = Tstat, nrank = 2, nsim = 20)   
pdf("pictures/c7_meg_T.pdf", height=4, width=6, bg = "white") 
    plot(meg_env_t)
dev.off()

save.image("ws/ws07.rws")


