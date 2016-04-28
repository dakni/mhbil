############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 10
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and   
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 060
## Date of last changes: 08.09.2015
## Data: Megaliths, Bronze Age Barrows and Medieval villages
## Author of data: gshdl
## Purpose: didactic
## Content: 1. empirical models, 2. theoretical models
## Description: applies some basic techniques of
##             interaction analysis
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #

# 0. Preparation ===============================
#  working direktory
wd <- "~/modproj_qaam"
wd <- "/home/fon/daten/analyse/modproj_qaam"
setwd(wd)
load("4ws/ws09.rws")

# 1.  empirical models  ===============================

# indicators
library(spdep)
voro_nb_del <- tri2nb(coordinates(fs_vd_spdf), row.names=row.names(as.data.frame(fs_vd_spdf))) 
moran.test(fs_vd_spdf@data$fs_vd, nb2listw(voro_nb_del, style="W"))


# dd1
spoints <- data.frame(cbind(x=(-seq(20:1)*1000) + 3571203, y=rep(6036796, 20)))
coordinates(spoints)=~x+y
proj4string(spoints)  <- CRS(as.character(crs1)) 
i_kde <- extract(raster(sgdf_meg_dens), spoints)
mdistance <- i_kde[1] - i_kde

pdf("6pictures/c10_dd1.pdf", height=4, width=6, bg = "white") 
    plot(mdistance, col="black", pch=16, xlab="spatial distance (km)", ylab= "density distance")
    lines(mdistance,lty=1,col="black")
dev.off() 


# dd3
dist <- seq(from=0, to=33000, by=1000)
cdist <- id <- 1:length(dist)
dresult <- cbind(id, dist, cdist[] <- 0)
samppt   <- spsample(sgdf_srtm, 500,  type="regular")
i_kde <- extract(raster(sgdf_meg_dens), samppt)
ref_kde <- extract(raster(sgdf_meg_dens), cbind(3571203,6036796))
dmat <- matrix(1:length(dist)*length(sgdf_meg_dens@data$v), nrow=length(dist), ncol=length(sgdf_meg_dens@data$v))
dmat[] <- as.double(dmat[] <- NA)
mmean <- function (x) {mean(x, na.rm=TRUE)}
edist <-function(x1,x2,y1,y2) {sqrt((x1 - x2)^2 + (y1 - y2)^2)}

for (i in seq_along(samppt@coords[,1])) {
    x1 <- 3571203; y1 <- 6036796
    x2 <- samppt@coords[i,1]
    y2 <- samppt@coords[i,2]
    sdist <-edist(x1,x2,y1,y2)
    dind <- floor(sdist/1000) + 1
    dmat[dind,i] <- abs(i_kde[i] - ref_kde)
    }
dresult[,3] <-  apply(dmat, 1, mmean)

pdf("6pictures/c10_dd3.pdf", height=4, width=6, bg = "white") 
    plot(x=dresult[,2], y=dresult[,3], col="black", pch=16, xlab="spatial distance (m)", ylab= "density distance")
    lines(x=dresult[,2], y=dresult[,3],lty=1,col="black")
dev.off() 


# dd2
dist <- seq(from=0, to=33000, by=1000)
cdist2 <- cdist <- id <- 1:length(dist)
dresult <- cbind(id, dist, cdist[] <- 0, cdist2[] <- 0)
samppt   <- spsample(sgdf_srtm, 500,  type="regular")
i_kde <- extract(raster(sgdf_meg_dens), samppt)
ref_kde <- extract(raster(sgdf_meg_dens), cbind(3571203,6036796))
dmat <- matrix(1:length(dist)*length(sgdf_meg_dens@data$v), nrow=length(dist), ncol=length(sgdf_meg_dens@data$v))
dmat[] <- as.double(dmat[] <- NA)
mmean <- function (x) {mean(x, na.rm=TRUE)}
edist <-function(x1,x2,y1,y2) {sqrt((x1 - x2)^2 + (y1 - y2)^2)}

for (i in seq_along(samppt@coords[,1])) {
    x1 <- 3571203; y1 <- 6036796
    x2 <- samppt@coords[i,1]; y2 <- samppt@coords[i,2]
    if (x1 > x2) {
    sdist <-edist(x1,x2,y1,y2)
    dind <- floor(sdist/1000) + 1
    dmat[dind,i] <- abs(i_kde[i] - ref_kde)
    }
}
dresult[,3] <-  apply(dmat, 1, mmean)

for (i in seq_along(samppt@coords[,1])) {
    x1 <- 3571203; y1 <- 6036796
    x2 <- samppt@coords[i,1]; y2 <- samppt@coords[i,2]
    if (x1 < x2) {
        sdist <-edist(x1,x2,y1,y2)
        dind <- floor(sdist/1000) + 1
        dmat[dind,i] <- abs(i_kde[i] - ref_kde)
    }
}
dresult[,4] <-  apply(dmat, 1, mmean)

pdf("6pictures/c10_dd2.pdf", height=4, width=6, bg = "white") 
    plot(x=dresult[,2], y=dresult[,4], col="grey", pch=16, xlab="spatial distance (m)", ylab= "density distance")
    lines(x=dresult[,2], y=dresult[,4],lty=1,col="grey")
    points(x=dresult[,2], y=dresult[,3], pch=16, col="black")
    lines(x=dresult[,2], y=dresult[,3],lty=1,col="black")
dev.off() 

# dd7
dist <- seq(from=0, to=33000, by=1000)
cdist <- id <- 1:length(dist)
dresult <- cbind(id, dist, cdist[] <- 0)
dresulta <- dresult
samppt   <- spsample(sgdf_srtm, 500,  type="regular")
i_kde <- extract(raster(sgdf_meg_dens), samppt)
#ref_kde <- extract(raster(sgdf_meg_dens), cbind(3571203,6036796))
dmat <- matrix(1:length(dist)*length(sgdf_meg_dens@data$v), nrow=length(dist), ncol=length(sgdf_meg_dens@data$v))
dmat[] <- as.double(dmat[] <- NA)
mmean <- function (x) {mean(x, na.rm=TRUE)}
edist <-function(x1,x2,y1,y2) {sqrt((x1 - x2)^2 + (y1 - y2)^2)}

for (j in seq_along(samppt@coords[,1])) {
    for (i in seq_along(samppt@coords[,1])) {
        x1 <- samppt@coords[j,1]; y1 <- samppt@coords[j,2]
        x2 <- samppt@coords[i,1]; y2 <- samppt@coords[i,2]
        if (x1 == x2 & y2 >= y1) {
            sdist <-edist(x1,x2,y1,y2)
            dind <- floor(sdist/1000) + 1
            dmat[dind,i] <- abs(i_kde[i] - i_kde[j])
        }
    }
    dresult[,3] <-  apply(dmat, 1, mmean)
    dresulta[,3] <- dresulta[,3] + dresult[,3]
    dresult <- cbind(id, dist, cdist[] <- 0)
    }

pdf("6pictures/c10_dd7.pdf", height=4, width=6, bg = "white") 
    plot(x=dresulta[,2], y=dresulta[,3], col="black", pch=16, xlab="spatial distance (m)", ylab= "density distance")
    lines(x=dresulta[,2], y=dresulta[,3],lty=1,col="black")
dev.off() 

# dd9
dist <- seq(from=0, to=33000, by=1000)
cdist <- id <- 1:length(dist)
dresult <- cbind(id, dist, cdist[] <- 0)
dresulta <- dresult
samppt   <- spsample(sgdf_srtm, 500,  type="regular")
i_kde <- extract(raster(sgdf_meg_dens), samppt)
#ref_kde <- extract(raster(sgdf_meg_dens), cbind(3571203,6036796))
dmat <- matrix(1:length(dist)*length(sgdf_meg_dens@data$v), nrow=length(dist), ncol=length(sgdf_meg_dens@data$v))
dmat[] <- as.double(dmat[] <- NA)
mmean <- function (x) {mean(x, na.rm=TRUE)}
edist <-function(x1,x2,y1,y2) {sqrt((x1 - x2)^2 + (y1 - y2)^2)}

for (j in seq_along(samppt@coords[,1])) {
    for (i in seq_along(samppt@coords[,1])) {
        x1 <- samppt@coords[j,1]; y1 <- samppt@coords[j,2]
        x2 <- samppt@coords[i,1]; y2 <- samppt@coords[i,2]
            sdist <-edist(x1,x2,y1,y2)
            dind <- floor(sdist/1000) + 1
            dmat[dind,i] <- abs(i_kde[i] - i_kde[j])
    }
    dresult[,3] <-  apply(dmat, 1, mmean)
    dresulta[,3] <- dresulta[,3] + dresult[,3]
    dresult <- cbind(id, dist, cdist[] <- 0)
}

pdf("6pictures/c10_dd9.pdf", height=4, width=6, bg = "white") 
    plot(x=dresulta[,2], y=dresulta[,3], col="black", pch=16, xlab="spatial distance (m)", ylab= "density distance")
    lines(x=dresulta[,2], y=dresulta[,3],lty=1,col="black")
dev.off() 

# 2.  theoretical models  ===============================
# fall off curves

ddecay1 <- function (d, k, j) {i <- k/d^j; return(i)}         # power
ddecay2 <- function (d, k, j) {i <- k*2.718282^(-j*d); return(i)}    # exponential 
ddecay3 <- function (d, k, j) {i <- (k/(d+k)^j); return(i)}   # pareto 
ddecay4 <- function (d, k, j) {i <- k*2.718282^(-j*(d^2); return(i)}    # gauss 
k <- 3
j <- 3
xval <- seq(0,1,0.01)

pdf("6pictures/c10_distancedecay.pdf", height=3, width=6, bg = "white") 
    par( mai = c(1, 1, 0.1, 0.1))
    matplot(xval, ddecay1(xval, k, j)*0.00005, ylim=c(0,1), type="l", xlab=expression(paste("distance")), ylab="interaction")
    lines(xval, ddecay2(xval, k, j)*0.3, type="l", lty=2)
    lines(xval, ddecay3(xval*10, k, j)*9, type="l", lty=3)
    lines(xval, ddecay4(xval*10, k, j)*0.2, type="l", lty=4)
    legend("topright", legend =c("power","exponential","Pareto","Gauss"),lty=c(1,2,3,4))
dev.off() 

# gravity models

library(tripack)
fsd <- tri.mesh(cent, duplicate = 'remove')   # ermittelt die delauny trinangulation
fsnn <- neighbours(fsd)     # Vektor mit Listen der Nachbarn 
LinesList <- list()         # leere LinesList Liste erstellen
sldf <- c();deldf_i <- c();deldf_x1 <- c();deldf_y1 <- c();deldf_k <- c();deldf_x2 <- 
    c();deldf_y2 <- c();deldf_name <- c();deldf_dens1 <- c();deldf_dens2 <- c()
for(i in seq(along=cent@coords[,1])) {         
    pid1 <- i                                  
    x1 <- cent@coords[i,1]
    y1 <- cent@coords[i,2]
    dens1 <- cent@data$meg[i]
    for(k in seq(along=(fsnn[i][[1]]))) {   
        pid2 <- fsnn[[i]][k]    
        if (pid2 > pid1) {     
            x2 <- cent@coords[pid2,1]
            y2 <- cent@coords[pid2,2] 
            dens2 <- cent@data$meg[pid2]
            m <- matrix(data = c(x1,x2,y1,y2), nrow=2, ncol=2)  
            L <- Line(m); LL <- list(L)
            name  <- paste("edge", "_", pid1,"_", pid2, sep="")
            LLL <- Lines(LL, ID = name)             
            LinesList[length(LinesList)+1] <- LLL   
            sldf[length(sldf)+1] <- name         
            j <- length(deldf_i) + 1
            deldf_i[j]    <- i
            deldf_x1[j]   <- x1; deldf_y1[j]   <- y1
            deldf_k[j]    <- pid2
            deldf_x2[j]   <- x2; deldf_y2[j]   <- y2
            deldf_name[j] <- name
            deldf_dens1[j] <- dens1
            deldf_dens2[j] <- dens2
        }
    }
}
deldf_c <- data.frame(deldf_i,deldf_x1,deldf_y1,deldf_k,deldf_x2,deldf_y2,deldf_name,deldf_dens1,deldf_dens2)
dist  <- sqrt((deldf_c[2] - deldf_c[5])^2 + (deldf_c[3] - deldf_c[6])^2 )
dimnames(dist) <- list(NULL,"dist")
deldf_c$dist   <- dist
inter <- deldf_c$deldf_dens1 * deldf_c$deldf_dens2 / deldf_c$dist^2
dimnames(inter) <- list(NULL, "inter")
deldf_c$inter   <- inter * 1e+22

library(igraph)
library(spdep)
co <- coordinates(cent)  
coords=as.matrix(coordinates(cent))
ids <- row.names(as.data.frame(cent))
meg_cent_nb_del <- tri2nb(coordinates(cent), row.names=ids) 
m <- nb2mat(meg_cent_nb_del)             
g <- graph.adjacency(m, mode="lower", weighted=T)
g <- set.vertex.attribute(g, "x", index=V(g), coordinates(cent)[,1])   
g <- set.vertex.attribute(g, "y", index=V(g), coordinates(cent)[,2]) 
g <- set.vertex.attribute(g, "dens", index=V(g), cent@data$meg) 
g <- set.edge.attribute(g, "distance2", index=E(g), deldf_c$dist)
g <- set.edge.attribute(g, "inter", index=E(g), deldf_c$inter)
E(g)$weight <- deldf_c$inter 
plot(g, edge.width=sqrt(deldf_c$inter), layout = coordinates(cent))

LinesList <- list()         # leere LinesList Liste erstellen
for(i in seq(along=deldf_c[,1])) {         
    pid1 <- i                                  
    x1 <- deldf_c$deldf_x1[i]; x1 <- deldf_c$deldf_y1[i]
    x2 <- deldf_c$deldf_x2[i]; y2 <- deldf_c$deldf_y2[i]
    m <- matrix(data = c(x1,x2,y1,y2), nrow=2, ncol=2)  
    L <- Line(m); LL <- list(L)
    LLL <- Lines(LL, ID = deldf_name[i])             
    LinesList[length(LinesList)+1] <- LLL   
    }
sl <- SpatialLines(LinesList, proj4string = CRS(as.character(crs1)))
sdf_inter_meg_cent <- SpatialLinesDataFrame(sl, deldf_c, match.ID = FALSE)

pdf("6pictures/c10_gravity_megcent.pdf", height=3.2, width=6, bg = "white") 
    par( mai = c(0, 0, 0, 0.2))
    plot(raster(sgdf_srtm), col = gray.colors(25, start =   0.97, end = 0.4))
    for (i in 1:length(sdf_inter_meg_cent@data$inter)){
    lines(sdf_inter_meg_cent[i,], lwd=sqrt(sdf_inter_meg_cent@data$inter[i]))}
dev.off() 



save.image("4ws/ws10.rws")

