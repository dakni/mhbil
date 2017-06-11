############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 8
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 06
## Date of last changes: 11.08.2015
## Data: Megaliths, Bronze Age Barrows and Medieval villages
## Author of data: gshdl
## Purpose: didactic
## Content: 1. territoriality, 2. Boundaries of cultural areas
##          3. empirical models, 4. theoretical models
## Description: applies some basic techniques of
##             boundary analysis
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)     
############################################## #

# 0. Preparation ===============================
#  working direktory
wd <- "~/modproj_qaam"
wd <- "/home/fon/daten/analyse/modproj_qaam"
setwd(wd)
load("4ws/ws07.rws")

# 3. empirical models ===============================

library(spatstat)
ch_meg <- convexhull(ppp_meg)
ch_tum <- convexhull(ppp_tum)
plot(ch_tum, border="grey", main="")
plot(ch_meg, add=TRUE)
points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)


buf_tum <- dilation(ppp_tum, 1000, polygonal=TRUE, tight=F) 
buf_meg <- dilation(ppp_meg, 1000, polygonal=TRUE, tight=F) 
plot(buf_tum, border="grey", main="")
plot(buf_meg, add=TRUE)
points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6,  col="grey")  
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)  

par(mfcol=c(1,2))
    plot(ch_tum, border="grey", main="")
    plot(ch_meg, add=TRUE)
    points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
    points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)
    
    plot(buf_tum, border="grey", main="")
    plot(buf_meg, add=TRUE)
    points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6,  col="grey")  
    points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4) 

library("classInt")
nb_meg <- classIntervals(dens_samp@data$meg, style =  "fisher", dataPrecision = NULL)  
nb_tum <- classIntervals(dens_samp@data$tum, style = "fisher", dataPrecision = NULL)  
#pdf("6pictures/c8_contour.pdf", height=4, width=6, bg = "white") 
par(mai = c(0, 0, 0, 0))
    contour(sgdf_meg_dens, add=F, method = "edge", levels = nb_meg$brks, drawlabels = F)    
    contour(sgdf_tum_dens, add=T, method = "edge", levels  = nb_tum$brks, drawlabels = F, col="grey")  
    points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
    points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)
#dev.off() 


ddif <- sgdf_meg_dens
ddif@data$v <- sgdf_meg_dens$v - sgdf_tum_dens$v

pdf("6pictures/c8_contourdiff.pdf", height=2.5, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    contour(sgdf_meg_dens, add=F, method = "edge", levels = nb_meg$brks, drawlabels = F)    
    contour(sgdf_tum_dens, add=T, method = "edge", levels  = nb_tum$brks, drawlabels = F, col="grey")  
    points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
    points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)
    
    contour(ddif, add=F, levels = c(0), drawlabels = F) 
    points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
    points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)
par(mfcol=c(1,1))
dev.off() 


dens_samp2 <- dens_samp[,1:3]
dens_samp2@data[,1] <- dens_samp2@data[,1]  / max(dens_samp2@data[,1])
dens_samp2@data[,2] <- dens_samp2@data[,2]  / max(dens_samp2@data[,2])
dens_samp2@data[,3] <- dens_samp2@data[,3]  / max(dens_samp2@data[,3])

distances <- dist(dens_samp2@data, method = "euclidean")
hc <- hclust(distances, method="centroid")
pdf("6pictures/c8_clusterhc.pdf", height=4, width=6, bg = "white") 
#par(mai = c(0, 0, 0, 0))
    plot(hc, labels=FALSE)
dev.off() 

library("cluster")
widthssum <- c(
    sum(pam(dens_samp2@data, 2, metric = "euclidean")$silinfo$clus.avg.widths),
    sum(pam(dens_samp2@data, 3, metric = "euclidean")$silinfo$clus.avg.widths),
    sum(pam(dens_samp2@data, 4, metric = "euclidean")$silinfo$clus.avg.widths),
    sum(pam(dens_samp2@data, 5, metric = "euclidean")$silinfo$clus.avg.widths),
    sum(pam(dens_samp2@data, 6, metric = "euclidean")$silinfo$clus.avg.widths),
    sum(pam(dens_samp2@data, 7, metric = "euclidean")$silinfo$clus.avg.widths),
    sum(pam(dens_samp2@data, 8, metric = "euclidean")$silinfo$clus.avg.widths),
    sum(pam(dens_samp2@data, 9, metric = "euclidean")$silinfo$clus.avg.widths))
pdf("6pictures/c8_clustersil.pdf", height=4, width=6, bg = "white") 
plot(widthssum, type ="b", pch=16)
dev.off() 

dens_samp_clus <- pam(dens_samp2@data, 4, metric = "euclidean")
dens_samp2@data <- cbind(dens_samp2@data,  dens_samp_clus$clustering)
names(dens_samp2)[names(dens_samp2) == 'dens_samp_clus$clustering'] <- 'clus'

pdf("6pictures/c8_cluster.pdf", height=4, width=6, bg = "white") 
par(mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(20, start = 0.8, end = 0.2))
    points(dens_samp, pch=dens_samp2@data$clus)
dev.off() 

clus1 <- c(mean(dens_samp@data[dens_samp_clus$clustering==1, 1]),
            mean(dens_samp@data[dens_samp_clus$clustering==1, 2]),
            mean(dens_samp@data[dens_samp_clus$clustering==1, 3]))
clus2 <- c(mean(dens_samp@data[dens_samp_clus$clustering==2, 1]),
           mean(dens_samp@data[dens_samp_clus$clustering==2, 2]),
           mean(dens_samp@data[dens_samp_clus$clustering==2, 3]))
clus1
clus2

dpd <- ppp_meg
dpd$d <- density(ppp_meg, 1000, edge=TRUE, at="points")          
k <- 1                                                         
s <- -20000000                                                 
dpd$id <- seq(along=dpd$x); dpd$v_id <- 0; dpd$v_x <- 0
dpd$v_y  <- 0
dpd$v_d  <- 0
dpd$dist <- 0; dpd$e <- 0 
dpd.d <- dist(cbind(dpd$x,dpd$y),upper = T)               
dpd.m <- as.matrix(dpd.d)                                   
maxd <- max(dpd$d)
maxm <- max(dpd.m)
dpd$d <- dpd$d * maxm / maxd                           
dpd.e <- dpd.m                               
for (i in seq(along=dpd$d))   {                     
   dpd.e[i,] <- (dpd$d - dpd$d[i]) - (k * dpd.m[i,])         
   w <-max(dpd.e[i,])                                        
   wi <- which(dpd.e[i,] == w)                                
   if (i == wi | w < s)    {dpd$v_id[i] <- 0}              
   else                    {dpd$v_id[i] <- wi}
    }
dpd$v_id  


for (i in seq(along=dpd$d))   {                           
    if (dpd$v_id[i] > 0)  {
        dpd$v_x[i] <- dpd$x[dpd$v_id[i]]
        dpd$v_y[i] <- dpd$y[dpd$v_id[i]]
        dpd$v_d[i] <- dpd$d[dpd$v_id[i]]
        dpd$dist[i]<- dpd.m[i,dpd$v_id[i]]
        dpd$e[i]   <- dpd.e[i,dpd$v_id[i]]
        }
     else {
        dpd$v_x[i] <- dpd$x[i]
        dpd$v_y[i] <- dpd$y[i]
        dpd$v_d[i] <- dpd$d[i]
        }
    }
dpd

LinesList <- list()                                   
for(i in seq(along=dpd$id)) {                   
    m <- matrix(data = c(dpd$x[i],dpd$v_x[i],dpd$y[i], dpd$v_y[i]), nrow=2, ncol=2)   
    L <- Line(m)
    LL <- list(L)
    name  <- paste("Zuordnung_", dpd$id,"_", dpd$v_id, sep="") 
    LLL <- Lines(LL, ID = name[i])                        
    LinesList[length(LinesList)+1] <- LLL  
    }
sl <- SpatialLines(LinesList, proj4string = CRS(crs1))

pdf("6pictures/c8_denscluster.pdf", height=4, width=6, bg = "white") 
par(mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(20, start = 0.8, end = 0.2) )
    lines(sl)
    points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.6, col="black")  
dev.off() 



# 4. theoretical models ===============================
library(spatstat)
bb   <- bbox(sgdf_srtm)    
win  <- owin(xrange=c(bb[1,1],bb[1,2]), yrange= c(bb[2,1],bb[2,2]), unitname="m")
tum_dens <- density(ppp_tum, kernel="gaussian", sigma=1000, dimyx=c(36,56), w=win,  edge=TRUE, at="pixels")
library(maptools)
sgdf_tum_dens     <- as.SpatialGridDataFrame.im(tum_dens)

ras <- sgdf_tum_dens
r   <- 3000      
ras@data$v[which(is.na(ras@data$v))] <- 0
m <- max(ras@data$v)
s <- m / 20                          
indmax  <- c()
indplan <- c()

for (i in seq(along=ras@data$v))  {      
    x <- coordinates(ras)[i,1]
    y <- coordinates(ras)[i,2]
    z <- ras@data$v[i]
    indx <- which((coordinates(ras)[,1] > x - r) &  (coordinates(ras)[,1] < x + r))
    indy <- which((coordinates(ras)[,2] > y - r) &  (coordinates(ras)[,2] < y + r))
    indxy <- intersect(indx,indy)
    if (max(ras@data[indxy,1]) == z & z > s)  {indmax[length(indmax)+1]   <- i}  
    if (sd(ras@data[indxy,1]) == 0)  {indplan[length(indplan)+1] <- i}
    rm(indx)
    rm(indy)
    rm(indxy)
    }

mn  <- length(indmax)      
mx  <- coordinates(ras)[indmax,1]
my  <- coordinates(ras)[indmax,2]
mx2 <- coordinates(ras)[indplan,1]
my2 <- coordinates(ras)[indplan,2]
mz  <- ras@data[indmax,1]
maxima <- data.frame(cbind(mx,my,mz))

pdf("6pictures/c8_denscent.pdf", height=4, width=6, bg = "white") 
par(mai = c(0, 0, 0, 0))
    image(ras, col = gray.colors(20, start = 0.8, end = 0.2))
    points(maxima$mx, maxima$my, pch=16, col="black")  
dev.off() 

library(deldir)
try <- deldir(maxima[,1],maxima[,2],plot=TRUE,wl='te')

cent <- data.frame(cbind(id=seq(1:length(maxima[,1])), x=maxima[,1],y=maxima[,2],meg=0, tum=0))
coordinates(cent)=~x+y
proj4string(cent)  <- CRS(as.character(crs1))
cent_meg  <- overlay(x=sgdf_meg_dens, y=cent)
cent_tum  <- overlay(x=sgdf_tum_dens, y=cent)
cent@data$meg <- cent_meg@data$v
cent@data$tum <- cent_tum@data$v
dens_samp@data <- cbind(dens_samp@data,cent=0)

for(i in seq(along=(dens_samp@data$cent))) {
     d1 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[1])^2 + (dens_samp@data$tum[i]  - cent@data$tum[1])^2)  
     d2 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[2])^2 + (dens_samp@data$tum[i]  - cent@data$tum[2])^2) 
     d3 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[3])^2 + (dens_samp@data$tum[i]  - cent@data$tum[3])^2) 
     d4 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[4])^2 + (dens_samp@data$tum[i]  - cent@data$tum[4])^2) 
     d5 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[5])^2 + (dens_samp@data$tum[i]  - cent@data$tum[5])^2) 
     d6 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[6])^2 + (dens_samp@data$tum[i]  - cent@data$tum[6])^2) 
     d7 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[7])^2 + (dens_samp@data$tum[i]  - cent@data$tum[7])^2) 
     d8 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[8])^2 + (dens_samp@data$tum[i]  - cent@data$tum[8])^2) 
     d9 <- sqrt((dens_samp@data$meg[i] - cent@data$meg[9])^2 + (dens_samp@data$tum[i]  - cent@data$tum[9])^2) 
     d <- c(d1,d2,d3,d4,d5,d6,d7,d8,d9)
     mindist <- min(d1,d2,d3,d4,d5,d6,d7,d8,d9)
     id  <- which(d == mindist)
     dens_samp@data$cent[i] <- id
     }

pdf("6pictures/c8_voro.pdf", height=2.5, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(20, start = 0.8, end = 0.2))
    plot(try, add=TRUE)

    image(sgdf_srtm, col = gray.colors(20, start = 0.8, end = 0.2))
    points(dens_samp, pch=dens_samp@data$cent)
par(mfcol=c(1,1))
dev.off() 

save.image("4ws/ws08.rws")




