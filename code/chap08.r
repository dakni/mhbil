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

spdf_meg <- df_meg <- read.table("../data/meg_dw.csv", sep=";", header=TRUE)
spdf_tum <- df_tum <- read.table("../data/tum_dw.csv", sep=";", header=TRUE)
spdf_vil <- df_vil <- read.table("../data/villages.csv", sep=",", header=TRUE)

library(sp)
crs1 <- CRS("+init=epsg:31467")
coordinates(spdf_meg)= ~x+y
coordinates(spdf_tum)= ~x+y
coordinates(spdf_vil)= ~x+y

proj4string(spdf_meg) <- CRS("+init=epsg:31467")
proj4string(spdf_tum) <- CRS("+init=epsg:31467")
proj4string(spdf_vil) <- CRS("+init=epsg:4326")
spdf_vil <- spTransform(x = spdf_vil,
                        CRSobj = CRS("+init=epsg:31467")
                        )

library(rgdal)
srtm <- readGDAL("../data/dw_gk3_50_ag.asc")
proj4string(srtm) <- CRS("+init=epsg:31467")
names(srtm@data) <- "srtm"

library(spatstat)
ppp_meg <- ppp(x = df_meg$x,
               y = df_meg$y,
               window = owin(xrange = srtm@bbox[1,],
                             yrange = srtm@bbox[2,]
                             )
               )
ppp_tum <- ppp(x = df_tum$x,
               y = df_tum$y,
               window = owin(xrange = srtm@bbox[1,],
                             yrange = srtm@bbox[2,]
                             )
               )
ppp_vil <- ppp(x = spdf_vil@coords[,1],
               y = spdf_vil@coords[,2],
               window = owin(xrange = srtm@bbox[1,],
                             yrange = srtm@bbox[2,]
                             )
               )
sdev <- 2 * (mean(nndist(ppp_meg)) + mean(nndist(ppp_tum)))
win  <- owin(xrange = c(bbox(srtm)[1,1],bbox(srtm)[1,2]),
             yrange = c(bbox(srtm)[2,1],bbox(srtm)[2,2]),
             unitname = "m"
             )

meg_dens <- density(ppp_meg, kernel="gaussian", sigma=sdev, dimyx=c(36,56), w=win,  edge=TRUE, at="pixels")
tum_dens <- density(ppp_tum, kernel="gaussian", sigma=sdev, dimyx=c(36,56), w=win,  edge=TRUE, at="pixels")
vil_dens <- density(ppp_vil, kernel="gaussian", sigma=sdev, dimyx=c(36,56), w=win,  edge=TRUE, at="pixels")
samppt   <- spsample(x = srtm, n = 500,  type = "random")

library(raster)
meg_dens_samp <- raster::extract(x = raster(meg_dens),  y = samppt)
tum_dens_samp <- raster::extract(x = raster(tum_dens),  y = samppt)
vil_dens_samp <- raster::extract(x = raster(vil_dens),  y = samppt)

dens_samp <- data.frame(meg = meg_dens_samp,
                        tum = tum_dens_samp,
                        vil = vil_dens_samp
                        )

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
nb_meg <- classIntervals(
  var = dens_samp$meg,
  style =  "fisher"
)
nb_tum <- classIntervals(
  var = dens_samp$tum,
  style = "fisher"
)  

contour(x = meg_dens, add = F,
        method = "edge",
        levels = nb_meg$brks,
        drawlabels = F, main = "")    
contour(x = tum_dens,
        add = T, method = "edge",
        levels  = nb_tum$brks,
        drawlabels = F, col="grey")  
points(x = ppp_tum$x, y = ppp_tum$y,
       pch = 17, cex = .6,  col = "grey")
points(x = ppp_meg$x, y = ppp_meg$y,
       pch = 16, cex = .4)



ddif <- meg_dens
ddif$v <- meg_dens$v - tum_dens$v

par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
contour(x = meg_dens, add = F,
        method = "edge",
        levels = nb_meg$brks,
        drawlabels = F, main = "")    
contour(x = tum_dens,
        add = T, method = "edge",
        levels  = nb_tum$brks,
        drawlabels = F, col="grey")  
points(x = ppp_tum$x, y = ppp_tum$y,
       pch = 17, cex = .6,  col = "grey")
points(x = ppp_meg$x, y = ppp_meg$y,
       pch = 16, cex = .4)
    
contour(ddif, add=F, levels = c(0), drawlabels = F) 
points(ppp_tum$x, ppp_tum$y, pch=17, cex=0.6, col="grey")  
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4)
par(mfcol=c(1,1))


dens_samp2 <- as.data.frame(apply(X = dens_samp,
                               MARGIN = 2,
                               FUN = function(x) { x / max(x)}
                               )
                            )

distances <- dist(dens_samp2, method = "euclidean")
hc <- hclust(distances, method="centroid")
plot(hc, labels=FALSE)

library("cluster")
widthssum <- c()

for ( i in 1:8 ) {
  widthssum[i] <-
    sum(
      pam(x = dens_samp2,
          k = i+1,
          metric = "euclidean"
          )$silinfo$clus.avg.widths
    )
}
plot(widthssum, type ="b", pch=16)

dens_samp2$clustering <- pam(x = dens_samp2, 
                             k = 4,
                             metric = "euclidean"
                             )$clustering
samppt <- SpatialPointsDataFrame(coords = samppt@coords,
                                 data = dens_samp2,
                                 proj4string = samppt@proj4string
                                 )

image(srtm,
      col = gray.colors(20, start = 0.8, end = 0.2)
      )
points(samppt,
       pch = samppt$clustering
       )

clus1 <- c(mean(dens_samp[dens_samp2$clustering==1, 1]),
           mean(dens_samp[dens_samp2$clustering==1, 2]),
           mean(dens_samp[dens_samp2$clustering==1, 3])
           )

clus2 <- c(mean(dens_samp[dens_samp2$clustering==2, 1]),
           mean(dens_samp[dens_samp2$clustering==2, 2]),
           mean(dens_samp[dens_samp2$clustering==2, 3])
           )
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

LLL <- list()                                   
for(i in seq(along=dpd$id)) {                   
    L <- Line(coords = rbind(c(dpd$x[i],dpd$y[i]),
                             c(dpd$v_x[i], dpd$v_y[i])
                             )
              )
    LL <- list(L)
    name  <- paste("Belonging_", dpd$id,"_", dpd$v_id, sep="") 
    LLL[[i]] <- Lines(LL, ID = name[i])    
    }
sl <- SpatialLines(LLL, proj4string = CRS("+init=epsg:31467"))

image(srtm, col = gray.colors(20, start = 0.8, end = 0.2) )
lines(sl)
points(x = ppp_meg$x, y = ppp_meg$y, pch = 16, cex = .6, col = "black")  



# 4. theoretical models ===============================
ras <- raster(tum_dens)

r   <- 3000      
m <- max(values(ras))
s <- m / 20                          
indmax  <- c()
indplan <- c()

for (i in seq(along=values(ras)))  {      
    x <- coordinates(ras)[i,1]
    y <- coordinates(ras)[i,2]
    z <- values(ras)[i]
    indx <- which((coordinates(ras)[,1] > x - r) &  (coordinates(ras)[,1] < x + r))
    indy <- which((coordinates(ras)[,2] > y - r) &  (coordinates(ras)[,2] < y + r))
    indxy <- intersect(indx,indy)
    if (max(values(ras)[indxy]) == z & z > s)  {
        indmax[length(indmax)+1]   <- i
    }  
}

mn  <- length(indmax)      
mx  <- coordinates(ras)[indmax,1]
my  <- coordinates(ras)[indmax,2]
mz  <- values(ras)[indmax]
maxima <- data.frame(cbind(mx,my,mz))

plot(ras, col = gray.colors(20, start = 0.8, end = 0.2))
points(maxima$mx, maxima$my, pch=16, col="black")   

library(deldir)
try <- deldir(maxima[,1],maxima[,2],plot=TRUE,wl='te')

cent <- data.frame(id = seq(1:length(maxima[,1])),
                   x = maxima[,1],
                   y = maxima[,2],
                   meg = 0,
                   tum = 0
                   )
coordinates(cent) <- ~x+y
proj4string(cent) <- CRS("+init=epsg:31467")

cent_meg <- raster::extract(x = raster(meg_dens),
                            y = cent)
cent_tum <- raster::extract(x = raster(tum_dens),
                            y = cent)
cent@data$meg <- cent_meg
cent@data$tum <- cent_tum
dens_samp <- cbind(dens_samp, cent=0)

dens_samp$cent <- 0
d <- c()

for(i in seq_along(dens_samp[,1])) {        
  for (j in seq_along(cent@data[,1])) {
    d[j] <- sqrt((dens_samp$meg[i] - cent@data$meg[j])^2 + (dens_samp$tum[i]  - cent@data$tum[j])^2)
  }
  mindist <- min(d)
  id  <- which(d == mindist)
  dens_samp$cent[i] <- id
}

samppt$cent <- dens_samp$cent

par(mfrow = c(1,2), mai = c(0, 0, 0, 0))
image(srtm, col = gray.colors(20, start = 0.8, end = 0.2) )
plot(try, add=TRUE)
image(srtm, col = gray.colors(20, start = 0.8, end = 0.2))
points(samppt, pch = samppt$cent)
