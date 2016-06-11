############################################## #
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
## Author of data: gshdl
## Purpose: didactic
## Content: 10. preparation, 1. random numbers, 2. differential 
##          equations, 3. point patterns, 4. cellular automaton, 
##          5. ABM
## Description: applies some basic simulation techniques
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
#load("4ws/ws12.rws")

# 1.  Random numbers ===========================

r1 <- runif(11, 10, 25)
r1
r2 <- runif(11, 10, 25)    
r2

set.seed(444)        
r3 <- runif(11, 10, 25)      
r3
set.seed(444) 
r4<- runif(11, 10, 25)    
r4

set.seed(444) 
r5 <- runif(11, 10, 25)     
r5 <- floor(r5)
r5

set.seed(444) 
r6 <- sample(10:25,11)       
r6


set.seed(444) 
r7<- sample(10:25, 11, replace=T)
r7

set.seed(444) 
vt <- c("a","b","c","d")
r8<- sample(1:4, 11, replace=T)
vt[r8]

set.seed(444) 
sample(vt, 11, replace=T)

set.seed(444) 
r9 <- runif(25, 10, 25)        
pdf("6pictures/c12_r9.pdf", height=4, width=6, bg = "white") 
    hist(r9, breaks=10:26, col="grey")
dev.off() 

set.seed(444) 
r10 <- runif(25000, 10, 25)     
pdf("6pictures/c12_r10.pdf", height=4, width=6, bg = "white") 
    hist(r10, breaks=10:26, col="grey")
dev.off() 

set.seed(444) 
r11 <- rnorm(25000, 10, 25)   
pdf("6pictures/c12_r11.pdf", height=4, width=6, bg = "white") 
    hist(r11[r11 > 10 & r11 < 25], breaks=10:25, col="grey", main="Histogram of r11", xlab="r11")
dev.off() 

pdf("6pictures/c12_r11.pdf", height=4, width=6, bg = "white") 
    hist(r11, col="grey", main="Histogram of r11", xlab="r11")
dev.off() 


pdf("6pictures/c12_r9-10.pdf", height=3, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(1, 1, 0.5, 0))
    hist(r9, breaks=10:26, col="grey")
par(mai = c(1, 0.5, 0.5, 0))
    hist(r10, breaks=10:26, col="grey")
par(mfcol=c(1,1))
dev.off() 

pdf("6pictures/c12_r11-11.pdf", height=3, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(1, 1, 0.5, 0))
    hist(r11[r11 > 10 & r11 < 25], breaks=10:25, col="grey", main="Histogram of r11", xlab="r11")
par(mai = c(1, 0.5, 0.5, 0))
    hist(r11, col="grey", main="Histogram of r11", xlab="r11")
par(mfcol=c(1,1))
dev.off() 

# 2.  Differential equations ===========================



# 3.  Point patterns  ===========================

x1 <- bbox(sgdf_srtm)[1,1]; x2 <- bbox(sgdf_srtm)[1,2]
y1 <- bbox(sgdf_srtm)[2,1]; y2 <- bbox(sgdf_srtm)[2,2]
set.seed(444) 
xs1 <-  runif(100, x1, x2) 
ys1 <-  runif(100, y1, y2) 
spdf_rpp1 <- SpatialPointsDataFrame(cbind(xs1,ys1), as.data.frame(cbind(xs1,ys1)), proj4string=CRS(as.character(crs1)))

set.seed(444) 
xs2 <-  rnorm(100, mean(c(x1, x2)), 4000)   
ys2 <-  rnorm(100, mean(c(y1, y2)), 2000)  
spdf_rpp2 <- SpatialPointsDataFrame(cbind(xs2,ys2), as.data.frame(cbind(xs2,ys2)), proj4string=CRS(as.character(crs1)))

pdf("6pictures/c12_rpp1-2.pdf", height=2.0, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
    points(spdf_rpp1, pch=16)
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
    points(spdf_rpp2, pch=16)
par(mfcol=c(1,1))
dev.off() 

xs3 <- c(); ys3 <- c()
rel_srtm <- sgdf_srtm
rel_srtm@data[is.na(rel_srtm@data)] <- 0
rel_srtm@data[,1] <- rel_srtm@data[,1] / max(rel_srtm@data[,1])
rel_srtm <- raster(rel_srtm)
i <- 0
set.seed(555) 
while (i < 101) {
    xs <-  runif(1, x1, x2)   
    ys <-  runif(1, y1, y2) 
    rn <-  runif(1,  0,  1) 
    rp_alt <- extract(rel_srtm, SpatialPoints(cbind(xs,ys)), proj4string=CRS(as.character(crs1)))
    if ((rp_alt + rn) > 1) {
        i <- i + 1
        xs3[i] <- xs
        ys3[i] <- ys
        }
    }
spdf_rpp3 <- SpatialPointsDataFrame(cbind(xs3,ys3), as.data.frame(cbind(xs3,ys3)), proj4string=CRS(as.character(crs1)))
image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
points(spdf_rpp3, pch=16)


pdf("6pictures/c12_dweight.pdf", height=3, width=6, bg = "white") 
plot(7*seq(0,5,0.01),df(seq(0,5,0.01), 20, 20),type="l", xlab="dist (km)", ylab="v")
dev.off() 

set.seed(333) 
xs4 <- c( runif(1, x1, x2)); ys4 <- c(runif(1, y1, y2))
edist <- function(a){sqrt((a[1] - a[3])^2 + (a[2] - a[4])^ 2)}
i <- 0
set.seed(555) 
while (i < 10) {
    xs <-  runif(1, x1, x2)   
    ys <-  runif(1, y1, y2)
    rn <-  runif(1,  0,  1) 
    rdf <- data.frame(xs4,ys4,xs,ys)
    dmin <- min(edist(rdf))
    dval <- df(dmin/7000, 20, 20)
    if ((dval + rn) > 1.4) {
        i <- i + 1
        xs4[i] <- xs
        ys4[i] <- ys
    }
}
spdf_rpp4 <- SpatialPointsDataFrame(cbind(xs4,ys4), as.data.frame(cbind(xs4,ys4)), proj4string=CRS(as.character(crs1)))
image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
points(spdf_rpp4, pch=16)

pdf("6pictures/c12_rpp3-4.pdf", height=2.0, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
    points(spdf_rpp3, pch=16)
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
    points(spdf_rpp4, pch=16)
par(mfcol=c(1,1))
dev.off() 

library(spatstat)
ppspec <- list(cif="strauss",par=list(beta=2,gamma=0.2,r=0.7), w=c(bb[1,1],bb[1,2],bb[2,1],bb[2,2]))
ppsim <- rmh(model=ppspec,start=list(n.start=200), control=list(nrep=10,nverb=5))
pdf("6pictures/c12_ppsim.pdf", height=4, width=6, bg = "white") 
    par(mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4)) 
    points(ppsim, pch=16)
dev.off() 

# 4.  Cellular automata ===========================
library(sp)
gt = GridTopology(cellcentre.offset=c(x1,y1),cellsize=c(2500,2500),cells=c(100, 64))
gt = SpatialGrid(gt, proj4string=CRS(as.character(crs1)))
df <- data.frame(coordinates(gt))
df[,1] <- 0
cgrid1 <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))
cgrid1@data[4555,1] <- 1
cgrid0 <- cgrid1
cgrid2 <- cgrid1
for (j in 1: 21){
    for (i in (1:length(coordinates(cgrid1)[,1]))) {
        gx1 <- coordinates(cgrid1)[i,1]
        gy1 <- coordinates(cgrid1)[i,2]
        rdf <- data.frame(coordinates(cgrid1)[,1],coordinates(cgrid1)[,2],gx1,gy1)
        ddist <- edist(rdf)
        nind  <- which(ddist < 3600)
        cgrid2@data[i,1] <- max(cgrid1@data[nind,1])
        }
    cgrid1 <- cgrid2
    if (j==5)  {cgrid5   <- cgrid2}
    if (j==10) {cgrid10  <- cgrid2}
    if (j==20) {cgrid20  <- cgrid2}
}

pdf("6pictures/c12_cellauto1.pdf", height=4.0, width=6, bg = "white") 
par(mfrow=c(2,2), mai = c(0, 0, 0, 0))
    image(cgrid0,  col=c("gray90", "gray45"))
    image(cgrid5,  col=c("gray90", "gray45"))
    image(cgrid10, col=c("gray90", "gray45"))
    image(cgrid20, col=c("gray90", "gray45"))
par(mfcol=c(1,1))
dev.off() 


gt = GridTopology(cellcentre.offset=c(x1,y1),cellsize=c(2500,2500),cells=c(100, 64))
gt = SpatialGrid(gt, proj4string=CRS(as.character(crs1)))
df <- data.frame(coordinates(gt))
df[,1] <- 0
cgrid1 <- SpatialGridDataFrame(gt, df, proj4string = CRS(as.character(crs1)))
cgridx <- cgrid1
cgrid1@data[4555,1] <- 1
cgridx@data[3550:3565,1] <- 2
image(cgrid1, col=c("grey", "black"))
cgrid0 <- cgrid1
for (j in 1: 21){
    for (i in (1:length(coordinates(cgrid1)[,1]))) {
        gx1 <- coordinates(cgrid1)[i,1]
        gy1 <- coordinates(cgrid1)[i,2]
        rdf <- data.frame(coordinates(cgrid1)[,1],coordinates(cgrid1)[,2],gx1,gy1)
        ddist <- edist(rdf)
        nind  <- which(ddist < 3600)
        cgrid2@data[i,1] <- max(cgrid1@data[nind,1])
    }
    cgrid2@data[3550:3565,1] <- 0
    cgrid1 <- cgrid2
    if (j==10) {cgrid10  <- cgrid2}
    if (j==15) {cgrid15  <- cgrid2}
    if (j==20) {cgrid20  <- cgrid2}
}
cgrid0@data[,1]  <- cgrid0@data[,1]  + cgridx@data[,1]
cgrid10@data[,1] <- cgrid10@data[,1] + cgridx@data[,1]
cgrid15@data[,1] <- cgrid15@data[,1] + cgridx@data[,1]
cgrid20@data[,1] <- cgrid20@data[,1] + cgridx@data[,1]

pdf("6pictures/c12_cellauto2.pdf", height=4.0, width=6, bg = "white") 
par(mfrow=c(2,2), mai = c(0, 0, 0, 0))
    image(cgrid0,  col= c("gray90", "gray45", "gray0"))
    image(cgrid10, col= c("gray90", "gray45", "gray0"))
    image(cgrid15, col= c("gray90", "gray45", "gray0"))
    image(cgrid20, col= c("gray90", "gray45", "gray0"))
par(mfcol=c(1,1))
dev.off() 


# 5.  ABM ===========================
n <- 100   # number of itterations
gt = GridTopology(cellcentre.offset=c(x1,y1),cellsize=c(500,500),cells=c(56, 36))
gt = SpatialGrid(gt, proj4string=CRS(as.character(crs1)))
abmpt <- SpatialPoints(gt)
proj4string(abmpt)  <- CRS(as.character(crs1))
cn <- length(abmpt@coords[,1])
abm_alt <- SpatialGridDataFrame(gt, data.frame(extract(rel_srtm, abmpt)), proj4string = CRS(as.character(crs1)))
abm_path <- SpatialGridDataFrame(gt, over(abmpt, as(path, 'SpatialGridDataFrame')), proj4string = CRS(as.character(crs1)))
abm_path@data[,1] <- abm_path@data[,1] / max(abm_path@data[,1])
abm_soil <- abm_path
abm_soil@data[,1] <- 1
abm <- data.frame(id=1:70, type=0, gc1=0, gc2=0, dist=0, dens1=0, dens2=0, densa1=0, densa2=0, alt1=0, alt2=0, soil1=0, soil2=0, path1=0, path2=0, eval1=0, eval2=0)
abm[,2] <- c(rep("a",30),rep("b",30),rep("c",10))
abm[,3] <- sample(1:cn, 70)     
# abm[,4] <- sample(1:cn, 70)  
# abm[,5] <-  apply(cbind(coordinates(abm_alt)[abm[,3],], coordinates(abm_alt)[abm[,4],]),1,  edist2) / 34000
# abm[,10] <- abm_alt@data[abm[,3],1]
# abm[,11] <- abm_alt@data[abm[,4],1]
# abm[,12] <- abm_soil@data[abm[,3],1]
# abm[,13] <- abm_soil@data[abm[,4],1]
# abm[,14] <- abm_path@data[abm[,3],1]
# abm[,15] <- abm_path@data[abm[,4],1]

#functions
edist2  <- function(x){sqrt((x[1] - x[3])^2 + ((x[2] - x[4])^2))}  # euklidische Distanz (x1,y1,x2,y2
gau1   <- function(x, sd){dnorm(edist2(x), mean=0, sd=sd)}       # euklidische Distanzen gewichtet mit normalverteilung des statischen kernels
kde <- function(pp,grid,sd1){    # pp: zwei spalten cbind(x,y); grid: spgdf; sd1: stdev für kde; f=skalierungsfaktor
    grid2 <- grid
    for (i in seq(along=grid2@data[,1])){       
        pdist <- cbind(coordinates(grid2)[i,1],coordinates(grid2)[i,2],pp[,1],pp[,2])### x,y,x1,y1
        grid2@data[i,1]<- sum(apply(pdist,1,gau1,sd=sd1))}
    grid2@data[,1] <- grid2@data[,1] / max(grid2@data[,1])
    return(grid2)}
fline <- function(x,x1,x2,y1,y2) {
    m <- (y2-y1) / (x2-x1)
    b <- (y2 - x2*(y2-y1) / (x2-x1))
    y <- m*x + b
    return(y)}
fdist_a <- function(x){
    if (x<0.2)  {y <- fline(x,0,0.2,0,1)}
    else        {y <- fline(x,0.2,1,1,0)}
    return(y)}
fdist_b <- function(x){
    if (x<0.5)  {y <- fline(x,0,0.5,0,1)}
    else        {y <- fline(x,0.5,1,1,0)}
    return(y)}
fdist_c <- function(x){
    if (x<0.9)  {y <- fline(x,0,0.9,0,1)}
    else        {y <- fline(x,0.9,1,1,0)}
    return(y)}
falt_a <- function(x){
    if (x<0.7)  {y <- fline(x,0,0.7,0.2,1)}
    else        {y <- fline(x,0.7,1,1,0.2)}
    return(y)}
falt_b <- function(x){
    if (x<0.3)  {y <- fline(x,0,0.3,0.2,1)}
    else        {y <- fline(x,0.3,1,1,0.2)}
    return(y)}
falt_c <- function(x){y=0.5
    return(y)}
fsoil_a <- function(x){
    y <- fline(x,0,1,0,1)
    return(y)}
fsoil_b <- function(x){y=0.5
    return(y)}
fsoil_c <- function(x){y=0.5
    return(y)}
fdens_a <- function(x){
    if (x<0.3)  {y <- fline(x,0,0.3,0,1)}
    else        {y <- fline(x,0.3,1,1,0)}
    return(y)}
fdens <- function(x){
    if (x<0.7)  {y <- fline(x,0,0.7,0,1)}
    else        {y <- fline(x,0.7,1,1,0)}
    return(y)}
fpath_c <- function(x){
    y <- fline(x,0,1,1,0)
    return(y)}
fpath <- function(x){
    y <- fline(x,0,1,1,0.5)
    return(y)}

# map of starting distribution
sp_a1 <- SpatialPoints(coordinates(abm_alt)[abm[which(abm[,2]=="a"),3],], proj4string = CRS(as.character(crs1)))
sp_b1 <- SpatialPoints(coordinates(abm_alt)[abm[which(abm[,2]=="b"),3],], proj4string = CRS(as.character(crs1)))
sp_c1 <- SpatialPoints(coordinates(abm_alt)[abm[which(abm[,2]=="c"),3],], proj4string = CRS(as.character(crs1)))

#loop
for (i in 1:n) {
    abm[,4] <- sample(1:cn, 70)  
    while (length(intersect(abm[,3],abm[,4])) > 0) {
        abm[,4] <- sample(1:cn, 70)  }
#density    
    kde_all <- kde(coordinates(abm_alt)[abm[,3],], abm_alt, 2000)
    kde_a <- kde(coordinates(abm_alt)[abm[which(abm[,2]=="a"),3],], abm_alt, 2000)
    kde_all2 <- kde(coordinates(abm_alt)[abm[,4],], abm_alt, 2000)
    kde_a2 <- kde(coordinates(abm_alt)[abm[which(abm[,2]=="a"),4],], abm_alt, 2000)
#sample
    abm[,5] <-  apply(cbind(coordinates(abm_alt)[abm[,3],], coordinates(abm_alt)[abm[,4],]),1,  edist2) / 34000
    abm[,6] <- kde_all@data[abm[,3],1]
    abm[,7] <- kde_all2@data[abm[,3],1]
    abm[,8] <- kde_a@data[abm[,3],1]
    abm[,9] <- kde_a2@data[abm[,3],1]
    abm[,10] <- abm_alt@data[abm[,3],1]
    abm[,11] <- abm_alt@data[abm[,4],1]
    abm[,12] <- abm_soil@data[abm[,3],1]
    abm[,13] <- abm_soil@data[abm[,4],1]
    abm[,14] <- abm_path@data[abm[,3],1]
    abm[,15] <- abm_path@data[abm[,4],1]
#weight
    abm[which(abm[,2]=="a"),5]  <- apply(data.frame(abm[which(abm[,2]=="a"),5]),  1, fdist_a)
    abm[which(abm[,2]=="a"),6]  <- apply(data.frame(abm[which(abm[,2]=="a"),6]),  1, fdens)
    abm[which(abm[,2]=="a"),7]  <- apply(data.frame(abm[which(abm[,2]=="a"),7]),  1, fdens)
    abm[which(abm[,2]=="a"),8]  <- apply(data.frame(abm[which(abm[,2]=="a"),8]),  1, fdens_a)
    abm[which(abm[,2]=="a"),9]  <- apply(data.frame(abm[which(abm[,2]=="a"),9]),  1, fdens_a)
    abm[which(abm[,2]=="a"),10] <- apply(data.frame(abm[which(abm[,2]=="a"),10]), 1, falt_a)
    abm[which(abm[,2]=="a"),11] <- apply(data.frame(abm[which(abm[,2]=="a"),11]), 1, falt_a)
    abm[which(abm[,2]=="a"),12] <- apply(data.frame(abm[which(abm[,2]=="a"),12]), 1, fsoil_a)
    abm[which(abm[,2]=="a"),13] <- apply(data.frame(abm[which(abm[,2]=="a"),13]), 1, fsoil_a)
    abm[which(abm[,2]=="a"),14] <- apply(data.frame(abm[which(abm[,2]=="a"),14]), 1, fpath)
    abm[which(abm[,2]=="a"),15] <- apply(data.frame(abm[which(abm[,2]=="a"),15]), 1, fpath)
    abm[which(abm[,2]=="b"),5]  <- apply(data.frame(abm[which(abm[,2]=="b"),5]),  1, fdist_b)
    abm[which(abm[,2]=="b"),6]  <- apply(data.frame(abm[which(abm[,2]=="b"),6]),  1, fdens)
    abm[which(abm[,2]=="b"),7]  <- apply(data.frame(abm[which(abm[,2]=="b"),7]),  1, fdens)
    abm[which(abm[,2]=="b"),8]  <- apply(data.frame(abm[which(abm[,2]=="b"),8]),  1, fdens)
    abm[which(abm[,2]=="b"),9]  <- apply(data.frame(abm[which(abm[,2]=="b"),9]),  1, fdens)
    abm[which(abm[,2]=="b"),10] <- apply(data.frame(abm[which(abm[,2]=="b"),10]), 1, falt_b)
    abm[which(abm[,2]=="b"),11] <- apply(data.frame(abm[which(abm[,2]=="b"),11]), 1, falt_b)
    abm[which(abm[,2]=="b"),12] <- apply(data.frame(abm[which(abm[,2]=="b"),12]), 1, fsoil_b)
    abm[which(abm[,2]=="b"),13] <- apply(data.frame(abm[which(abm[,2]=="b"),13]), 1, fsoil_b)
    abm[which(abm[,2]=="b"),14] <- apply(data.frame(abm[which(abm[,2]=="b"),14]), 1, fpath)
    abm[which(abm[,2]=="b"),15] <- apply(data.frame(abm[which(abm[,2]=="b"),15]), 1, fpath)
    abm[which(abm[,2]=="c"),5]  <- apply(data.frame(abm[which(abm[,2]=="c"),5]),  1, fdist_c)
    abm[which(abm[,2]=="c"),6]  <- apply(data.frame(abm[which(abm[,2]=="c"),6]),  1, fdens)
    abm[which(abm[,2]=="c"),7]  <- apply(data.frame(abm[which(abm[,2]=="c"),7]),  1, fdens)
    abm[which(abm[,2]=="c"),8]  <- apply(data.frame(abm[which(abm[,2]=="c"),8]),  1, fdens)
    abm[which(abm[,2]=="c"),9]  <- apply(data.frame(abm[which(abm[,2]=="c"),9]),  1, fdens)
    abm[which(abm[,2]=="c"),10] <- apply(data.frame(abm[which(abm[,2]=="c"),10]), 1, falt_c)
    abm[which(abm[,2]=="c"),11] <- apply(data.frame(abm[which(abm[,2]=="c"),11]), 1, falt_c)
    abm[which(abm[,2]=="c"),12] <- apply(data.frame(abm[which(abm[,2]=="c"),12]), 1, fsoil_c)
    abm[which(abm[,2]=="c"),13] <- apply(data.frame(abm[which(abm[,2]=="c"),13]), 1, fsoil_c)
    abm[which(abm[,2]=="c"),14] <- apply(data.frame(abm[which(abm[,2]=="c"),14]), 1, fpath_c)
    abm[which(abm[,2]=="c"),15] <- apply(data.frame(abm[which(abm[,2]=="c"),15]), 1, fpath_c)
#valid
    abm[,16] <- apply(data.frame(abm[,6],abm[,10],abm[,12],abm[,14]), 1, sum)
    abm[,17] <- apply(data.frame(abm[,5],abm[,7],abm[,11],abm[,13],abm[,15]), 1, sum)
    abm[which(abm[,2]=="a"),16]  <- abm[which(abm[,2]=="a"),16] + abm[which(abm[,2]=="a"), 8]
    abm[which(abm[,2]=="a"),17]  <- abm[which(abm[,2]=="a"),16] + abm[which(abm[,2]=="a"), 9]
#soil
    abm_soil@data[abm[,3],1] <- abm_soil@data[abm[,3],1] * 0.9
# move
    for (j in 1:length(abm[,3])) {if(abm[j,17]<abm[j,16]) {abm[j,3] <- abm[j,4]}}
}
#maps
sp_a1 <- SpatialPoints(coordinates(abm_alt)[abm[which(abm[,2]=="a"),3],], proj4string = CRS(as.character(crs1)))
sp_b1 <- SpatialPoints(coordinates(abm_alt)[abm[which(abm[,2]=="b"),3],], proj4string = CRS(as.character(crs1)))
sp_c1 <- SpatialPoints(coordinates(abm_alt)[abm[which(abm[,2]=="c"),3],], proj4string = CRS(as.character(crs1)))

pdf("6pictures/c12_abm.pdf", height=2.0, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
    points(sp_a, pch=15, col="black", cex=0.5)
    points(sp_b, pch=19, col="black", cex=0.5)
    points(sp_c, pch=17, col="black", cex=0.5)
    image(sgdf_srtm, col = gray.colors(25, start = 0.97, end = 0.4))   
    points(sp_a1, pch=15, col="black", cex=0.5)
    points(sp_b1, pch=19, col="black", cex=0.5)
    points(sp_c1, pch=17, col="black", cex=0.5)
par(mfcol=c(1,1))
dev.off() 


save.image("4ws/ws12.rws")
