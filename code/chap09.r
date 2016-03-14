############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 9
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 06
## Date of last changes: 12.08.2015
## Data: Megaliths, Bronze Age Barrows and Medieval villages
## Author of data: gshdl
## Purpose: didactic
## Content: 1. network and transportation, 2. supra regional level 
##          3. regional level, 4. local level
##          5. characterising elements and networks
## Description: applies some basic techniques of
##             network analysis
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #

# 0. Preparation ===============================
#  working direktory
wd <- "~/modproj_qaam"
wd <- "/home/fon/daten/analyse/modproj_qaam"
setwd(wd)
load("4ws/ws08.rws")

# 3.  regional level  ===============================

library("spdep")
co <- coordinates(cent)  
coords <- as.matrix(coordinates(cent))
ids <- row.names(as.data.frame(cent))
wts <- fs[,1]; wts[] <- 1
fs_nb_del <- tri2nb(co, row.names=ids)       
del <- nb2lines(fs_nb_del, wts=wts, coords=coords,  proj4string =  CRS(as.character(crs1)))
fs_nb_soi <- graph2nb(soi.graph(fs_nb_del, co),   row.names=ids)
soi <- nb2lines(fs_nb_soi, wts=wts, coords=coords,  proj4string =  CRS(as.character(crs1)))
fs_nb_gabriel <- graph2nb(gabrielneigh(co),  row.names=ids)  
gabriel <- nb2lines(fs_nb_gabriel, wts=wts,   coords=coords, proj4string =  CRS(as.character(crs1)))
fs_nb_relative <- graph2nb(relativeneigh(co),  row.names=ids) 
relative <- nb2lines(fs_nb_relative, wts=wts,  coords=coords, proj4string =  CRS(as.character(crs1)))

pdf("6pictures/c9_ngraph.pdf", height=6, width=6, bg = "white") 
    par(mfrow=c(2,2), mai = c(0, 0, 1, 0))
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
    lines(del)
    title(main="Delaunay", col.main="black")
    image(sgdf_srtm, col = gray.colors(25, start =  0.97, end = 0.4))  
    lines(soi)
    title(main="SOI", col.main="black")
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
    lines(gabriel)
    title(main="Gabriel", col.main="black")
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
    lines(relative)
    title(main="Relative Neighbour", col.main="black")
    par(mfrow=c(1,1))
dev.off() 


# 4.  local level  ===============================
library(rgdal)
library(tripack)
fsd <- tri.mesh(spdf_meg, duplicate = 'remove')   # ermittelt die delauny trinangulation
fsnn <- neighbours(fsd)     # Vektor mit Listen der Nachbarn 
LinesList <- list()         # leere LinesList Liste erstellen
sldf <- c();deldf_i <- c();deldf_x1 <- c();deldf_y1 <- c();deldf_k <- c();deldf_x2 <- 
    c();deldf_y2 <- c();deldf_name <- c()
for(i in seq(along=spdf_meg@coords[,1])) {         
    pid1 <- i                                  
    x1 <- spdf_meg@coords[i,1]
    y1 <- spdf_meg@coords[i,2]
    for(k in seq(along=(fsnn[i][[1]]))) {   
        pid2 <- fsnn[[i]][k]    
        if (pid2 > pid1) {     
            x2 <- spdf_meg@coords[pid2,1]
            y2 <- spdf_meg@coords[pid2,2] 
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
        }
    }
}
deldf <- data.frame(deldf_i,deldf_x1,deldf_y1,deldf_k,deldf_x2,deldf_y2,deldf_name)
sldf2 <- data.frame(sldf)
sl <- SpatialLines(LinesList, proj4string = CRS(as.character(crs1)))
sdf <- SpatialLinesDataFrame(sl, sldf2, match.ID = FALSE)
writeOGR(sdf, "./5result", "delaunay_meg", driver="ESRI Shapefile", overwrite_layer=TRUE)
write(rbind(deldf_i,deldf_x1,deldf_y1,deldf_k,deldf_x2,deldf_y2,deldf_name), file 
      = "./5result/delaunay_meg.csv", ncolumns = 7, sep = ";")

coords <- as.matrix(coordinates(cent))
ids <- row.names(as.data.frame(1:length(spdf_meg@coords[,1])))
wts <- spdf_meg@coords[,1]; wts[] <- 1
meg_nb_del <- tri2nb(spdf_meg@coords, row.names=ids)       
meg_del <- nb2lines(meg_nb_del, wts=wts, coords=spdf_meg@coords,  proj4string =  CRS(as.character(crs1)))
meg_nb_soi <- graph2nb(soi.graph(meg_nb_del, spdf_meg@coords),   row.names=ids)
meg_soi <- nb2lines(meg_nb_soi, wts=wts, coords=spdf_meg@coords,  proj4string =  CRS(as.character(crs1)))
meg_nb_gabriel <- graph2nb(gabrielneigh(spdf_meg@coords),  row.names=ids)  
meg_gabriel <- nb2lines(meg_nb_gabriel, wts=wts,   coords=spdf_meg@coords, proj4string =  CRS(as.character(crs1)))
meg_nb_relativ <- graph2nb(relativeneigh(spdf_meg@coords),  row.names=ids) 
meg_relativ <- nb2lines(meg_nb_relativ, wts=wts,  coords=spdf_meg@coords, proj4string =  CRS(as.character(crs1)))

pdf("6pictures/c9_ngraphmeg.pdf", height=6, width=6, bg = "white") 
par(mfrow=c(2,2), mai = c(0, 0, 1, 0))
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
    lines(meg_del)
    title(main="Delaunay", col.main="black")
    image(sgdf_srtm, col = gray.colors(25, start =  0.97, end = 0.4))  
    lines(meg_soi)
    title(main="SOI", col.main="black")
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
    lines(meg_gabriel)
    title(main="Gabriel", col.main="black")
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
    lines(meg_relativ)
    title(main="Relative Neighbour", col.main="black")
par(mfrow=c(1,1))
dev.off() 

nb_k1 <- knn2nb(knearneigh(spdf_meg@coords, k=1),row.names=ids) 
nb_k2 <- knn2nb(knearneigh(spdf_meg@coords, k=2),row.names=ids)   
nb_k3 <- knn2nb(knearneigh(spdf_meg@coords, k=3),row.names=ids)    
nb_k4 <- knn2nb(knearneigh(spdf_meg@coords, k=4),row.names=ids)  

pdf("6pictures/c9_knn.pdf", height=6, width=6, bg = "white") 
par(mfrow=c(2,2), mai = c(0, 0, 1, 0))
image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))
    plot.nb(nb_k1, spdf_meg@coords, pch= 16, col="black", points=TRUE, add=TRUE, arrows=FALSE, cex=0.4)
    title("k = 1")
image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))
    plot.nb(nb_k2, spdf_meg@coords,  pch= 16, col="black", points=TRUE, add=TRUE, arrows=FALSE, cex=0.4)
    title("k = 2")
image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))
    plot.nb(nb_k3, spdf_meg@coords,  pch= 16, col="black", points=TRUE, add=TRUE, arrows=FALSE, cex=0.4)
    title("k = 3")
image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))
    plot.nb(nb_k4, spdf_meg@coords,  pch= 16, col="black", points=TRUE, add=TRUE, arrows=FALSE, cex=0.4)
    title("k = 4")
par(mfrow=c(1,1))
dev.off() 

ras_ridges <- sgdf_tum_dens
ras        <- sgdf_tum_dens
cols <- sgdf_tum_dens@grid@cells.dim[1]
ras_ridges@data$v <- 1        # initial value of 1 for ridges
ras@data$v[is.na(ras@data$v)] <- 10000000000

for (i in 1:(length(ras@data$v)-(2*cols)-2))  { 
    ind <- c(i,i+1,i+2,i+cols,i+cols+1,i+cols+2,i+(2*cols),i+(2*cols)+1,i+(2*cols)+2)
    ind_min1 <- which(ras@data$v[ind]==min(ras@data$v[ind]))
    ind_min2 <-ind[ind_min1]
    ras_ridges@data$v[ind_min2] <- 0
}

pdf("6pictures/c9_dridge_barrows.pdf", height=4, width=6, bg = "white") 
par( mai = c(0, 0, 0, 0))
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
    image(ras_ridges, col = gray.colors(2, start = 0.3, end = 0.7), zlim=c(1,1), add=TRUE)   
dev.off() 

library(raster)
library(gdistance)
ras <- raster(sgdf_srtm)
nz <- 8   # neigbourhood number: 4,8,16
projection(ras) <- crs1
ras <- focal(ras, w=matrix(1/9,nrow=3,ncol=3), NAonly=TRUE)
plot(ras, col = gray.colors(25, start = 0.97, end = 0.4))
starts <- cbind(cent@coords[,1],cent@coords[,2])

# cost functions
# Tobler1993 velocity
tobler1993a <- function(s){6 * exp(-3.5 * abs(s + 0.05))}      # km/h
tobler1993b <- function(s){0.36 * exp(-3.5 * abs(s + 0.05))}   # m/min
# Minetti2002 Metabolische Kosten in J/(kg*m) für laufen
minetti2002w <- function(s){(280.5 * s^5 - 58.7 * s^4 - 76.8 * s^3 + 51.9 * s^2 + 19.6  * s + 2.5)}
minetti2002wi <- function(s){1/(280.5 * s^5 - 58.7 * s^4 - 76.8 * s^3 + 51.9 * s^2 + 19.6 * s + 2.5)}
# Herzog2012 (after data from Minetti) walking. metabolic costs J/(kg*m)
herzog2012_w <- function(s){(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 * s^3  + 93.419 * s^2 + 19.825 * s + 1.64)}
herzog2012_wi <- function(s){1/(1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 *   s^3 + 93.419 * s^2 + 19.825 * s + 1.64)}

pdf("6pictures/c9_cfunc.pdf", height=3, width=6, bg = "white") 
par( mai = c(1, 1, 0.1, 0.1))
    matplot(seq(-0.5,0.5,0.01), herzog2012_wi(seq(-0.5,0.5,0.01)), type="l", xlab=expression(paste("slope (", Delta, "h/", Delta, "d)")), ylab="metabolic costs (kg m/J)")
    lines(seq(-0.5,0.5,0.01), minetti2002wi(seq(-0.5,0.5,0.01)), type="l", lty=2)
dev.off() 

# auxilliary function
hdiff <- function(x){x[2]-x[1]}

# transitional object 
hd <- transition(ras,hdiff,nz,symm=TRUE)    
slope <- geoCorrection(hd,scl=FALSE)             # transitional object with slope
adj <- adjacent(x=ras, cells=1:ncell(ras), direction=nz) 
cost <- slope       
cost[adj] <- herzog2012_wi(slope[adj]) 
conduct <- geoCorrection(cost, scl=FALSE) # conductivity=cost/dist; time=1/conductivity


cost_surface1 <- accCost(conduct, starts[1,])
cost_surface2 <- accCost(conduct, starts[2,])
cost_surface3 <- accCost(conduct, starts[3,])
cost_surface4 <- accCost(conduct, starts[4,])
cost_surface5 <- accCost(conduct, starts[5,])
cost_surface6 <- accCost(conduct, starts[6,])
cost_surface7 <- accCost(conduct, starts[7,])
cost_surface8 <- accCost(conduct, starts[8,])
cost_surface9 <- accCost(conduct, starts[9,])

#  voronoi in economic space
csm1 <- cost_surface1 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
csm2 <- cost_surface2 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,2))
csm2 <- subs(x=csm2, y=df)
csm3 <- cost_surface3 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,3))
csm3 <- subs(csm3, df)
csm4 <- cost_surface4 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,4))
csm4 <- subs(csm4, df)
csm5 <- cost_surface5 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,5))
csm5 <- subs(csm5, df)
csm6 <- cost_surface6 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,6))
csm6 <- subs(csm6, df)
csm7 <- cost_surface7 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,7))
csm7 <- subs(csm7, df)
csm8 <- cost_surface8 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,8))
csm8 <- subs(csm8, df)
csm9 <- cost_surface9 == min(cost_surface1,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,9))
csm9 <- subs(csm9, df)
csm.a <- csm1 + csm2 + csm3 + csm4 + csm5 + csm6 + csm7 + csm8 + csm9

#  weighted voronoi in economic space
cost_surface1b <- cost_surface1 * 2.5
cost_surface6b <- cost_surface6 * 0.7
csm1 <- cost_surface1b == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
csm2 <- cost_surface2 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,2))
csm2 <- subs(x=csm2, y=df)
csm3 <- cost_surface3 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,3))
csm3 <- subs(csm3, df)
csm4 <- cost_surface4 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,4))
csm4 <- subs(csm4, df)
csm5 <- cost_surface5 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,5))
csm5 <- subs(csm5, df)
csm6 <- cost_surface6b == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,6))
csm6 <- subs(csm6, df)
csm7 <- cost_surface7 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,7))
csm7 <- subs(csm7, df)
csm8 <- cost_surface8 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,8))
csm8 <- subs(csm8, df)
csm9 <- cost_surface9 == min(cost_surface1b,cost_surface2,cost_surface3,cost_surface4,cost_surface5,cost_surface6b,cost_surface7,cost_surface8,cost_surface9)
df <- data.frame(id=c(0,1), v=c(0,9))
csm9 <- subs(csm9, df)
csm.b <- csm1 + csm2 + csm3 + csm4 + csm5 + csm6 + csm7 + csm8 + csm9

pdf("6pictures/c09_voroeco.pdf", height=2.0, width=6, bg = "white") 
par(mfcol=c(1,2), mai = c(0, 0, 0, 0))
image(as(csm.a, "SpatialGridDataFrame"), col = gray.colors(25, start = 0.97, end = 0.4))
points(starts, pch=16, cex=0.8) 
image(as(csm.b, "SpatialGridDataFrame"), col = gray.colors(25, start = 0.97, end = 0.4))
points(starts, pch=16, cex=0.8) 
par(mfcol=c(1,1))
dev.off() 


# Delaunay start points
library(deldir)
try <- deldir(starts[,1], starts[,2], plot=TRUE,wl='tr')         
LinesList <- list()      
id <- seq(along=try$delsgs[,1])
for(i in seq(along=try$delsgs[,1])) {    
    m <- matrix(data = c(try$delsgs[i,1],try$delsgs[i,3],try$delsgs[i,2],
                         try$delsgs[i,4]), nrow=2, ncol=2)   
    L <- Line(m)
    LL <- list(L)
    name  <- paste("edge", "_", try$delsgs[i,5],"_", try$delsgs[i,6], sep="")
    LLL <- Lines(LL, ID = name)               
    LinesList[length(LinesList)+1] <- LLL     
}
deldf2 <- data.frame(try$delsgs[,5], try$delsgs[,1], try$delsgs[,2], try$delsgs[,6], 
                    try$delsgs[,3], try$delsgs[,4], paste("edge", "_", try$delsgs[,5],"_"), 
                    try$delsgs[,6], sep="")
cols <- c("a","b","c","d","e","f","g","h","i")
colnames(deldf2) <- cols
sl <- SpatialLines(LinesList, proj4string = CRS(as.character(crs1))
starts2 <- SpatialLinesDataFrame(sl, deldf2, match.ID = FALSE)

# least cost path
for(i in 1:length(starts2[,1])){
    i1 <- i*2-1
    i2 <- i*2
    s <- SpatialPoints(cbind(starts2@data$b[i],starts2@data$c[i]))
    z <- SpatialPoints(cbind(starts2@data$e[i],starts2@data$f[i]))
    sz <- shortestPath(conduct, s, z, output="SpatialLines")
    zs <- shortestPath(conduct, z, s, output="SpatialLines")
    sz@lines[[1]]@ID <- as.character(i1)    
    zs@lines[[1]]@ID <- as.character(i2)           
    if(i==1){sdf <-rbind(sz,zs)}
    if(i>1){sdf <- rbind(sdf,sz,zs, makeUniqueIDs = TRUE)} 
    if(i==1){df <- cbind(c(1,2), c("sz","zs"), c(starts2@data$g[i]))}
    if(i>1){df <- cbind(c(df[,1],i1,i2), c(df[,2],"sz","zs"), c(df[,3],
                                                                starts2@data$g[i],starts2@data$g[i]))}
}
lcp_df <- as.data.frame(df)
lcp_sldf <- SpatialLinesDataFrame(sdf,lcp_df, match.ID = FALSE) 

pdf("6pictures/c9_lcp.pdf", height=3.2, width=6, bg = "white") 
par( mai = c(0, 0, 0, 0.2))
    #image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4))  
plot(raster(sgdf_srtm), col = gray.colors(25, start =   0.97, end = 0.4))
    lines(lcp_sldf)
    points(starts, pch=16, cex=0.8)  
dev.off() 

# least cost with random walk
path <- ras
path[] <- 0
for(i in 1:length(starts2[,1])){
    s <- SpatialPoints(cbind(starts2@data$b[i],starts2@data$c[i]))
    z <- SpatialPoints(cbind(starts2@data$e[i],starts2@data$f[i]))
    path <- max(path, passage(conduct, s, z, theta=0.0001, totalNet="total"))
}

pdf("6pictures/c9_lcp2.pdf", height=3.2, width=6, bg = "white") 
par( mai = c(0, 0, 0, 0.2))
    plot(path, col=gray.colors(80, start=0.9, end=0.1, gamma=0.2))
    points(SpatialPoints(cbind(starts2@data$b,starts2@data$c)), pch=16, cex=0.8)  
dev.off() 

                   
pdf("6pictures/c9_lcp1-2.pdf", height=2, width=6, bg = "white")               
par(mfrow=c(1,2), mai = c(0, 0, 0, 0.0), mar = c(0, 0, 0, 0.0))  
   image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4)) 
   #plot(raster(sgdf_srtm), col = gray.colors(25, start =   0.97, end = 0.4))
   lines(lcp_sldf)
   points(starts, pch=16, cex=0.8)     
   image(as(path, 'SpatialGridDataFrame') , col = gray.colors(80, start=0.9, end=0.1, gamma=0.2)) 
   #plot(path, col=gray.colors(80, start=0.9, end=0.1, gamma=0.2))
   points(SpatialPoints(cbind(starts2@data$b,starts2@data$c)), pch=16, cex=0.8)  
par(mfrow=c(1,1), mai = c(0, 0, 0, 0.0))           
dev.off() 
                   
                   
                   
                   
                   
                   
# the same with megalith density values
ras2 <- raster(sgdf_meg_dens)
projection(ras2) <- crs1
ras2 <- focal(ras2, w=matrix(1/9,nrow=3,ncol=3), NAonly=TRUE)
hd2 <- transition(ras2,mean,nz,symm=TRUE)    
slope2 <- geoCorrection(hd2,scl=FALSE)             # transitional object with slope
adj2  <- adjacent(x=ras2, cells=1:ncell(ras2), direction=nz) 
cost2 <- slope2       
cost2[adj2] <- herzog2012_wi(slope2[adj2]) 
conduct2 <- geoCorrection(cost2, scl=FALSE) # conductivity=cost/dist; time=1/conductivity
path2 <- ras2
path2[] <- 0
for(i in 1:length(starts2[,1])){
    s <- SpatialPoints(cbind(starts2@data$b[i],starts2@data$c[i]))
    z <- SpatialPoints(cbind(starts2@data$e[i],starts2@data$f[i]))
    path2 <- max(path2, passage(conduct2, s, z, theta=0.0001, totalNet="total"))
}
plot(path2)

# sampling
pathn_topo <- extract(path, spdf_meg)
pathn_tumdens <- extract(path2, spdf_meg)
mean(pathn_topo)
mean(pathn_tumdens)

# 5.  characterising  ===============================
library(igraph)
library(spdep)
#co <- coordinates(spdf_meg)  
#coords=as.matrix(coordinates(spdf_meg))
ids <- row.names(as.data.frame(spdf_meg))
meg_nb_del <- tri2nb(coordinates(spdf_meg), row.names=ids) 
m <- nb2mat(meg_nb_del)             
g <- graph.adjacency(m, mode="lower", weighted=T)
g <- set.vertex.attribute(g, "x", index=V(g), coordinates(spdf_meg)[,1])   
g <- set.vertex.attribute(g, "y", index=V(g), coordinates(spdf_meg)[,2]) 

dpd        <- deldf
dist  <- sqrt((dpd[2] - dpd[5])^2 + (dpd[3] - dpd[6])^2 )
dist2 <- dist^2
dimnames(dist) <- list(NULL,"dist")
dimnames(dist2) <- list(NULL,"dist2")
dpd$dist   <- dist
dpd$dist2  <- dist2
g          <- set.edge.attribute(g, "distance2", index=E(g), dist2)
E(g)$weight <- dist

c.degree <- degree(g, v=V(g))
c.closness <- closeness(g, v=V(g))
c.betweenness <- betweenness(g, v=V(g))
c.bonacich.power <- bonpow(g, nodes=V(g))

ctab <- data.frame(cbind(id = V(g)+1, x = fs[,1], y = fs[,2], degree=c.degree, closness=c.closness, betweenness=c.betweenness, bpower=c.bonacich.power))
coordinates(ctab)=~x+y  
proj4string(ctab)  <- CRS(as.character(crs1)) 

pdf("6pictures/c9_cent.pdf", height=4, width=6, bg = "white") 
par( mai = c(0, 0, 0, 0))
    #plot(raster(sgdf_srtm), col = gray.colors(25, start = 0.97, end = 0.4))   
    image(sgdf_srtm, col = gray.colors(25, start =   0.97, end = 0.4)) 
    points(ctab, pch=16, cex=sqrt(ctab$betweenness)/50)
dev.off() 



# ?????????????????????????????????????????????????
spoints <- data.frame(cbind(x=(seq(1:44)*500) + 3552000, y=rep(6031270, 44)))
coordinates(spoints)=~x+y
proj4string(spoints)  <- CRS(as.character(crs1)) 

i_kde <- extract(raster(sgdf_meg_dens), spoints)
#i_kde <- overlay(x=sgdf_meg_dens, y=spoints)
#i_kde_tum <- overlay(x=sgdf_tum_dens, y=spoints)
i_kde_tum <- extract(raster(sgdf_tum_dens), spoints)
#names(i_kde)[names(i_kde) == 'v'] <- 'meg'
#i_kde   <- cbind(i_kde@data,i_kde_tum@data$v)
i_kde   <- cbind(i_kde,i_kde_tum)
#names(i_kde)[names(i_kde) == 'i_kde_tum@data$v'] <- 'tum'

cdistance <- sqrt((i_kde[1,1] - i_kde[,1])^2 + (i_kde[,2] - i_kde[,2])^2)

plot(cdistance, col="black", pch=16)
lines(cdistance,lty=2,col="black")

save.image("4ws/ws09.rws")
