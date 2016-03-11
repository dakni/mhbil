############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 3
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 06
## Date of last changes: 03.08.2015
## Data: villages
## Author of data: gshdl
## Purpose: didactic
## Content: 1. one dimension, 2. two dimensions
## Description: applies some basic techniques of
##              regression and interpolation
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #

# 0. Preparation ===============================
#  working direktory
wd <- "~/modproj_qaam"
wd <- "/home/fon/daten/analyse/modproj_qaam"
setwd(wd)

# 2.3. R-Script ===============================

# the lines with basic function will not appear in this script 

file_meg    <- "1data/meg_dw.csv" 
file_tum    <- "1data/tum_dw.csv" 
file_coast1 <- "2geodata/coast_gk3.shp" 
file_coast  <- "coast_gk3"  
file_srtm   <- "2geodata/dw_gk3_50_ag.asc"
file_vil    <- "1data/villages.xls" 

crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs"
crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

install.packages('sp') 
library(sp)      
library(proj4)           
library(rgdal)         
library(spatstat)     
library(RSQLite)    
library(gdata) 

df_meg <- read.table(file_meg, sep=';', header=TRUE)
df_tum <- read.table(file_tum, sep=';', header=TRUE)

spdf_meg <- read.table(file_meg, sep=';', header=TRUE)
spdf_tum <- read.table(file_tum, sep=';', header=TRUE)
coordinates(spdf_meg)=~x+y  
coordinates(spdf_tum)=~x+y  

avs <- paste(wd,"/2geodata",sep="")
coast <- readOGR(avs, p4s=NULL, file_coast) 

df_vil_wgs84 <- read.xls(file_vil, 1)
spdf_vil_wgs84 <- df_vil_wgs84
coordinates(spdf_vil_wgs84)=~x+y  

sgdf_srtm <- read.asciigrid(file_srtm) 

spdf_vil_wgs84@coords[,1]        
spdf_vil_wgs84@coords[,"x"]  
spdf_vil_wgs84@coords[,2]  

df_vil_coord <- project(cbind(spdf_vil_wgs84@coords[,1], spdf_vil_wgs84@coords[,2]), crs1)
df_vil_k <- cbind(x=df_vil_coord[,1],  y=df_vil_coord[,2])
df_vil <- data.frame(id=df_vil_wgs84[,1], village=as.character(df_vil_wgs84[,2]),  AD=df_vil_wgs84[,3])
spdf_vil <- SpatialPointsDataFrame(df_vil_k, as.data.frame(df_vil), proj4string=CRS(as.character(crs1)))

proj4string(spdf_meg)  <- CRS(as.character(crs1)) 
proj4string(spdf_tum)  <- CRS(as.character(crs1)) 
proj4string(sgdf_srtm) <- CRS(as.character(crs1)) 

bb = bbox(sgdf_srtm)    
win <- owin(xrange=c(bb[1,1],bb[1,2]), yrange=c( bb[2,1],bb[2,2]), unitname="m")
spdf_meg <-remove.duplicates(spdf_meg, zero=0, remove.second=TRUE)
spdf_tum <-remove.duplicates(spdf_tum, zero=0, remove.second=TRUE)
spdf_vil <-remove.duplicates(spdf_vil, zero=0,remove.second=TRUE)
ppp_meg <- ppp(spdf_meg@coords[,1], spdf_meg@coords[,2], window=win)
ppp_tum <- ppp(spdf_tum@coords[,1], spdf_tum@coords[,2], window=win)
ppp_vil <- ppp(spdf_vil@coords[,1], spdf_vil@coords[,2], window=win)




library(raster)
srtm_slope  <- terrain(raster(sgdf_srtm), opt='slope')
srtm_aspect <- terrain(raster(sgdf_srtm), opt='aspect')
srtm_shade  <- hillShade(srtm_slope, srtm_aspect, 40, 270)
top.colors = colorRampPalette(c("#618CB5", "#23B0EE", "#81BB7C", "#E5CE98", "#B89E83", "#9E5D4C"), bias=1.2) 


# pdf("6pictures/c03_map.pdf", height=3.28, width=6, bg = "white") 
#     par(mai = c(0, 0, 0, 0.2), mar = c(0, 0, 0, 0.2))
#     plot(raster(sgdf_srtm), col = top.colors(25))
#     plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
#     contour(raster(sgdf_srtm), levels=c(0.2), labcex=0.001, cex=1, add=TRUE)
#     points(spdf_tum,pch = 19, cex = 0.6, col = "#D0043A")
#     points(spdf_meg,pch = 17, cex = 0.8, col = "#074A9D")
#     points(spdf_vil,pch = 15, cex = 1, col = "#D08204")
#     legend("bottomright", cex=0.7, legend=c("Bronze Age Barrows","Neolithic Megaliths","Mediaeval Villages"),pch=c(19,17,15), col=c("#D0043A","#074A9D","#D08204"))
#     scalebar(d = 5000, cex=0.7, divs = 2, below="m", type = "bar", xy=c(3573000, 6026500), label= c("0","2.5","5"), adj = c(.5,-1.3))
#     text(3581800, 6025000, "Altitude (m)", cex=0.8)
# dev.off() 


png("6pictures/c03_map.png", height=3.28, width=6, units="in", res=600, bg = "white") 
    par(mai = c(0, 0, 0, 0.2), mar = c(0, 0, 0, 0.2))
    plot(raster(sgdf_srtm), col = top.colors(25))
    plot(srtm_shade, col=grey(seq(from=0,to=1,by=0.02), alpha=0.60), legend=FALSE, add=TRUE,  cex=0.8)
    contour(raster(sgdf_srtm), levels=c(0.2), labcex=0.001, cex=1, add=TRUE)
    points(spdf_tum,pch = 19, cex = 0.4, col = "#D0043A")
    points(spdf_meg,pch = 17, cex = 0.5, col = "#074A9D")
    points(spdf_vil,pch = 15, cex = 0.8, col = "black")
    legend("bottomright", cex=0.7, legend=c("Bronze Age Barrows","Neolithic Megaliths","Mediaeval Villages"),pch=c(19,17,15), col=c("#D0043A","#074A9D","black"))
    scalebar(d = 5000, cex=0.7, divs = 2, below="km", type = "bar", xy=c(3573500, 6026500), label= c("0","2.5","5"), adj = c(.5,-1.3))
    text(3581600, 6025000, "Altitude (m)", cex=0.8)
    text(3565000, 6039600, "Baltic Sea", cex=0.8,  font=3)
    text(3562000, 6030000, "Dänischer Wohld", cex=0.8,  font=3)
dev.off() 



library(OpenStreetMap)
openmap(bbox(sgdf_srtm)[2,2],bbox(sgdf_srtm)[1,1],9,'maprequest')



# 2.4. Helpful functions in R ===============================

#loops
a <- 5; b <- 0; c <- 100; i=0
while (i < c) {
    b <- b + 1
    i <- i + a
    }
b    


a <- 6; b <- 0; c <- 100; i <- 0
repeat{
    b <- b + 1
    i <- i + a
    if (i==c | i>(c-a)) {break}
    }
b

library(foreach)
foreach(a=c(0.05,0.24,0.3), b=1:3) %do% (sin(a) / b)

foreach(a=c(0.05,0.24,0.3), b=1:3, .combine='c') %do% (sin(a) / b)

library(doMC)
registerDoMC(cores=4)
foreach(a=c(0.05,0.24,0.3), .combine='c') %dopar% (sin(a))

a1 <- c(2,5,4,7,5,9,7)
a2 <- c(1,2,3,4,5,6,7)
a3 <- c(8,4,6,2,3,7,8)
a  <- cbind(a1,a2,a3)
apply(a,1,sum)

#tidy
id <- c(1,2,3,4,5,6)
diameter <- c(3,6,4,4,2,9)
lengtth <- c(23,32,12,22,16,77)
colour <- c("red","red","blue","red","blue","green")
finds <- rbind(id,diameter,lengtth,colour)
finds

finds <- t(finds)
finds

colnames(finds)[3]  <- "length"

finds <- data.frame(id=as.numeric(finds[,1]), diameter=as.numeric(finds[,2]), length=as.numeric(finds[,3]), colour=factor(finds[,4]))
str(finds)


library(reshape)
finds_m <- melt(finds, id.vars="id")
finds_m



finds_m_f <- subset(finds_m, variable=="colour")
finds_m_n <- subset(finds_m, variable=="diameter" | variable=="length")

str(finds_m_n)

finds_m_n$value <- as.numeric(finds_m_n$value)

finds_w_n <-cast(finds_m_n)
finds_w_n
finds_w_f <-cast(finds_m_f)

finds_w <- merge(finds_w_n, finds_w_f)
finds_w

#sqldf
library(sqldf)
sqldf('select diameter, length from finds_w where colour is "red"')

sqldf('select id, diameter, length, colour from finds_w_n natural join finds_w_f where colour is "blue"')


library(data.table)

dt_finds_w_n <- data.table(finds_w_n)
setkey(dt_finds_w_n, id)
dt_finds_w_f <- data.table(finds_w_f)
setkey(dt_finds_w_f, id)




save.image("4ws/ws03.rws")







