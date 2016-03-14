############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 6
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 08
## Date of last changes: Fr 21. Aug 15:29:08 CEST 2015
## 
## Data:
## - Megaliths, Bronze Age Barrows and Medieval villages
## - SRTM
## Author of data: gshdl, NASA?!
##
## Purpose: didactic
## Content: 
## Description: Locational Analysis
##
## Licence data: -
## Licence Script: GPL 
##    (http://www.gnu.org/licenses/gpl-3.0.html)
############################################## #


## R code for Chapter 6 "Location and Characterisation"
## ==================================================


## Prerequisites
## ==================================================

setwd("PATH_TO_YOUR_WORKING_DIRECTORY")

## define variables
## --------------------------------------------------
## We presume that you store the data according to
## their type, i.e. tables are in a "1data" folder
## and geodata are in a "2geodata" folder
## later we will store results in a folder "3results"
file_meg <- "1data/meg.dw2.csv"
file_tum <- "1data/tum.dw2.csv"
file_vil <- "1data/villages.xls"
sgdf_srtm <- "2geodata/dw_gk3_50_ag.asc"

crs1 <- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1+x_0=3500000 +y_0=0 +ellps=WGS84 +units=m +no_defs"
crs2 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

## load library
library(sp)
library(gdata)
library(rgdal)
library(proj4)
library(raster)
#library(rasterVis)
library(dplyr) # masks a lot of other functions...make sure to select the right function; if unsure use library::function
library(dismo)

## load data
df_meg <- read.table(file_meg, sep=";", header=TRUE)
df_tum <- read.table(file_tum, sep=";", header=TRUE)

spdf_meg <- read.table(file_meg, sep=";", header=TRUE)
spdf_tum <- read.table(file_tum, sep=";", header=TRUE)

coordinates(spdf_meg)= ~x+y # requires sp
coordinates(spdf_tum)= ~x+y

df_vil_wgs84 <- read.xls(file_vil, 1) # requires gdata
spdf_vil_wgs84 <- df_vil_wgs84
coordinates(spdf_vil_wgs84)= ~x+y

## load SRTM
sgdf_srtm <- readGDAL(sgdf_srtm) # requires rgdal
sgdf_srtm@proj4string@projargs <- crs1
names(sgdf_srtm@data) <- "srtm" # change the name from "band1" to "srtm"

## assing coordinate reference systems (CRS)
proj4string(spdf_meg) <- CRS(as.character(crs1))
proj4string(spdf_tum) <- CRS(as.character(crs1))
proj4string(sgdf_srtm) <- CRS(as.character(crs1))

## reproject villages.xls --> requires proj4string
df_vil_coord <- proj4::project(cbind(spdf_vil_wgs84@coords[,1],
                                     spdf_vil_wgs84@coords[,2]),
                               crs1)

df_vil_k <- cbind(x=df_vil_coord[,1]+3500000,
                  y=df_vil_coord[,2]) # +3500000 due to GK zoning

df_vil <- data.frame(id=df_vil_wgs84[,1],
                     village=as.character(df_vil_wgs84[,2]),
                     AD=df_vil_wgs84[,3])

spdf_vil <- SpatialPointsDataFrame(df_vil_k,
                                   as.data.frame(df_vil),
                                   proj4string=CRS(as.character(crs1)))



## make sure that all variables have the same CRS
spdf_vil@proj4string@projargs == spdf_meg@proj4string@projargs &&
    spdf_meg@proj4string@projargs == spdf_tum@proj4string@projargs &&
        spdf_meg@proj4string@projargs == sgdf_srtm@proj4string@projargs


## compare coordinates
head(spdf_tum@coords)
head(spdf_meg@coords)
head(spdf_vil@coords)
sgdf_srtm@bbox

## From this point on we are back in the book
## ==================================================
sgdf_srtm <- readGDAL(sgdf_srtm) # requires rgdal
names(sgdf_srtm@data) <- "srtm" # change the name from "band1" to "srtm"
is.projected(sgdf_srtm)
sgdf_srtm@proj4string@projargs <- crs1
is.projected(sgdf_srtm)
spdf_meg@proj4string@projargs == sgdf_srtm@proj4string@projargs

## create raster object from load SRTM
srtm <- raster(sgdf_srtm)
srtm
summary(srtm)

pdf("./3results/ch6_srtm.pdf", height=8, width=12, bg = "white")
plot(srtm,
     col = gray.colors(n = 25, start = 0.1, end = .9),
     legend.lab="Altitude (m)",
     horizontal = FALSE,
     legend.width = 1,
     cex.axis=.9,
     tcl=-.3,
     mgp = c(3,.2,0)
     )
points(spdf_tum,pch = 19, cex = .8, col = "black",bg="black")
points(spdf_meg,pch = 22, cex = 1.5, col = "black", bg="white")
points(spdf_vil,pch = 23, cex = 2.5, col = "black", bg = "white")
legend("bottomright",legend=c("Bronze Age Barrows","Megaliths","Villages"),pch=c(19,22,23))
scalebar(d = 5000, divs = 2, below="Meter",type = "bar",xy=c(3573800,6025500),adj = c(.5,-1.3))
dev.off()

## extract characteristics from locations
spdf_meg@data$srtm <- extract(x = srtm, y = spdf_meg, buffer = 200, fun = median)
spdf_tum@data$srtm <- extract(x = srtm, y = spdf_tum, buffer = 200, fun = median)
spdf_vil@data$srtm <- extract(x = srtm, y = spdf_vil, buffer = 200, fun = median)

## t.test
t.test(x=spdf_meg@data$srtm, y=spdf_tum@data$srtm)
t.test(x=spdf_tum@data$srtm, y=spdf_vil@data$srtm)
t.test(x=spdf_meg@data$srtm, y=spdf_vil@data$srtm)


library("KernSmooth")
ks_meg <- bkde(spdf_meg@data$srtm,
               kernel="normal",
               bandwidth=3,
               gridsize=201,
               range.x = c(summary(srtm)[1],summary(srtm)[5]))

ks_tum <- bkde(spdf_tum@data$srtm,
               kernel="normal",
               bandwidth=3,
               gridsize=201,
               range.x = c(summary(srtm)[1],summary(srtm)[5]))

ks_vil2 <- bkde(spdf_vil@data$srtm,
                kernel="normal",
                bandwidth=3,
                gridsize=201,
                range.x = c(summary(srtm)[1],summary(srtm)[5]))

srtm_char <- data.frame(ks_meg,ks_tum$y,ks_vil2$y)
colnames(srtm_char) <- c("altitude","meg","tum","vil")

library(reshape)
srtm_char2<- melt(srtm_char, id.vars = "altitude")
head(srtm_char2,3)
tail(srtm_char2,3)

library(ggplot2)

pdf("./3results/ch6_srtm_char.pdf", height=4, width=6, bg = "white")
srtm_char_plot <- ggplot(srtm_char2, aes(x=altitude, y=value)) +
    geom_line(aes(linetype = variable)) +
    labs(x="Altitude (m)",y="density",legend="") +
    theme_bw(base_size = 12) +
    theme(legend.position="bottom",
          legend.title=element_blank()
          ) +
        scale_linetype_discrete(labels = c("Megaliths","Bronze Age Barrows", "Villages"))
srtm_char_plot
dev.off()
#ggsave("../ssmod_v06/figures/ch6_srtm_char.pdf", srtm_char_plot, scale= 1.2,height=4, width=6, bg = "white")

## NOT IN THE TEXT
## interesting...now there is no differnces any more...
t.test(x=srtm_char$meg, y=srtm_char$vil)
t.test(x=srtm_char$tum, y=srtm_char$meg)
t.test(x=srtm_char$vil, y=srtm_char$tum)

## calculate terrain paramerts
ter.par <- terrain(x = srtm, opt = c("slope","aspect","tpi"), neighbors = 4)

# TPI for different neighborhood size:
tpi.win <- function(x, w=5) {
    m <- matrix(1/(w^2-1), nc=w, nr=w)
    m[ceiling(0.5 * length(m))] <- 0
    f <- focal(x, m)
    x - f
}
tpi.large <- tpi.win(x = srtm, w=15)
ter.par <- brick(x = c(tpi.large,ter.par))
colnames(ter.par@data@values)[1] <- "tpi_15"
colnames(ter.par@data@values)[2] <- "tpi_5"
names(ter.par)[1] <- "tpi_15"
names(ter.par)[2] <- "tpi_5"

summary(ter.par)

srtm.shade <- hillShade(slope = ter.par$slope,
                        aspect = ter.par$aspect,
                        angle = 15,
                        direction = 200,
                        normalize = TRUE)

pdf("./3results/ch6_terr_par.pdf", height=10, width=12, bg = "white")
par(mfrow = c(2,2))
plot(ter.par, 1, col = grey(c(8:0/8)),horizontal = FALSE,legend.width = 1,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
plot(srtm.shade, col = grey(0:200/200,alpha = .3), legend = FALSE,add=TRUE)
points(spdf_tum,pch = 19, cex = .5, col = "black",bg="black")
points(spdf_meg,pch = 19, cex = .5, col = "black", bg="white")
points(spdf_vil,pch = 19, cex = .5, col = "black", bg = "white")
plot(ter.par, 2, col = grey(c(0:4/4)),horizontal = FALSE,legend.width = 1,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
plot(srtm.shade, col = grey(0:200/200,alpha = .5), legend = FALSE,add=TRUE)
points(spdf_tum,pch = 19, cex = .5, col = "black",bg="black")
points(spdf_meg,pch = 19, cex = .5, col = "black", bg="white")
points(spdf_vil,pch = 19, cex = .5, col = "black", bg = "white")
plot(ter.par, 3, col = grey(c(10:0/10)),horizontal = FALSE,legend.width = 1,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
plot(srtm.shade, col = grey(0:200/200,alpha = .3), legend = FALSE,add=TRUE)
points(spdf_tum,pch = 19, cex = .5, col = "black",bg="black")
points(spdf_meg,pch = 19, cex = .5, col = "black", bg="white")
points(spdf_vil,pch = 19, cex = .5, col = "black", bg = "white")
plot(ter.par, 4, col = grey(c(8:0/10)),horizontal = FALSE,legend.width = 1,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
points(spdf_tum,pch = 19, cex = .5, col = "black",bg="black")
points(spdf_meg,pch = 19, cex = .5, col = "black", bg="white")
points(spdf_vil,pch = 19, cex = .5, col = "black", bg = "white")
dev.off()


## extract values of terrain parameters
tp_meg <- extract(x = ter.par, y = spdf_meg, buffer = 200, fun = median,df=TRUE)[,-1]
tp_tum <- extract(x = ter.par, y = spdf_tum, buffer = 200, fun = median,df=TRUE)[,-1]
tp_vil <- extract(x = ter.par, y = spdf_vil, buffer = 200, fun = median,df=TRUE)[,-1]

library(foreach)
ks_tp_meg <- foreach (i=1:4) %do% lapply(X = tp_meg[i],
                                         FUN = function(x){bkde(x[i],
                                                                kernel="normal",
                                                                bandwidth=3,
                                                                gridsize=201,
                                                                range.x = c(summary(ter.par[[i]])[1],summary(ter.par[[i]])[5])
                                                                )
                                         })

ks_tp_tum <- foreach (i=1:4) %do% lapply(X = tp_tum[i],
                                         FUN = function(x){bkde(x[i],
                                                                kernel="normal",
                                                                bandwidth=3,
                                                                gridsize=201,
                                                                range.x = c(summary(ter.par[[i]])[1],summary(ter.par[[i]])[5])
                                                                )
                                         })

ks_tp_vil <- foreach (i=1:4) %do% lapply(X = tp_vil[i],
                                         FUN = function(x){bkde(x[i],
                                                                kernel="normal",
                                                                bandwidth=3,
                                                                gridsize=201,
                                                                range.x = c(summary(ter.par[[i]])[1],summary(ter.par[[i]])[5])
                                                                )
                                         })

ter_char_tpi15 <- data.frame(ks_tp_meg[[1]][[1]],
                             ks_tp_tum[[1]][[1]]$y,
                             ks_tp_vil[[1]][[1]]$y)

colnames(ter_char_tpi15) <- c("tpi_15","meg","tum","vil")

library(reshape)
ter_char_tpi15_2 <- melt(ter_char_tpi15, id.vars = "tpi_15")

library(ggplot2)
pdf("./3results/ch6_ter_char_tpi15.pdf", height=4, width=6, bg = "white")
ter_char_tpi15_plot <- ggplot(ter_char_tpi15_2, aes(x=tpi_15, y=value)) +
    geom_line(aes(linetype = variable, color = variable)) +
    labs(x="TPI",y="density",legend="") +
    theme_bw(base_size = 12) +
    theme(legend.position="bottom",
          legend.title=element_blank()
          ) +
        scale_linetype_manual(values=c("solid", "dashed", "dotted"),
                              labels = c("Megaliths","Bronze Age Barrows", "Villages")) +
            scale_color_manual(values=c('#999999','#000000','#000000'),
                               labels = c("Megaliths","Bronze Age Barrows", "Villages"))
ter_char_tpi15_plot
dev.off()

ter_char_slope <- data.frame(ks_tp_meg[[3]][[1]],
                             ks_tp_tum[[3]][[1]]$y,
                             ks_tp_vil[[3]][[1]]$y)

colnames(ter_char_slope) <- c("slope","meg","tum","vil")

library(reshape)
ter_char_slope_2 <- melt(ter_char_slope, id.vars = "slope")

library(ggplot2)
pdf("./3results/ch6_ter_char_slope.pdf", height=4, width=6, bg = "white")
ter_char_slope_plot <- ggplot(ter_char_slope_2, aes(x=slope, y=value)) +
    geom_line(aes(linetype = variable, color = variable)) +
    labs(x="Slope (radians)",y="density",legend="") +
    theme_bw(base_size = 12) +
    theme(legend.position="bottom",
          legend.title=element_blank()
          ) +
        scale_linetype_manual(values=c("solid", "dashed", "dotted"),
                              labels = c("Megaliths","Bronze Age Barrows", "Villages")) +
            scale_color_manual(values=c('#999999','#000000','#000000'),
                               labels = c("Megaliths","Bronze Age Barrows", "Villages"))
ter_char_slope_plot
dev.off()


#### PREDICTIVE MODELING ###

#### Inductive Approach
#### --------------------------------------------------

ter.par2 <- ter.par
srtm[srtm<=0] <- NA
ter.par2$tpi_15[is.na(srtm)] <- NA
ter.par2$tpi_5[is.na(srtm)] <- NA
ter.par2$slope[is.na(srtm)] <- NA
ter.par2$aspect[is.na(srtm)] <- NA


## define the size of test and training data; in this case 30 vs 70 %
trainSize <- round(nrow(spdf_meg@data) * 0.7)
testSize <- nrow(spdf_meg@data) - trainSize

## select the training and test cases
## do this in a way that the data is reproducable

set.seed(333)
training_indices <- sample(seq_len(nrow(tp_meg)),size=trainSize)
trainSet <- tp_meg[training_indices,]
testSet <- tp_meg[-training_indices,]

## binary addition

## reclassify the geomorphometric parameter raster based on the trainSet characteristics
q.tpi15 <- quantile(tp_meg$tpi_15, probs = c(.1,.9))
q.tpi5 <- quantile(tp_meg$tpi_5, probs = c(.1,.9))
q.sl <- quantile(tp_meg$slope, probs = c(.1,.9))
q.as <- quantile(tp_meg$aspect, probs = c(.1,.9)) 

rcl.tpi15 <- c(-Inf, q.tpi15[1], 0,
            q.tpi15[1],q.tpi15[2], 1,
            q.tpi15[2],+Inf, 0
            )
rcl.tpi15 <- matrix(rcl.tpi15, ncol = 3, byrow = TRUE)
ter.par2$tpi_15 <- reclassify(x = ter.par$tpi_15,rcl = rcl.tpi15)

rcl.tpi5 <- c(-Inf, q.tpi5[1], 0,
            q.tpi5[1],q.tpi5[2], 1,
            q.tpi5[2],+Inf, 0
            )
rcl.tpi5 <- matrix(rcl.tpi5, ncol = 3, byrow = TRUE)
ter.par2$tpi_5 <- reclassify(x = ter.par$tpi_5,rcl = rcl.tpi5)

rcl.sl <- c(-Inf, q.sl[1], 0,
            q.sl[1],q.sl[2], 1,
            q.sl[2],+Inf, 0
            )
rcl.sl <- matrix(rcl.sl, ncol = 3, byrow = TRUE)
ter.par2$slope <- reclassify(x = ter.par$slope,rcl = rcl.sl)

rcl.as <- c(-Inf, q.as[1], 0,
            q.as[1],q.as[2], 1,
            q.as[2],+Inf, 0
            )
rcl.as <- matrix(rcl.as, ncol = 3, byrow = TRUE)
ter.par2$aspect <- reclassify(x = ter.par$aspect,rcl = rcl.as)

ba <- overlay(ter.par2, fun=function(w,x,y,z){return(w+x+y+z)}, unstack = TRUE)

pdf("./3results/ch6_pm_ba.pdf", height=4, width=6, bg = "white")
plot(ba,col = grey(c(7:3/8)),legend.width = 1,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
points(spdf_meg[rownames(testSet),],pch = 19, cex = .3)
dev.off()

## get the predicted values for the testSet
testSet.ba <- extract(x = ba, y = spdf_meg[rownames(testSet),])
table(testSet.ba)


## weighted binary addition

## calculating weights
weights <- apply(X = tp_meg, MARGIN = 2, FUN = function(x){sqrt(1/sd(x)*range(x)[2]-range(x)[1])})

wba <- overlay(ter.par2, fun=function(w,x,y,z){return((w*weights[1])+(x*weights[2])+(y*weights[3])+(z*weights[4]))}, unstack = TRUE)

rcl.wba <- c(-Inf, quantile(wba)[1], 0,
             quantile(wba)[1], quantile(wba)[2], 1,
             quantile(wba)[2], quantile(wba)[3], 2,
             quantile(wba)[3], quantile(wba)[4], 3,
             quantile(wba)[4], quantile(wba)[5], 4
            )

rcl.wba <- matrix(rcl.wba, ncol = 3, byrow = TRUE)
wba.rc <- reclassify(x = wba, rcl = rcl.wba)

pdf("./3results/ch6_pm_wba.pdf", height=4, width=6, bg = "white")
plot(wba.rc,col = grey(c(7:3/8)),legend.width = 1,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
points(spdf_meg[rownames(testSet),],pch = 19, cex = .3)
dev.off()


## get the predicted values for the testSet
#testSet.wba <- extract(x = wba, y = spdf_meg[rownames(testSet),])
testSet.wba.r <- extract(x = wba.rc, y = spdf_meg[rownames(testSet),])
#table(testSet.wba)
table(testSet.wba.r)

## ## reclass to four classes based on quantiles
## testSet.wba.r <- data.frame("0"=0,"1"=0,"2"=0,"3"=0,"4"=0)
## testSet.wba.r$X0 <- length(testSet.wba[testSet.wba>quantile(wba)[1] & testSet.wba<quantile(wba)[2]])
## testSet.wba.r$X1 <- length(testSet.wba[testSet.wba>quantile(wba)[1] & testSet.wba<quantile(wba)[2]])
## testSet.wba.r$X2 <- length(testSet.wba[testSet.wba>quantile(wba)[2] & testSet.wba<quantile(wba)[3]])
## testSet.wba.r$X3 <- length(testSet.wba[testSet.wba>quantile(wba)[3] & testSet.wba<quantile(wba)[4]])
## testSet.wba.r$X4 <- length(testSet.wba[testSet.wba>quantile(wba)[4]])


## gain statistics
## 1-(%Area/%Sites)

## -- in progress -- ## 
## ## gain function for raster objects
## gain <- function(area,sites){
##     # percent area
##     prop.table(table(area@data@values))*100
##     # percent sites
## 
## }

## %Area
pa.ba <- prop.table(table(ba@data@values))*100
pa.wba <- prop.table(table(wba.rc@data@values))*100

## %sites
ps.ba <- prop.table(table(testSet.ba))*100
ps.wba <- prop.table(testSet.wba.r)*100

## gain
1-(pa.ba[5]/ps.ba[4])
1-(pa.wba[5]/ps.wba[5])


## LOGISTIC REGRESSION
## extract geomorphometric values 
## --------------------------------------------------
set.seed(123)

## at random locations [n = amount of megatlihs], i.e. absence data
## first create random points [seed already defined -> for reproducability]
rand_points <- randomPoints(mask = srtm,
                            p = spdf_meg,
                            n = length(spdf_meg)) # requires dismo
## extract data
rand_points <- extract(x = ter.par, y = rand_points, buffer = 200, fun = median,df=TRUE)[,-1]

## bring empiricial and artifical points together
emp_ran <- c(rep(1, nrow(tp_meg)), rep(0, nrow(rand_points)))
geom_data <- data.frame(cbind(as.factor(emp_ran),rbind(tp_meg, rand_points)))
names(geom_data)[1] <- "emp_ran"
geom_data <- geom_data[complete.cases(geom_data)==TRUE,]
head(geom_data,2)
tail(geom_data,2)

glm1 <- glm(emp_ran ~ tpi_15+tpi_5+slope+aspect, data=geom_data, family = binomial(link=logit))
summary(glm1)
glm1 

## do not use aspect - scale of measurement

## check for correlation in the TPIs
cor(geom_data$tpi_15,geom_data$tpi_5)

glm2 <- glm(emp_ran ~ slope+tpi_15, data=geom_data, family = binomial(link=logit))
summary(glm2)
glm2

ge1 <- evaluate(p = tp_meg, a = rand_points, model = glm1)
ge2 <- evaluate(p = tp_meg, a = rand_points, model = glm2)

tr1 <- threshold(ge1, "spec_sens")
pdf("./3results/ch6_pm_lr_1.pdf", height=4, width=6, bg = "white")
plot(pg1 > tr1, main="glm1 presence (white) and absence (gray)",col = grey(c(1:2/2)),legend = FALSE,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
points(spdf_meg,pch=19,cex=.3)
dev.off()

tr2 <- threshold(ge2, "spec_sens")
pdf("./3results/ch6_pm_lr_2.pdf", height=4, width=6, bg = "white")
plot(pg2 > tr2, main="glm2 presence (white) and absence (gray)",col = grey(c(1:2/2)),legend = FALSE,cex.axis=.9,tcl=-.3,mgp = c(3,.2,0))
points(spdf_meg,pch=19,cex=.3)
dev.off()

##### GAIN ####

## %Area
pa.ba <- prop.table(table(ba@data@values))*100
pa.wba <- prop.table(table(wba.rc@data@values))*100

## %sites
ps.ba <- prop.table(table(testSet.ba))*100
ps.wba <- prop.table(table(testSet.wba.r))*100

## gain
g.ba <- 1-(pa.ba[5]/ps.ba[4])
g.wba <- 1-(pa.wba[5]/ps.wba[5])

## gain model glm1
## percent area positive
pa.glm1 <- (100*length(pg1[pg1 > tr1]))/length(pg1[!is.na(pg1)])
## percent sites on positive areas
meg.logit <- extract(x = pg1, y = spdf_meg)
ps.glm1 <- (100*length(meg.logit[meg.logit > tr1]))/length(meg.logit)
##gain
g.glm1 <- 1-pa/ps

## gain model glm2
## percent area positive
pa.glm2 <- (100*length(pg2[pg2 > tr2]))/length(pg2[!is.na(pg2)])
## percent sites on positive areas
meg.logit <- extract(x = pg2, y = spdf_meg)
ps.glm2 <- (100*length(meg.logit[meg.logit > tr2]))/length(meg.logit)
##gain
g.glm2 <- 1-pa/ps

gain <- data.frame(BA=c(pa.ba[5],ps.ba[4],g.ba),
                WBA=c(pa.wba[5],ps.wba[5],g.wba),
                GLM1=c(pa.glm1,ps.glm1,g.glm1),
                GLM2=c(pa.glm2,ps.glm2,g.glm2),
                row.names = c("p_a","p_s","gain")
                )
library(knitr)
kable(gain, format = "latex")
