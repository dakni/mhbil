############################################## #
## R-Script QAAM 1, Modelling ... Landscapes
##                  chapter 5
## =========================================== =
## Project: Modelling Human Behaviour in
##          Landscapes - Basic concepts and 
##          modelling elements
##          Springer Series:
##          Quantitative Archaeology and
##          Archaeological modelling
## Authors: Oliver Nakoinz & Daniel Knitter
## Version: 06
## Date of last changes: 07.08.2015
## Data: Megaliths, Bronze Age Barrows and Medieval villages
## Author of data: gshdl
## Purpose: didactic
## Content: 1. regression, 2. interpolation
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
load("4ws/ws04.rws")

# 1. Regression ================================
library(gdata) 

file_vil    <- "1data/villages.xls" 
df_vil <- read.xls(file_vil)
# Correlation
X <- df_vil[,3]
E_vil <- df_vil[,1]

cor(cbind(X,E_vil), method="pearson")
# Regression T1, T2
pdf("6pictures/c5_T1T2vil.pdf", height=4, width=6, bg = "white") 
    plot(X,E_vil,col="black", pch=16)
    lines(c(df_vil[1,3],df_vil[13,3]),c(df_vil[1,1], df_vil[13,1]) ) 
    years <- 1259:1350
    y <- 1 + 0.003* (years-1259)^2
    lines (years,y)
dev.off() 

E <-  E_vil
T1_vil <- 1 + (12/91)*(X - 1259)
T2_vil <- 1 + 0.003* (X-1259)^2
sum(((E-T1_vil)^2)/T1_vil)
sum(((E-T2_vil)^2)/T2_vil)

yr <- X - 1259 
T3_vil <- lm(E ~ yr)
coef(T3_vil)
plot(X- 1259,E_vil,col="black", pch=16)
abline(T3_vil)


pdf("6pictures/c5_T3vil.pdf", height=4, width=6, bg = "white") 
    plot(X- 1259,E_vil,col="black", pch=16)
    abline(T3_vil)
dev.off() 

fT3_vil <- fitted(T3_vil)
sum(((E-fT3_vil)^2)/fT3_vil)

summary(T3_vil)$r.squared

T4_vil <- lm(E ~ I(yr^2))
coef(T4_vil)
fT4_vil <- fitted(T4_vil)
sum(((E-fT4_vil)^2)/fT4_vil)
summary(T4_vil)$r.squared

T5_vil <- lm(E ~ poly(yr, 4, raw=TRUE))
coef(T5_vil)
fT5_vil <- fitted(T5_vil)
sum(((E-fT5_vil)^2)/fT5_vil)
summary(T5_vil)$r.squared

pdf("6pictures/c5_T5vil.pdf", height=4, width=6, bg = "white") 
    plot(X- 1259,E_vil,col="black", pch=16)
    lines(yr,fT5_vil)
dev.off() 

chi <- 1:20
for(i in seq(1:20)){
    T6_vil <- lm(E ~ poly(yr, i, raw=TRUE))
    fT6_vil <- fitted(T6_vil)
    chi[i] <- sum(((E-fT6_vil)^2)/fT6_vil)
}
chi
pdf("6pictures/c5_chi.pdf", height=4, width=6, bg = "white") 
    plot(chi,col="black", pch=16)
    lines(chi)
dev.off() 

r <- 1:20
for(i in seq(1:20)){
    T7_vil <- lm(E ~ poly(yr, i, raw=TRUE))
    fT7_vil <- fitted(T6_vil)
    r[i] <- summary(T7_vil)$r.squared
}
r
pdf("6pictures/c5_r.pdf", height=4, width=6, bg = "white") 
    plot(r,col="black", pch=16)
    lines(r)
dev.off() 

# bootstrapping
data <- data.frame(x=yr, y=E)
its <- 100
dd <- 6
res <-data.frame()
for (d in 1:dd) {
    ei <- numeric(d)
    eo <- numeric(d)
    for (it in 1:its)     {
        sel    <- sample(x=1:13, size=8,replace=F)
        train  <- data[sel,]
        seli   <- !(1:13 %in% sel)
        test   <- data[seli,]
        T7_vil <- lm(y~poly(x,d),data=train)
        ei[it] <- mean(T7_vil$residuals^2)
        pred   <- predict(T7_vil, newdata=data.frame(data[seli,]))
        eo[it] <- mean((test-pred)^2)
    }
    eiv <- mean(ei)
    eov <- mean(eo)
    res <- rbind(res,data.frame(d=d,eiv=eiv,eov=eov))
}
pdf("6pictures/c5_boot.pdf", height=4, width=6, bg = "white") 
    plot(res$d,res$eiv,col="grey",log="y",ylim=c(min(res$eiv),max(res$eov)), pch=16, type="b")
    points(res$d,res$eov, type ="b", pch=16, col="black")
dev.off() 

# decomposition
T3_vil <- lm(E ~ yr)
fT3_vil <- fitted(T3_vil)
pred3 <- predict(T3_vil)
E3 <- E-pred3

T8_vil <- lm(E3 ~ poly(yr, 4, raw=TRUE))
fT8_vil <- fitted(T8_vil)
pred8 <- predict(T8_vil)
E8 <- E3-pred8

T9_vil <- lm(E8 ~ poly(yr, 12, raw=TRUE))
fT9_vil <- fitted(T9_vil)

pdf("6pictures/c5_decomp.pdf", height=4, width=6, bg = "white") 
    plot(X- 1259,E_vil,col="white", pch=16, ylim=c(-2,12))
    lines(yr,fT3_vil, col="grey")
    curve(coef(T8_vil)[1] + coef(T8_vil)[2]*x + coef(T8_vil)[3]*x^2 + coef(T8_vil)[4]*x^3 + coef(T8_vil)[5]*x^4, add = TRUE)
    lines(yr,fT9_vil, col="grey", lty=2)
dev.off() 



# 2. Interpolation ===============================

library(gstat)
library(sp)

meg_idw     <- idw(fs_vd_spdf@data$fs_vd ~ 1, fs_vd_spdf, sgdf) ###################################################################################
image(meg_idw, col = gray.colors(25, start = 0.4, end = 0.97))
contour(meg_idw, add=T) 
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4) 

library(maptools)
library(spatstat)
rw <- 500  
sd <- 2000
dens_r <- density(ppp_meg, sd, eps=rw, edge=TRUE, at="pixels") 
#plot(dens_r, col = gray.colors(25, start = 0.97, end = 0.4))
df_dens_r <- as.SpatialGridDataFrame.im(dens_r)
meg_kde_samppoints <- spsample(df_dens_r, 1000,  type="random")
meg_kde_samp <- overlay(x=df_dens_r, y= meg_kde_samppoints)
##### meg_kde_samp <- over(meg_kde_samppoints, df_dens_r)
vt2    <- variogram(meg_kde_samp@data$v ~ 1, meg_kde_samp)
v.fit2 <- fit.variogram(vt2, vgm(1, "Gau", 5000, 1), fit.sills = TRUE, fit.ranges = TRUE, fit.method = 7)
plot(vt2,v.fit2)
k2     <- krige(meg_kde_samp@data$v ~ 1, meg_kde_samp, df_dens_r, v.fit2, nmin = 3, maxdist = 10000, nmax = 8)
image(k2, col = gray.colors(25, start = 0.97,  end = 0.4))
contour(k2, add=T)    
points(ppp_meg$x, ppp_meg$y, pch=16, cex=0.4) 

save.image("4ws/ws05.rws")






