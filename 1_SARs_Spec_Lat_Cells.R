# SAR models Speciation-Latitude (See Figure 1 of manuscript)
# By Adrian Garcia
# Last update: 16.3.2025
rm(list=ls())
gc()

#Packages
library(spdep)          
library(ncf)
library(ggpubr)

setwd("C:/Users/Administrator/Dropbox/PANDEMIC PRIORITIES/5.LGD_Speciation/")
setwd("C:/Users/garciaa/Dropbox/MANUSCRIPTS/2025/LGSpeciation/")

load("THIRD TRY/Docs_Submission_NEE/Codes/SARs_Rich_Lat_CELLS.RData")

###DATA
Data_all<- read.csv("Codes/Github/Data/XY_RichSpec_all.csv")
Data_all<- Data_all[, 2:5]

####OLS
ols<-lm(Data_all$meanDR ~ abs(Data_all$y))
plot( abs(Data_all$y),Data_all$meanDR, pch=19, cex= 1.5)
abline(ols, lwd=2, col="brown2")
summary(ols)
res.ols <- residuals(ols)

nb5k <- spdep::dnearneigh(x=as.matrix(Data_all[,1:2]), d1 = 0, d2 = 5, longlat = FALSE)
str(nb5k, list.len=5, give.attr = F)
lw5k <-nb2listw(nb5k, style = "W", zero.policy = TRUE)
str(lw5k$weights, list.len=5, give.attr = F)

#####MORAN TEST (to detect spatial autocorrelation)
MORAN<-lm.morantest(ols, lw5k, zero.policy = T)
MORAN$p.value

#Correlogram for OLS model residuals
cor.ols.res<-correlog(as.numeric(Data_all$y), as.numeric(Data_all$x), z=residuals(ols), na.rm=T, increment=1, resamp=1)

#Set plotting options to plot correlogram
par(mar=c(5,5,2,0.1), mfrow=c(1,1))
dev.off()
#Plot correlogram
PlotOLS_all<-plot(cor.ols.res$correlation, type="l", pch=19, cex=1, col="gray20", lwd=1,
                  ylim=c(-1.5,.5 ), xlab="distance", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) +
  abline(h=0, lwd=2, col="red") +
  title(main="OLS residuals Richness~SpeciationRate_ALL TAXA", cex=5)

# annotate
legend(x=200, y=1.5, legend=c("Residuals OLS"), pch=c(1), bty="n", cex=1, ncol = 5)

###### SAR ANALYSIS FOR ALL SPECIES
# 1.1) Global connectivity (neighbourhood) matrix
# Define connectivity matrix (0/1)
nbdist<-dnearneigh(x=as.matrix(Data_all[,1:2]), d1=0, d2=5) 

# Compute the Euclidean distance between neighbouring sites
neigh.dist<-nbdists(nbdist,
                    as.matrix(Data_all[,1:2]), longlat=F)

# Compute the inverse distance weigthed matrix
inverse<-lapply(neigh.dist, function(x) (1/(x^2)))

# Coding style W = row standardised
nlw <- nb2listw(neighbours=nbdist, 
                glist=inverse,style="W", zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets

# 1.2 ) Global SAR model
sar_ALL <- spatialreg::errorsarlm(meanDR ~ abs(y),data = Data_all,
                                  listw = nlw, zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets

summary(sar_ALL, Nagelkerke=TRUE)

###PLOT SAR
library(ggplot2)
RegresALL<-ggplot(data = Data_all, aes(x =abs(y) , y = meanDR)) + 
  geom_point(color='gray75', size=2, alpha=.3) +
  ggtitle("ALL AMPHIBIANS")+
  #geom_smooth(data =GlobCorr, method = "lm",color="red",se=FALSE)+
  geom_abline(slope=coef(sar_ALL)[[3]], intercept = coef(sar_ALL)[[2]], linewidth = 1, col="red", linetype="dashed", alpha=.5) +
  theme_classic()
RegresALL

#######ANALYSES BY ORDER##############
###DATA
Data_ANURA<- read.csv("Codes/Github/Data/XY_RichSpec_anura.csv")
Data_ANURA<- Data_ANURA[, 2:5]

####OLS
ols_ANURA<-lm(Data_ANURA$mean.speciationAnura~abs(Data_ANURA$Lat))
summary(ols_ANURA)
res.ols_ANURA <- residuals(ols_ANURA)

#Correlogram for OLS model residuals
cor.ols.res_ANURA<-correlog(as.numeric(Data_ANURA$Lat), as.numeric(Data_ANURA$Lon), z=residuals(ols_ANURA), na.rm=T, increment=1, resamp=1)

#Set plotting options to plot correlogram
par(mar=c(5,5,2,0.1), mfrow=c(1,1))
dev.off()

#Plot correlogram
PlotOLS_ANURA<-plot(cor.ols.res_ANURA$correlation, type="l", pch=19, cex=1, col="gray40", lwd=1,
                    ylim=c(-1.1,1.1 ), xlab="distance", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) +
  abline(h=0, lwd=2, col="red") +
  title(main="OLS residuals SpeciationRate_ANURA~Latitude", cex=5)

# annotate
legend(x=200, y=1.5, legend=c("Residuals OLS_ANURA"), pch=c(1), bty="n", cex=1, ncol = 5)

# 1.1) Global connectivity (neighbourhood) matrix

# Define connectivity matrix (0/1)
nbdist_ANURA<-dnearneigh(x=as.matrix(Data_ANURA[,1:2]), d1=0, d2=5)

# Compute the Euclidean distance between neighbouring sites
neigh.dist_ANURA<-nbdists(nbdist_ANURA,
                          as.matrix(Data_ANURA[,1:2]), longlat=F)

# Compute the inverse distance weigthed matrix
inverse_ANURA<-lapply(neigh.dist_ANURA, function(x) (1/(x^2)))

# Coding style W = row standardised
nlw_ANURA <- nb2listw(neighbours=nbdist_ANURA, 
                      glist=inverse_ANURA,
                      style="W", 
                      zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets

# 1.2 ) Global SAR model
sar__ANURA <- spatialreg::errorsarlm(mean.speciationAnura ~ abs(Lat),
                                     data = Data_ANURA,
                                     listw = nlw_ANURA, 
                                     zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets
summary(sar__ANURA, Nagelkerke=TRUE)

###PLOT SAR
library(ggplot2)
Regres_ANURA<-ggplot(data = Data_ANURA, aes(x = abs(Lat)  , y = mean.speciationAnura )) + 
  geom_point(color='gray75', size=2, alpha=.3) +
  ggtitle("FROGS and TOADS")+
  #geom_smooth(data =GlobCorr, method = "lm",color="red",se=FALSE)+
  geom_abline(slope=coef(sar__ANURA)[[3]], intercept = abs(coef(sar__ANURA)[[2]]), size = 1, col="red", linetype="dashed", alpha=.5)+
  theme_classic()
Regres_ANURA

###CAUDATA
Data_CAUDATA<- read.csv("Codes/Github/Data/XY_RichSpec_caudata.csv")
Data_CAUDATA<- Data_CAUDATA[, 2:5]

####OLS
ols_CAUDATA<-lm(Data_CAUDATA$mean.speciationCaudata~abs(Data_CAUDATA$Lat))
summary(ols_CAUDATA)
res.ols_CAUDATA <- residuals(ols_CAUDATA)

#Correlogram for OLS model residuals
cor.ols.res_CAUDATA<-correlog(as.numeric(Data_CAUDATA$Lat), as.numeric(Data_CAUDATA$Lon), z=residuals(ols_CAUDATA), na.rm=T, increment=1, resamp=1)

#Set plotting options to plot correlogram
par(mar=c(5,5,2,0.1), mfrow=c(1,1))
dev.off()

#Plot correlogram
PlotOLS_CAUDATA<-plot(cor.ols.res_CAUDATA$correlation[2:240], type="l", pch=19, cex=1, col="gray40", lwd=1,
                      ylim=c(-2,2 ), xlab="distance", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) +
  abline(h=0, lwd=2, col="red") +
  title(main="OLS residuals Richness~SpeciationRate_CAUDATA", cex=5)

# annotate
legend(x=200, y=1.5, legend=c("Residuals OLS_CAUDATA"), pch=c(1), bty="n", cex=1, ncol = 5)

# Define connectivity matrix (0/1)
nbdist_CAUDATA<-dnearneigh(x=as.matrix(Data_CAUDATA[,1:2]), d1=0, d2=5) 

# Compute the Euclidean distance between neighbouring sites
neigh.dist_CAUDATA<-nbdists(nbdist_CAUDATA,
                            as.matrix(Data_CAUDATA[,1:2]), longlat=F)

# Compute the inverse distance weigthed matrix
inverse_CAUDATA<-lapply(neigh.dist_CAUDATA, function(x) (1/(x^2)))

# Coding style W = row standardised
nlw_CAUDATA <- nb2listw(neighbours=nbdist_CAUDATA, 
                        glist=inverse_CAUDATA,
                        style="W", 
                        zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets

# 1.2 ) Global SAR model
sar__CAUDATA <- spatialreg::errorsarlm(mean.speciationCaudata ~ abs(Lat),
                                       data = Data_CAUDATA,
                                       listw = nlw_CAUDATA, 
                                       zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets
summary(sar__CAUDATA, Nagelkerke=TRUE)

###PLOT SAR
Regres_CAUDATA<-ggplot(data = Data_CAUDATA, aes(x = abs(Lat) , y = mean.speciationCaudata)) + 
  geom_point(color='gray75', size=2, alpha=.3) +
  ggtitle("SALAMANDERS")+
  #geom_smooth(data =GlobCorr, method = "lm",color="red",se=FALSE)+
  geom_abline(slope=coef(sar__CAUDATA)[[3]], intercept = coef(sar__CAUDATA)[[2]], size = 1, col="red", linetype="dashed", alpha=.5)+
  theme_classic()
Regres_CAUDATA

###GYMNOPHIONA
Data_GYMNO<- read.csv("Codes/Github/Data/XY_RichSpec_gymno.csv")
Data_GYMNO<- Data_GYMNO[, 2:5]

####OLS
ols_GYMNO<-lm(Data_GYMNO$mean.speciationGymno~abs(Data_GYMNO$Lat))
summary(ols_GYMNO)
res.ols_GYMNO <- residuals(ols_GYMNO)

#Correlogram for OLS model residuals
cor.ols.res_GYMNO<-correlog(as.numeric(Data_GYMNO$Lat), as.numeric(Data_GYMNO$Lon), z=residuals(ols_GYMNO), na.rm=T, increment=1, resamp=1)

#Set plotting options to plot correlogram
par(mar=c(5,5,2,0.1), mfrow=c(1,1))
dev.off()

#Plot correlogram
PlotOLS_GYMNO<-plot(cor.ols.res_GYMNO$correlation[2:232], type="l", pch=19, cex=1, col="gray40", lwd=1,
                    ylim=c(-1.5,1 ), xlab="distance", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) +
  abline(h=0, lwd=2, col="red") +
  title(main="OLS residuals Richness~SpeciationRate_GYMNO", cex=5)

# annotate
legend(x=200, y=1.5, legend=c("Residuals OLS_GYMNO"), pch=c(1), bty="n", cex=1, ncol = 5)

# Define connectivity matrix (0/1)
nbdist_GYMNO<-dnearneigh(x=as.matrix(Data_GYMNO[,1:2]), d1=0, d2=5) 

# Compute the Euclidean distance between neighbouring sites
neigh.dist_GYMNO<-nbdists(nbdist_GYMNO,
                          as.matrix(Data_GYMNO[,1:2]), longlat=F)

# Compute the inverse distance weigthed matrix
inverse_GYMNO<-lapply(neigh.dist_GYMNO, function(x) (1/(x^2)))

# Coding style W = row standardised
nlw_GYMNO <- nb2listw(neighbours=nbdist_GYMNO, 
                      glist=inverse_GYMNO,
                      style="W", 
                      zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets

# 1.2 ) Global SAR model
sar__GYMNO <- spatialreg::errorsarlm(mean.speciationGymno ~ abs(Lat),
                                     data = Data_GYMNO,
                                     listw = nlw_GYMNO, 
                                     zero.policy=TRUE) # Use zero policy to avoid error for any empty neighbour sets
summary(sar__GYMNO, Nagelkerke=TRUE)

###PLOT SAR
Regres_GYMNO<-ggplot(data = Data_GYMNO, aes(x = abs(Lat) , y = mean.speciationGymno)) + 
  geom_point(color='gray75', size=2, alpha=.3) +
  ggtitle("CAECILIANS")+
  #geom_smooth(data =GlobCorr, method = "lm",color="red",se=FALSE)+
  geom_abline(slope=coef(sar__GYMNO)[[3]], intercept = coef(sar__GYMNO)[[2]], size = 1, col="red", linetype="dashed", alpha=.5)+
  theme_classic()
Regres_GYMNO

plots<-list(RegresALL,Regres_ANURA,Regres_CAUDATA,Regres_GYMNO)

ggarrange(plotlist = plots, ncol = 1, nrow = 4) #Plot saved as PDF (10x3 Portrait)

save.image("C:/Users/Administrator/Dropbox/MANUSCRIPTS/2025/LGSpeciation/Codes/SARs_Rich_Lat_CELLS.RData")

