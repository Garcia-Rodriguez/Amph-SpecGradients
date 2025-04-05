# Mapping speciation rates
# Author: Adrian Garcia-Rodriguez
# Last update: 16.3.2025
rm(list=ls())
gc()

#Load Packages
pcks<-c("viridis","rasterVis","maps","letsR","sp","rgdal","raster","Rfast", "ggplot2", "sf")
sapply(pcks, require, character.only = TRUE)
#sapply(pcks, install.packages, character.only = TRUE)

#setwd("C:/Users/garciaa/Dropbox/MANUSCRIPTS/2025/LGSpeciation")
setwd("~/LGSpeciation")

#Load Functions
source("~/DFtoRASTER.R")
source("~/addRates.to.PAM.R")
source("~/PAM_mean.sd.cv_func.R")

####SPECIATION RATES (DR)
ratesFull<-read.csv("~/DR_100trees.csv")

####PAM
load("~/amphi_pam.Rdata")
amphi.pam<- amphi_pam$Presence_and_Absence_Matrix

###################
#   ALL SPECIES   #
###################
PAM_specrateFULL<-addRates.to.PAM(PAM = amphi.pam, rates = ratesFull) ##This function crops the PAM based on the species with calculated speciation rate and replace these values for the 1's (presences) in the PAM 

## Estimate mean, sd and cv speciation per grid cell
PAMmean.sd.cv_FULL<- PAM_mean.sd.cv(df = PAM_specrateFULL) # 6443 spp   # This function uses the PAM_specrate object created above to estimate statistics across rows

## Create Raster maps for each statistic
mean.sp.rastFULL <- DFtoRaster(PAMmean.sd.cv_FULL,d = 3)
sd.sp.rastFULL <- DFtoRaster(PAMmean.sd.cv_FULL,d = 4)
cv.sp.rastFULL <- DFtoRaster(PAMmean.sd.cv_FULL,d = 5)
rich.rastFULL <- DFtoRaster(PAMmean.sd.cv_FULL,d = 6)
max.rastFULL <- DFtoRaster(PAMmean.sd.cv_FULL,d = 7)
min.rastFULL <- DFtoRaster(PAMmean.sd.cv_FULL,d = 8)
range.rastFULL <- DFtoRaster(PAMmean.sd.cv_FULL,d = 9)

## Stack rasters
stack.FULL<- stack(rich.rastFULL, mean.sp.rastFULL, sd.sp.rastFULL, cv.sp.rastFULL,  min.rastFULL, max.rastFULL,range.rastFULL)

## Add proper names
names(stack.FULL) <- c("Spp richness ALL" ,"mean DR ALL", "sd DR ALL", "cv DR ALL", "min DR ALL", "max DR ALL", "range DR ALL" )

##Plots to check (Choose number 1:"Spp.richness.ALL" 2:"mean.DR.ALL"  3:"sd.DR.ALL"  4:"cv.DR.ALL"  5:"min.DR.ALL" 6:max.DR.ALL", 7: "range.DR.ALL")    
x=c(1,2,4,6)  
plot(stack.FULL[[x]], col= rev(plasma(256)), main= names(stack.FULL[[x]]))
  names(stack.FULL)

###################
#     ANURANS     #
###################

## Add rates to PAM
ratesAnurans <- ratesFull[ratesFull$ORDER=="ANURA",]
PAM_specrateANURA <- addRates.to.PAM(PAM = amphi.pam, rates = ratesAnurans) #5681 spp

###Estimate mean, sd and cv speciation per grid cell
PAMmean.sd.cv_ANURA<- PAM_mean.sd.cv(df = PAM_specrateANURA)

## Stack rasters
stack.ANURA <- stack(rich.rastANURA, mean.sp.rastANURA, sd.sp.rastANURA, cv.sp.rastANURA,  min.sp.rastANURA, max.sp.rastANURA,range.sp.rastANURA)

## Add proper names
names(stack.ANURA) <- c("Spp richness ANURA" ,"mean DR ANURA", "sd DR ANURA", "cv DR ANURA", "min DR ANURA", "max DR ANURA", "range DR ANURA" )

##Plots to check (Choose number 1:"Spp.richness.ANURA" 2:"mean.DR.ANURA"  3:"sd.DR.ANURA"  4:"cv.DR.ANURA"  5:"min.DR.ANURA" 6:max.DR.ANURA", 7: "range.DR.ANURA")    
x=c(1,2,4,6)  
plot(stack.ANURA[[x]], col= rev(plasma(256)), main= names(stack.ANURA[[x]]))
names(stack.ANURA)

###################  
#     CAUDATA     #
###################

###Add rates to PAM
ratesCaudata <- ratesFull[ratesFull$ORDER=="CAUDATA",]
PAM_specrateCAUDATA <- addRates.to.PAM(PAM = amphi.pam, rates = ratesCaudata)

###Estimate mean, sd and cv speciation per grid cell
PAMmean.sd.cv_CAUDATA<- PAM_mean.sd.cv(df = PAM_specrateCAUDATA)

###Create Rasters
mean.sp.rastCAUDATA <- DFtoRaster(PAMmean.sd.cv_CAUDATA,d = 3) 
sd.sp.rastCAUDATA <- DFtoRaster(PAMmean.sd.cv_CAUDATA,d = 4) 
cv.sp.rastCAUDATA <- DFtoRaster(PAMmean.sd.cv_CAUDATA,d = 5)
rich.rastCAUDATA <- DFtoRaster(PAMmean.sd.cv_CAUDATA,d = 6)
max.sp.rastCAUDATA <- DFtoRaster(PAMmean.sd.cv_CAUDATA,d = 7) 
min.sp.rastCAUDATA <- DFtoRaster(PAMmean.sd.cv_CAUDATA,d = 8)
range.sp.rastCAUDATA <- DFtoRaster(PAMmean.sd.cv_CAUDATA,d = 9)

## Stack rasters
stack.CAUDATA <- stack(rich.rastCAUDATA, mean.sp.rastCAUDATA, sd.sp.rastCAUDATA, cv.sp.rastCAUDATA,  min.sp.rastCAUDATA, max.sp.rastCAUDATA,range.sp.rastCAUDATA)

## Add proper names
names(stack.CAUDATA) <- c("Spp richness CAUDATA" ,"mean DR CAUDATA", "sd DR CAUDATA", "cv DR CAUDATA", "min DR CAUDATA", "max DR CAUDATA", "range DR CAUDATA")

##Plots to check (Choose number 1:"Spp.richness.CAUDATA" 2:"mean.DR.CAUDATA"  3:"sd.DR.CAUDATA"  4:"cv.DR.CAUDATA"  5:"min.DR.CAUDATA" 6:max.DR.CAUDATA", 7: "range.DR.CAUDATA")    
x=c(1,2,4,6)  
plot(stack.CAUDATA[[x]], col= rev(plasma(256)), main= names(stack.CAUDATA[[x]]))
names(stack.CAUDATA)

###########################
#       GYMNOPHIONA       #
###########################

###Add rates to PAM
ratesGymno <- ratesFull[ratesFull$ORDER=="GYMNOPHIONA",]
PAM_specrateGYMNO <- addRates.to.PAM(PAM = amphi.pam, rates = ratesGymno)

###Estimate speciation statistics per grid cell
PAMmean.sd.cv_GYMNO<- PAM_mean.sd.cv(df = PAM_specrateGYMNO)

##Plots
par(mfrow=c(2,2))
plot(mean.sp.rastGYMNO, col= rev(plasma(256)), main= paste0("Mean speciation DR GYMNOPHIONA"," (",ncol(PAM_specrateGYMNO)-2 ," species)"))

plot(sd.sp.rastGYMNO, col= rev(plasma(256)), main= paste0("SD DR GYMNOPHIONA"," (",ncol(PAM_specrateGYMNO)-2 ," species)"))

plot(cv.sp.rastGYMNO, col= rev(inferno(256)), main= paste0("CV DR GYMNOPHIONA"," (",ncol(PAM_specrateGYMNO)-2 ," species)"))

plot(rich.rastGYMNO, col= rev(inferno(256)), main= paste0("CV DR GYMNOPHIONA"," (",ncol(PAM_specrateGYMNO)-2 ," species)"))


