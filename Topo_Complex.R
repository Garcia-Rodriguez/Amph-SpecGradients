####PREDICTORS
###TOPOGRAPHIC COMPLEXITY
alt<- raster("~/DATA/alt.bil")
plot(alt)

###METRICS SD_alt, Roughness, TRI, TPI
#1 degree 100x100 elevation variability (SD_alt)
TCI_1dg<-aggregate(alt, fact=100, fun=sd, na.rm=TRUE)
TCI_1dg<-resample(TCI_1dg, biom_realm)
names(TCI_1dg)<-"SD_Alt"
plot(TCI_1dg)

###Roughness: the difference between the maximum and the minimum value of a cell and its 8 surrounding cells.
Rough<-terrain(alt, opt='roughness',8)
Rough_1dg<-resample(Rough, biom_realm)
plot(Rough_1dg)

#TRI (Terrain Ruggedness Index):  mean of the absolute differences between the value of a cell and the value of its 8 surrounding cells
TRI<-terrain(alt, opt='TRI', neighbors=8)
TRI_1dg<-resample(TRI, biom_realm)
plot(TRI_1dg)

##Topographic Position Index:  difference between the value of a cell and the mean value of its 8 surrounding cells
TPI<-terrain(alt, opt='TPI', neighbors=8)
TPI_1dg<-resample(TPI, biom_realm)
plot(TPI_1dg)

TOPOMETRICS<-stack(TCI_1dg, Rough_1dg, TRI_1dg, TPI_1dg)
plot(TOPOMETRICS)

##EXTRACT DATA
XY<-cbind(PAM_REGS$`Longitude(x)`, PAM_REGS$`Latitude(y)`)
TopoComplexCell<- as.data.frame(extract(TOPOMETRICS,XY))
TopoComplexRegs<-as.data.frame(cbind(as.character(PAM_REGS$BioRegion),XY,TopoComplexCell))
names(TopoComplexRegs)<-c("BioRegion", "Lon", "Lat", "SD_Alt", "Rough", "TRI", "TPI")

###MEANS PER REGION
library(dplyr)

MeanTopoComplexReg<-TopoComplexRegs %>%
  group_by(BioRegion) %>%
  summarise_at(vars(-c(Lon,Lat)), funs(mean(., na.rm=TRUE)))
