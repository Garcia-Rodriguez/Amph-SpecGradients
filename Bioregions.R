### Empty environment
rm(list=ls());gc()

####Load libraries
library(raster)
library(rgdal)
library(maps)
library(mapdata)

###SPECIATION RATE LAYER to define resolution of the layer to rasterize BIOMES AND REALMS
Spec.Rate<-raster("~/DATA/Speciation.asc")
mod.r<-raster(ncol=360,nrow=126,xmn=-180, xmx=180, ymn=-55, ymx=71)

####OLSON TERRESTRIAL BIOMES (2001)
olson<-rgdal::readOGR(dsn ="~/DATA/REALMS", layer="REALMS_OLSON")
unique(olson$BIOME)

####REALMS AS JETZ&FINE
JetzFine<-rgdal::readOGR(dsn ="~/DATA/REALMS", layer="REALMS_JETZ&FINE")

###RASTERIZE BIOMES
biome.r <- rasterize(olson, mod.r, field = 'BIOME')
###RASTERIZE REALMS
realm.r <- rasterize(JetzFine, mod.r, field = 'WWF_REALM2')

##STACK layers
biom_realm<-stack(biome.r, realm.r)
names(biom_realm)<-c("biom", "realm")

####GET BIOME REALM FOR EACH CELL ON PAM
load("~DATA/amphi_pam.Rdata")
PAM<-amphi_pam$Presence_and_Absence_Matrix
XY<-PAM[ ,1:2]
cell_Bioregs<-as.data.frame(extract(biom_realm, XY))
PAM_Bioregs<-cbind(cell_Bioregs,PAM)

###CHANGE CODES FOR NAMES
plot(biome.r==7)
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==1]<-"MOIST"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==2]<-"DRY"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==3]<-"MOIST"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==4]<-"TEMP"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==5]<-"TEMP"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==6]<-"BOR"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==7]<-"DRY"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==8]<-"GRASS"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==9]<-"GRASS"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==10]<-NA
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==11]<-"TUNDRA"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==12]<-"MED"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==13]<-"DESERT"
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==14]<-NA
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==98]<-NA
PAM_Bioregs["biom"][PAM_Bioregs["biom"]==99]<-NA

##REALMS
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==1]<-"Afrotropics"
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==2]<-NA
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==3]<-"Australia"
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==4]<-"IndoMalay"
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==5]<-"Madagascar"
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==6]<-"North America"
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==7]<-"South America"
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==8]<- NA
PAM_Bioregs["realm"][PAM_Bioregs["realm"]==9]<-"Eurasia"

BioRegion<-paste(PAM_Bioregs$biom, PAM_Bioregs$realm, sep="_")
unique(BioRegion) 

###ADD REGS TO PAM
PAM_REGS<-cbind(BioRegion, PAM_Bioregs)
#write.csv(PAM_REGS, "~/RESULTS/PAM+BIOREGS.csv")

#REMOVE ROWS OF NOT ANALIZED BIOREGS
PAM_REGS<-PAM_REGS[!(PAM_REGS$BioRegion=="TEMP_IndoMalay"|
PAM_REGS$BioRegion=="GRASS_IndoMalay"|
PAM_REGS$BioRegion=="GRASS_Afrotropics"|
PAM_REGS$BioRegion=="DESERT_IndoMalay"| 
PAM_REGS$BioRegion=="DRY_South America"|
PAM_REGS$BioRegion=="MOIST_Eurasia"| 
PAM_REGS$BioRegion=="MOIST_North America"|
PAM_REGS$BioRegion=="MOIST_Australia" |
PAM_REGS$BioRegion=="TUNDRA_NA"|
PAM_REGS$BioRegion=="NA_NA"|
PAM_REGS$BioRegion=="NA_Eurasia"|
PAM_REGS$BioRegion=="TEMP_NA"| 
PAM_REGS$BioRegion=="MED_NA"|
PAM_REGS$BioRegion=="NA_South America"| 
PAM_REGS$BioRegion=="NA_IndoMalay"|
PAM_REGS$BioRegion=="NA_Afrotropics"|
PAM_REGS$BioRegion=="DRY_NA"|
PAM_REGS$BioRegion=="NA_North America"|
PAM_REGS$BioRegion=="DRY_NA"|
PAM_REGS$BioRegion=="MOIST_NA"| 
PAM_REGS$BioRegion=="NA_Madagascar"|
PAM_REGS$BioRegion=="DESERT_NA") ,] 

##NAMES STUDIED BIOREGIONS
Bioregs_names<-unique(PAM_REGS$BioRegion)# MUST BE TOTAL 32 BIOREGS
