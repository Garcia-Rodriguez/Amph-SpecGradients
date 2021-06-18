# Amph-SpecGradients
Compilation of analyses included in "Causes and consequences of the latitudinal variation in amphibian speciation rates"

Includes:

Bioregions: Definition of studied bioregions using layers from Olson et al (2001, BioScience) for terrestrial ecoregions
and the reclassification used by Jetz and Fine (2015, PlosOne)

Topo.Complex: Global estimation of metrics of terrain complexity (SDaltitude, Roughness, Terrain Ruggedness Index and Topographic Position Index)  
based on 30' elevation layer available at www.worldclim.org

Clim.Stab: For this predictor we estimated two metrics, Mean Climate Change Velocity and Mean Climatic Euclidian distances. In both cases we used 11 periods with
a maximum temporal depth of 3.3 Mya. Paleoclimatic data was obtained from Brown et al (2018 SciData)

NRI: Estimation of Mean Net Relatedness Index for each studied bioregion. For this we used a set of 10 random topologies from the posterior probability 
distribution of trees from Jetz and Pyron (2018, Nature Ecology and Evolution). For each bioregion we first defined species pools, then prunned the trees 
based on such species composition and estimated NRI.

OLS_SARs: Regressions between mean Speciation and Absolute Latitude using an assemblage based approach. We conducted this analyses for the entire amphibian radiation and
also at the Order level
