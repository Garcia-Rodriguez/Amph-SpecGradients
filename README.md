# Amph-SpecGradients
Compilation of analyses included in "The latitudinal variation in amphibian speciation rates revisited"

Includes:

Study Units
===========
Bioregions: Definition of studied bioregions using layers from Olson et al (2001, BioScience) for terrestrial ecoregions and the reclassification used by Jetz and Fine (2012, PlosOne)

Predictors
==========
Topo.Complex: Global estimation of metrics of terrain complexity (SDaltitude, Roughness, Terrain Ruggedness Index and Topographic Position Index) based on 30' elevation layer available at www.worldclim.org

Clim.Stab: For this predictor we estimated two metrics, Mean Climate Change Velocity and Mean Climatic Euclidian distances. In both cases we used 11 periods with
a maximum temporal depth of 3.3 Mya. Paleoclimatic data was obtained from Brown et al (2018 SciData)

NRI: Estimation of Mean Net Relatedness Index for each studied bioregion. For this we used a set of 10 random topologies from the posterior probability 
distribution of trees from Jetz and Pyron (2018, Nature Ecology and Evolution). For each bioregion we first defined species pools, then prunned the trees 
based on such species composition and estimated NRI.

Estimations for Time Integrated Area, Area-Productivity, and Temperature were obtained from the Supplementary Material in Jetz W, Fine PVA (2012) Global Gradients in Vertebrate Diversity Predicted by Historical Area-Productivity Dynamics and Contemporary Environment. PLoS Biol 10(3): e1001292. https://doi.org/10.1371/journal.pbio.1001292 
