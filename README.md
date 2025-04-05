# Amph-SpecGradients
Codes and data used for the analyses included in "The latitudinal variation in amphibian speciation rates revisited"  by Garcia Rodriguez et al.

This repository includes:
1. The code for the Spatial Autoregressive model implemented to test for the relatioship between grid-cell mean speciation rates and absolute latitudes
   
2. The code the spatialization and mapping of estimated speciation rates for the overall Amphibian radiation and for each taxonomic Order separately. For this we developed functions (see files ending in "_fun.R") that should be loaded when running the code.
   
4. The code for the Spatial Autoregressive model testing for the relationship between the six predictor variables consider an mean speciation rates at the bioregion level.

   For this analysis we considered bioregions as study unit, using layers from Olson et al (2001, BioScience: https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2) for terrestrial ecoregions and the        reclassification used by Jetz and Fine (2012, PlosOne: https://doi.org/10.1371/journal.pbio.1001292)

As predictors we tested:
   Topo.Complex: Global estimation of metrics of terrain complexity (Terrain Ruggedness Index) based on 30' elevation layer available at www.worldclim.org
   
   Clim.Stab: For this predictor we estimated two metrics, Mean Climate Change Velocity and Mean Climatic Euclidian distances. In both cases we used 11 periods with
   a maximum temporal depth of 3.3 Mya. Paleoclimatic data was obtained from Brown et al (2018, SciData: https://doi.org/10.1038/sdata.2018.254).
   
   NRI: Estimation of Mean Net Relatedness Index for each studied bioregion. For this we used a set of 10 random topologies from the posterior probability 
   distribution of trees from Jetz and Pyron (2018, Nature Ecology and Evolution: https://doi.org/10.1038/s41559-018-0515-5). For each bioregion we first defined species pools, then prunned the trees 
   based on such species composition and estimated NRI.
   
   Estimations for Time Integrated Area, Productivity, and Temperature were obtained from the Supplementary Material in Jetz W, Fine PVA (2012) Global Gradients in Vertebrate Diversity Predicted by Historical 
   Area-Productivity Dynamics and Contemporary Environment. PLoS Biol 10(3): e1001292. https://doi.org/10.1371/journal.pbio.1001292 
