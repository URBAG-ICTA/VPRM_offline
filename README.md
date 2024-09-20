#VPRM_offline
Vegetation Photosynthesis and Respiration Model offline written in Python based on the VPRM model 
described at Mahadevan et al. (2008) and the forecasting version described in Chen et al. (2020). This code has been used in the study "Heat and drought events alter biogenic capacity to balance CO2 budget in south-western Europe" (DOI: 10.22541/essoar.171136782.20335714/v1).

The code is composed of two programs:
Offline_VPRM_Domain.py - for estimating the biogenic CO2 fluxes for a 2D domain.

Offline_VPRM_point.py - for estimating the biogenic CO2 fluxes for a station point. 

The input data is composed of:
-Meteorological data stored in the data/ directory providing from the ECMWF ERA5 (data/ERA5) database or from WRF simulations (data/WRF_9km) of the same domain. 
-MODIS indexes (EVI and LSWI) in the data/MODIS directory alongside to vegetation cover
-VPRM parameters optimized (data/VPRMparameters).


References:
Chen, J., Gerbig, C., Marshall, J., & Uwe Totsche, K. (2020). Short-Term forecasting of regional biospheric CO2 fluxes in Europe using a light-use-efficiency model (VPRM, MPI-BGC version 1.2). Geoscientific Model Development, 13(9), 4091–4106. https://doi.org/10.5194/gmd-13-4091-2020

Mahadevan, P., Wofsy, S. C., Matross, D. M., Xiao, X., Dunn, A. L., Lin, J. C., … Gottlieb, E. W. (2008). A satellite-based biosphere parameterization for net ecosystem CO2 exchange: Vegetation Photosynthesis and Respiration Model (VPRM). Global Biogeochemical Cycles, 22(2). https://doi.org/10.1029/2006GB002735
