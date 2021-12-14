 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VPRM offline for EU sites to compare to FLUXNET observations.
Input: ECMWF met data: temperature, downward shortwave radiation
       or WRF met data: tempeature, downward shortwave radiation
       MODIS data     : EVI, LSWI, Land cover type.
       Flux obs.: flux measurement from FLUXNET2015
"""


from netCDF4 import Dataset
import numpy as np
from OfflineVPRM import julian
from scipy import interpolate
import datetime







def plot_swvl1(sitename, year, lat, lon, tile, input_origin = 'ERA5', level = 1):
    """
    This function prepares the inputs for the vprm_station_for_morris function.
    Input:
        - sitename, str code of the station
        - year, int year to simulate
        - lat, latitude of the station
        - lon, longitude of the station
        - tile, list with the tile indexes where the station data is 
        - input_origin, input option of the meterology data 
    Outputs:
        - evi, numpy array of EVI indexes
        - lswi, numpy array of LSWI indexes
        - evimax, numpy array of max EVI values in the year
        - evimin, numpy array of min EVI values in the year
        - lswimax, numpy array of max LSWI values in the year
        - lswimin, numpy array of min LSWI values in the year
        - temp, numpy array of temperatures
        - rad, numpy array of surface solar radiation downwards    
    """


    if level == 1:
        levs = ['1']
        
        coefs = [1]
        
    
    elif level == 2:
        levs = ['1','2']
        coefs = [7/25, 18/25]

    elif level == 3:
        levs = ['1','2','3']
        coefs = [7/70, 21/70, 42/70]
        
    elif level == 4:
        levs = ['1','2','3','4']
        coefs = [7/150, 21/150, 72/150, 50/150]
    

    
    first_day = datetime.datetime(year,1,1)
    last_day = datetime.datetime(year+1,1,1)
    num_days = (last_day - first_day).days
    time_steps = num_days *48

    input_origin = 'ERA5' #options are 'ERA5' or 'WRF' or 'OBS' in case there are observations otherwise 'ERA5'



    ### Information of input and output
    workpath = './' #Set your work path
    if input_origin == 'ERA5' or 'OBS':
        Metpath = workpath + 'data/ERA5/'
    elif input_origin == 'WRF':
        Metpath = workpath + 'data/WRF_9km/'
    MODISpath = workpath + 'data/MODIS/'



    
    
    fjul_out = (julian(1,1, year)) + np.arange(0,num_days, step = 1/48)

    

    """
    2. INITIALIZATION OF METEOROLOGY
    """
    print('getting met data at ' + sitename + ' station')
            
    SM = np.array([])
    for month in range(12):
        met_nc = Dataset(Metpath + 'ERA5_SM_' + str(month+1).zfill(2) + '_'+ str(year) + '.nc', 'r')
        if month == 0:
            lat_era5 = np.array(met_nc.variables['latitude'])
            lat_era5 = lat_era5[::-1] #Define latitudes from lower values to higher values
            lon_era5 = np.array(met_nc.variables['longitude'])
        time_era5 = np.array(met_nc.variables['time'])
        swvl_era5 = np.array(met_nc.variables['swvl'+levs[0]])
        swvl_era5 = swvl_era5[:,::-1,:]

                
                
        swvl_out = interpolate.interpn((time_era5,lat_era5, lon_era5), swvl_era5, (time_era5, lat, lon), method='linear')
        if swvl_out[0] < 0 :
            swvl_out = interpolate.interpn((time_era5,lat_era5, lon_era5), swvl_era5, (time_era5, lat, lon), method='nearest')
        SM = np.concatenate((SM, swvl_out))
            
        met_nc.close()
    
    

    fjul = (julian(1,1, year)) + np.arange(0,num_days, step = 1/24)



    f = interpolate.interp1d(fjul, SM, fill_value = 'extrapolate')
    Sm = f(fjul_out) #Interpolated to flux time steps
    
    

    
    
    return  Sm
    

    
        

        