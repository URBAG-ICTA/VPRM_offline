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
from src.get_modis_point import get_modis_point
from OfflineVPRM import julian
from scipy import interpolate
import datetime



def flatten_list_2d(l):
    """
      Flattens a 2D list to 1D list

      l = [[1, 2], [3, 4], [5]]

      flatten_list_2d(l)

      [1, 2, 3, 4, 5]
    """
    flat = [x for sublist in l for x in sublist]

    return flat

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]



def vprm_station_for_morris(sitename, year, iveg, params, EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad, Sm):
    """
    This function estimates the mean GPP, RESP and NEE annually and seasonally for a year and a 
    selected ICOS station using the VPRM model.
    Input:
        - sitename, str code of the station
        - year, int year to simulate
        - params, list of floats with the VPRM parameters [[lambda1, par01, alpha1, beta1, Tmin1, Tmax1, Topt1],[lambda2, par02, alpha2, beta2, Tmin2, Tmax2, Topt2],...]
        - evi, numpy array of EVI indexes
        - lswi, numpy array of LSWI indexes
        - evimax, float max EVI value in the year
        - evimin, float min EVI value in the year
        - lswimax, float max LSWI value in the year
        - lswimin, float min LSWI value in the year
        - temp, numpy array of temperatures
        - rad, numpy array of surface solar radiation downwards
    Output:
        - List with GPPan, GPPs1, GPPs2, GPPs3, GPPs4, RESPan, RESPs1, RESPs2, RESPs3, RESPs4, NEEan, NEEs1, NEEs2, NEEs3, NEEs4, RMSEan, RMSEs1, RMSEs2, RMSEs3, RMSEs4
        where an means: annually mean
        s1: first season mean Dec-Feb
        s2: second season mean Mar-May
        s3: third season mean Jun-Aug
        s4: fourth season mean Sep-Nov
    """
    iveg = iveg -1
    TMIN = params[4]
    TMAX = params[5]
    TOPT = params[6]
    T_low = [4,0,2,3,0,0,0,-999,0]
    tlow = T_low[iveg]
    
    Tscale = ((Temp - TMIN)*(Temp-TMAX))/(((Temp-TMIN)*(Temp-TMAX))-((Temp-TOPT)*(Temp-TOPT)))
    Tscale[Tscale < 0] = 0
    #modification for so-called "xeric systems", comprising shrublands and grasslands
    #these have different dependencies on ground water.
    
    if iveg in [3, 6]:
        Wscale = (LSWI - LSWImin)/(LSWImax - LSWImin)
    else:
        Wscale = (1 + LSWI)/(1 + LSWImax)

    Wscale[Wscale < 0] = 0
    Pscale = (1 + LSWI)/2

    if iveg == 0:
        Pscale[:] = 1
        
    if iveg in [1, 2, 3, 5, 7, 8]:
        threshmark = 0.55
        evithresh = EVImin + (threshmark*(EVImax-EVImin))
        phenologyselect = np.where(EVI[:] > evithresh)
        Pscale[phenologyselect] = 1
    #by default, grasslands and savannas never have pScale=1
    Pscale[Pscale < 0] = 0

    lambdaGPP = params[0]
    radZero = params[1]
    
    q = params[7]
    sm_thres = params[8]
    
    betascale = np.zeros(shape=np.shape(Sm))
    betascale[:] = 1
    
    stressSelect = np.where(Sm <= sm_thres)
    if len(stressSelect[0]) > 0:
        betascale[stressSelect] = q*(Sm[stressSelect]-sm_thres)**2 + np.ones(shape=np.shape(Sm[stressSelect]))
    
    betascale[betascale<0] = 0
    #plt.plot(betascale)
    GEE = betascale*lambdaGPP*Tscale*Wscale*Pscale*EVI*Rad/(1 + (Rad/radZero))*(-1)
    
    GEE[GEE > 0] = 0
    GEE = GEE *3600
    
    alpha = params[2]
    beta = params[3]
    Temp[Temp<tlow] = tlow
    
    RSP = Temp*alpha + beta
    
    RSP = RSP *3600
    
    NEE = GEE + RSP
    
    return GEE, RSP, NEE





def preprocess_vprm_for_morris(sitename, year, lat, lon, tile, input_origin = 'ERA5', level = 1):
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


    """
    1. ESTIMATION OF HALFHOURLY EVI/LSWI VARIATIONS
    """
    print('getting MODIS for ', sitename, ' station ', year)
    data = get_modis_point(year=year, lat = lat, lon = lon, tile = tile, MODISpath = MODISpath)
    
    
    fjul = (julian(1,1, year)-1) + data[0]
    
    
    fjul_out = (julian(1,1, year)) + np.arange(0,num_days, step = 1/48)
    data = [fjul, data[1], data[2]]
    
    if np.isnan(data[1]).all():
        data[1][:] = 0.5
        print("all EVI missing for",sitename)
    if np.isnan(data[2]).all():
        data[2][:] = 0.5
        print("all LSWI missing for",sitename)
    EVI = np.empty(shape=(time_steps))
    LSWI = np.empty(shape=(time_steps))
    
    
    f = interpolate.interp1d(fjul, data[1], fill_value = 'extrapolate')
    EVI[:] = f(fjul_out)
    f = interpolate.interp1d(fjul, data[2], fill_value = 'extrapolate')
    LSWI[:] = f(fjul_out)
    

    """
    2. INITIALIZATION OF METEOROLOGY
    """
    print('getting met data at ' + sitename + ' station')
    
  
    
    TEMP = np.array([])
    RAD = np.array([])
 
    for month in range(12):
        met_nc = Dataset(Metpath + 'ERA5_' + str(month+1).zfill(2) + '_'+ str(year) + '.nc', 'r')
        if month == 0:
            lat_era5 = np.array(met_nc.variables['latitude'])
            lat_era5 = lat_era5[::-1] #Define latitudes from lower values to higher values
            lon_era5 = np.array(met_nc.variables['longitude'])
        time_era5 = np.array(met_nc.variables['time'])
        temp_era5 = np.array(met_nc.variables['t2m']) - 273.15
        temp_era5 = temp_era5[:,::-1,:]
        temp_out = interpolate.interpn((time_era5,lat_era5, lon_era5), temp_era5, (time_era5, lat, lon))
        TEMP = np.concatenate((TEMP, temp_out))

        ssrd_era5 = np.array(met_nc.variables['ssrd'])
        ssrd_era5 = ssrd_era5[:,::-1,:]/3600
        ssrd_out = interpolate.interpn((time_era5,lat_era5, lon_era5), ssrd_era5, (time_era5, lat, lon))
        RAD = np.concatenate((RAD, ssrd_out))
        met_nc.close()
            
    SM = np.array([])
    for month in range(12):
        met_nc = Dataset(Metpath + 'ERA5_SM_' + str(month+1).zfill(2) + '_'+ str(year) + '.nc', 'r')
        if month == 0:
            lat_era5 = np.array(met_nc.variables['latitude'])
            lat_era5 = lat_era5[::-1] #Define latitudes from lower values to higher values
            lon_era5 = np.array(met_nc.variables['longitude'])
        time_era5 = np.array(met_nc.variables['time'])
        swvl_era5 = np.array(met_nc.variables['swvl'+levs[0]])
        swvl_era5 = swvl_era5[:,::-1,:]*coefs[0]
                
        if level >1:
            for k in range(1,len(levs)):
                ext_swvl_era5 = np.array(met_nc.variables['swvl'+levs[k]])
                ext_swvl_era5 = ext_swvl_era5[:,::-1,:]*coefs[k]
                swvl_era5 = swvl_era5 + ext_swvl_era5
                
                
        swvl_out = interpolate.interpn((time_era5,lat_era5, lon_era5), swvl_era5, (time_era5, lat, lon), method='linear')
        if swvl_out[0] < 0 :
            swvl_out = interpolate.interpn((time_era5,lat_era5, lon_era5), swvl_era5, (time_era5, lat, lon), method='nearest')
        SM = np.concatenate((SM, swvl_out))
            
        met_nc.close()
    
    

    fjul = (julian(1,1, year)) + np.arange(0,num_days, step = 1/24)

    f = interpolate.interp1d(fjul, TEMP, fill_value = 'extrapolate')
    Temp = f(fjul_out) #Interpolated to flux time steps


    f = interpolate.interp1d(fjul, RAD, fill_value = 'extrapolate')
    Rad = f(fjul_out) #Interpolated to flux time steps

    f = interpolate.interp1d(fjul, SM, fill_value = 'extrapolate')
    Sm = f(fjul_out) #Interpolated to flux time steps
    
    
    """
    3. ESTIMATION OF MAX/MIN OF EVI/LSWI VARIATIONS
    """
    EVImax = np.max(EVI)
    EVImin = np.min(EVI)
    LSWImax = np.max(LSWI)
    LSWImin = np.min(LSWI)
    
    
    return EVI, EVImax, EVImin, LSWI, LSWImax, LSWImin, Temp, Rad, Sm
    

    
        

        