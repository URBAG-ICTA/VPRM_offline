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
import pandas as pd
from os import listdir
from calendar import monthrange
#from datetime import datetime, timedelta
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



def vprm_station_for_morris(sitename, year, iveg, params, EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad):
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
    GEE = lambdaGPP*Tscale*Wscale*Pscale*EVI*Rad/(1 + (Rad/radZero))*(-1)
    
    GEE[GEE > 0] = 0
    GEE = GEE *3600
    
    alpha = params[2]
    beta = params[3]
    Temp[Temp<tlow] = tlow
    
    RSP = Temp*alpha + beta
    
    RSP = RSP *3600
    
    NEE = GEE + RSP
    
    return GEE, RSP, NEE





def preprocess_vprm_for_morris(sitename, year, lat, lon, tile, input_origin = 'ERA5'):
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
    
    first_day = datetime.datetime(year,1,1)
    last_day = datetime.datetime(year+1,1,1)
    num_days = (last_day - first_day).days
    time_steps = num_days *48

    input_origin = 'ERA5' #options are 'ERA5' or 'WRF' or 'OBS' in case there are observations otherwise 'ERA5'
    wrf_domain = 1 #Only needed in case input_origin = 'WRF'


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
    
    OBS_TA = False
    OBS_SW_IN = False
    if input_origin == 'OBS':
        fls = listdir(StationDataPath)
        fls = [x for x, y in zip(fls, [(sitename in file) for file in fls]) if y == True]
        fls = [x for x, y in zip(fls, [(str(year) in file) for file in fls]) if y == True]
        df_obs = pd.read_table(StationDataPath+fls[0], sep=',')
    
        if 'SW_IN' in df_obs.columns:
            df_obs.loc[df_obs['SW_IN'] < -9990,'SW_IN'] = np.nan
            if df_obs['SW_IN'].isna().sum() < 1000:
                RAD = df_obs['SW_IN'].values
                OBS_SW_IN = True

        if 'TA' in df_obs.columns:
            df_obs.loc[df_obs['TA'] < -9990,'TA'] = np.nan
            if df_obs['TA'].isna().sum() < 1000:
                TEMP = df_obs['TA'].values
                OBS_TA = True 

        print(OBS_TA, OBS_SW_IN)
        if OBS_TA:
            nans, x= nan_helper(TEMP)
            TEMP[nans]= np.interp(x(nans), x(~nans), TEMP[~nans])
        if OBS_SW_IN:
            nans, x= nan_helper(RAD)
            RAD[nans]= np.interp(x(nans), x(~nans), RAD[~nans])
    
    

    if input_origin == 'ERA5' or input_origin == 'WRF' or (input_origin == 'OBS' and (not OBS_SW_IN  or not OBS_TA)):
        if not OBS_TA:
            TEMP = np.array([])
        if not OBS_SW_IN:
            RAD = np.array([])
        for month in range(12):
            if input_origin == 'ERA5' or (input_origin == 'OBS' and (not OBS_SW_IN  or not OBS_TA)):
                met_nc = Dataset(Metpath + 'ERA5_' + str(month+1).zfill(2) + '_'+ str(year) + '.nc', 'r')
                if month == 0:
                    lat_era5 = np.array(met_nc.variables['latitude'])
                    lat_era5 = lat_era5[::-1] #Define latitudes from lower values to higher values
                    lon_era5 = np.array(met_nc.variables['longitude'])

                    dlat = abs(lat-lat_era5)
                    dlon = abs(lon -lon_era5)
                    sela = np.where(dlat == np.min(dlat))[0][0]
                    selo = np.where(dlon == np.min(dlon))[0][0]
                    if lat_era5[sela] >= lat:
                        if lon_era5[selo] >= lon:
                            ISW = selo -1
                            JSW = sela -1   
                        else:
                            ISW = selo
                            JSW = sela -1
                    else:
                        if lon_era5[selo] >= lon:
                            ISW = selo -1
                            JSW = sela 
                        else:
                            ISW = selo
                            JSW = sela 
                    factorNE = ((lat - lat_era5[JSW])/(lat_era5[JSW+1] - lat_era5[JSW]))*((lon - lon_era5[ISW])/(lon_era5[ISW+1] - lon_era5[ISW]))
                    factorSE = ((lat_era5[JSW + 1] - lat)/(lat_era5[JSW+1] - lat_era5[JSW]))*((lon - lon_era5[ISW])/(lon_era5[ISW+1] - lon_era5[ISW]))
                    factorSW = ((lat_era5[JSW + 1] - lat)/(lat_era5[JSW+1] - lat_era5[JSW]))*((lon_era5[ISW + 1] - lon)/(lon_era5[ISW+1] - lon_era5[ISW]))
                    factorNW = ((lat - lat_era5[JSW])/(lat_era5[JSW+1] - lat_era5[JSW]))*((lon_era5[ISW + 1] - lon)/(lon_era5[ISW+1] - lon_era5[ISW]))
                    
                if not OBS_TA:
                    temp_era5 = np.array(met_nc.variables['t2m']) - 273.15
                    temp_era5 = temp_era5[:,::-1,:]
                    temp_out = factorNE*temp_era5[:,JSW + 1, ISW + 1] + factorNW*temp_era5[:,JSW + 1, ISW] + factorSE*temp_era5[:,JSW, ISW + 1] + factorSW*temp_era5[:,JSW, ISW]
                    TEMP = np.concatenate((TEMP, temp_out))
                if not OBS_SW_IN:
                    ssrd_era5 = np.array(met_nc.variables['ssrd'])
                    ssrd_era5 = ssrd_era5[:,::-1,:]/3600
                    ssrd_out = factorNE*ssrd_era5[:,JSW + 1, ISW + 1] + factorNW*ssrd_era5[:,JSW + 1, ISW] + factorSE*ssrd_era5[:,JSW, ISW + 1] + factorSW*ssrd_era5[:,JSW, ISW]
                    RAD = np.concatenate((RAD, ssrd_out))
                met_nc.close()
            elif input_origin == 'WRF':
                dds = monthrange(year, month)[1]
                for dd in range(1, dds + 1):
                    met_nc = Dataset(Metpath+'wrfout_d'+str(wrf_domain).zfill(2)+'_'+str(year)+'-'+str(month).zfill(2) + '-' + str(dd).zfill(2) + '_00:00:00','r') 
                    if month == 0 and dd ==1:
                        lat_wrf = np.array(met_nc.variables['XLAT'])
                        lon_wrf = np.array(met_nc.variables['XLONG'])
                        dist = abs(lat-lat_wrf) + abs(lon-lon_wrf)
                        res = np.where(dist == np.min(dist))
                        sela = res[0][0]
                        selo = res[1][0]
                        if lat_wrf[sela] >= lat:
                            if lon_wrf[selo] >= lon:
                                ISW = selo -1
                                JSW = sela -1   
                            else:
                                ISW = selo
                                JSW = sela -1
                        else:
                            if lon_wrf[selo] >= lon:
                                ISW = selo -1
                                JSW = sela 
                            else:
                                ISW = selo
                                JSW = sela 
                        factorNE = ((lat - lat_wrf[JSW])/(lat_wrf[JSW+1] - lat_wrf[JSW]))*((lon - lon_wrf[ISW])/(lon_wrf[ISW+1] - lon_wrf[ISW]))
                        factorSE = ((lat_wrf[JSW + 1] - lat)/(lat_wrf[JSW+1] - lat_wrf[JSW]))*((lon - lon_wrf[ISW])/(lon_wrf[ISW+1] - lon_wrf[ISW]))
                        factorSW = ((lat_wrf[JSW + 1] - lat)/(lat_wrf[JSW+1] - lat_wrf[JSW]))*((lon_wrf[ISW + 1] - lon)/(lon_wrf[ISW+1] - lon_wrf[ISW]))
                        factorNW = ((lat - lat_wrf[JSW])/(lat_wrf[JSW+1] - lat_wrf[JSW]))*((lon_wrf[ISW + 1] - lon)/(lon_wrf[ISW+1] - lon_wrf[ISW]))
                           
                    temp_wrf = np.array(met_nc.variables['T2']) - 273.15
                    temp_out = factorNE*temp_wrf[:,JSW + 1, ISW + 1] + factorNW*temp_wrf[:,JSW + 1, ISW] + factorSE*temp_wrf[:,JSW, ISW + 1] + factorSW*temp_wrf[:,JSW, ISW]
                    TEMP = np.concatenate((TEMP, temp_out))

                    ssrd_wrf = np.array(met_nc.variables['SWDOWN'])
                    ssrd_out = factorNE*ssrd_wrf[:,JSW + 1, ISW + 1] + factorNW*ssrd_wrf[:,JSW + 1, ISW] + factorSE*ssrd_wrf[:,JSW, ISW + 1] + factorSW*ssrd_wrf[:,JSW, ISW]
                    RAD = np.concatenate((RAD, ssrd_out))
                    met_nc.close()
    
    

    fjul = (julian(1,1, year)) + np.arange(0,num_days, step = 1/24)
    if not OBS_TA:
        f = interpolate.interp1d(fjul, TEMP, fill_value = 'extrapolate')
        Temp = f(fjul_out) #Interpolated to flux time steps
    else:
        Temp = TEMP
    if not OBS_SW_IN:
        f = interpolate.interp1d(fjul, RAD, fill_value = 'extrapolate')
        Rad = f(fjul_out) #Interpolated to flux time steps
    else:
        Rad = RAD
    
    
    """
    3. ESTIMATION OF MAX/MIN OF EVI/LSWI VARIATIONS
    """
    EVImax = np.max(EVI)
    EVImin = np.min(EVI)
    LSWImax = np.max(LSWI)
    LSWImin = np.min(LSWI)
    
    
    return EVI, EVImax, EVImin, LSWI, LSWImax, LSWImin, Temp, Rad
    
    
def separate_variable_seasonally(variable, year):
    """
    This function takes a numpy array of semihourly values and separates it in four numpy asrrays one for each season:
    Inputs:
        - variable, numpy array with semihourly values for an entire year
        - year, year to study
    Outputs:
        - variable1, numpy array with semihourly values from Dec-Feb
        - variable2, numpy array with semihourly values from Mar-May
        - variable3, numpy array with semihourly values from Jun-Aug
        - variable4, numpy array with semihourly values from Sep-Oct
    """
    
    s1_i_date = datetime.datetime(year,1,1)
    s2_i_date = datetime.datetime(year,3,1)
    s3_i_date = datetime.datetime(year,6,1)
    s4_i_date = datetime.datetime(year,9,1)
    s5_i_date = datetime.datetime(year,12,1)
    
    start_s1 = 0
    start_s2 = (s2_i_date - s1_i_date).days*48
    start_s3 = (s3_i_date - s1_i_date).days*48
    start_s4 = (s4_i_date - s1_i_date).days*48
    start_s5 = (s5_i_date - s1_i_date).days*48
    
    split_variable = np.split(variable, [start_s2, start_s3, start_s4, start_s5])
    
        

        