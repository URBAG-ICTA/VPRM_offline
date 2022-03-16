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
from src.get_sentinel_point import get_sentinel_point
import src.WriteVPRMConstants as WriteVPRMConstants
from src.OfflineVPRM import julian
from scipy import interpolate
import matplotlib.pyplot as plt
import datetime

"""
0. Initialization
"""

year = 2018

snames = ['FR-Bil', 'FR-Aur'] #Stations to simulate

first_day = datetime.datetime(year,1,1)
last_day = datetime.datetime(year+1,1,1)
num_days = (last_day - first_day).days
time_steps = num_days *48

input_origin = 'ERA5' #options are 'ERA5' or 'WRF' or 'OBS' in case there are observations otherwise 'ERA5'
wrf_domain = 1 #Only needed in case input_origin = 'WRF'
satellite_origin = 'SENTINEL2'

tag = 'SENTINEL2_150m' #tag for simulation output

### Information of input and output
workpath = './' #Set your work path
outpath = workpath + 'VPRMoutput/' #Path to output file
stations_file = '/data/co2flux/common/rsegura/DATA/FLUXNET/Stations_info.csv'
StationDataPath = '/data/co2flux/common/rsegura/DATA/FLUXNET/'

if input_origin == 'ERA5' or 'OBS':
    Metpath = workpath + 'data/ERA5/'
elif input_origin == 'WRF':
    Metpath = workpath + 'data/WRF_9km/'

if satellite_origin == 'MODIS':
    MODISpath = workpath + 'data/MODIS/'
elif satellite_origin == 'SENTINEL2':
    SENTINELpath = workpath + 'data/SENTINEL2/'
    res = '_150m'

###Other settings
vprm_par_name = 'vprmopt.EU2007.local.par.csv'
parapath = workpath + 'data/VPRMparameters/'

#get parameters, reduce to 8 vegetation classes (no 4 evergreens needed, one enough)
VegClass = 8
vprmConstants = WriteVPRMConstants.WriteVPRMConstants(outdir = outpath, nveg = 8)
#### 8 vegetation classes: evergreen,deciduous,mixed forest,shrubland, savanna, cropland, grassland, others######
T_low = [4,0,2,3,0,0,0,-999]

igbp_dict = {'ENF':0,'EBF':0, 'DNF':1, 'DBF':1, 'MF':2, 'CSH':3, 'OSH':3, 'WS':4, 'SAV':4, 'GRA':6, 'CRO':5, }


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


"""
1. READ INPUT TOWER DATA
"""
stations_df = pd.read_csv(stations_file, sep = ',')

#snames = stations_df['Station']
vegfra_types = ['Evergreen','Decid', 'Mixfrst','Shrub','Savan','Crop','Grass','Other']
vprm_class = ["Evergreen Forest","Deciduous Forest","Mixed Forest","Shrubland","Savannas","Cropland","Grassland","Others"]

stations_df.set_index(stations_df['Station'], inplace=True)

"""
for station in stations_df['Station'].unique():
    fls = listdir(StationDataPath)
    fls = [x for x, y in zip(fls, [(station in file) for file in fls]) if y == True]
    fls = [x for x, y in zip(fls, [(str(year) in file) for file in fls]) if y == True]
    if len(fls) == 0:
        stations = stations.drop(stations[stations['Station'] == station].index)


snames = stations['Station'].unique()
print(snames)
stations.set_index(stations['Station'], inplace=True)
nsites = len(stations)

#snames = np.delete(snames, 6)
"""

output_df = pd.DataFrame()
### Loop over sites
for sitename in snames:
    print('Start processing at ' + sitename + ' station')
    
    lat = stations_df.loc[sitename, 'Latitude']
    lon = stations_df.loc[sitename, 'Longitude']
    tile = [stations_df.loc[sitename, 'tile_h'], stations_df.loc[sitename, 'tile_v']]
    veg_type = stations_df.loc[sitename, 'IGBP']
    iveg = igbp_dict[veg_type]
    #iveg = veg_type - 1
    """
    2. ESTIMATION OF HOURLY EVI/LSWI VARIATIONS
    """
    if satellite_origin == 'MODIS':
        print('getting MODIS for ', sitename, ' station ', year)
        data = get_modis_point(year=year, lat = lat, lon = lon, tile = tile, MODISpath = MODISpath, sitename=sitename)
    elif satellite_origin == 'SENTINEL2':
        data = get_sentinel_point(year=year, SENTINELpath= SENTINELpath, res = res, sitename= sitename)
    
    fjul = (julian(1,1, year)-1) + data[0]
    
    
    fjul_out = (julian(1,1, year)) + np.arange(0,365, step = 1/48)
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
    
    
    #LSWI = LSWI*(2)
    """
    3. INITIALIZATION OF METEOROLOGY
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
                met_nc = Dataset(Metpath+'ERA5_'+str(month+1).zfill(2)+'_' + str(year)+'.nc', 'r')
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
    
    

    fjul = (julian(1,1, year)) + np.arange(0,365, step = 1/24)
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
    fig, ax = plt.subplots(figsize = (10,10))
    if OBS_SW_IN:
        ax.plot(RAD_OBS, 'k-', label='OBS')
    ax.plot(Rad, 'b-', label = 'ERA5')
    plt.legend()
    plt.show()
    plt.close()
    """
    
    """
    4. ESTIMATION OF MAX/MIN OF EVI/LSWI VARIATIONS
    """
    EVImax = np.max(EVI)
    EVImin = np.min(EVI)
    LSWImax = np.max(LSWI)
    LSWImin = np.min(LSWI)
    
    """
    5. ESTIMATION OF SCALAR EFFECTS ON GPP PRODUCTS
    """
    TMIN = vprmConstants.loc[iveg, 'tempMin']
    TMAX = vprmConstants.loc[iveg, 'tempMax']
    TOPT = vprmConstants.loc[iveg, 'tempOpt']
    
    Tscale = ((Temp - TMIN)*(Temp-TMAX))/(((Temp-TMIN)*(Temp-TMAX))-((Temp-TOPT)*(Temp-TOPT)))
    Tscale[Tscale < 0] = 0
    #modification for so-called "xeric systems", comprising shrublands and grasslands
    #these have different dependencies on ground water.

    
    
    if iveg in [3, 6]:
        Wscale = (LSWI - LSWImin)/(LSWImax - LSWImin)
    else:
        Wscale = (1 + LSWI)/(1 + LSWImax)

    Wscale[Wscale <0] = 0
    Pscale = (1 + LSWI)/2

    if iveg == 0:
        Pscale[:] = 1
        
    if iveg in [1, 2, 3, 5, 7]:
        threshmark = 0.55
        evithresh = EVImin + (threshmark*(EVImax-EVImin))
        phenologyselect = np.where(EVI[:] > evithresh)
        Pscale[phenologyselect] = 1
    Pscale[Pscale < 0] = 0
    
    #Pscale[:] = 1
        #by default, grasslands and savannas never have pScale=1
    """
    6. HOURLY GEE (mol/km2/hr) ESTIMATIONS
    """
    lambdaGPP = vprmConstants.loc[iveg, 'lambdaGPP.sw']
    radZero = vprmConstants.loc[iveg, 'swradZero']
    GEE = lambdaGPP*Tscale*Wscale*Pscale*EVI*Rad/(1 + (Rad/radZero))*(-1)
    
    GEE[GEE > 0] = 0
    GEE = GEE *3600
    
    """
    7. HOURLY RSP (mol/km2/hr)
    """
    
    alpha = vprmConstants.loc[iveg, 'alphaResp']
    beta = vprmConstants.loc[iveg, 'intResp']
    
    tlow = T_low[iveg]
    Temp[Temp< tlow] = tlow
    RSP = Temp*alpha + beta
    
    RSP = RSP *3600
    
    NEE = GEE + RSP
    
    sep = np.modf(fjul_out)
    base = datetime.datetime(year=1960, month=1, day=1)
    date_list = [base + datetime.timedelta(days=x) +datetime.timedelta(hours=y) for x, y in zip(sep[1], sep[0]*24)]
    
    if sitename == snames[0]:
        output_df['Times'] = date_list
        
    output_df[sitename + '_GEE'] = GEE
    output_df[sitename + '_RSP'] = RSP
    output_df[sitename + '_NEE'] = NEE
    output_df[sitename + '_EVI'] = EVI
    output_df[sitename + '_LSWI'] = LSWI


output_df.to_csv(outpath+'VPRM.'+tag+'_'+str(year)+'.csv', index = False, header=True)
    
    