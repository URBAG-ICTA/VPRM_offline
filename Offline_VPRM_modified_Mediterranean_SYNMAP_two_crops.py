#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Offline VPRM model for WRF domains

The spatial resolution is the same that in WRF domain.
Outputs daily fluxes, hourly with times from 0:23, meaning time
averaged fluxes (1:00 means fluxes between 1:00 and 2:00)
"""

"""
1. Initialization
"""

from netCDF4 import Dataset
import numpy as np
from os import listdir
from calendar import monthrange
#from datetime import datetime, timedelta
import datetime
import xesmf as xe
from src.OfflineVPRM_modified_SYNMAP_two_crops import *
from src.convert_to_hours_since_1900 import convert_to_hours_since_1990


#Set the period to study
year = 2022
month1 = 1
month2 = 12

mms = np.arange(month1, month2+1)

#Information about the input and output

input_origin = 'ERA5-Land' #options are 'ERA5' or 'WRF' 
wrf_domain = 1 

workpath = './'
tag = 'SYNMAP_mod_two_crops'
datp = workpath + 'VPRMoutput/'
metpath = workpath + 'data/' + input_origin + '/' #Path to the met files
evilswipath = workpath + 'data/MODIS_SYN/'

T_low = [4,4,4,0,2,3,0,0,0,-999]


"""
2. Read MODIS indices and veg cover
"""

fls = listdir(evilswipath)
fls = [x for x, y in zip(fls, [(str(year) in file) for file in fls]) if y == True] 
fls = [x for x, y in zip(fls, [('d0' + str(wrf_domain) in file) for file in fls]) if y == True] 


for fl in fls:
    vn = fl[11:(len(fl)-12)].lower()
    print('Reading ', vn, ' from NetCDF file.' )
    nc = Dataset(evilswipath+fl, 'r')
    if vn == 'evi':
        sday = nc.variables['start_day_of_year']
        lat_out = np.array(nc.variables['lat'])
        ny = len(lat_out)
        lon_out = np.array(nc.variables['lon'])
        nx = len(lat_out[0])
    if vn == 'veg_fra':
        dat = np.array(nc.variables['vegetation_fraction_map'])
        VEG_FRA = dat
    elif vn =='xlat_c':
        dat = np.array(nc.variables['XLAT_C'])
        lat_out_b = dat
    elif vn =='xlong_c':
        dat = np.array(nc.variables['XLONG_C'])
        lon_out_b = dat
    else:
        dat = np.array(nc.variables[vn])
    
    if vn == 'evi':
        EVI = dat
    elif vn == 'evi_max':
        EVI_MAX = dat
    elif vn == 'evi_min':
        EVI_MIN = dat
    elif vn == 'lswi':
        LSWI = dat
    elif vn == 'lswi_max':
        LSWI_MAX = dat
    elif vn == 'lswi_min':
        LSWI_MIN = dat
    
    print(vn, dat.shape)
for iveg in range(len(VEG_FRA)):
    EVI_MAX[iveg] = np.max(EVI[iveg], axis=0)
    LSWI_MAX[iveg] = np.max(LSWI[iveg], axis=0)
    EVI_MIN[iveg] = np.min(EVI[iveg], axis=0)
    LSWI_MIN[iveg] = np.min(LSWI[iveg], axis=0)

"""
3. Read T and DSWF fields form ECMWF
"""

soiltyp = Dataset('data/ERA5-Land/soil_type.nc','r')

field_cap_era5 = np.array(soiltyp.variables['field_cap'][0])
pwp_era5 = np.array(soiltyp.variables['pwp'][0])

soiltyp.close()
field_cap_era5[field_cap_era5 ==0] = 1

firstTF = True
for mm in mms:
    if input_origin == 'ERA5':
        met_nc = Dataset(metpath+'ERA5_'+str(mm).zfill(2)+'_'+str(year)+'.nc', 'r')
        lat_era5 = np.array(met_nc.variables['latitude'])
        lat_era5 = lat_era5[::-1] #Define latitudes from lower values to higher values
        nlat = len(lat_era5)
        lon_era5 = np.array(met_nc.variables['longitude'])
        nlon = len(lon_era5)
        time_era5 = np.array(met_nc.variables['time'])
        temp_era5 = np.array(met_nc.variables['t2m'])
        ssrd_era5 = np.array(met_nc.variables['ssrd'])
        ssrd_era5[ssrd_era5 < -30000] = np.nan
        temp_era5[temp_era5 < -30000] = np.nan
        swvl1_era5 = np.array(met_nc.variables['swvl1'])
        swvl2_era5 = np.array(met_nc.variables['swvl2'])
        

        lon_era5_b = np.linspace(lon_era5[0]-(0.25/2), lon_era5[-1]+(0.25/2), len(lon_era5)+1)
        lat_era5_b = np.linspace(lat_era5[0]-(0.25/2), lat_era5[-1]+(0.25/2), len(lat_era5)+1)
        met_nc.close()
    elif input_origin == 'ERA5-Land':
        met_nc = Dataset(metpath+'ERA5_T2M_'+str(mm).zfill(2)+'_'+str(year)+'.nc', 'r')
        lat_era5 = np.array(met_nc.variables['latitude'])
        lat_era5 = lat_era5[::-1] #Define latitudes from lower values to higher values
        nlat = len(lat_era5)
        lon_era5 = np.array(met_nc.variables['longitude'])
        nlon = len(lon_era5)
        time_era5 = np.array(met_nc.variables['time'])
        temp_era5 = np.array(met_nc.variables['t2m'])
        met_nc.close()
        met_nc = Dataset(metpath+'ERA5_SSRD_'+str(mm).zfill(2)+'_'+str(year)+'.nc', 'r')
        
        ssrd_era5 = np.array(met_nc.variables['ssrd'])
        met_nc.close()
        met_nc = Dataset(metpath+'ERA5_SM2_'+str(mm).zfill(2)+'_'+str(year)+'.nc', 'r')

        swvl2_era5 = np.array(met_nc.variables['swvl2'])

        swvl2_era5[swvl2_era5 < -30000] = np.nan

        #swvl_era5 = 0.5*(swvl1_era5+swvl2_era5)
        swvl_era5 = swvl2_era5
        theta_era5 = (swvl_era5-pwp_era5)/(field_cap_era5-pwp_era5)
        theta_era5[theta_era5 > 1] = 1
        theta_era5[theta_era5 < 0] = 0
        
        lon_era5_b = np.linspace(lon_era5[0]-(0.1/2), lon_era5[-1]+(0.1/2), len(lon_era5)+1)
        lat_era5_b = np.linspace(lat_era5[0]-(0.1/2), lat_era5[-1]+(0.1/2), len(lat_era5)+1)
        met_nc.close()
    
    #Prepare output file in netcdf format
    output = Dataset(datp+'VPRM.'+tag+'_'+str(mm).zfill(2)+'_'+str(year)+'.nc', 'w', format='NETCDF4')
    lat_dim = output.createDimension('lat', ny)     # latitude axis
    lon_dim = output.createDimension('lon', nx)    # longitude axis
    time_dim = output.createDimension('time', None) # unlimited axis (can be appended to).
    vprm_class_dim = output.createDimension('vprm_classes', 8)
    
    lat = output.createVariable('lat', np.float32, ('lat', 'lon'))
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lat[:] = lat_out
    lon = output.createVariable('lon', np.float32, ('lat','lon'))
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    lon[:] = lon_out
    time_out = output.createVariable('time', np.float64, ('time',))
    time_out.units = 'hours since 1900-01-01'
    time_out.long_name = 'time'
    # Define a 3D variable to hold the data
    
    GEE_output = output.createVariable('EBIO_GEE',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
    GEE_output.units = 'mol/km2/h' # degrees Kelvin
    GEE_output.standard_name = 'gross_ecosystem_exchange' # this is a CF standard name
    
    R_output = output.createVariable('EBIO_RES',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
    R_output.units = 'mol/km2/h' # degrees Kelvin
    R_output.standard_name = 'respiration' # this is a CF standard name

    NEE_output = output.createVariable('EBIO_NEE',np.float64,('time','lat','lon')) # note: unlimited dimension is leftmost
    NEE_output.units = 'mol/km2/h' # degrees Kelvin
    NEE_output.standard_name = 'net_ecosystem_exchange' # this is a CF standard name
    """
    GEE_output = output.createVariable('EBIO_GEE',np.float64,('time','vprm_classes','lat','lon')) # note: unlimited dimension is leftmost
    GEE_output.units = 'mol/km2/h' # degrees Kelvin
    GEE_output.standard_name = 'gross_ecosystem_exchange' # this is a CF standard name

    R_output = output.createVariable('EBIO_RES',np.float64,('time','vprm_classes','lat','lon')) # note: unlimited dimension is leftmost
    R_output.units = 'mol/km2/h' # degrees Kelvin
    R_output.standard_name = 'respiration' # this is a CF standard name

    NEE_output = output.createVariable('EBIO_NEE',np.float64,('time','vprm_classes','lat','lon')) # note: unlimited dimension is leftmost
    NEE_output.units = 'mol/km2/h' # degrees Kelvin
    NEE_output.standard_name = 'net_ecosystem_exchange' # this is a CF standard name 
    """
    VEG_FRAC_output = output.createVariable('VEG_FRAC',np.float64,('vprm_classes','lat','lon')) # note: unlimited dimension is leftmost
    VEG_FRAC_output.units = 'unitless' # degrees Kelvin
    VEG_FRAC_output.standard_name = 'vprm_vegetation_fractions'
    VEG_FRAC_output[:] = VEG_FRA
    
    
    dds = monthrange(year, mm)[1]
    for dd in range(1, dds + 1):
        
    
        start_day = datetime.datetime(year= year, month = mm, day = dd, hour = 0)        
        ti = (dd -1)*24
        tf = dd *24
        if input_origin == 'ERA5':
            temp = temp_era5[ti:tf, ::-1]
            ssrd = ssrd_era5[ti:tf, ::-1]
            ssrd = ssrd / 3600
            time = time_era5[ti:tf]
            
            """
            Regrid to VPRM scale
            """ 
            
            grid_in = {'lon': lon_era5, 'lat': lat_era5,
                       'lon_b': lon_era5_b, 'lat_b': lat_era5_b}
            
            grid_out = {'lon': lon_out, 'lat': lat_out,
                        'lon_b': lon_out_b, 'lat_b': lat_out_b}

            regridder = xe.Regridder(grid_in, grid_out, 'bilinear') #Regridding using bilinear interpolation

            temp_out = regridder(temp)
            ssrd_out = regridder(ssrd)
        if input_origin == 'ERA5-Land':
            temp = temp_era5[ti:tf, ::-1]
            ssrd = ssrd_era5[ti:tf, ::-1]
            ssrd[0] = 0
            ssrd[1:] = ssrd[1:] - ssrd[:-1]
            ssrd = ssrd / 3600
            time = time_era5[ti:tf]
            theta = theta_era5[ti:tf,::-1]
            """
            Regrid to VPRM scale
            """

            grid_in = {'lon': lon_era5, 'lat': lat_era5,
                       'lon_b': lon_era5_b, 'lat_b': lat_era5_b}

            grid_out = {'lon': lon_out, 'lat': lat_out,
                        'lon_b': lon_out_b, 'lat_b': lat_out_b}

            regridder = xe.Regridder(grid_in, grid_out, 'nearest_s2d') #Regridding using bilinear interpolation

            temp_out = regridder(temp)
            ssrd_out = regridder(ssrd)
            theta_out = regridder(theta)
        
             
        #create hourly temperatures
        
        """
        Calculate fluxes using pre-optimized parameters
        """
        start_day = datetime.datetime(year= year, month = mm, day = dd)
        start_hrs = 0
        delt = 1 #1 hourly temp and dswf (interpolated), so want hourly fluxes
        evi_times = np.arange(0,np.shape(EVI)[1])*8+year*1000

        if mm == 1 and dd ==1:
            initial_day = True
            GEE_previous = np.empty((24,10,np.shape(temp)[1], np.shape(temp)[2]))
            GEE_previous[:] = np.nan
        else:
            initial_day = False
        flxs, GEE_past_day = offlineVPRM(Temp = temp_out, Rad = ssrd_out, Sm = theta_out, start_mdy = start_day,  
                           start_hrs = start_hrs, evi = EVI, lswi = LSWI, vegFracMap = VEG_FRA,  
                           evi_max = EVI_MAX, evi_min = EVI_MIN, lswi_max = LSWI_MAX, lswi_min = LSWI_MIN, datp = datp,
                           usepar = False, delt = delt, evi_times = evi_times, 
                           tlow = T_low, initial_day = initial_day, GPP_previous = GEE_previous)
      
        
        GEE_previous = np.copy(GEE_past_day)
                    
        GEE_output[ti:tf,:,:] = flxs['gee']
        R_output[ti:tf,:,:] = flxs['resp']
        NEE_output[ti:tf,:,:] = flxs['nee']
        """
        GEE_output[ti:tf,:,:] = flxs['GEEpft']
        R_output[ti:tf,:,:] = flxs['RESpft']
        NEE_output[ti:tf,:,:] = flxs['NEEpft']
        """
        time_out[ti:tf] = time


    output.close()
            
            
        

    


