#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#script to read modis indices at specific location (e.g. flux site)
#from preprocessed, lowess fltered data at 1 km resolution
"""
import pyreadr
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

def get_modis_point(year, lat, lon, tile, MODISpath, sitename):
    #year: desired year
    #lat, lon: cartesian coordinates for point extraction
    #configfile: config.r file (absolute path) used for VPRM preprocessing

    #returns:
    #doy:  day of year (start of 8-day period during which MODIS indices were acquired)
    #evi,lswi: MODIS indices
    
    #  source(configfile)
    
    #find tile id 
    #  tile <- system(paste(ldope,"bin/tile_id -proj=SIN -lon=",lon," -lat=",lat,sep=""),intern=TRUE)
    #above gives e.g. "Tile ID (h, v): 11 9" "Processing done ! "  
    #  tile <- as.numeric(strsplit(tile,split=" ")[[1]][5:6]) # tile[1] = h; tile[2]= v
    
    
    
    tilenm = 'h' + str(tile[0]).zfill(2) + 'v' + str(tile[1]).zfill(2) + '.' + str(year) + '.ier1km'
    
    
    
    TLUT = pyreadr.read_r(MODISpath+'TLUT.RData')

    TLUT = TLUT['TLUT']
    TLUT.columns = ['H', 'V', 'LAT.LL', 'LON.LL', 'LAT.UR', 'LON.UR']
    tile_id = TLUT.loc[(TLUT['H'] == tile[0]) & (TLUT['V'] == tile[1] )]
    
    lat_syn = np.linspace(24,56, num=3841) 
    lat_syn = lat_syn[:-1] #projinfo$lat_syn obtained from projinfo e.g. is lower left corner coordinates of 1 km pixels
    lon_syn = np.linspace(-21, 32, num = 6361)
    lon_syn = lon_syn[:-1]
    

    lat_id = lat_syn[lat_syn -tile_id['LAT.LL'].values[0] + 0.001 >= 0] #add small offset to account for rounding error
    
    lon_id = lon_syn[lon_syn -tile_id['LON.LL'].values[0] + 0.001 >= 0] #add small offset to account for rounding error

    #now lat_id and lon_id are lat-lon starting with lower left corner of tile
    #now get tile pointers
    dlat = abs(lat - lat_id)
    dlon = abs(lon - lon_id)
    #select pixel for which lower left corner has minimum positive distance from desired location
    #dlat[dlat< 0] = 100
    #dlon[dlon < 0] = 100
    sela = np.where(dlat == np.min(dlat))
    selo = np.where(dlon == np.min(dlon))
    
    if sitename == 'IT-Ro2':
        sela = 287
        selo = 1430
    print(lat_id[sela], lon_id[selo])
    #extract point data from ncdf files for given tile
    
    evif = Dataset(MODISpath+'evi.smooth.'+tilenm+'.nc')
    evi = np.array(evif.variables['indices'])
    evi = evi[:,selo,sela]
    evi = np.reshape(evi, (len(evi)))
    evif.close()
    evi[evi < -99999] = np.nan
    
    lswif = Dataset(MODISpath+'lswi.smooth.'+tilenm+'.nc')
    lswi = np.array(lswif.variables['indices'])
    lswi = lswi[:,selo,sela]
    lswi = np.reshape(lswi, (len(lswi)))
    lswif.close()
    lswi[lswi < -99999] = np.nan
    
    doy = np.arange(0,len(evi))*8 +1
    
    return [doy, evi, lswi]
    
    
