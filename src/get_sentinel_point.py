#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#script to read modis indices at specific location (e.g. flux site)
#from preprocessed, lowess fltered data at 1 km resolution
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_sentinel_point(year, SENTINELpath, sitename, res = '_150m'):
    #year: desired year
    #SENTINELpath: path to sentinel data
    #res: label to indicate the correct buffer zone around the station
    

    #returns:
    #doy:  day of year 
    #evi,lswi: MODIS indices
    

    df = pd.read_csv(SENTINELpath+sitename+res+'.csv', sep = ',')
    df['Date'] = pd.to_datetime(df['Date'], format='%Y-%m-%d')
    df = df.sort_values(by='Date')
    df.set_index(df['Date'], inplace=True)
    df = df.loc[df.index.year == year]

    evi = df['EVI_mean'].values
    lswi = df['LSWI_mean'].values
    
    doy = np.arange(0,len(evi)) +1
    
    return [doy, evi, lswi]
    
    