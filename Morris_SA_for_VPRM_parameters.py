#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform a sensitivity Analysis using the Morris Method to obtain the VPRM parameters for a station
"""

from SALib.sample.morris import sample
from SALib.analyze.morris import analyze
import numpy as np
from Offline_VPRM_for_Morris import preprocess_vprm_for_morris
from Offline_VPRM_for_Morris import vprm_station_for_morris
from Offline_VPRM_for_Morris import flatten_list_2d
from sys import exit
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt
import datetime
"""
1. Introduction
"""

sitename = 'FRFBn' #code of the station
StationDataPath = '/home/users/rsegura/Stations_data/' #path to observed data
iveg = 1 #PFT index of the station
year = 2015 #Year to perform the analysis
input_origin = 'ERA5' #Input origin for the meteorological data
pathout = '/home/users/rsegura/Scripts/plots/'
#PFT for the Iberian Peninsula - Stations
#1. Evergreen needleleaf Forest - ESES1,FRFBn,FRLBr
#2. Broadleaf deciduous Forest - PTEsp, 
#3. Mixed Forest
#4. Mixed shrubland/grassland
#5. Savanna
#6. Dryland cropland and pasture
#7. Grassland
#8. Others
#9. Cropland/woodland mosaic


"""
2. Observational data
"""

fls = listdir(StationDataPath)
fls = [x for x, y in zip(fls, [(sitename in file) for file in fls]) if y == True]
fls = [x for x, y in zip(fls, [(str(year) in file) for file in fls]) if y == True]

if len(fls) == 0:
    exit('There is not file for this station and this year.')

Flux_file = fls[0]

df = pd.read_csv(StationDataPath+'Stations_info.csv', sep=',')
lat = df.loc[df['Station'] == sitename, 'Latitude'].values[0]
lon = df.loc[df['Station'] == sitename, 'Longitude'].values[0]
tile_h = df.loc[df['Station'] == sitename, 'tile_h'].values[0]
tile_v = df.loc[df['Station'] == sitename, 'tile_v'].values[0]
tile = [tile_h, tile_v]

df_obs = pd.read_table(StationDataPath+Flux_file, sep=',')
df_obs['TIMESTAMP_START']= pd.to_datetime(df_obs['TIMESTAMP_START'], format='%Y%m%d%H%M')
if 'NEE_PI' in df_obs.columns:
    label = 'NEE_PI'
else:
    label = 'FC'
df_obs.loc[df_obs[label] < -9990, label] = np.nan
df_obs[label] = df_obs[label]*3600
df_obs.set_index('TIMESTAMP_START', inplace=True)
df_obs = df_obs[[label]]

NEE_obs = df_obs[label].values

"""
3. Run the preprocess_vprm_for_morris function to obtain the MODIS indices and the meteorology
"""

EVI, EVImax, EVImin, LSWI, LSWImax, LSWImin, Temp, Rad = preprocess_vprm_for_morris(sitename = sitename, year = year, lat = lat, lon = lon, tile = tile, input_origin = 'ERA5')

"""
4. Declare problem of the sensitivity analysis.
"""

problem = {
    'num_vars': 7,
    'names': ['lambdaGPP', 'radZero', 'alpha', 'beta', 'Tmin', 'Tmax', 'Topt'],
    'bounds': [[0.2, 0.45],
               [100, 250],
               [0.1, 0.2],
               [0.0, 1],
               [-1, 1],
               [39, 41],
               [19, 21]]
}

X = sample(problem, 1000, num_levels=6)

"""
5. Run the model for each set of parameters
"""

Y = []
RMSE = []

for i in range(len(X)):
    GEE, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, X[i], EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad)
    df_case = df_obs.copy(deep=True)
    df_case['Simulated GEE'] = GEE
    df_case['Simulated RSP'] = RSP
    df_case['Simulated NEE'] = NEE
    mean_NEE = df_case['Simulated NEE'].mean()
    Y.append(mean_NEE)
    df_case = df_case.dropna(axis=0)
    RMSE.append(np.mean((df_case['Simulated NEE'] - df_case[label]) ** 2) ** .5)
    
"""
6. Analyze  Morris indices
"""
Y = np.array(Y)

Si = analyze(problem, X, Y, conf_level=0.95,
             print_to_console=True, num_levels=6)


"""
7. Check best configuration output
"""
ind = np.where(RMSE == np.min(RMSE))[0][0]
print(X[ind])
GEE, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, X[ind], EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad)
df_opt = df_obs.copy(deep=True)
wrf_convert = 24*44/1000000
df_opt['NEE'] = NEE*wrf_convert
df_opt[label] = df_opt[label]*wrf_convert
unit = '($\mathrm{g_{CO_2} m^{-2} day^{-1}}$)'
"""
df_obs = df_obs[[label]]
df_obs[label] = NEE_obs*wrf_convert

df_opt = df_opt[['NEE']]
"""
"""
NEE = NEE*wrf_convert
NEE_obs = NEE_obs*wrf_convert
GEE = GEE*wrf_convert
RSP = RSP*wrf_convert
unit = '($\mathrm{g_{CO_2} m^{-2} day^{-1}}$)'
plt.plot(df_opt.index, NEE_obs, 'k-', label='NEE obs')
plt.plot(df_opt.index, GEE, 'g-', label='GEE')
plt.plot(df_opt.index, NEE, 'b-', label='NEE')
plt.plot(df_opt.index, RSP, 'r-', label='RSP')
plt.ylabel('Flux '+unit, fontsize=10)
plt.legend(loc='best', ncol= 4)
plt.title(sitename)
plt.savefig(pathout+sitename+'_optimized.png', dpi=100)
plt.close()
#plt.plot(df_final.index, NEE_obs, 'k-', label='NEE obs')
#plt.plot(df_final.index, NEE_obs, 'k-', label='NEE obs')
"""


params = [0.3084, 270.2, 0.1797, 0.8800, 0, 40, 20]

GEE, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, params, EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad)
#df_nonopt = df_obs.copy()
df_opt['NEE_prior'] = NEE*wrf_convert
#df_nonopt = df_nonopt[['NEE_prior']]


begin = [datetime.datetime(year, 1, 1, 0, 0, 0), datetime.datetime(year, 2, 1, 0, 0, 0), datetime.datetime(year, 3, 1, 0, 0, 0),
             datetime.datetime(year, 4, 1, 0, 0, 0), datetime.datetime(year, 5, 1, 0, 0, 0), datetime.datetime(year, 6, 1, 0, 0, 0),
             datetime.datetime(year, 7, 1, 0, 0, 0), datetime.datetime(year, 8, 1, 0, 0, 0), datetime.datetime(year, 9, 1, 0, 0, 0),
             datetime.datetime(year, 10, 1, 0, 0, 0), datetime.datetime(year, 11, 1, 0, 0, 0), datetime.datetime(year, 12, 1, 0, 0, 0)]
end = [datetime.datetime(year, 1, 31, 23, 0, 0), datetime.datetime(year, 2, 28, 23, 0, 0), datetime.datetime(year, 3, 31, 23, 0, 0),
           datetime.datetime(year, 4, 30, 23, 0, 0), datetime.datetime(year, 5, 31, 23, 0, 0), datetime.datetime(year, 6, 30, 23, 0, 0),
           datetime.datetime(year, 7, 31, 23, 0, 0), datetime.datetime(year, 8, 31, 23, 0, 0), datetime.datetime(year, 9, 30, 23, 0, 0),
           datetime.datetime(year, 10, 31, 23, 0, 0), datetime.datetime(year, 11, 30, 23, 0, 0), datetime.datetime(year, 12, 31, 23, 0, 0)]

sim_nee = []
sim_nonopt_nee = []
obs_nee = []
for month in range(len(begin)):

    df_opt_m = df_opt[begin[month]:end[month]]
    df_opt_mean = df_opt_m.groupby([df_opt_m.index.hour]).mean()
    obs_nee.append(df_opt_mean[label])

    sim_nee.append(df_opt_mean.NEE)

    sim_nonopt_nee.append(df_opt_mean.NEE_prior)
    

sim_nee = flatten_list_2d(sim_nee)
obs_nee = flatten_list_2d(obs_nee)
sim_nonopt_nee = flatten_list_2d(sim_nonopt_nee)

time_day = np.arange(0,288)
fig,ax = plt.subplots(figsize=(6.5,3))
plt.subplots_adjust(left=0.13, right=0.81, top=0.9, bottom=0.12)
colors = ['b','m','c']
ax.plot(time_day, obs_nee, linewidth=1, color='k', label='OBS')
ax.plot(time_day, sim_nee, linewidth=1, color = 'b', label='VPRM_opt')
ax.plot(time_day, sim_nonopt_nee, linewidth=1, color = 'r', label='VPRM_prior')

legend=ax.legend(loc='center left', shadow=False, fontsize=9, ncol=1, handletextpad=0.5, bbox_to_anchor=(1., 0.512))
ax.set_xlim(0, 288)
ax.set_ylim(-80, 40)
for tt in range(24,288,24):
    ax.axvline(tt, color='y', linewidth=0.8)
ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
plt.ylabel('Flux '+unit, fontsize=10)
major_ticks = np.arange(0, 288, 12)
minor_ticks = np.arange(0, 288, 2)
ax.xaxis.set_ticks(major_ticks)
#ax.set_xticklabels([0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12])
ax.set_xticklabels(['','Jan\n2015','','Feb','','Mar','','Apr','','May','','Jun','','Jul','','Aug','','Sep','','Oct','','Nov','','Dec'])
ax.xaxis.set_ticks(minor_ticks, minor = True)
ax.xaxis.set_tick_params(which='major', labelsize=8)
ax.set_title(sitename, fontsize=11)
plt.savefig(pathout+sitename+'_diurnal_'+str(year)+'.png',dpi=300)
