#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Perform a sensitivity Analysis using the Morris Method to obtain the VPRM parameters for a station
"""

import numpy as np
from src.Offline_VPRM_with_soil_threshold_averaged import preprocess_vprm_for_morris
from src.Offline_VPRM_with_soil_threshold_averaged import vprm_station_for_morris
from src.Offline_VPRM_with_soil_threshold_averaged import flatten_list_2d
from sys import exit
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt
import datetime

#
#1. Introduction
#

sitename = 'IT-SR2' #code of the station
name_fluxnet = 'IT-SR2'
StationDataPath = '/home/satellites19/rsegura/Stations_data/FLUXNET/' #path to observed data
iveg =  1 #PFT index of the station
year = 2015 #Year to perform the analysis
input_origin = 'ERA5' #Input origin for the meteorological data
pathout = '/home/users/rsegura/Scripts/plots/'
saturations = [0.403, 0.439, 0.430, 0.520, 0.614, 0.766]
field_capacities = [0.244, 0.347, 0.383, 0.448, 0.541, 0.663]

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


soil_levels = [4,4,4,3,3,3,3,1]

level = soil_levels[iveg-1]


#
#2. Observational data
#


fls = listdir(StationDataPath)
path = [x for x, y in zip(fls, [(name_fluxnet in file) for file in fls]) if y == True]


fls = listdir(StationDataPath+path[0])

fls = [x for x, y in zip(fls, [('FULLSET_HH' in file) for file in fls]) if y == True]


if len(fls) == 0:
    exit('There is not file for this station and this year.')

Flux_file = fls[0]

df = pd.read_csv(StationDataPath+'Stations_info.csv', sep=',')
lat = df.loc[df['Station'] == sitename, 'Latitude'].values[0]
lon = df.loc[df['Station'] == sitename, 'Longitude'].values[0]
tile_h = df.loc[df['Station'] == sitename, 'tile_h'].values[0]
tile_v = df.loc[df['Station'] == sitename, 'tile_v'].values[0]
tile = [tile_h, tile_v]

soiltype = df.loc[df['Station'] == sitename, 'ESDAC'].values[0]
sat = saturations[soiltype-1]


df_obs = pd.read_table(StationDataPath+path[0]+'/'+Flux_file, sep=',')
df_obs['TIMESTAMP_START']= pd.to_datetime(df_obs['TIMESTAMP_START'], format='%Y%m%d%H%M')
df_obs.set_index(df_obs['TIMESTAMP_START'],inplace=True)
df_obs = df_obs.loc[df_obs.index.year == year]

df_obs['GPP_NT_VUT_REF'] = df_obs['GPP_NT_VUT_REF']*-1

df_obs = df_obs[['NEE_VUT_REF','GPP_NT_VUT_REF', 'RECO_NT_VUT_REF' ]]

NEE_obs = df_obs['NEE_VUT_REF'].values
GEE_obs = df_obs['GPP_NT_VUT_REF'].values

#
#3. Run the preprocess_vprm_for_morris function to obtain the MODIS indices and the meteorology
#

EVI, EVImax, EVImin, LSWI, LSWImax, LSWImin, Temp, Rad, Sm = preprocess_vprm_for_morris(sitename = sitename, year = year, lat = lat, lon = lon, tile = tile, input_origin = 'ERA5', level = level)

Sm = Sm/sat
plt.plot(df_obs.index, Sm)

#
#4. Declare problem of the sensitivity analysis.
#





#VPRM parameters from WRF
#PAR0       275.4595, 254.4188, 446.0888, 70.3829, 682.0, 1132.2, 527.9303, 0.00 &
#lambda     0.22577703, 0.21489270, 0.16293380, 0.29311134, 0.1141, 0.08626603, 0.11930965, 0.00, &
#alpha      0.28773167, 0.18056630, 0.24447911, 0.05464646, 0.0049, 0.09231632, 0.1245603, 0.00, &
#beta       -1.09316696, 0.83641734, -0.48669162, -0.12080592, 0.0000, 0.28788863, 0.01743361, 0 /
#Tmin       0.,0.,0.,2.,2.,5.,2.,0.
#Tmax       40,40,40,40,40,40,40,40,
#Topt       20.,20.,20.,20.,20.,22.,18.,0.     

#New parameters  Stocker et al. 2019
#q         -0.5???
#sm_thres   0.1???




lambdaGPP = np.linspace(0.2, 0.45, 11)
radZero =  np.linspace(100, 450, 11)

lambda_prior = 0.22577703
radZero_prior = 275.4595
alpha = 0.28773167
beta = -1.09316696
Tmin = 0
Tmax = 40
Topt = 20
q = np.linspace(-20, -0., 11)
sm_thres = np.linspace(0.4, 0.9, 11)
q_prior = 0
sm_thres_prior = 0.



X = []

for i in range(len(lambdaGPP)):
    for j in range(len(radZero)):
        for k in range(len(q)):
            for l in range(len(sm_thres)):
                X.append([lambdaGPP[i], radZero[j], alpha, beta, Tmin, Tmax, Topt, q[k], sm_thres[l]])

X = np.array(X)


#
#5. Run the model for each set of parameters
#

print('Running VPRM for '+sitename +' station')
Y = []
RMSE = []

for i in range(len(X)):
    GEE, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, X[i], EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad, Sm)
    df_case = df_obs.copy(deep=True)
    df_case['Simulated GEE'] = GEE/3600
    df_case['Simulated RSP'] = RSP/3600
    df_case['Simulated NEE'] = NEE/3600
    mean_NEE = df_case['Simulated NEE'].mean()
    Y.append(mean_NEE)
    df_case = df_case.dropna(axis=0)
    RMSE.append(np.mean((df_case['Simulated GEE'] - df_case['GPP_NT_VUT_REF']) ** 2) ** .5)
    

#
#7. Check best configuration output
#
ind = np.where(RMSE == np.min(RMSE))[0][0]
print(X[ind])
GEE, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, X[ind], EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad, Sm)
df_opt = df_obs.copy(deep=True)

df_opt['NEE_opt'] = NEE/3600
df_opt['GEE_opt'] = GEE/3600


unit = '($\mathrm{\mu mol CO_2/m^2 s}$)'


#
#8. Check prior configuration
#
params = [lambda_prior, radZero_prior, alpha, beta, Tmin, Tmax, Topt, q_prior, sm_thres_prior]

GEE, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, params, EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad, Sm)
#df_nonopt = df_obs.copy()
df_opt['GEE_prior'] = GEE/3600
#df_nonopt = df_nonopt[['NEE_prior']]


#
#9.Plot results
#
begin = [datetime.datetime(year, 1, 1, 0, 0, 0), datetime.datetime(year, 2, 1, 0, 0, 0), datetime.datetime(year, 3, 1, 0, 0, 0),
             datetime.datetime(year, 4, 1, 0, 0, 0), datetime.datetime(year, 5, 1, 0, 0, 0), datetime.datetime(year, 6, 1, 0, 0, 0),
             datetime.datetime(year, 7, 1, 0, 0, 0), datetime.datetime(year, 8, 1, 0, 0, 0), datetime.datetime(year, 9, 1, 0, 0, 0),
             datetime.datetime(year, 10, 1, 0, 0, 0), datetime.datetime(year, 11, 1, 0, 0, 0), datetime.datetime(year, 12, 1, 0, 0, 0)]
end = [datetime.datetime(year, 2, 1, 0, 0, 0), datetime.datetime(year, 3, 1, 0, 0, 0), datetime.datetime(year, 4, 1, 0, 0, 0),
             datetime.datetime(year, 5, 1, 0, 0, 0), datetime.datetime(year, 6, 1, 0, 0, 0), datetime.datetime(year, 7, 1, 0, 0, 0),
             datetime.datetime(year, 8, 1, 0, 0, 0), datetime.datetime(year, 9, 1, 0, 0, 0), datetime.datetime(year, 10, 1, 0, 0, 0),
             datetime.datetime(year, 11, 1, 0, 0, 0), datetime.datetime(year, 12, 1, 0, 0, 0), datetime.datetime(year+1, 1, 1, 0, 0, 0)]
sim_nee = []
sim_gee = []
sim_nonopt_gee = []
obs_nee = []
obs_gee = []
for month in range(len(begin)):

    df_opt_m = df_opt[begin[month]:end[month]]
    df_opt_mean = df_opt_m.groupby([df_opt_m.index.hour]).mean()
    obs_nee.append(df_opt_mean['NEE_VUT_REF'])
    obs_gee.append(df_opt_mean['GPP_NT_VUT_REF'])
    

    sim_nee.append(df_opt_mean.NEE_opt)
    sim_gee.append(df_opt_mean.GEE_opt)

    sim_nonopt_gee.append(df_opt_mean.GEE_prior)
    

sim_nee = flatten_list_2d(sim_nee)
sim_gee = flatten_list_2d(sim_gee)
obs_nee = flatten_list_2d(obs_nee)
obs_gee = flatten_list_2d(obs_gee)
sim_nonopt_gee = flatten_list_2d(sim_nonopt_gee)

time_day = np.arange(0,288)
fig,ax = plt.subplots(figsize=(6.5,3))
plt.subplots_adjust(left=0.13, right=0.81, top=0.9, bottom=0.12)
colors = ['b','m','c']
#ax.plot(time_day, obs_nee, 'k-', linewidth=1, label='NEE OBS')
ax.plot(time_day, obs_gee, 'k-', linewidth=1,label='GEE obs')

#ax.plot(time_day, sim_nee, 'k--', linewidth=1,label='NEE opt')
ax.plot(time_day, sim_gee, 'g--', linewidth=1,label='GEE opt')

ax.plot(time_day, sim_nonopt_gee, 'r--', linewidth=1,label='GEE prior')

legend=ax.legend(loc='center left', shadow=False, fontsize=9, ncol=1, handletextpad=0.5, bbox_to_anchor=(1., 0.512))
ax.set_xlim(0, 288)
ax.set_ylim(-30, 5)
for tt in range(24,288,24):
    ax.axvline(tt, color='y', linewidth=0.8)
ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
plt.ylabel('Flux '+unit, fontsize=10)
major_ticks = np.arange(0, 288, 12)
minor_ticks = np.arange(0, 288, 2)
ax.xaxis.set_ticks(major_ticks)
#ax.set_xticklabels([0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12,0,12])
ax.set_xticklabels(['','Jan\n'+str(year),'','Feb','','Mar','','Apr','','May','','Jun','','Jul','','Aug','','Sep','','Oct','','Nov','','Dec'])
ax.xaxis.set_ticks(minor_ticks, minor = True)
ax.xaxis.set_tick_params(which='major', labelsize=8)
ax.set_title(sitename, fontsize=11)
plt.savefig(pathout+sitename+'_diurnal_'+str(year)+'.png',dpi=300)
plt.show()
plt.close()

df_opt = df_opt[['GPP_NT_VUT_REF', 'GEE_opt', 'GEE_prior']]


initial_date = datetime.datetime(year, 1,1)
final_date = datetime.datetime(year+1, 1,1)

delta = final_date- initial_date       # as timedelta
dates = []
GPP_obs = []
GPP_opt = []
GPP_prior = []
for i in range(delta.days):
    day = initial_date + datetime.timedelta(days=i)
    dates.append(day)
    short = df_opt.loc[(df_opt.index >= day) & (df_opt.index < day +datetime.timedelta(days=1))]
    GPP_obs.append(short['GPP_NT_VUT_REF'].mean())
    GPP_opt.append(short['GEE_opt'].mean())
    GPP_prior.append(short['GEE_prior'].mean())

fig,ax = plt.subplots(figsize=(6.5,3))
plt.subplots_adjust(left=0.13, right=0.81, top=0.9, bottom=0.12)
colors = ['b','m','c']
#ax.plot(time_day, obs_nee, 'k-', linewidth=1, label='NEE OBS')
ax.plot(dates, GPP_obs, 'k-', linewidth=1,label='GEE obs')


ax.plot(dates, GPP_opt, 'g--', linewidth=1,label='GEE opt')

ax.plot(dates, GPP_prior, 'r--', linewidth=1,label='GEE prior')

legend=ax.legend(loc='center left', shadow=False, fontsize=9, ncol=1, handletextpad=0.5, bbox_to_anchor=(1., 0.512))
ax.set_ylim(-30, 5)

plt.ylabel('Flux '+unit, fontsize=10)

ax.set_title(sitename, fontsize=11)
plt.savefig(pathout+sitename+'_daily_'+str(year)+'.png',dpi=300)
plt.show()
plt.close()

params_dict = {'lambdaGPP':lambdaGPP, 'radZero':radZero,  'q':q, 'sm_thres':sm_thres}

df_params = pd.DataFrame()

df_params['lambdaGPP'] =  X[:,0]
df_params['radZero'] = X[:,1]
df_params['q'] = X[:,7]
df_params['sm_thres'] = X[:,8]
df_params['RMSE'] = RMSE

plt.plot(np.arange(len(RMSE)), RMSE,'bo', np.where(RMSE==np.min(RMSE))[0][0], np.min(RMSE), 'ro' , markersize=1)
plt.show()
plt.close()




for param in params_dict:
    RMSE_for_param = []
    for i in range(len(params_dict[param])):
        val = params_dict[param][i]
        RMSE_for_param.append(df_params.loc[df_params[param] == val, 'RMSE'].mean())
        
    
    plt.plot(params_dict[param], RMSE_for_param)
    plt.xlabel(param)
    plt.ylabel('RMSE '+unit)
    plt.show()
    plt.close()
    
    



