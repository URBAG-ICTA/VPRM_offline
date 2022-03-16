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

iveg = 5
method = 'NT'

#stations = ['FR-Bil', 'IT-Lav', 'IT-Ren', 'IT-SR2', 'IT-SRo'] #ENF
#train_sites = ['FR-Bil', 'FR-Bil','FR-Bil', 'IT-Lav','IT-Lav','IT-Lav', 'IT-Ren','IT-Ren','IT-Ren',  'IT-SR2','IT-SR2','IT-SRo'] #ENF
#train_years = [2015, 2016, 2017, 2011, 2012, 2014, 2011, 2012, 2013, 2013,  2017, 2011] #ENF
#stations = ['FR-Pue', 'IT-Cp2'] #EBF
#stations = ['IT-CA1', 'IT-CA3', 'IT-Col', 'IT-Isp', 'IT-Ro2'] #DBF
stations = ['ES-LJu', 'ES-Amo', 'IT-Noe'] #OSH
#years = [2011, 2012, 2013, 2014, 2015, 2016, 2017]

stations =[ 'ES-Abr', 'ES-LM1']

"""
train_sites = ['FR-Pue', 'FR-Pue', 'FR-Pue', 'IT-Cp2'] #EBF
train_years = [2011, 2013, 2014, 2013] #EBF
"""
train_sites = ['ES-Abr', 'ES-LM1', 'ES-LM1']
train_years = [2016, 2015, 2017]

#train_sites = ['IT-Noe', 'IT-Noe', 'ES-LJu', 'ES-LJu', 'ES-LJu', 'ES-Amo', 'ES-Amo', 'ES-Amo']
#train_years = [2012, 2013, 2009, 2011, 2012, 2010, 2011, 2012]

#train_sites = [ 'IT-Col', 'IT-Isp', 'IT-Ro2']
#train_years = [ 2011, 2014, 2012]




StationDataPath = '/data/co2flux/common/rsegura/DATA/FLUXNET/' #path to observed data
input_origin = 'ERA5' #Input origin for the meteorological data
pathout = '/data/co2flux/common/rsegura/Scripts/plots/VPRMopt/'

wilting_points = [0.059, 0.151, 0.133, 0.279, 0.335, 0.267]
saturations = [0.403, 0.439, 0.430, 0.520, 0.614, 0.766]
field_capacities  = [0.244, 0.347, 0.383, 0.448, 0.541, 0.663]
soil_levels = [4,4,4,3,3,3,3,1]
#level = soil_levels[iveg-1]
#levels = [200, 200, 200, 200, 200] #ENF
#levels = [200,200] #EBF
#levels = [200, 200, 200, 200, 200] #DBF
#levels = [100, 100, 100] #OSH
levels = [100, 100] #SAV
#saturations = [0.3952, 0.4073, 0.3795, 0.43, 0.43] #ENF
#saturations = [0.43, 0.4366] #EBF
#saturations = [0.439, 0.439, 0.439, 0.46, 0.439] #DBF
#saturations = [0.439, 0.439, 0.439] #OSH
saturations = [0.439, 0.439]

#PFT for the Iberian Peninsula - Stations
#1. Evergreen needleleaf Forest
#2. Broadleaf deciduous Forest
#3. Mixed Forest
#4. Shrubland
#5. Savannah
#6. Cropland
#7. Grassland
#8. Others
#9. Evergreen broadleaf forest

#VPRM parameters from WRF
#PAR0       275.4595, 254.4188, 446.0888, 70.3829, 682.0, 1132.2, 527.9303, 0.00 &
#lambda     0.22577703, 0.21489270, 0.16293380, 0.29311134, 0.1141, 0.08626603, 0.11930965, 0.00, &
#alpha      0.28773167, 0.18056630, 0.24447911, 0.05464646, 0.0049, 0.09231632, 0.1245603, 0.00, &
#beta       -1.09316696, 0.83641734, -0.48669162, -0.12080592, 0.0000, 0.28788863, 0.01743361, 0 /
#Tmin       0.,0.,0.,2.,2.,5.,2.,0.
#Tmax       40,40,40,40,40,40,40,40,
#Topt       20.,20.,20.,20.,20.,22.,18.,0.     

#New parameters  Stocker et al. 2019
#q         -5 
#sm_thres   0.8


lambdaGPP = np.linspace(0.2, 0.35, 11)
lambda_prior = [0.22577703, 0.21489270, 0.16293380, 0.29311134, 0.1141, 0.08626603, 0.11930965, 0.00, 0.22577703]
lambda_prior = lambda_prior[iveg-1]

radZero =  np.linspace(400, 800, 11)
radZero_prior = [275.4595, 254.4188, 446.0888, 70.3829, 682.0, 1132.2, 527.9303, 0.00, 275.4595]
radZero_prior = radZero_prior[iveg-1]

alpha = [0.28773167, 0.18056630, 0.24447911, 0.05464646, 0.0049, 0.09231632, 0.1245603, 0.00, 0.28773167]
alpha = alpha[iveg-1]

beta = [-1.09316696, 0.83641734, -0.48669162, -0.12080592, 0.0000, 0.28788863, 0.01743361, 0, 1.09316696]
beta = beta[iveg-1]

Tmin = [0.,0.,0.,2.,2.,5.,2.,0.,0.]
Tmin = Tmin[iveg-1]

Tmax = [40,40,40,40,40,40,40,40,40]
Tmax = Tmax[iveg-1]

Topt = [20.,20.,20.,20.,20.,22.,18.,0.,20.]
Topt = Topt[iveg-1]

q = np.linspace(-30, -10., 11)
q_prior = 0

sm_thres = np.linspace(0.3, 0.8, 11)

sm_thres_prior = 0.

unit = '($\mathrm{\mu mol CO_2/m^2 s}$)'

X = []
lambdas = []
radZeros = []
qs = []
sm_thress = [] 

for i in range(len(lambdaGPP)):
    for j in range(len(radZero)):
        for k in range(len(q)):
            for l in range(len(sm_thres)):
                X.append([lambdaGPP[i], radZero[j], alpha, beta, Tmin, Tmax, Topt, q[k], sm_thres[l]])
                lambdas.append(lambdaGPP[i])
                radZeros.append(radZero[j])
                qs.append(q[k])
                sm_thress.append(sm_thres[l])

X = np.array(X)


df = pd.read_csv(StationDataPath+'Stations_info.csv', sep=',')




siteyear_df = pd.DataFrame()
siteyear_df['lambdaGPP'] = lambdas
siteyear_df['radZero'] = radZeros
siteyear_df['q'] = qs
siteyear_df['sm_thres'] = sm_thress
siteyear_labels = []

points = 0
for j, sitename in enumerate(train_sites):
    print(sitename)
    lat = df.loc[df['Station'] == sitename, 'Latitude'].values[0]
    lon = df.loc[df['Station'] == sitename, 'Longitude'].values[0]
    tile_h = df.loc[df['Station'] == sitename, 'tile_h'].values[0]
    tile_v = df.loc[df['Station'] == sitename, 'tile_v'].values[0]
    tile = [tile_h, tile_v]
    soiltype = df.loc[df['Station'] == sitename, 'ESDAC'].values[0]
    k = stations.index(sitename)
    sat = saturations[k]
    #wp = wilting_points[soiltype-1]
    level = levels[k]
    fls = listdir(StationDataPath)
    path = [x for x, y in zip(fls, [(sitename in file) for file in fls]) if y == True]

    fls = listdir(StationDataPath+path[0])
    fls = [x for x, y in zip(fls, [('FULLSET_HH' in file) for file in fls]) if y == True]
    Flux_file = fls[0]
    
    df_obs = pd.read_table(StationDataPath+path[0]+'/'+Flux_file, sep=',')
    df_obs['TIMESTAMP_START']= pd.to_datetime(df_obs['TIMESTAMP_START'], format='%Y%m%d%H%M')
    df_obs.set_index(df_obs['TIMESTAMP_START'],inplace=True)
    df_obs.loc[df_obs['GPP_'+method+'_VUT_REF'] < -9990, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs.loc[df_obs['NEE_VUT_REF_QC'] == 3, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs.loc[df_obs['NEE_VUT_REF_QC'] == 2, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs['GPP_'+method+'_VUT_REF'] = df_obs['GPP_'+method+'_VUT_REF']*-1
    df_obs.loc[df_obs['GPP_'+method+'_VUT_REF'] > 0, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs = df_obs[['NEE_VUT_REF','GPP_'+method+'_VUT_REF', 'RECO_'+method+'_VUT_REF' ]]
    
    year = train_years[j]
    df_year = df_obs.loc[df_obs.index.year == year]
        

        
    print(year)
        
    EVI, EVImax, EVImin, LSWI, LSWImax, LSWImin, Temp, Rad, Sm = preprocess_vprm_for_morris(sitename = sitename, year = year, lat = lat, lon = lon, tile = tile, input_origin = 'ERA5', level = level)

    Sm = (Sm)/sat
    plt.plot(df_year.index, Sm)
    plt.ylabel('Soil moisture fraction')
    plt.title(sitename)
    plt.show()
    plt.close()
    
    print('VPRM running for '+sitename+' '+str(year))
    SSE = []
    for i in range(len(X)):
        GEE, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, X[i], EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad, Sm)
        df_case = df_year.copy(deep=True)
        df_case['Simulated GEE'] = GEE/3600
        df_case = df_case.dropna(axis=0)
        SSE.append(np.sum((df_case['Simulated GEE'] - df_case['GPP_'+method+'_VUT_REF']) ** 2))
    points += len(df_case)
    siteyear_df[sitename+'_'+str(year)] = np.array(SSE)
    siteyear_labels.append(sitename+'_'+str(year))
     
        
siteyear_df['SSE'] =  siteyear_df[siteyear_labels].sum(axis=1)
siteyear_df['MSE'] = siteyear_df['SSE']/points
siteyear_df['RMSE'] = siteyear_df['MSE']**0.5

plt.plot(siteyear_df['RMSE'])
plt.show()
plt.close()
    
minimum_RMSE = siteyear_df[siteyear_df['RMSE'] == siteyear_df['RMSE'].min()]
print('RMSE optimized ', minimum_RMSE)
lambda_opt = minimum_RMSE['lambdaGPP'].values[0]
radZero_opt = minimum_RMSE['radZero'].values[0]
q_opt = minimum_RMSE['q'].values[0]
sm_thres_opt = minimum_RMSE['sm_thres'].values[0]

"""
lambda_opt = 0.4
radZero_opt = 175
q_opt = -30
sm_thres_opt = 0.4
"""
SSE = 0

params_opt = [lambda_opt, radZero_opt, alpha, beta, Tmin, Tmax, Topt, q_opt, sm_thres_opt]

print(params_opt)
params_prior = [lambda_prior, radZero_prior, alpha, beta, Tmin, Tmax, Topt, q_prior, sm_thres_prior]

print("Simulating with optimized configuration")


for j, sitename in enumerate(train_sites):
    print(sitename)
    lat = df.loc[df['Station'] == sitename, 'Latitude'].values[0]
    lon = df.loc[df['Station'] == sitename, 'Longitude'].values[0]
    tile_h = df.loc[df['Station'] == sitename, 'tile_h'].values[0]
    tile_v = df.loc[df['Station'] == sitename, 'tile_v'].values[0]
    tile = [tile_h, tile_v]
    soiltype = df.loc[df['Station'] == sitename, 'ESDAC'].values[0]
    k = stations.index(sitename)
    sat = saturations[k]
    #wp = wilting_points[soiltype-1]
    level = levels[k]
    fls = listdir(StationDataPath)
    path = [x for x, y in zip(fls, [(sitename in file) for file in fls]) if y == True]

    fls = listdir(StationDataPath+path[0])
    fls = [x for x, y in zip(fls, [('FULLSET_HH' in file) for file in fls]) if y == True]
    Flux_file = fls[0]
    
    df_obs = pd.read_table(StationDataPath+path[0]+'/'+Flux_file, sep=',')
    df_obs['TIMESTAMP_START']= pd.to_datetime(df_obs['TIMESTAMP_START'], format='%Y%m%d%H%M')
    df_obs.set_index(df_obs['TIMESTAMP_START'],inplace=True)
    df_obs.loc[df_obs['GPP_'+method+'_VUT_REF'] < -9990, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs.loc[df_obs['NEE_VUT_REF_QC'] == 3, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs.loc[df_obs['NEE_VUT_REF_QC'] == 2, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs['GPP_'+method+'_VUT_REF'] = df_obs['GPP_'+method+'_VUT_REF']*-1
    df_obs.loc[df_obs['GPP_'+method+'_VUT_REF'] > 0, 'GPP_'+method+'_VUT_REF'] = np.nan
    df_obs = df_obs[['NEE_VUT_REF','GPP_'+method+'_VUT_REF', 'RECO_'+method+'_VUT_REF' ]]
    
    year = train_years[j]
    df_year = df_obs.loc[df_obs.index.year == year]
        
    if len(df_year) == 0:
        continue
    if sitename == 'FR-Bil' and year == 2014:
        continue
        
    print(year)
        
    EVI, EVImax, EVImin, LSWI, LSWImax, LSWImin, Temp, Rad, Sm = preprocess_vprm_for_morris(sitename = sitename, year = year, lat = lat, lon = lon, tile = tile, input_origin = 'ERA5', level = level)

    Sm = (Sm)/sat
        
        

    GEE_opt, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, params_opt, EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad, Sm)
    GEE_prior, RSP, NEE = vprm_station_for_morris(sitename, year, iveg, params_prior, EVI, LSWI, EVImax, EVImin, LSWImax, LSWImin, Temp, Rad, Sm)
    
    df_case = df_year.copy(deep=True)
    df_case['Simulated GEE'] = GEE_prior/3600
    df_case = df_case.dropna(axis=0)
    SSE += (np.sum((df_case['Simulated GEE'] - df_case['GPP_'+method+'_VUT_REF']) ** 2))
    
    
    GEE_opt = GEE_opt/3600
    GEE_prior = GEE_prior/3600

    df_opt = df_year.copy(deep=True)
    df_opt['GEE_obs'] = df_opt['GPP_'+method+'_VUT_REF']
    df_opt['GEE_opt'] = GEE_opt
    df_opt['GEE_prior'] = GEE_prior
        
    begin = [datetime.datetime(year, 1, 1, 0, 0, 0), datetime.datetime(year, 2, 1, 0, 0, 0), datetime.datetime(year, 3, 1, 0, 0, 0),
             datetime.datetime(year, 4, 1, 0, 0, 0), datetime.datetime(year, 5, 1, 0, 0, 0), datetime.datetime(year, 6, 1, 0, 0, 0),
             datetime.datetime(year, 7, 1, 0, 0, 0), datetime.datetime(year, 8, 1, 0, 0, 0), datetime.datetime(year, 9, 1, 0, 0, 0),
             datetime.datetime(year, 10, 1, 0, 0, 0), datetime.datetime(year, 11, 1, 0, 0, 0), datetime.datetime(year, 12, 1, 0, 0, 0)]
    end = [datetime.datetime(year, 2, 1, 0, 0, 0), datetime.datetime(year, 3, 1, 0, 0, 0), datetime.datetime(year, 4, 1, 0, 0, 0),
           datetime.datetime(year, 5, 1, 0, 0, 0), datetime.datetime(year, 6, 1, 0, 0, 0), datetime.datetime(year, 7, 1, 0, 0, 0),
           datetime.datetime(year, 8, 1, 0, 0, 0), datetime.datetime(year, 9, 1, 0, 0, 0), datetime.datetime(year, 10, 1, 0, 0, 0),
           datetime.datetime(year, 11, 1, 0, 0, 0), datetime.datetime(year, 12, 1, 0, 0, 0), datetime.datetime(year+1, 1, 1, 0, 0, 0)]


        
    opt_gee = []
    obs_gee = []
    prior_gee = []
        
    for month in range(len(begin)):
        df_opt_m = df_opt[begin[month]:end[month]]
        df_opt_mean = df_opt_m.groupby([df_opt_m.index.hour]).mean()
        obs_gee.append(df_opt_mean['GEE_obs'])
        opt_gee.append(df_opt_mean['GEE_opt'])
        prior_gee.append(df_opt_mean['GEE_prior'])
            
    opt_gee = flatten_list_2d(opt_gee)
    obs_gee = flatten_list_2d(obs_gee)
    prior_gee = flatten_list_2d(prior_gee)


    time_day = np.arange(0,288)
    fig,ax = plt.subplots(figsize=(6.5,3))
    plt.subplots_adjust(left=0.13, right=0.81, top=0.9, bottom=0.12)
    colors = ['b','m','c']
    ax.plot(time_day, obs_gee, 'k-', linewidth=1,label='GEE obs')
    ax.plot(time_day, opt_gee, 'g--', linewidth=1,label='GEE opt')
    ax.plot(time_day, prior_gee, 'r--', linewidth=1,label='GEE prior')
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

    df_opt = df_opt[['GEE_obs', 'GEE_opt', 'GEE_prior']]


    initial_date = datetime.datetime(year, 1,1)
    final_date = datetime.datetime(year+1, 1,1)
        
    delta = final_date- initial_date       # as timedelta
    dates = []
    GEE_obs = []
    GEE_opt = []
    GEE_prior = []
    for i in range(delta.days):
        day = initial_date + datetime.timedelta(days=i)
        dates.append(day)
        short = df_opt.loc[(df_opt.index >= day) & (df_opt.index < day +datetime.timedelta(days=1))]
        GEE_obs.append(short['GEE_obs'].mean())
        GEE_opt.append(short['GEE_opt'].mean())
        GEE_prior.append(short['GEE_prior'].mean())

    fig,ax = plt.subplots(figsize=(6.5,3))
    plt.subplots_adjust(left=0.13, right=0.81, top=0.9, bottom=0.12)
    ax.plot(dates, GEE_obs, 'k-', linewidth=1,label='GEE obs')
    ax.plot(dates, GEE_opt, 'g--', linewidth=1,label='GEE opt')
    ax.plot(dates, GEE_prior, 'r--', linewidth=1,label='GEE prior')
    ax2 = ax.twinx()
    ax2.plot(df_year.index, Sm, color ='blue' , linewidth=1)
    ax2.set_ylabel("Fraction of soil moisture to field capacity",color="blue",fontsize=10)
        
    legend=ax.legend(loc='center left', shadow=False, fontsize=9, ncol=3, handletextpad=0.5, bbox_to_anchor=(0.3, 1.05))
    ax.set_ylim(-20, 5)
    ax.set_ylabel('Flux '+unit, fontsize=10)
    ax.set_title(sitename, fontsize=11, loc='left' )
    plt.savefig(pathout+sitename+'_daily_'+str(year)+'.png',dpi=300)
    plt.show()
    plt.close()
        
    plt.plot(df_year.index, EVI)
    plt.title(sitename)
    plt.ylabel('EVI')
    plt.show()
    plt.close()
        
    plt.plot(df_year.index, LSWI)
    plt.title(sitename)
    plt.ylabel('LSWI')
    plt.show()
    plt.close()
        
    plt.plot(df_year.index, Temp)
    plt.title(sitename)
    plt.ylabel('Temp')
    plt.show()
    plt.close()

MSE = SSE/points
RMSE = MSE**0.5
print('RMSE prior  ', RMSE)