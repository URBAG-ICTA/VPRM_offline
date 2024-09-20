#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VPRM code to obtain the R, GEE and NEE fluxes
"""
import numpy as np
import datetime 
import src.WriteVPRMConstants_modified_SYNMAP_two_crops as WriteVPRMConstants
from sys import exit
import copy
from scipy.signal import convolve

def julian(m, d, y):
    df = datetime.datetime(year = y, month = m, day = d)
    di = datetime.datetime(year = 1960, month = 1, day = 1)
    diff = df - di
    return diff.days


def extract_wrf_times_from_evi(evi, start_mdy, start_hrs, nhrs, evi_times, delt):
    #to extract hourly interpolated evi and lswi
    
    jday0 = julian(start_mdy.month, start_mdy.day, start_mdy.year)
    doyuse = evi_times
    doyuse = doyuse - start_mdy.year*1000
    fjdays = jday0 - julian(1,1, start_mdy.year) + (start_hrs/24) + np.arange(0, ((nhrs-1)/delt)+1)/(24/delt)
    
    evih = np.empty(( np.shape(evi)[0], nhrs, np.shape(evi)[2], np.shape(evi)[3]))
    for i in range(len(fjdays)):
        fjday = fjdays[i]
        evi_id1 = np.where(abs(fjday-doyuse) == min(abs(fjday-doyuse)))[0][0]
        if fjday > doyuse[evi_id1]:
            evi_id2 = evi_id1 + 1
        else:
            evi_id2 = evi_id1 - 1
        if evi_id2 < 1: 
            evi_id2 = evi_id1
        if evi_id2 > len(doyuse) -1: 
            evi_id2 = evi_id1
        evi1 = evi[:,evi_id1,:,:]
        evi2 = evi[:,evi_id2,:,:]
        if evi_id2 != evi_id1:
            evih[:,i,:,:] = (evi1*(doyuse[evi_id2] - fjday) + evi2*(fjday-doyuse[evi_id1])) /(doyuse[evi_id2]-doyuse[evi_id1])
        else:
            evih[:,i,:,:] = evi1
    
    return evih

def compute_daily_GPP(gpp, gpp_previous):
    conv = np.ones(25, dtype=gpp.dtype)
    conv[0] = 0
    gpp_together = np.concatenate([gpp_previous, gpp]) 
    GPP_new = np.apply_along_axis(lambda m: np.convolve(m, conv, mode='valid'), axis=0, arr=gpp_together)
    #GPP_new = convolve(gpp, conv, mode='valid')
    #initial = [GPP_new[0] for x in range(24)]
    #GPP_new = np.concatenate((np.array(initial), GPP_new))
    return GPP_new*(-1)/86400    

def compute_daily_GPP_initial(gpp):
    GPP_new = np.sum(gpp, axis= 0)
    GPP_new = GPP_new[np.newaxis,:,:]
    GPP_new = np.repeat(GPP_new, 24, axis=0)
    return GPP_new*(-1)/86400




def offlineVPRM(Temp, Rad, Sm, start_mdy, start_hrs,  evi, lswi, vegFracMap,
                      evi_max, evi_min, lswi_max, lswi_min, datp, usepar=False,
                      delt=1, evi_times=np.linspace(2005041,2005305,8),
                      tlow= None,returnJackList=False,
                      initial_day = False, GPP_previous = None):
    
    nhrs = np.shape(Temp)[0]
    evi = extract_wrf_times_from_evi(evi = evi, start_mdy = start_mdy, start_hrs = start_hrs, nhrs = nhrs, evi_times = evi_times, delt = delt)
    lswi = extract_wrf_times_from_evi(evi = lswi, start_mdy = start_mdy, start_hrs = start_hrs, nhrs = nhrs, evi_times = evi_times, delt = delt)
    #Temp = np.transpose(Temp, axes=(1,2,0))
    #Rad = np.transpose(Rad, axes=(1,2,0))

    
    month = start_mdy.month
    nVegClass = 8
    
    if month in [6,7,8]:
        vprmConstants = WriteVPRMConstants.WriteVPRMConstants_summer(outdir = datp, nveg = 8)
    else:
        vprmConstants = WriteVPRMConstants.WriteVPRMConstants(outdir = datp, nveg = 8)
    """
    vprmConstants = WriteVPRMConstants.WriteVPRMConstants(outdir = datp, nveg = 10)
    """
    checkdim = nVegClass

    
    if (Temp.size != Rad.size) or (Temp.size != evi.size/checkdim) or (Temp.size != lswi.size/checkdim):
        exit('1. Temperature and radiation must have same dimensions.')
    
    Temp = Temp - 273.15
    
    gppTot = np.zeros(np.shape(Temp))
    respTot = np.zeros(np.shape(Temp))
    Tscale = np.zeros(np.shape(evi))
    Wscale = np.zeros(np.shape(evi))    
    SMscale = np.zeros(np.shape(evi))
    NEEpft = np.zeros(np.shape(evi))
    GEEpft = np.zeros(np.shape(evi))
    RESpft = np.zeros(np.shape(evi))
    R0_out = np.zeros(np.shape(evi))
    RT_out = np.zeros(np.shape(evi))
    SMstress = np.zeros(np.shape(evi))
    #print(np.shape(np.repeat(evi_max[np.newaxis], np.shape(evi)[-1], axis=0)))

    evim = evi
    lswim = lswi
    evi_max = np.repeat(evi_max[:,np.newaxis,:,:], nhrs, axis=1)
    evi_min = np.repeat(evi_min[:,np.newaxis,:,:], nhrs, axis=1)
    lswi_max = np.repeat(lswi_max[:,np.newaxis,:,:], nhrs, axis=1)
    lswi_min = np.repeat(lswi_min[:,np.newaxis,:,:], nhrs, axis=1)
    evi_maxm = evi_max
    evi_minm = evi_min
    lswi_maxm = lswi_max
    lswi_minm = lswi_min
    
        
    GPP_next_day = []

    for k in range(nVegClass-1):
        #Assume last reclass category is water/ice/concrete, which doesn't contribute to first order flux
        #Determine total influence for all vegetation from this uberclass


            
        vegFrack = vegFracMap[k]
        vegFrac = np.repeat(vegFrack[np.newaxis,:,:], nhrs, axis=0)
        evi = evim[k]
        lswi = lswim[k]
        evi_max = evi_maxm[k]
        evi_min = evi_minm[k]
        lswi_max = lswi_maxm[k]
        lswi_min = lswi_minm[k]
        
        #get EVI parameters

        lambdaGPP = vprmConstants.loc[k, 'lambdaGPP.sw']
        radZero = vprmConstants.loc[k, 'swradZero']
        Q = vprmConstants.loc[k, 'Q']
        theta_thres = vprmConstants.loc[k, 'theta_thres']
        E0 = vprmConstants.loc[k, 'E0']
        R0 = vprmConstants.loc[k, 'R0']
        k2 = vprmConstants.loc[k, 'k2']
        k1 = vprmConstants.loc[k, 'k1']
        gamma = vprmConstants.loc[k, 'gamma'] 
        
        #temperature parameters, assume calculations in degrees C
        tempOpt = vprmConstants.loc[k, 'tempOpt']
        tempMax = vprmConstants.loc[k, 'tempMax']
        tempMin = vprmConstants.loc[k, 'tempMin']
        
        tempScalar = ((Temp - tempMin)*(Temp-tempMax))/(((Temp-tempMin)*(Temp-tempMax))-((Temp-tempOpt)*(Temp-tempOpt)))
        tempScalar[np.where(Temp > tempMax)] = 0
        tempScalar[np.where(Temp < tempMin)] = 0
        
        #SIMPLIFICATION: Have EVI max/min, LSWI max values ONLY by vegetation map, Pscale is done by comparison of EVI to EVI max
        #Check on eviMax and eviMin and lswiMax--these should be the limits
        #lswi = lswi + 0.5*lswi #Increase 50% lswi
        
        lswi[np.where(lswi > lswi_max)] = lswi_max[np.where(lswi > lswi_max)]
        lswi[np.where(lswi < lswi_min)] = lswi_min[np.where(lswi < lswi_min)]
        
        evi[np.where(evi > evi_max)] = evi_max[np.where(evi > evi_max)]
        evi[np.where(evi < evi_min)] = evi_min[np.where(evi < evi_min)]
        
        #modification for so-called "xeric systems", comprising shrublands and grasslands
        #these have different dependencies on ground water.
        
        #if (k in [3,6]):
        #    wScalar = (lswi - lswi_min)/(lswi_max - lswi_min)
        #else:
        wScalar = (1 + lswi)/(1 + lswi_max)
        
        #wScalar[:] = 1 #Cancel water stress factor
        pScalar = (1 + lswi)/2
        
        if  (k ==0):
            pScalar[:] = 1
        
        if (k in [1,2,3,5]):  #if decid, mixed, shrub, crop, or other
            threshmark = 0.55
            evithresh = evi_min + (threshmark*(evi_max - evi_min))
            phenologyselect = np.where(evi > evithresh)
            pScalar[phenologyselect] = 1
        #pScalar[:] = 1 #Cancel phenology stress factor
        
        #by default, grasslands and savannas never have pScale=1      
  
    
        #get PAR, conversion PPFD ~= 4.6 (umol/m2/s per W/m2) * SWRF {Integration by Licor, ref McGee, 1972};
        #But also MUST account for cosine factor, easiest way to do that is to fit to data from Fitzjarrald at HF
        #yields results of 1.9 for the factor.
        #Update: No longer a factor devan has fit to sw radiation AND par, simply use the appropriate numbers
        #####PAR CONVERSION ISSUES#####
        

        #Water stress factor
        smScalar = np.ones(shape=np.shape(Sm))
        smScalar[np.where(Sm < theta_thres)] = smScalar[np.where(Sm < theta_thres)]+Q*(Sm[np.where(Sm < theta_thres)] -theta_thres)
        smScalar[np.where(smScalar < 0)] = 0
        smScalar[np.where(smScalar > 1)] = 1


        
        radScalar = 1/(1 + (Rad/radZero))
        
        #Error Check evi,lswi,tempScalar
        evi[np.isnan(evi)] = 0
        lswi[np.isnan(lswi)] = 0
        tempScalar[np.isnan(tempScalar)] = 0
        tempScalar[np.where(tempScalar < 0)] = 0
        wScalar[np.isnan(wScalar)] = 0
        pScalar[np.isnan(pScalar)] = 0
        wScalar[np.where(wScalar < 0)] = 0 
        pScalar[np.where(pScalar < 0)] = 0
        Tscale[k] = tempScalar
        Wscale[k] = wScalar
        SMscale[k] = smScalar
        #Determine GPP
        #NOTE UNITS--VPRM outputs GPP and Respiration in umol/m2/s (conveniently, what is needed here); when multiplied by
        #    influence (ppm/(umol/m2/s)) get ppm 
        #xiaoGPP<-lambdaGPP*eta.zero*tempScalar*wScalar*pScalar*evi*radScalar*swrad*vegFrac*-1
        #NOTE: eliminated eta.zero--now rolled into lambda
        
        
        
        GPP = lambdaGPP*tempScalar*wScalar*pScalar*smScalar*evi*radScalar*Rad*(-1)
        #want vegetative uptake to be negative with respect to the atmosphere, so multiply by negative one 
        #for symmetry    
        gpp = GPP*3600
        GEEpft[k] = gpp*vegFrac
        ###If it is the initial day then there is no previous gpp data to compute the initial day respiration. So the daily GPP for the 
        ### resp is approximated from the initial day GPP.
        if initial_day:
            GPP_daily = compute_daily_GPP_initial(gpp)
        else: 
            GPP_daily = compute_daily_GPP(gpp, GPP_previous[k])
        ###Storing the hourly GPP for the next day.
        GPP_next_day.append(gpp)
        
        #Determine respiration
        Temp0 = copy.deepcopy(Temp)
        if tlow == None:
            Temp0[np.where(Temp0 < 0)] = 0 #set below zero temperature to zero
        else:
            Temp0[np.where(Temp0 < tlow[k])] = tlow[k] #set below "tlow" temperature to "tlow" (veg specific)
        
        Tref = 15
        T0 = -46.02
        tempArrhen = np.exp(E0*((1/(Tref-T0))-(1/(Temp0-T0))))
        tempArrhen[np.isnan(tempArrhen)] = 0
        
        smStress = np.tanh(2*np.pi*k1*Sm + np.arctanh(gamma))
        smStress[np.isnan(smStress)] = 0
        #smStress = Sm/(K1+Sm)
        devanResp = (R0+k2*GPP_daily)*tempArrhen*smStress
        resp = devanResp*3600
        R0_out[k] = (R0+k2*GPP_daily)
        RT_out[k] = tempArrhen
        SMstress[k] = smStress
        RESpft[k] = resp*vegFrac
        NEEpft[k] = (resp+gpp)*vegFrac
        gppTot = gppTot + gpp*vegFrac
        respTot = respTot + resp*vegFrac
        
    
    neeTot = gppTot + respTot
    
    resultLst = {'gee':gppTot, 'resp':respTot, 'nee':neeTot, 'Rad':Rad, 'Temp': Temp, 'evi':np.transpose(evim,(1,0,2,3)), 'lswi':np.transpose(lswim,(1,0,2,3)), 'Tscale':np.transpose(Tscale,(1,0,2,3)), 'Wscale':np.transpose(Wscale,(1,0,2,3)), 'SMscale':np.transpose(SMscale,(1,0,2,3)), 'GEEpft':np.transpose(GEEpft, (1,0,2,3)), 'NEEpft':np.transpose(NEEpft, (1,0,2,3)), 'RESpft':np.transpose(RESpft, (1,0,2,3)), 'R0':np.transpose(R0_out,(1,0,2,3)), 'RT':np.transpose(RT_out,(1,0,2,3)), 'SMstress':np.transpose(SMstress,(1,0,2,3)) }
    
    
    return resultLst, GPP_next_day
        
