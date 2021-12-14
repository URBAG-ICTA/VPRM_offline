#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VPRM code to obtain the R, GEE and NEE fluxes
"""
import numpy as np
import datetime 
import src.WriteVPRMConstants as WriteVPRMConstants
from sys import exit
import copy

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

    

def offlineVPRM(Temp, Rad, start_mdy, start_hrs,  evi, lswi, vegFracMap,
                      evi_max, evi_min, lswi_max, lswi_min, datp, usepar=False, returnNEEOnly=False,
                      returnLongList=True, delt=1, evi_times=np.linspace(2005041,2005305,8),
                      vprm_par_name="vprmoptNONCERES.meas.mod.local.par", tlow= None,returnJackList=False,parapath='./data/VPRMpara/'):
    
    nhrs = np.shape(Temp)[0]
    evi = extract_wrf_times_from_evi(evi = evi, start_mdy = start_mdy, start_hrs = start_hrs, nhrs = nhrs, evi_times = evi_times, delt = delt)
    lswi = extract_wrf_times_from_evi(evi = lswi, start_mdy = start_mdy, start_hrs = start_hrs, nhrs = nhrs, evi_times = evi_times, delt = delt)
    #Temp = np.transpose(Temp, axes=(1,2,0))
    #Rad = np.transpose(Rad, axes=(1,2,0))

    

    nVegClass = 8
    vprmConstants = WriteVPRMConstants.WriteVPRMConstants(outdir = datp, nveg = 8)


    checkdim = nVegClass

    
    if (Temp.size != Rad.size) or (Temp.size != evi.size/checkdim) or (Temp.size != lswi.size/checkdim):
        exit('1. Temperature and radiation must have same dimensions.')
    
    Temp = Temp - 273.15
    if returnLongList:
        vegFracLst = {0:None, 1:None, 2:None, 3:None, 4:None, 5:None, 6:None}
        geeLst = {0:None, 1:None, 2:None, 3:None, 4:None, 5:None, 6:None}
        respLst = {0:None, 1:None, 2:None, 3:None, 4:None, 5:None, 6:None}
        wscalarLst = {0:None, 1:None, 2:None, 3:None, 4:None, 5:None, 6:None}
        pscalarLst = {0:None, 1:None, 2:None, 3:None, 4:None, 5:None, 6:None}
        radscalarLst = {0:None, 1:None, 2:None, 3:None, 4:None, 5:None, 6:None}
        paramLst = vprmConstants
    
    if returnJackList:
        neeLst = {0:None, 1:None, 2:None, 3:None, 4:None, 5:None, 6:None}
        paramLst = vprmConstants
    
    geeTot = np.zeros(np.shape(Temp))
    respTot = np.zeros(np.shape(Temp))
    
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
        if usepar:
            lambdaGPP = vprmConstants.loc[k, 'lambdaGPP.par']
        else:
            lambdaGPP = vprmConstants.loc[k, 'lambdaGPP.sw']
        alphaResp = vprmConstants.loc[k, 'alphaResp']
        intResp = vprmConstants.loc[k, 'intResp']
        
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
        
        if (k in [3,6]):
            wScalar = (lswi - lswi_min)/(lswi_max - lswi_min)
        else:
            wScalar = (1 + lswi)/(1 + lswi_max)
        
        #wScalar[:] = 1 #Cancel water stress factor
        pScalar = (1 + lswi)/2
        
        if  (k == 0):
            phenologyselect = np.arange(0, evi.size)
            pScalar[:] = 1
        
        if (k in [1,2,3,5,7]):  #if decid, mixed, shrub, crop, or other
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
        
        if usepar:
            radZero = vprmConstants.loc[k, 'parZero']
        else:
            radZero = vprmConstants.loc[k, 'swradZero']
        
        radScalar = 1/(1 + (Rad/radZero))
        
        #Error Check evi,lswi,tempScalar
        evi[np.isnan(evi)] = 0
        lswi[np.isnan(lswi)] = 0
        tempScalar[np.isnan(tempScalar)] = 0
        tempScalar[np.where(tempScalar < 0)] = 0
        wScalar[np.isnan(wScalar)] = 0
        pScalar[np.isnan(pScalar)] = 0
        wScalar[np.where(wScalar < 0)] = 0 #sometimes negative at borders 
        pScalar[np.where(pScalar < 0)] = 0
        
        
        #Determine GPP
        #NOTE UNITS--VPRM outputs GPP and Respiration in umol/m2/s (conveniently, what is needed here); when multiplied by
        #    influence (ppm/(umol/m2/s)) get ppm 
        #xiaoGPP<-lambdaGPP*eta.zero*tempScalar*wScalar*pScalar*evi*radScalar*swrad*vegFrac*-1
        #NOTE: eliminated eta.zero--now rolled into lambda
        
        
        
        GPP = lambdaGPP*tempScalar*wScalar*pScalar*evi*radScalar*Rad*vegFrac*(-1)
        #want vegetative uptake to be negative with respect to the atmosphere, so multiply by negative one 
        #for symmetry    
        gee = GPP*3600
        
        #Determine respiration
        Temp0 = copy.deepcopy(Temp)
        if tlow == None:
            Temp0[np.where(Temp0 < 0)] = 0 #set below zero temperature to zero
        else:
            Temp0[np.where(Temp0 < tlow[k])] = tlow[k] #set below "tlow" temperature to "tlow" (veg specific)
        
        devanResp = (Temp0*alphaResp + intResp)*vegFrac
        resp = devanResp*3600
        
        geeTot = geeTot + gee
        respTot = respTot + resp
        
        if returnLongList:
            vegFracLst[k] = vegFrac
            geeLst[k] = gee
            respLst[k] = resp
            wscalarLst[k] = wScalar
            pscalarLst[k] = pScalar
            radscalarLst[k] = radScalar
        
        if returnJackList:
            neeLst[k] = gee + resp
    
    neeTot = geeTot + respTot
    if returnNEEOnly:
        return neeTot
    
    resultLst = {'gee':geeTot, 'resp':respTot, 'nee':neeTot, 'Rad':Rad, 'Temp': Temp, 'evi':evi, 'lswi':lswi }
    
    if returnLongList:
        debugLst = {'lambdaGPP':lambdaGPP, 'tempScalar':tempScalar, 'wScalar':wScalar, 'pScalar':pScalar, 'radScalar':radScalar}
        return [resultLst, geeLst, respLst, vegFracLst, paramLst, debugLst]
    
    if returnJackList:
        return [neeLst, paramLst]
    
    return resultLst
        
