#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write VPRM constants into a variable a dataFrame and a text file.
"""

import pandas as pd
import numpy as np
import datetime

lambdaGPP_factor = 0.7
shrubfactor = 0.3
grassfactor = 0.3

#
#                 VPRM Class    STILT-VPRM class   Validation Site
#Evergreen A		1A		1       	NOBS
#Evergreen B		1B		2		Niwot
#Evergreen C		1C		3		Oregon
#Evergreen D		1D		4		Donaldson
#Deciduous		2		5		Havard
#Mixed forest		3		6		Howland
#Shrubland		4		7		Lucky-Hill
#Savanna		5		8		Tonzi
#Cropland-Soy		6A		9		Mead-S2
#Cropland-Maize		6B		9		Mead-S2
#Grassland		7		10		Vaira
#Others			8		11		--

vnames = ["Evergreen A","Evergreen B","Evergreen C","Evergreen D","Deciduous","Mixed forest","Shrubland","Savanna","Cropland","Grassland","Others"]

#lambdaGPP now includes the 0.044 eta.zero factor, so that is removed from
#processing--no longer necessary


def WriteVPRMConstants(outdir='./', nveg = 8):
    
    lambdaGPP_sw = lambdaGPP_factor*np.array([0.499,
                                              0.282,
                                              0.3084,
                                              0.271,
                                              0.1955,
                                              0.2856,
                                              0.0874,
                                              0.1141,
                                              0.1350,
                                              0.1748,
                                              0.00])
    
    lambdaGPP_par = lambdaGPP_factor*np.array([0.263,
                                              0.148,
                                              0.187,
                                              0.142,
                                              0.130,
                                              0.192,
                                              0.2187*shrubfactor,
                                              0.0859,
                                              0.115,
                                              0.334*grassfactor,
                                              0.00])
    
    alphaResp = np.array([0.2672,
                          0.2668,
                          0.1797,
                          0.1878,
                          0.1495,
                          0.2258,
                          0.0239,
                          0.0049,
                          0.1699,
                          0.0881,
                          0.00])
    
    intResp = np.array([0,
                       0,
                       0.8800,
                       0,
                       0.8233,
                       0.4321,
                       0,
                       0,
                       -0.0144,
                       0.05843,
                       0])
    
    tempMin = np.array([0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        2,
                        2,
                        5,
                        2,
                        0])
    
    tempMax = np.full((11), 40)
    
    tempOpt = np.array([20,
                        20,
                        20,
                        20,
                        20,
                        20,
                        20,
                        20,
                        22,
                        18,
                        0])
    
    parZero = np.array([237,
                        400,
                        496,
                        577,
                        616,
                        392,
                        690,
                        1297,
                        1439,
                        300,
                        0.00])
    
    swradZero = np.array([124,
                       210,
                       270.2,
                       303,
                       271.4,
                       236.6,
                       363,
                       682,
                       690,
                       229.1,
                       0.00])
    
    eviMax = np.array([np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.00])
    
    lswiMax = np.array([np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.00])
    
    eviMin = np.array([np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    0.00])

    outLst = pd.DataFrame()
    outLst['Class'] = vnames
    outLst['lambdaGPP.sw'] = lambdaGPP_sw
    outLst['lambdaGPP.par'] = lambdaGPP_par
    outLst['alphaResp'] = alphaResp
    outLst['intResp'] = intResp
    outLst['tempMin'] = tempMin
    outLst['tempMax'] = tempMax
    outLst['tempOpt'] = tempOpt
    outLst['parZero'] = parZero
    outLst['swradZero'] = swradZero
    outLst['eviMin'] = eviMin
    outLst['eviMax'] = eviMax
    outLst['lswiMax'] = lswiMax
    
    dateLn = datetime.datetime.now()
    outfile = outdir + 'vprmConstants.' + str(dateLn.day) + str(dateLn.month) + str(dateLn.year) + '.txt'
    
    outLst.to_csv(outfile, header=True, index=None, sep=' ')
    outLst11 = outLst

    outLst = outLst.iloc[[2, 4,5,6,7,8,9,10]]
    outLst = outLst.reset_index(drop=True)
    outLst8 = outLst
    outfile = outdir + 'vprmConstants8.' + str(dateLn.day) + str(dateLn.month) + str(dateLn.year) + '.txt'
    
    outLst.to_csv(outfile, header=True, index=None, sep=' ')
    if nveg == 11:
        return outLst11
    elif nveg == 8:
        return outLst8
    else:
        return None