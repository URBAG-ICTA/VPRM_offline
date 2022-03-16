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
#                 VPRM Class
#Evergreen 		1
#Deciduous		2
#Mixed forest		3
#Shrubland		4
#Savanna		5
#Cropland		6
#Grassland		7
#Others			8

vnames = ["Evergreen","Deciduous","Mixed forest","Shrubland","Savanna","Cropland","Grassland","Others"]

#lambdaGPP now includes the 0.044 eta.zero factor, so that is removed from
#processing--no longer necessary


def WriteVPRMConstants(outdir='./', nveg = 8):
    
    lambdaGPP_sw =                  np.array([0.22577703,
                                              0.21489270,
                                              0.16293380,
                                              0.29311134,
                                              0.1141,
                                              0.08626603,
                                              0.11930965,
                                              0.00])
    
    print(lambdaGPP_sw)
    lambdaGPP_par = lambdaGPP_factor*np.array([0.187,
                                              0.130,
                                              0.192,
                                              0.2187*shrubfactor,
                                              0.0859,
                                              0.115,
                                              0.334*grassfactor,
                                              0.00])
    
    alphaResp = np.array([0.28773167,
                          0.18056630,
                          0.24447911,
                          0.05464646,
                          0.0049,
                          0.09231632,
                          0.1245603,
                          0.00])
    
    intResp = np.array([-1.09316696,
                       0.83641734,
                       -0.48669162,
                       -0.12080592,
                       0.0000,
                       0.28788863,
                       0.01743361,
                       0])
    
    tempMin = np.array([0,
                        0,
                        0,
                        2,
                        2,
                        5,
                        2,
                        0])
    
    tempMax = np.full((8), 40)
    
    tempOpt = np.array([20,
                        20,
                        20,
                        20,
                        20,
                        22,
                        18,
                        0])
    
    parZero = np.array([496,
                        616,
                        392,
                        690,
                        1297,
                        1439,
                        300,
                        0.00])
    
    swradZero = np.array([275.4595,
                       254.4188,
                       446.0888,
                       70.3829,
                       682.0,
                       1132.2,
                       527.9303,
                       0.00])
    
    eviMax = np.array([np.nan,
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
                    0.00])
    
    eviMin = np.array([np.nan,
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

    outLst8 = outLst
    outfile = outdir + 'vprmConstants8.' + str(dateLn.day) + str(dateLn.month) + str(dateLn.year) + '.txt'
    
    outLst.to_csv(outfile, header=True, index=None, sep=' ')
    if nveg == 8:
        return outLst8
    else:
        return None
