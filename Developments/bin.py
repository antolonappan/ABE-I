#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 11:41:53 2017

@author: anto
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

dat = np.loadtxt('MvsL.dat')

Mas = dat[:,0]
Lum = dat[:,1]

Log10Mass = np.log10(Mas)

minMass = math.floor(np.min(Log10Mass))
maxMass = math.ceil(np.max(Log10Mass))

StepSize = 0.2

NoOfbins = (maxMass - minMass)/StepSize

if not NoOfbins.is_integer():
        raise Exception('No of Bins is not a whole no.')
else:
        NoOfbins = int(NoOfbins)
Mass=[]
Luminosity=[]
SigmaLum=[]
for i in range(NoOfbins):
        start = minMass+(StepSize*i)
        stop = minMass+(StepSize*(i+1))
        select = np.where((Log10Mass>start) & (Log10Mass<stop))
        if not len(select[0]) == 0:
                Mass.append(np.average(Log10Mass[select]))
                Luminosity.append(np.average(Lum[select]))
                SigmaLum.append(1./len(select[0]))

PC = pearsonr(Mass,Luminosity)
def fit_func(x, m, c):
    return m*x + c
params = curve_fit(fit_func, Mass, Luminosity)
param = params[0]

plt.scatter(Mass,Luminosity)
plt.text(8.2,-34, 'PPMCC : %f'%PC[0])
plt.text(8.2,-34.5, 'Slope : %f'%param[0])
plt.text(8.2,-35, 'Intercept : %f'%param[1])
