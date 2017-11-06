#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:18:48 2017

@author: anto
"""
import numpy as np
from mapper import map
import _pickle as pl
import matplotlib.pyplot as plt

with open('dataUpd.p', 'rb') as fr:
        data = pl.load(fr)

with open('mass.p', 'rb') as fr:
        mass = pl.load(fr)

if not mass.keys() == data.keys():
        raise Exception('FileError: Mass data is not Matching with data set')

bh_id = list(data.keys())
mass_x = []
luminosity_y = []
i = 0
le = len(bh_id)
for bh in bh_id:
    M = map(data[bh])
    mass_x.append(mass[bh]*10e10)
    luminosity_y.append(M.avg_luminosity)
    i += 1
    print('%d/%d' % (i, le))

da = np.array([mass_x,luminosity_y])
np.savetxt('MvsL.dat', da.T)
