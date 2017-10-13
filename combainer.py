#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 10:20:35 2017

@author: anto
"""
import numpy as np
import _pickle as pl
import glob
import os

i = 0
bh_cat_files = glob.glob('catalogue/bh*')
for file in bh_cat_files:
        with open(file, 'rb') as reader:
                globals()['data%d' % i] = pl.load(reader)
        i += 1

for j in range(i):
        if j == 0:
                data = np.concatenate((globals()['data%d' % j],
                                       globals()['data%d' % (j+1)]), axis=0)
        elif j == 1:
                pass
        else:
                data = np.concatenate((data, globals()['data%d' % j]), axis=0)
with open('catalogue/bh_cat.p', 'wb') as writer:
        pl.dump(data, writer)

for file in bh_cat_files:
        os.remove(file)

snap_cat = {}
snap_cat_files = glob.glob('catalogue/snap*')
for file in snap_cat_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        snap_cat.update(dic)
with open('catalogue/snap_cat.p', 'wb') as writer:
        pl.dump(snap_cat, writer)

for file in snap_cat_files:
        os.remove(file)

data = {}
data_files = glob.glob('data/data*')
for file in data_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        data.update(dic)
with open('data/data.p', 'wb') as writer:
        pl.dump(data, writer)

for file in data_files:
        os.remove(file)

mass = {}
mass_files = glob.glob('data/mass*')
for file in mass_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        mass.update(dic)
with open('data/mass.p', 'wb') as writer:
        pl.dump(mass, writer)

for file in mass_files:
        os.remove(file)
