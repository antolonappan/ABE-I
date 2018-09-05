#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Style:PEP-8
"""
Created on Thu Sep  7 10:20:35 2017

@author: anto
"""
__author__ = "Anto I Lonappan"
__copyright__ = "Copyright 2017, MB-II analysis Project, Presidency University"
__date__ = "7-Sep-2017"
__version__ = "0.4"
__email__ = "antoidicherian@gmail.com"
__status__ = "Completed"

# PYTHON LIBRARIES
import os
import sys
import glob
import numpy as np
import _pickle as pl
from configparser import ConfigParser as cnf


# CONFIGURATION
conf = cnf()
conf.read('abe.ini')
output_dir = conf.get('live', 'cat_mak_out')
output_dir_lum = conf.get('live', 'dat_anl_out')
delete = conf.get('misc', 'delete_dump')
last_program = conf.get('live', 'last_program')
mode_run = conf.get('misc', 'mode_run')

if 'lum' in sys.argv:
    print('Combaining LUM')
    lum = {}
    lum_files = glob.glob('%s/lum*' % output_dir_lum)
    for file in lum_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        lum.update(dic)
    with open('%s/lum.p' % output_dir_lum, 'wb') as writer:
        pl.dump(lum, writer)
    
    if delete is 'T':
        for file in lum_files:
            os.remove(file)
    sys.exit()

elif mode_run == 'automated':
    if not last_program == 'CatalogMaker_MPI.py': 
          print("""
                I found Last program that you have run was {}
                Before running this code please run DataEditor.py.
                Also see abe.ini
                """.format(last_program))
          raise Exception('PreviousRunError')
    else:
        pass
elif mode_run == 'individual':
    pass
############################### bh_cat #######################################
print('Combaining bh_cat ')
i = 0
bh_cat_files = glob.glob('%s/bh*' % output_dir)
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
with open('%s/bh_cat.p' % output_dir, 'wb') as writer:
        pl.dump(data, writer)

if delete is 'T':
    for file in bh_cat_files:
        os.remove(file)
############################### snap_cat ######################################
print('Combaining snap_cat')
snap_cat = {}
snap_cat_files = glob.glob('%s/snap*' % output_dir)
for file in snap_cat_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        snap_cat.update(dic)
with open('%s/snap_cat.p' % output_dir, 'wb') as writer:
        pl.dump(snap_cat, writer)

if delete is 'T':
    for file in snap_cat_files:
        os.remove(file)

################################# data #######################################
print('Combaining Data')
data = {}
data_files = glob.glob('%s/data*' % output_dir)
for file in data_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        data.update(dic)
with open('%s/data.p' % output_dir, 'wb') as writer:
        pl.dump(data, writer)

if delete is 'T':
    for file in data_files:
        os.remove(file)

################################# mass #######################################
print('Combaining Mass')
mass = {}
mass_files = glob.glob('%s/mass*' % output_dir)
for file in mass_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        mass.update(dic)
with open('%s/mass.p' % output_dir, 'wb') as writer:
        pl.dump(mass, writer)
if delete is 'T':
    for file in mass_files:
        os.remove(file)
################################ Acc ##########################################
print('Combaining Acc')
acc = {}
acc_files = glob.glob('%s/acc*' % output_dir)
for file in acc_files:
        with open(file, 'rb') as reader:
                dic = pl.load(reader)
        acc.update(dic)
with open('%s/acc.p' % output_dir, 'wb') as writer:
        pl.dump(acc, writer)

if delete is 'T':
    for file in acc_files:
        os.remove(file)
###############################################################################
conf.set('live', 'last_program', 'Combiner.py')
with open('abe.ini', 'w') as fw:
    conf.write(fw)
    

