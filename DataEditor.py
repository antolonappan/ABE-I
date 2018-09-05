#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Style:PEP-8
"""
Created on Thu Sep  7 10:20:35 2017

@author: anto
"""
__author__ = "Anto I Lonappan"
__copyright__ = "Copyright 2017, MB-II analysis Project, Presidency University"
__date__ = "14-Sep-2017"
__version__ = "0.0"
__email__ = "antoidicherian@gmail.com"
__status__ = "Completed"


# PYTHON LIBRARIES
import numpy as np
import _pickle as pl
from configparser import ConfigParser as cnf

# CONFIGURATION
conf = cnf()
conf.read('abe.ini')
output_dir = conf.get('live', 'cat_anl_out')
last_program = conf.get('live', 'last_program')
mode_run = conf.get('misc', 'mode_run')

if mode_run == 'automated':
    if not last_program == 'CatalogAnalyser_MPI.py': 
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

XH = 0.76

with open('%s/dataUpd.p' % output_dir, 'rb') as fr:
        data = pl.load(fr)
Real_data = data.copy()

bhs = list(data.keys())


def temperature(internal_E, ElectronAb):
        return 120.99 * (2./3.) * internal_E/(ElectronAb *
                                              XH + (1.-XH) * 0.25 + XH)


for bh in bhs:
        dummy = data[bh]
        temp = temperature(dummy[:, 0], dummy[:, -1])
        dummy1 = np.array([temp, dummy[:, 1], dummy[:, 2], dummy[:, 3],
                           dummy[:, 4], dummy[:, 5]])
        data[bh] = dummy1.T

with open('%s/dataEdt.p' % output_dir,'wb') as fw:
        pl.dump(data, fw)
print('Edited')

conf.set('live', 'last_program', 'DataEditor.py')
with open('abe.ini', 'w') as fw:
       conf.write(fw)
    
