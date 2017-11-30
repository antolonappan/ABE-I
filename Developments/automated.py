#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Style:PEP-8


import sys
import glob
sys.path.insert(0, '/home/anto/Desktop/test0/modules')
import logging
import numpy as np
import gad_snapread as snr
#import mapper
from numba import jit


__author__ = "Anto I Lonappan"
__copyright__ = "Copyright 2017, MB-II analysis Project, Presidency University"
__credit__ = ["Alankar Dutta, Rudrani K C"]
__date__ = "17-Aug-2017"
__version__ = "0.0"
__email__ = "antoidicherian@gmail.com"
__status__ = "Developing"

box = 25
loc = ""
folder = "snapshot"

formatter = logging.Formatter('%(asctime)s %(message)s')
def setup_logger(name, level=logging.INFO):
        handler = logging.FileHandler('log/'+name)
        handler.setFormatter(formatter)
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger
SNAP = setup_logger('snapshot.log')
BH = setup_logger('bh.log')
SER = setup_logger('search.log')
def snapread(filename):
        return snr.readsnap(filename,suppress=True)


def selection(gas_position,bh_position,i):
        bh_box_array = np.where(((
                gas_position[:, 0] >= bh_position[i, 0]-box) & (
                gas_position[:, 0] <= bh_position[i, 0]+box)) & ((
                gas_position[:, 1] >= bh_position[i, 1]-box) & (
                gas_position[:, 1] <= bh_position[i, 1]+box)) & ((
                gas_position[:, 2] >= bh_position[i, 2]-box) & (
                gas_position[:, 2] <= bh_position[i, 2]+box
                       )))[0]
        return bh_box_array


if'-v' in sys.argv:
        verbose = True
else:
        verbose = False
vprint = print if verbose else lambda *a, **k: None


#files.sort()


if __name__ == "__main__":
#   try:
    files = glob.glob('snapshot/snap*')
    for file in files:
      snap = snapread(file)
      bh_id = snap.read_id('BH')
      l = len(bh_id)
      bh_position = snap.read_pos('BH')
      gas_position = snap.read_pos('GAS')
      gas_density = snap.read_density()
      gas_interalE = snap.read_U()
      gas_smoothingL= snap.read_hsml()
      vprint('%d BHs found in %s'%(l, file))
      SNAP.info('%d BHs found in %s'%(l, file[9:]))
      files_t=files[:]
      files_t.remove(file)


      for i in range(l):
        vprint('     Analysing BH %d/%d in %s'%(i+1,l,file[9:]))
        BH.info('Analysing BH %d/%d in %s'%(i+1,l,file[9:]))
        bh_box_array = selection(gas_position, bh_position, i)
        bh_loc_gas_position = gas_position[bh_box_array]
        bh_loc_gas_density = gas_density[bh_box_array]
        bh_loc_gas_interalE = gas_interalE[bh_box_array]
        bh_loc_gas_smoothingL= gas_smoothingL[bh_box_array]
#for i in range(len(bh_position[:, 0])):
#        print (i)

        for j in files_t:
                vprint('          Searching for ranges in %s'%j[9:])
                SER.info('Searching for ranges in %s for %d/%d BH of %s'%(j[9:],i+1,l,file[9:]))
                snap_t = snapread(j)
                gas_position_t = snap_t.read_pos('GAS')
                bh_position_t = snap.read_pos('BH')
                #print (gas_position_t)
                bh_box_array_t = selection(gas_position_t, bh_position, i)
                #print (bh_box_array_t)
                if len(bh_box_array_t) is not 0:
                        gas_density_t = snap_t.read_density()
                        gas_interalE_t = snap_t.read_U()
                        gas_smoothingL_t= snap_t.read_hsml()
                        bh_loc_gas_position=np.concatenate((bh_loc_gas_position,gas_position_t[bh_box_array_t]))
                        bh_loc_gas_density=np.concatenate((bh_loc_gas_density,gas_density_t[bh_box_array_t]))
                        bh_loc_gas_interalE = np.concatenate((bh_loc_gas_interalE,gas_interalE_t[bh_box_array_t]))
                        bh_loc_gas_smoothingL=np.concatenate((bh_loc_gas_smoothingL,gas_smoothingL_t[bh_box_array_t]))
                        vprint('                    %s Particles found in %s'%(len(bh_box_array_t),j))
                else:
                        vprint('                    Not found in %s'%j)

        dat = np.array([bh_loc_gas_interalE,bh_loc_gas_density,
                        bh_loc_gas_position[:,0],bh_loc_gas_position[:,1],
                        bh_loc_gas_smoothingL])
        np.savetxt('data/%d.dat'%bh_id[i], dat.T)
#   except:
#        print(sys.exc_info()[])
