#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Style:PEP-8
"""

              ****************
              * No. of files *
              ****************
                    /|\
                   / | \
                  /  |  \
        ###############################
        # fork files to N No.of Nodes #
        ###############################
               |    |    |    |
  ex:N = 4 ==> |    |    |    |
             __|____|____|____|___
             |     N-Node        |
             |     ------        |
             |1.Read a snapfile  |
             |2.Find particles   |
             |  around every     |
             |  BHs in that file |
             |3.Update data dict |
             |  with key as BH ID|
             |  and correspoding |
             |  to its data array|
             |4.Read the nextfile|
             |  and go to step 2 |
             |___________________|
"""
__author__ = "Anto I Lonappan"
__copyright__ = "Copyright 2017, MB-II analysis Project, Presidency University"
__date__ = "14-Sep-2017"
__version__ = "0.4"
__email__ = "antoidicherian@gmail.com"
__status__ = "Developing"


# PYTHON LIBRARIES
import sys
import glob
import logging
import numpy as np
import _pickle as pl
from numba import jit
from mpi4py import MPI
# *WRITTEN LIBRARIES
sys.path.insert(0, '/media/rc/0a4532c4-f23a-43f0-bfa5-1e4579468872/modules')  # PATH OF MODULES
import gad_snapread as snr

# MPI INITIALIZATIION
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

box = 25  # SIZE OF SELECTION BOX IN KPC


# VERBOSE FOR PRINTING
verbose = False
if'-v' in sys.argv:
        verbose = True
vprint = print if verbose else lambda *a, **k: None

 MASS CUT-OFF
mass_sorting = False
if '-M' in sys.argv:
        mass_sorting = True
        M_i = sys.argv.index('-M')


def InitMass():
        MassBH = 1e7
        if (len(sys.argv)-1) > M_i:
                try:
                        MassBH = float(sys.argv[M_i+1])
                except ValueError:
                        pass
        return MassBH


# FUNCTION FOR LOGGER
formatter = logging.Formatter('%(asctime)s %(message)s')


def setup_logger(name, level=logging.INFO):
        handler = logging.FileHandler('logs/'+name)
        handler.setFormatter(formatter)
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger


# SNAPFILE READER
def snapread(filename):
        return snr.readsnap(filename, suppress=True)


# FUNCTION FOR FINDING GAS PARTICLES AROUND THE BH
@jit
def selection(gas_position, bh_position):
        bh_box_array = np.where(((
                gas_position[:, 0] >= bh_position[0]-box) & (
                gas_position[:, 0] <= bh_position[0]+box)) & ((
                   gas_position[:, 1] >= bh_position[1]-box) & (
                   gas_position[:, 1] <= bh_position[1]+box)) & ((
                       gas_position[:, 2] >= bh_position[2]-box) & (
                       gas_position[:, 2] <= bh_position[2]+box
                       )))[0]
        return bh_box_array

Acc_dict = {}
Mass_dict = {}
Data_dict = {}
snap_dict = {}
bh_id_cat, bh_x_cat, bh_y_cat, bh_z_cat = [], [], [], []

SNAP = setup_logger('snap%d.log'%rank)
# MAIN FUNCTION FOR SEARCHING AND SAVING DATA
def localread(file):
        snap = snapread(file)
        bh_id = snap.read_id('BH')
        bh_mass = snap.read_mass('BH')
        bh_position = snap.read_pos('BH')
        bh_acc = snap.read_mbhdot()
        ini_len = len(bh_id)
        if '-M' in sys.argv:
                cutoff_mass = InitMass()
                cutoff_index = np.where(bh_mass*1e10 > cutoff_mass)
                bh_id = bh_id[cutoff_index]
                bh_position = bh_position[cutoff_index]
                bh_mass = bh_mass[cutoff_index]
                bh_acc = bh_acc[cutoff_index]
        fin_len = len(bh_id)
        SNAP.info('%s:: %d/%d' % (file,ini_len, fin_len))
        temp_snap_dict = {file: bh_id}
        snap_dict.update(temp_snap_dict)
        gas_position = snap.read_pos('GAS')
        gas_density = snap.read_density()
        gas_interalE = snap.read_U()
        gas_smoothingL = snap.read_hsml()
        vprint('%d/%d BHs selected from %s' % (fin_len, ini_len, file[11:]))


        for i in range(fin_len):
                vprint('     Analysing BH %d/%d in %s' % (i+1, fin_len,
                                                          file[11:]))

                bh_id_cat.append(str(bh_id[i]))
                bh_x_cat.append(bh_position[i, 0])
                bh_y_cat.append(bh_position[i, 1])
                bh_z_cat.append(bh_position[i, 2])
                bh_box_array = selection(gas_position, bh_position[i])
                bh_loc_gas_position = gas_position[bh_box_array]
                bh_loc_gas_density = gas_density[bh_box_array]
                bh_loc_gas_interalE = gas_interalE[bh_box_array]
                bh_loc_gas_smoothingL = gas_smoothingL[bh_box_array]

                dat = np.array([bh_loc_gas_interalE, bh_loc_gas_density,
                               bh_loc_gas_position[:, 0],
                               bh_loc_gas_position[:, 1],
                               bh_loc_gas_position[:, 2],
                               bh_loc_gas_smoothingL])
                temp_data = {'%d' % bh_id[i]: dat.T}
                Data_dict.update(temp_data)
                temp_mass_info = {'%d' % bh_id[i]: bh_mass[i]}
                Mass_dict.update(temp_mass_info)
                temp_acc_info = {'%d' % bh_id[i]: bh_acc[i]}
                Acc_dict.update(temp_acc_info)
               


if __name__ == "__main__":
    files = glob.glob('snapdir_068/snap*')
    
    local_n = int(len(files) / size)
    local_i = rank*int(local_n)
    local_f = local_i + local_n
    local_file = files[local_i:local_f]

    for file in local_file:
            localread(file)
            
    
    with open('data/data%d.p' % rank, 'wb') as fp:
                pl.dump(Data_dict, fp)
    with open('data/mass%d.p' % rank, 'wb') as fp:
                pl.dump(Mass_dict, fp)
    with open('data/acc%d.p' % rank, 'wb') as fp:
                pl.dump(Acc_dict, fp)
    
    bh_cat = np.array([bh_id_cat, bh_x_cat, bh_y_cat, bh_z_cat])

    with open('catalogue/bh_cat%d.p' % rank, 'wb') as fp:
                pl.dump(bh_cat.T, fp)

    with open('catalogue/snap_cat%d.p' % rank, 'wb') as fp:
                pl.dump(snap_dict, fp)
    
