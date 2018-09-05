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
__version__ = "0.5"
__email__ = "antoidicherian@gmail.com"
__status__ = "Completed"


# PYTHON LIBRARIES
import os
import sys
import glob
import logging
import numpy as np
import _pickle as pl
from numba import jit
from mpi4py import MPI
from configparser import ConfigParser as cnf
# *WRITTEN LIBRARIES
pwd = os.getcwd()
sys.path.insert(0, pwd + '/Modules')  # PATH OF MODULES
import gad_snapread as snr


# CONFIGURATION
conf = cnf()
conf.read('abe.ini')
snapshot_dir = conf.get('inputs', 'snapshot')
output_dir = conf.get('live', 'cat_mak_out')
log_dir = conf.get('live', 'cat_mak_log')
box = float(conf.get('misc', 'box'))
mass_cutoff = conf.get('misc', 'mass_cutoff')
lower_cutoff_value = float(conf.get('misc', 'lower_cutoff_value'))
upper_cutoff_value = float(conf.get('misc', 'upper_cutoff_value'))
last_program = conf.get('live', 'last_program')
mode_run = conf.get('misc', 'mode_run')

# MPI INITIALIZATIION
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# CHECK PREVIOUS PROGRAM
if rank == 0:
    if mode_run == 'automated':
        if not last_program == 'initial.py': 
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


# FUNCTION FOR LOGGER
formatter = logging.Formatter('%(asctime)s %(message)s')


def setup_logger(name, level=logging.INFO):
        handler = logging.FileHandler('%s/%s' % (log_dir, name))
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
        if mass_cutoff is 'T':
                cutoff_index = np.where((bh_mass*1e10 > lower_cutoff_value) &
                                         (bh_mass*1e10 < upper_cutoff_value))
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
        gas_electronAb = snap.read_ne()
        

        for i in range(fin_len):
                
                bh_id_cat.append(str(bh_id[i]))
                bh_x_cat.append(bh_position[i, 0])
                bh_y_cat.append(bh_position[i, 1])
                bh_z_cat.append(bh_position[i, 2])
                bh_box_array = selection(gas_position, bh_position[i])
                bh_loc_gas_position = gas_position[bh_box_array]
                bh_loc_gas_density = gas_density[bh_box_array]
                bh_loc_gas_interalE = gas_interalE[bh_box_array]
                bh_loc_gas_smoothingL = gas_smoothingL[bh_box_array]
                bh_loc_gas_electronAb = gas_electronAb[bh_box_array]

                dat = np.array([bh_loc_gas_interalE, bh_loc_gas_density,
                               bh_loc_gas_position[:, 0],
                               bh_loc_gas_position[:, 1],
                               bh_loc_gas_position[:, 2],
                               bh_loc_gas_smoothingL, bh_loc_gas_electronAb])
                temp_data = {'%d' % bh_id[i]: dat.T}
                Data_dict.update(temp_data)
                temp_mass_info = {'%d' % bh_id[i]: bh_mass[i]}
                Mass_dict.update(temp_mass_info)
                temp_acc_info = {'%d' % bh_id[i]: bh_acc[i]}
                Acc_dict.update(temp_acc_info)
               


if __name__ == "__main__":
    files = glob.glob(snapshot_dir + '/snap*')
    
    local_n = int(len(files) / size)
    local_i = rank*int(local_n)
    local_f = local_i + local_n
    local_file = files[local_i:local_f]

    for file in local_file:
            localread(file)
            
    
    with open(output_dir + '/data%d.p' % rank, 'wb') as fp:
                pl.dump(Data_dict, fp)
    with open(output_dir + '/mass%d.p' % rank, 'wb') as fp:
                pl.dump(Mass_dict, fp)
    with open(output_dir + '/acc%d.p' % rank, 'wb') as fp:
                pl.dump(Acc_dict, fp)
    
    bh_cat = np.array([bh_id_cat, bh_x_cat, bh_y_cat, bh_z_cat])

    with open(output_dir + '/bh_cat%d.p' % rank, 'wb') as fp:
                pl.dump(bh_cat.T, fp)

    with open(output_dir + '/snap_cat%d.p' % rank, 'wb') as fp:
                pl.dump(snap_dict, fp)

    if rank == 0:
        conf.set('live', 'last_program', 'CatalogMaker_MPI.py')
        with open('abe.ini', 'w') as fw:
               conf.write(fw)
    
