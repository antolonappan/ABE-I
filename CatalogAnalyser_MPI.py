#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Style : PEP-8
"""
Created on Thu Sep  7 12:10:18 2017
@author: anto

MPI PROCESS ILLUSTRATION:
************************


  /-->---->---\FILES/                                        <-----<-----<--\
  |            | 2 |                                        /               |
  |            | 1 |                                 No.of files            |
  |   _________\_0_/____________              _____________|_____________   |
  |   |         ROOT           |              |         NODES           |   |
  ^   |         ----           |              |         -----           |   ^
  |   |1.Read BHs catalogue    |              |                         |   |
  |   |2.Read gas position     |              |                         |   |
  |   |3.Sending gas position  |->->->->->->->|1.Receiving gas position |   |
  |   |4.Read other details    |              |2.Search for positions   |   |
  ^   |        from file       |              |3.If found: update a     |   ^
  |   |5.Search for positions  |              |  dictionary with array  |   |
  |   |6.If found: update data |              |  details and bh as key  |   |
  |   |7.Wait for other nodes  |          <-<-|4.send this dictionary   |   |
  |   |  to finish their search|         /    |  to root                |   |
  ^   |8.Receive found particle|<-<-<-<-<     |5.Wait for next iteration|   ^
  |   |  details from nodes    |              |  to start at ROOT       |   |
  |   |9.Update data           |              |                         |   |
  |   |________________________|              |_________________________|   |
  |                |                                       |                |
  ^-<-----<-----<--/                                       \-->------>----->-

TO RUN:
******
REQUIRMENTS:1. MPICH or MPICH2 and mpi4py.
            2. NUMBA from conda

PYTHON VERSION: Python 3.6.1.
                ==> For lower versions change '_pickle' to 'cPickle'

CMD:
    $mpiexec -n <no.of nodes> python catalogue.py
     ==> MPI load balancing is implemented, so that no.of nodes can also be odd
     ==> mpirun will also work, but mpiexec is recommended
"""
__author__ = "Anto I Lonappan"
__copyright__ = "Copyright 2017, MB-II analysis Project, Presidency University"
__date__ = "14-Sep-2017"
__version__ = "1.0"
__email__ = "antoidicherian@gmail.com"
__url__ = "https://github.com/antolonappan/Automated-Black-Hole-Explorer-from-GADGET-2"
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
#  *WRITTEN LIBRARIES
pwd = os.getcwd()
sys.path.insert(0, pwd + '/Modules')
import gad_snapread as snr

# CONFIGURATION
conf = cnf()
conf.read('abe.ini')
snapshot_dir = conf.get('inputs', 'snapshot')
box = float(conf.get('misc', 'box'))
input_dir = conf.get('live', 'cat_mak_out')
output_dir = conf.get('live', 'cat_anl_out')
log_dir = conf.get('live', 'cat_anl_log')
last_program = conf.get('live', 'last_program')
mode_run = conf.get('misc', 'mode_run')


# MPI DETAILS
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# CHECK PREVIOUS PROGRAM
if rank == 0:
    if mode_run == 'automated':
        if not last_program == 'Combiner.py': 
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

# MPI OBJECT INTIALIZATION
length = np.zeros(1)
selection_dict = {}

# COUSTUMIZED LOGGER FUNCTION
formatter = logging.Formatter('%(asctime)s %(message)s')


def setup_logger(name, level=logging.INFO):
        handler = logging.FileHandler('%s/%s' % (log_dir, name))  # PATH + NAME
        handler.setFormatter(formatter)
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger


# GADGET FILE READER
def snapread(filename):
        return snr.readsnap(filename, suppress=True)


# FUNCTION FOR SELECTING GAS PARTICLES 
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


# GETTING BLACK HOLES INFO
with open('%s/bh_cat.p' % input_dir, 'rb') as fr:
        bh_cat = pl.load(fr)
with open('%s/snap_cat.p' % input_dir, 'rb') as fr:
        snap_cat = pl.load(fr)
bh_id = bh_cat[:, 0]
bh_x = bh_cat[:, 1].astype(np.float)
bh_y = bh_cat[:, 2].astype(np.float)
bh_z = bh_cat[:, 3].astype(np.float)
bh_pos_t = np.array([bh_x, bh_y, bh_z])
bh_pos = bh_pos_t.T
bh_l = len(bh_pos)

# MPI: ROOT AND LOCAL ITERATIONS [START: STOP]
local_n = int(bh_l/(size))
local_i = (size-rank-1)*int(local_n)
local_f = local_i + local_n
local_bh_pos = bh_pos[local_i:local_f]
local_bh_id = bh_id[local_i:local_f]
local_slice = len(local_bh_id)
root_f = local_i + local_n + bh_l % size
root_bh_pos = bh_pos[local_i:root_f]
root_bh_id = bh_id[local_i:root_f]
root_slice = len(root_bh_id)


files = glob.glob('%s/snap*' % snapshot_dir)
LOG = setup_logger('process%d.log' % rank)
SNP = setup_logger('snapshot.log')


# ROOT NODE
if rank == 0:

        with open('%s/data.p' % input_dir, 'rb') as fr:
                data = pl.load(fr)

        for filename in files:
                LOG.info('@%d: Reading %s' % (rank, filename[9:]))
                snap = snapread(filename)
                gas_position = snap.read_pos('GAS')
                length[0] = len(gas_position)

                filen = {0: filename}
                for i in range(1, size):
                        comm.send(filen, dest=i, tag=i)

                #  SENDING LENGTH OF GAS POSITION ARRAY
                for i in range(1, size):
                        comm.Send(length, dest=i, tag=i)
                LOG.info('@%d: Details of gas_postion array has been sent' % (
                                                                         rank))

                #  SENDING GAS POSTIONS TO OTHER NODES
                for i in range(1, size):
                        comm.Send(gas_position, dest=i, tag=i)
                LOG.info('@%d: gas_postion array has been sent' % (rank))

                #  READING ALL OTHER REQUIRED DETAILS FROM FILE : LINE-95
                gas_density = snap.read_density()
                gas_interalE = snap.read_U()
                gas_smoothingL = snap.read_hsml()
                gas_electronAb = snap.read_ne()
                LOG.info('@%d: Read gas properties from %s' % (rank,
                                                               filename[9:]))

                # SEARCHING INSIDE THE FILE FOR ROOT SLICE OF BHs: LINE 84
                for i in range(root_slice):

                    if int(root_bh_id[i]) not in snap_cat[filename]:

                        LOG.info('@%d: searching gas positions of bh %s #%d/%d'
                                 % (rank, root_bh_id[i], i+1, root_slice))
                        dummy = selection(gas_position, root_bh_pos[i])

                        if dummy.size is not 0:  # CONDITON OF PARTICLE FOUND
                                LOG.info('@%d:  ^found %d particles'
                                         % (rank, len(dummy)))
                                selected_bh = root_bh_id[i]
                                selected_bh_data = data[selected_bh]

                                bh_loc_gas_position = gas_position[dummy]
                                bh_loc_gas_density = gas_density[dummy]
                                bh_loc_gas_interalE = gas_interalE[dummy]
                                bh_loc_gas_smoothingL = gas_smoothingL[dummy]
                                bh_loc_gas_electronAb = gas_electronAb[dummy]
                                dat = np.array([bh_loc_gas_interalE,
                                                bh_loc_gas_density,
                                                bh_loc_gas_position[:, 0],
                                                bh_loc_gas_position[:, 1],
                                                bh_loc_gas_position[:, 2],
                                                bh_loc_gas_smoothingL,
                                                bh_loc_gas_electronAb])
                                dat = dat.T

                                # UPDATING DATASET
                                bh_data_updater = np.concatenate((
                                                          selected_bh_data, dat
                                                                  ), axis=0)
                                data[selected_bh] = bh_data_updater
                                LOG.info('@%d: Data updated for %s' % (
                                                          rank, root_bh_id[i]))
                    else:
                            LOG.info('skipping %d' % int(root_bh_id[i]))

                # COLLECTING ALL DETAILS OF FOUND PARTICLES FROM OTHER NODES
                for i in range(1, size):
                        selection_dict = comm.recv(source=i, tag=i)
                        LOG.info('@%d:Selection dictionary received fom node%d'
                                 % (rank, i))
                        for bh in selection_dict.keys():
                                selected_bh_data = data[bh]
                                dummy = selection_dict[bh]
                                bh_loc_gas_position = gas_position[dummy]
                                bh_loc_gas_density = gas_density[dummy]
                                bh_loc_gas_interalE = gas_interalE[dummy]
                                bh_loc_gas_smoothingL = gas_smoothingL[dummy]
                                bh_loc_gas_electronAb = gas_electronAb[dummy]
                                dat = np.array([bh_loc_gas_interalE,
                                                bh_loc_gas_density,
                                                bh_loc_gas_position[:, 0],
                                                bh_loc_gas_position[:, 1],
                                                bh_loc_gas_position[:, 2],
                                                bh_loc_gas_smoothingL,
                                                bh_loc_gas_electronAb])
                                dat = dat.T

                                bh_data_updater = np.concatenate((
                                                          selected_bh_data, dat
                                                                  ), axis=0)

                                data[bh] = bh_data_updater
                                LOG.info('@%d: Data updated for %s' % (
                                                          rank, bh))
                SNP.info('%s' % filename)
        with open('%s/dataUpd.p' % output_dir, 'wb') as fw:  # SAVING THE FINAL DATASET
                pl.dump(data, fw)

        conf.set('live', 'last_program', 'CatalogAnalyser_MPI.py')
        with open('abe.ini', 'w') as fw:
               conf.write(fw)
# LOCAL NODES
else:
        for i in range(len(files)):

                # RECEIVING LENTH OF NEXT COMING DATA OF GAS
                filen = comm.recv(source=0, tag=rank)
                filename = filen[0]
                comm.Recv(length, source=0, tag=rank)
                LOG.info('@%d: Details of gas_postion array has been recieved'
                         % (rank))
                gas_position = np.empty((int(length[0]), 3))  # EMPTY ARRAY:LINE 100

                # RECEIVING GAS POSITION ARRAY
                comm.Recv(gas_position, source=0, tag=rank)
                LOG.info('@%d: gas_postion array has been recieved' % (rank))

                # SEARCHING INSIDE THE FILE FOR LOCAL SLICE OF BHs: LINE 84
                for i in range(local_slice):

                    if int(local_bh_id[i]) not in snap_cat[filename]:

                        LOG.info('@%d: searching gas positions of bh %s #%d/%d'
                                 % (rank, local_bh_id[i], i+1, local_slice))
                        dummy = selection(gas_position, local_bh_pos[i])
                        if dummy.size is not 0:
                                LOG.info('@%d:^ found %d particles'
                                         % (rank, len(dummy)))
                                temp_dict = {local_bh_id[i]: dummy}
                                selection_dict.update(temp_dict)
                    else:
                            LOG.info('skipping %d' % int(local_bh_id[i]))
                # SENDING SELECTION DICTIONARY TO ROOT
                comm.send(selection_dict, dest=0, tag=rank)
                LOG.info('@%d: Selection dictionary has been sent to root' %
                         rank)
                selection_dict = {}
