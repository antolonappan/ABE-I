#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Style:PEP-8

import _pickle as pl
import sys
import glob
import logging
import numpy as np
from mpi4py import MPI
sys.path.insert(0, '/media/anto/Seagate Expansion Drive/mb2/snapshots/modules')
import gad_snapread as snr
formatter = logging.Formatter('%(asctime)s %(message)s')
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

__author__ = "Anto I Lonappan"
__copyright__ = "Copyright 2017, MB-II analysis Project, Presidency University"
__credit__ = ["Alankar Dutta, Rudrani K C"]
__date__ = "31-Aug-2017"
__version__ = "0.3"
__email__ = "antoidicherian@gmail.com"
__status__ = "Developing"

box = 25

verbose = False
if'-v' in sys.argv:
        verbose = True
vprint = print if verbose else lambda *a, **k: None


mass_sorting = False
if '-M' in sys.argv:
        mass_sorting = True
        M_i = sys.argv.index('-M')
def InitMass():
        MassBH = 10e7
        if (len(sys.argv)-1) > M_i:
                try:
                        MassBH = float(sys.argv[M_i+1])
                except ValueError:
                        pass
        return MassBH


def setup_logger(name, level=logging.INFO):
        handler = logging.FileHandler('log/'+name)
        handler.setFormatter(formatter)
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger

val = 'st'
guessing = False
limit = {'sf': (0.66, 1.5), 'st': (.42, 2.378), 'et': (.25, 4.)}
if '-g' in sys.argv:
        guessing = True
        G_i = sys.argv.index('-g')
        try:
              if sys.argv[G_i+1] in limit.keys():
                     val = sys.argv[G_i+1]
        except IndexError:
              pass
def guesser(mean,dat):
        lower = limit[val][0]
        upper = limit[val][1]
        # TODO : ZeroDivision Error
        try:
            lx=len(np.where(dat[:,0]<mean[0])[0])/len(np.where(dat[:,0]>mean[0])[0])
            ly=len(np.where(dat[:,1]<mean[1])[0])/len(np.where(dat[:,1]>mean[1])[0])
            lz=len(np.where(dat[:,2]<mean[2])[0])/len(np.where(dat[:,2]>mean[2])[0])
            if not (lower < lx < upper) or not (lower < ly < upper) or not (lower < lz < upper):
                    return True
            else:
                    return False
        except ZeroDivisionError:
            return True


def snapread(filename):
        return snr.readsnap(filename, suppress=True)


def selection(gas_position, bh_position, i):
        bh_box_array = np.where(((
                gas_position[:, 0] >= bh_position[0]-box) & (
                gas_position[:, 0] <= bh_position[0]+box)) & ((
                gas_position[:, 1] >= bh_position[1]-box) & (
                gas_position[:, 1] <= bh_position[1]+box)) & ((
                gas_position[:, 2] >= bh_position[2]-box) & (
                gas_position[:, 2] <= bh_position[2]+box
                       )))[0]
        return bh_box_array


Mass_dict = {}
Data_dict = {}
def localread(file):
        snap = snapread(file)
        bh_id = snap.read_id('BH')
        bh_mass = snap.read_mass('BH')
        bh_position = snap.read_pos('BH')

        ini_len = len(bh_id)
        if '-M' in sys.argv:
                cutoff_mass = InitMass()
                cutoff_index = np.where(bh_mass*10e10 > cutoff_mass)
                bh_id = bh_id[cutoff_index]
                bh_position = bh_position[cutoff_index]
                bh_mass = bh_mass[cutoff_index]
        fin_len = len(bh_id)
        SELECT.info('%d/%d BHs selected from %s'%(fin_len, ini_len, file[11:]))

        gas_position = snap.read_pos('GAS')
        gas_density = snap.read_density()
        gas_interalE = snap.read_U()
        gas_smoothingL = snap.read_hsml()
        vprint('%d/%d BHs selected from %s'%(fin_len, ini_len, file[11:]))

        for i in range(fin_len):
                vprint('     Analysing BH %d/%d in %s' % (i+1, fin_len, file[11:]))
                LOG.info('STARTED: *BH: %d #file: %s $BH_no: %d/%d' % (bh_id[i], file[11:], i+1, fin_len))
                bh_box_array = selection(gas_position, bh_position[i], i)
                bh_loc_gas_position = gas_position[bh_box_array]
                bh_loc_gas_density = gas_density[bh_box_array]
                bh_loc_gas_interalE = gas_interalE[bh_box_array]
                bh_loc_gas_smoothingL = gas_smoothingL[bh_box_array]
                if guessing is True:
                        if guesser(bh_position[i], bh_loc_gas_position):
                                GUESS.info('*BH: %d #Snapshot: %s'%(bh_id[i],file[11:]))

                LOG.info('FINISHED: #file:%s $BH:%d/%d' % (file[11:], i+1, fin_len))

                dat = np.array([bh_loc_gas_interalE, bh_loc_gas_density,
                        bh_loc_gas_position[:,0], bh_loc_gas_position[:,1], bh_loc_gas_position[:,2],
                        bh_loc_gas_smoothingL])
                temp_data = {'%d'%bh_id[i]: dat.T}
                Data_dict.update(temp_data)
#                np.save('data/%d' % (bh_id[i]), dat.T)
                temp_mass_info = {'%d'%bh_id[i] : bh_mass[i]}
                Mass_dict.update(temp_mass_info)




if __name__ == "__main__":
    files = glob.glob('snapdir_068/snap*')

    local_n = int(len(files) / size)
    local_i = rank*int(local_n)
    local_f = local_i + local_n
    local_file = files[local_i:local_f]
    LOG = setup_logger('process%d.log' % rank)
    SNAP = setup_logger('snapshot%d.log' % rank)
    GUESS = setup_logger('guess.log')
    SELECT = setup_logger('selection.log')
    for file in local_file:
            localread(file)
            SNAP.info('%s' % file[11:])

    with open('data/data%d.p' % rank, 'wb') as fp:
                pl.dump(Data_dict, fp)
    with open('data/mass%d.p' % rank, 'wb') as fp:
                pl.dump(Mass_dict, fp)
