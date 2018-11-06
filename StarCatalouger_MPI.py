#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import glob
import logging
import numpy as np
import _pickle as pl
from mpi4py import MPI
from configparser import ConfigParser as cnf

pwd = os.getcwd()
sys.path.insert(0, pwd + '/Modules')
import gad_snapread as snr

conf = cnf()
conf.read('abe.ini')
snapshot_dir = conf.get('inputs', 'snapshot')
delete = conf.get('misc', 'delete_dump')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    if size%2 ==0:
        pass
    else:
        print("""Assuming the no of snapshot files is 1024,
                 the no of processor for MPI job is not a 
                 total divisor of 1024, better to run with 
                 8 or 4""")
        exit(2)

star_mass_arr, star_pos_x, star_pos_y, star_pos_z = [], [], [], []


def localread(file):
    snap = snapread(file)
    star_mass = snap.read_mass('STAR')
    star_position = snap.read_pos('STAR')
    for i in range(len(star_mass)):
        star_mass_arr.append(star_mass[i])
        star_pos_x.append(star_position[i,0])
        star_pos_y.append(star_position[i,1])
        star_pos_z.append(star_position[i,2])

files = glob.glob(snapshot_dir + '/snap*')

local_n = int(len(files)/size)
local_i = rank*int(local_n)
local_f = local_i + local_n
local_file = files[local_i:local_f]

for file in local_file:
    localread(file)

star_catalog = np.array([star_pos_x, star_pos_y, star_pos_z, star_mass_arr])
with open('star_cat%d.p' % rank, 'wb') as fp:
    pl.dump(star_catalog.T, fp)

if rank == 0:
    i = 0
    star_cat_files = glob.glob('star_cat*')
    for file in star_cat_files:
        with open(file, 'rb') as reader:
            globals()['data%d' % i] = pl.load(reader)
        i += 1
    for j in range(i):
        if j == 0:
            data = np.concatenate((globals()['data%d' % j], globals()['data%d' % (j+1)]), axis=0)
        elif j == 1:
            pass
        else:
            data = np.concatenate((data, globals()['data%d' % j]), axis=0)

    np.savetxt('star_cat.txt', data)

    if delete is 'T':
        for file in star_cat_files:
            os.remove(file)
