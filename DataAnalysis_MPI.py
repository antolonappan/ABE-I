#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:42:34 2017

@author: anto

MPI IMPLEMENTATION
******************


      __________________________              ___________________________ 
      |         ROOT           |              |         NODES           |
      |         ----           |              |         -----           |
      | ********************** |              |                         |
      | ** Mass Population *** |              |                         |
      | ********************** |              |                         |
      |1.Set the bin size      |              |                         |
      |2.Check all data set    |              |                         |
      |3.Bin the BHs           |              |                         |
      |4.Save the bin infos    |              |                         |
      |5.Send no of bins       |>->->->->->->-|1.Receive No. of bins    |
      |                        |              |                         |
  /-->---->----\Bins/          |              |              <------<------<\
  |   |        | 2 |           |              |             /           |   |
  |   |        | 1 |           |              |        No Of Bins       |   |
  |   |        \_0_/           |              |                         |   |
  ^   |          |             |        ->->->|1.Receive BHs            |   ^
  |   |   ******************   |       /      |2.Find the average lum   |   |
  |   |   *** Luminosity ***   |      /       |3.Stack the Lum Matrix   |   |
  |   |   ******************   |     /    -<-<|4.Send the Avg lum array |   |
  ^   |1.Select a bin          |    /    /  <-|5.Send the Stacked matrix|   ^
  |   |2.Divide the BHs into   |   /    /  /  |            |            |   |
  |   |  n No.of Nodes & send  |>-     /  /   |            |            |   |
  |   |3.Find the Avg Lum      |      /  /    |            |            |   |
  ^   |4.Stack the Lum Matrix  |     /  /     |            |            |   ^
  |   |5.Receive Avglum        |<-<-<  /      |            |            |   |
  |   |6.Receive Stacked matrix|<-<-<-<       |            |            |   |
  |   |7.Find and save the     |              |            |            |   |
  ^   |      Avg Lum and stack |              |            |            |   ^
  |   |____________|___________|              |____________|____________|   |
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
__date__ = "23-Nov-2017"
__version__ = "1.0"
__email__ = "antoidicherian@gmail.com"
__url__ = "https://github.com/antolonappan/ABE-I"
__status__ = "Completed"

# PYTHON LIBRARIES
import os
import sys
import math
import logging
import numpy as np
import _pickle as pl
from mpi4py import MPI
from configparser import ConfigParser as cnf
#  *WRITTEN LIBRARIES
pwd = os.getcwd()
sys.path.insert(0, pwd + '/Modules')
from mapper_v2 import map


# ######################## MPI DETAILS ########################################
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# CONFIGURATION
conf = cnf()
conf.read('abe.ini')
bin_size = float(conf.get('misc', 'bin'))
input_dir_mak = conf.get('live', 'cat_mak_out')
input_dir_anl = conf.get('live', 'cat_anl_out')
output_dir = conf.get('live', 'dat_anl_out')
log_dir = conf.get('live', 'dat_anl_log')
last_program = conf.get('live', 'last_program')

# CHECK PREVIOUS PROGRAM
if rank == 0:
    if not last_program == 'DataEditor.py': 
          print("""
                I found Last program that you have run was {}
                Before running this code please run DataEditor.py.
                Also see abe.ini
                """.format(last_program))
          raise Exception('PreviousRunError')
else:
    pass

# ################## COUSTUMIZED LOGGER FUNCTION ##############################
formatter = logging.Formatter('%(asctime)s %(message)s')


def setup_logger(name, level=logging.INFO):
        handler = logging.FileHandler('%s/%s' % (log_dir, name))  # PATH + NAME
        handler.setFormatter(formatter)
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger


PRO = setup_logger('process.log')
MP = setup_logger('MassPopulation.log')
LUM = setup_logger('Luminosity%d.log' % rank)

# ######################### DATA READS ########################################
with open('%s/dataEdt.p' % input_dir_anl, 'rb') as fr:
        data = pl.load(fr)

with open('%s/mass.p' % input_dir_mak, 'rb') as fr:
        mass = pl.load(fr)

with open('%s/acc.p' % input_dir_mak, 'rb') as fr:
        acc = pl.load(fr)
# ################# MPI COMMUNICATING VARIABLES ###############################
lumsel = {}
noofbin = np.zeros(1)
stack = np.zeros([100, 100])
agn_avg_c = np.zeros(1)
lum_mass_acc_dict = {}
stack_dict = {}
lumino = {}
# ############################ ROOT ###########################################
if rank == 0:
        PRO.info('Mass population binning initiated from root: See MassPopulation.log')
        # ********************* Mass Population *******************************
        # CHECK ALL DATA HAVE SAME BLACKHOLES
        if not mass.keys() == data.keys() == acc.keys():
                raise Exception("""MatchError:
                                  Check Keys in dataUpd.p, mass.p and acc.p""")
        else:
                MP.info('Data match test passed')
        # FINDING BIN DETAILS          
        Log10Mass = np.log10(np.array(list(mass.values()))*1e10)
        minMass = math.floor(np.min(Log10Mass))  # LOWER BIN VALUE
        maxMass = math.ceil(np.max(Log10Mass))  # UPPER BIN VALUE

        StepSize = bin_size # BIN SIZE
        NoOfbins = (maxMass - minMass)/StepSize
        
        # CHECK VALIDITY OF BIN SIZE
        if not NoOfbins.is_integer():
                raise Exception('No of Bins is not a whole no.')
        else:
                NoOfbins = int(NoOfbins)
        MP.info('%d No of bin found' % NoOfbins)
        bins = {}
        for i in range(NoOfbins):
                locals()['bin%d' % i] = []
        statfile = open("%s/bin.stat" % output_dir, "w")
        statfile.write("Bin Scale : Log10(mass code unit * 1e10)\n")
        statfile.write("Step Size : %f\n" % StepSize)
        statfile.write("No. of Bins : %d\n\n" % NoOfbins)
        statfile.write("| Bin |   Range    | No.of BHs |     Avg Mass      |\n")
        for i in range(NoOfbins):
                # BINNING 
                start = minMass+(StepSize*i)
                stop = minMass+(StepSize*(i+1))
                massavg = []
                for key, value in mass.items():
                        if start < np.log10(value*1e10) < stop:
                                locals()['bin%d' % i].append(key)
                                massavg.append(value*1e10)
                # SELECTING NON EMPTY BINS
                if not len(massavg) == 0:
                        statfile.write("|{:2d}   |({}, {:4}) |    {:3}    | {:18f}|\n".format(i+1,start, stop,len(massavg), np.average(massavg)))
                        temp_bin = {'bin%d' % i: locals()['bin%d' % i]}
                        bins.update(temp_bin)
                # REMOVING EMPTY BINS
                else:
                        MP.info('Excluding bin%d having range(%f-%f): No BHs found in this range' % (i, start, stop))
        statfile.close()
        # SAVING BIN INFO
        with open('%s/bin.p' % output_dir, 'wb') as fw:
                pl.dump(bins, fw)
        MP.info('Bin data saved to a file `bin.p`')
        MP.info('Mass population binning finished.')

        # *************************** Luminosity ******************************
        PRO.info('Luminosity computing initiated. See: Luminosity.log')
        # SEND NO OF BINS TO OTHER NODES
        for i in range(1, size):
                noofbin[0] = len(bins)
                comm.Send(noofbin, dest=i, tag=i)
                LUM.info('No of bins pickle has been sent to Node%d' % i)
        # ANALYSING EACH BIN
        for b in bins.keys():
                # DIVIDING JOBS FOR NODES
                bh_ids = list(bins[b])
                mass_p = [mass[mid] for mid in bh_ids]
                acc_p = [acc[aid] for aid in bh_ids]
                LUM.info('Mass Acc average calculated for %s' % b)
                lenbh = len(bh_ids)
                njobs = int(lenbh/size)
                mpilum = lenbh % size
                rootsel = njobs+mpilum
                LUM.info('%s is Selected' % b)
                # SENDING JOBS(BHs) TO OTHE NODES
                for i in range(1, size):
                        startl = (njobs*i)+mpilum
                        stopl = (njobs*(i+1))+mpilum
                        lumsel = {0: bh_ids[startl: stopl]}
                        comm.send(lumsel, dest=i, tag=i)
                        LUM.info('Selection dict containing BH ids has been sent to node %d' % i)
                value_agns = np.zeros([100, 100])
                avg_agns = 0
                dis = 0
                disl = len(bh_ids[:rootsel])
                # COMPUTING AVG LUM AND STACK FOR ROOT SELECTION
                for id in bh_ids[:rootsel]:
                    if len(data[id]) != 0:
                        M = map(data[id])
                        value_agn = M.mapmaker()
                        avg_agn = M.avg_luminosity(value_agn)
                        temp_lumino = {id: avg_agn}
                        lumino.update(temp_lumino)
                        value_agns += value_agn
                        avg_agns += avg_agn
                        dis += 1
                        LUM.info('Finished %s  %d/%d' % (id, dis, disl))
                stackr = value_agns
                avg_lumr = avg_agns
                # RECEIVE AVG LUM
                for pro in range(1, size):
                        comm.Recv(stack, source=pro, tag=pro)
                        LUM.info('Received stack from Node %d' % pro)
                        stackr += stack
                stack_avg = stackr/lenbh
                # RECEIVE STACK 
                for pro in range(1, size):
                        comm.Recv(agn_avg_c, source=pro, tag=pro)
                        LUM.info('Received Avg Luminosity from Node %d' % pro)
                        avg_lumr += agn_avg_c[0]
                lum_avg = avg_lumr/lenbh
                mass_avg = sum(mass_p)/lenbh
                acc_avg = sum(acc_p)/lenbh
                temp_lum_mass_acc = {b: [lum_avg, mass_avg, acc_avg]}
                temp_stack = {b: stack_avg}
                lum_mass_acc_dict.update(temp_lum_mass_acc)
                stack_dict.update(temp_stack)
                LUM.info('Updated stack and table')
        # SAVE DATA
        with open('%s/stack.p' % output_dir, 'wb') as fw:
                pl.dump(stack_dict, fw)
        with open('%s/tab.p' % output_dir, 'wb') as fw:
                pl.dump(lum_mass_acc_dict, fw)
        with open('%s/lum%d.p' % (output_dir, rank), 'wb') as fwt:
                pl.dump(lumino, fwt)
        LUM.info('Saved stack.p and tab.p')
# ########################### NODES ###########################################
else:
        # ************************* Luminosity *******************************
        comm.Recv(noofbin, source=0, tag=rank)
        LUM.info('Received no of bins %d' % noofbin[0])
        for i in range(int(noofbin[0])):
                lumsel = comm.recv(source=0, tag=rank)
                LUM.info('Selection dict received from root')
                value_agns = np.zeros([100, 100])
                avg_agns = 0
                dis = 0
                disl = len(lumsel[0])
                for id in lumsel[0]:
                   if len(data[id]) != 0:
                        M = map(data[id])
                        value_agn = M.mapmaker()
                        avg_agn = M.avg_luminosity(value_agn)
                        temp_lumino = {id: avg_agn}
                        lumino.update(temp_lumino)
                        avg_agns += avg_agn
                        value_agns += value_agn
                        dis += 1
                        LUM.info('Finished %s  %d/%d' % (id, dis, disl))
                stack = value_agns
                agn_avg_c[0] = avg_agns
                comm.Send(stack, dest=0, tag=rank)
                comm.Send(agn_avg_c, dest=0, tag=rank)
        with open('%s/lum%d.p' % (output_dir, rank), 'wb') as fwt:
                pl.dump(lumino, fwt)
