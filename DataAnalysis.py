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
  ^   |1.Select a bin          |>->->    /  <-|5.Send the Stacked matrix|   ^
  |   |2.Divide the BHs into   |        /  /  |            |            |   |
  |   |           n No.of Nodes|       /  /   |            |            |   |
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
__date__ = "21-Nov-2017"
__version__ = "1.0"
__email__ = "antoidicherian@gmail.com"
__url__ = "https://github.com/antolonappan/Automated-Black-Hole-Explorer-from-GADGET-2"
__status__ = "Completed"

import math
import logging
import numpy as np
from mapper import map
import _pickle as pl
import matplotlib.pyplot as plt
from mpi4py import MPI


# ######################## MPI DETAILS ########################################
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# ################## COUSTUMIZED LOGGER FUNCTION ##############################
formatter = logging.Formatter('%(asctime)s %(message)s')


def setup_logger(name, level=logging.INFO):
        handler = logging.FileHandler('logs/'+name)  # PATH + NAME
        handler.setFormatter(formatter)
        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger


PRO = setup_logger('process.log')
MP = setup_logger('MP/MassPopulation.log')
LUM = setup_logger('LUM/Luminosity%d.log' % rank)
PLT = setup_logger('PLT/Pots%d.log' % rank)

# ######################### DATA READS ########################################
with open('data/dataUpd.p', 'rb') as fr:
        data = pl.load(fr)

with open('data/mass.p', 'rb') as fr:
        mass = pl.load(fr)

with open('data/acc.p', 'rb') as fr:
        acc = pl.load(fr)

lumsel = {}
noofbin = np.zeros(1)
stack = np.zeros([100, 100])
agn_avg_c = np.zeros(1)
lum_mass_acc_dict = {}
stack_dict = {}
# ############################ ROOT ###########################################
if rank == 0:
        PRO.info('Mass population binning initiated from root: See MassPopulation.log')
        # ********************* Mass Population *******************************
        if not mass.keys() == data.keys() == acc.keys():
                raise Exception("""MatchError:
                                  Check Keys in dataUpd.p, mass.p and acc.p""")
        else:
                MP.info('Data match test passed')
        Log10Mass = np.log10(np.array(list(mass.values()))*1e10)
        minMass = math.floor(np.min(Log10Mass))
        maxMass = math.ceil(np.max(Log10Mass))

        StepSize = 0.2
        NoOfbins = (maxMass - minMass)/StepSize

        if not NoOfbins.is_integer():
                raise Exception('No of Bins is not a whole no.')
        else:
                NoOfbins = int(NoOfbins)
        MP.info('%d No of bin found' % NoOfbins)
        bins = {}
        for i in range(NoOfbins):
                locals()['bin%d' % i] = []
        statfile = open("bin.stat", "w")
        statfile.write("Bin Scale : Log10(mass code unit * 1e10)\n")
        statfile.write("Step Size : %f\n" % StepSize)
        statfile.write("No. of Bins : %d\n\n" % NoOfbins)
        statfile.write("| Bin |   Range    | No.of BHs |     Avg Mass      |\n")
        for i in range(NoOfbins):
                start = minMass+(StepSize*i)
                stop = minMass+(StepSize*(i+1))
                massavg = []
                for key, value in mass.items():
                        if start < np.log10(value*1e10) < stop:
                                locals()['bin%d' % i].append(key)
                                massavg.append(value*1e10)
                if not len(massavg) == 0:
                        statfile.write("|{:2d}   |({}, {:4}) |    {:3}    | {:18f}|\n".format(i+1,start, stop,len(massavg), np.average(massavg)))
                        temp_bin = {'bin%d' % i: locals()['bin%d' % i]}
                        bins.update(temp_bin)
                else:
                        MP.info('Excluding bin%d having range(%f-%f): No BHs found in this range' % (i, start, stop))
        statfile.close()
        with open('.bin.p', 'wb') as fw:
                pl.dump(bins, fw)
        MP.info('Bin data saved to a hidden file `.bin.p`')
        MP.info('Mass population binning finished.')

        # *************************** Luminosity ******************************
        PRO.info('Luminosity computing initiated. See: Luminosity.log')
        for i in range(1, size):
                noofbin[0] = len(bins)
                comm.Send(noofbin, dest=i, tag=i)
                LUM.info('No of bins pickle has been sent to Node%d' % i)
        for b in bins.keys():
                bh_ids = list(bins[b])
                mass_p = [mass[mid] for mid in bh_ids]
                acc_p = [acc[aid] for aid in bh_ids]
                LUM.info('Mass Acc average calculated for %s' % b)
                lenbh = len(bh_ids)
                njobs = int(lenbh/size)
                mpilum = lenbh % size
                rootsel = njobs+mpilum
                LUM.info('%s is Selected' % b)
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
                for id in bh_ids[:rootsel]:
                        M = map(data[id])
                        value_agn = M.mapmaker()
                        avg_agn = M.avg_luminosity(value_agn)
                        value_agns += value_agn
                        avg_agns += avg_agn
                        dis += 1
                        LUM.info('Finished %s  %d/%d' % (id, dis, disl))
                stackr = value_agns
                avg_lumr = avg_agns

                for pro in range(1, size):
                        comm.Recv(stack, source=pro, tag=pro)
                        LUM.info('Received stack from Node %d' % pro)
                        stackr += stack
                stack_avg = stackr/lenbh

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
        with open('stack.p', 'wb') as fw:
                pl.dump(stack_dict, fw)
        with open('tab.p', 'wb') as fw:
                pl.dump(lum_mass_acc_dict, fw)
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
                avg_agn = 0
                dis = 0
                disl = len(lumsel[0])
                for id in lumsel[0]:
                        M = map(data[id])
                        value_agn = M.mapmaker()
                        avg_agn += M.avg_luminosity(value_agn)
                        value_agns += value_agn
                        dis += 1
                        LUM.info('Finished %s  %d/%d' % (id, dis, disl))
                stack = value_agns
                agn_avg_c[0] = avg_agn
                comm.Send(stack, dest=0, tag=rank)
                comm.Send(agn_avg_c, dest=0, tag=rank)
