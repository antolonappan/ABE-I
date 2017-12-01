#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 01:51:41 2017

@author: anto
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.rcParams['text.usetex']=True

# 'BIN.DAT' IS A FILE I CREATED MANUALLY 
# FROM 'BIN.STATS', YOU CAN DO THE SAME 
# JUST REMOVE THE TABULAR FORM OF 'BIN.STATS'
# AND SAVE IT WITH .DAT EXTENSION
# OR YOU CAN SIMPLY TWEAK THE CODE 'DATAANALYSIS.PY'
dat = np.loadtxt('bin.dat')

x = dat[:,0]
rangel = dat[:,1]
rangeu = dat[:,2]
no = dat[:,3]
mass = dat[:,4]

ticks=[]
for i in range(len(x)):
    ticks.append(r'$bin%d_{%.1f$-$%.1f}$' % (i+1,rangel[i],rangeu[i]))
# COLOR = ? ASSUMES 13 BINS, ADD OR REMOVE COLORS IF THERE IS ANY CHANGE 
color=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'gold','olive','pink','wheat','purple','orange']
plt.bar(x,no,color=color)
plt.xticks(x,ticks,rotation=90)
for i in range(len(color)):
    globals()['p%d'%i] = mpatches.Patch(color=color[i], label='%.2e$M_\odot$'%mass[i])

plt.legend(handles=[globals()['p%d'%i] for i in range(len(color))])
plt.xlabel('$Log_{10}(flux)$')
plt.ylabel('$No.\; of\; BHs$')
plt.title('$Mass\;Population$')
plt.savefig('out/mass_pop.pdf',bbox_inches='tight')
