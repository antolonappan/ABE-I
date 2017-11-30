#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 13:37:29 2017

@author: anto
"""
import _pickle as pl
import numpy as np
import matplotlib.pyplot as plt
from mapper import map
plt.rcParams['text.usetex']=True

with open('stack.p', 'rb') as f:
    stack = pl.load(f)
    
fig, axs = plt.subplots(4, 3)

keys=list(stack.keys())

p = 0
for i in range(4):
    for j in range(3):
        axs[i, j].pcolormesh(np.log10(stack[keys[p]]+.0000000000000000000000000001), cmap='plasma')
        axs[i, j].set_title('$Bin\; %d$'%(p+1))
        axs[i, j].set_xticklabels([])
        axs[i, j].set_yticklabels([])
        if i>2:
            axs[i, j].set_xticklabels(['0','25','50'])
            axs[i, j].set_xlabel('$Kpc$')
        if j<1:
            axs[i, j].set_yticklabels(['0','50'])
            axs[i, j].set_ylabel('$Kpc$')
        p+=1

fig.subplots_adjust( hspace=1, wspace=0.5)
fig.savefig('out/stack_map.pdf',bbox_inches='tight')

p = 0
for i in range(4):
    for j in range(3):
        m=map(stack[keys[p]])
        x,y=m.profile(stack[keys[p]])
        axs[i, j].plot(x,y)
        axs[i, j].set_title('$Bin\; %d$'%(p+1))
        
        axs[i, j].set_xticklabels([])
#        axs[i, j].set_yticklabels([])
        if i>2:
            axs[i, j].set_xticklabels(['0','250','500'])

        
        p+=1

fig.subplots_adjust( hspace=1, wspace=0.5)
fig.savefig('out/stack_profile.pdf',bbox_inches='tight')
