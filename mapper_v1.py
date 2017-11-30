#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# The code is orginally  written in Matlab by, 
# Dr.Suchetana Chaterjee.
# Prof, Department of Physics,
# Presidency University, Kolkata.
# email: suchetana.physics@presiuniv.ac.in 

import numpy as np
import matplotlib.pyplot as plt
from numba import jit


__author__ = "Anto I Lonappan"
__copyright__ = "Copyright 2017, MB-II analysis Project"
__date__ = "17-Aug-2017"
__version__ = "0.0"
__email__ = "antoidicherian@gmail.com"
__status__ = "Developing"


class map(object):

   def __init__(self, A):
                self.A = A
                self.xpixels, self.ypixels = 100, 100
                self.xmax, self.xmin = np.max(self.A[:, 2]), np.min(self.A[:, 2])
                self.ymax, self.ymin = np.max(self.A[:, 3]), np.min(self.A[:, 3])
                self.lengthx, self.lengthy = self.xmax-self.xmin,self. ymax-self.ymin

   @jit
   def mapmaker(self):
    pixelsizex, pixelsizey = self.lengthx/self.xpixels, self.lengthy/self.ypixels
    numpart = np.size(self.A[:,0])
    #print(numpart)
    value_agn = np.zeros((self.xpixels,self.ypixels), dtype=np.float64)

    hmin, hmax = 0., 0.
    if (pixelsizex < pixelsizey):
        hmin=pixelsizex/2
        hmax=64*pixelsizex
    else:
        hmin=pixelsizey/2
        hmax=64*pixelsizey
    for n in range(0,numpart):
        pos0 = self.A[n,2]-self.xmin
        pos1 = self.A[n,3]-self.ymin
        h1= self.A[n,5]
        summ = 0
        if (h1< hmin):
            h=hmin
        elif(h1 > hmax):
            h = hmax
        else:
            h = h1
        nx=(h/(pixelsizex))+.5
        ny=(h/(pixelsizey))+.5
        dx=-nx
        while (dx <= nx):
            xx = dx * pixelsizex
            dy = -ny
            while (dy <= ny):
                yy = dy * pixelsizey
                r2 = xx * xx + yy * yy
                r = np.sqrt(r2)
                u = r/h
                if (r <= h):
                    if (u <= .5):
                        wk= 2.546479089470+15.278874536822*(u-1)*u*u
                    else:
                        wk = 5.092958178941*(1-u)*(1-u)*(1-u)
                else:
                    wk = 0.
                summ = summ + wk
                dy += 1
            dx += 1

        dx1=-nx
        while (dx1 <= nx):
            xxx = pos0 + dx1 * pixelsizex
            dy1=-ny
            while (dy1 <= ny):
                yyy = pos1 + dy1 * pixelsizey
                if (xxx >= 0 and yyy >= 0):
                    num = 0
                    k= int(np.ceil(xxx / pixelsizex))
                    l= int(np.ceil(yyy/pixelsizey))
                    if (k < self.xpixels):
                        if (l < self.ypixels):
                            xxxx = dx1 * pixelsizex
                            yyyy = dy1 * pixelsizey
                            r21 = xxxx*xxxx +yyyy*yyyy
                            r1 = np.sqrt(r21)
                            u1=r1/h
                            if (r1 <= h):
                                if (u1 <= 0.5):
                                    wk1=2.546479089470+15.278874536822*(u1-1)*u1*u1
                                else:
                                    wk1=5.092958178941*(1-u1)*(1-u1)*(1-u1)
                            else:
                                wk1=0.
                            num += wk1
                            #print (k,l,n)
                            #print(l)
                            value_agn [k,l] = value_agn[k,l] + ((num/summ) * ( self.A[n,1]**2.0) * (self.A[n,0]**0.5))
                dy1 += 1
            dx1 += 1
    #print(value_agn.size)
    return value_agn

   def profile(self,i):
           value_agn = self.mapmaker()
           comovevalue = 8.041
           pixdim = 1000/self.xpixels
           dr = 5
           ang, len_agn = [], []
           finalvalue = self.xpixels/2-1
           #print(finalvalue)

           for l in range(0,int(finalvalue-5)):
                   a_new = value_agn[int(self.xpixels/2):int((self.xpixels/2)+dr),int((self.ypixels/2)):int((self.ypixels/2)+dr)]
                   sa_new = np.size(a_new)
                   b_new= value_agn[int((self.xpixels/2)):int((self.xpixels/2)+(dr+1)),int((self.ypixels/2)):int((self.ypixels/2)+(dr+1))]
                   sb_new = np.size(b_new)
                   totval = (np.sum(np.sum(b_new))-np.sum(np.sum(a_new)))/(sb_new-sa_new)
                   ang.append(dr*pixdim)
                   len_agn.append(totval)
                   dr += 1

           TBL1 = [ang,len_agn]
           plt.plot(ang,len_agn)
           plt.savefig('%s.png'%i)
           plt.clf()

   def mapper(self, i):
           value_agn = self.mapmaker()
           fig = plt.figure()
           ax = fig.add_subplot(111)
           #print(value_agn)
           ax.pcolormesh(np.log10(value_agn+.0000000000000000000000000001), cmap='plasma')#[150:850,150:850]
           plt.colorbar(ax.pcolormesh(np.log10(value_agn+.0000000000000000000000000001), cmap='plasma'),label=r'$\log_{10}$(agn)')
           plt.xlabel('kpc',fontsize=15,fontname='times',fontweight='bold')
           plt.ylabel('kpc',fontsize=15,fontname='times',fontweight='bold')
           #plt.xticks(np.linspace(0,200,9),np.linspace(0,200,9))
           #plt.yticks(np.linspace(0,200,9),np.linspace(0,200,9))
           #ax.set_ytickslabel([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
           #plt.title(r'$\log_{10}$(agn)',fontsize=15,fontname='times',fontweight='bold')
           plt.savefig('%s.png'%i)
           plt.clf()

   @property
   def pixels(self):
           value_agn = self.mapmaker()
           return np.log10(value_agn+.0000000000000000000000000001)

   @property
   def avg_luminosity(self):
           p = self.pixels
           return np.average(p[0]+p[1])
