#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import numpy as np

class header(object):
    def __init__(self,filename):
        self.filename = filename
        self.f = open(self.filename,'rb')

    def consistent(self,start_byte,off_pos):    #off_pos is is the blocksize in bytes
        self.f.seek(off_pos+start_byte)
        end_byte = np.fromfile(self.f,dtype=np.uint32,count=1)[0]
        #print(end_byte)
        if start_byte == end_byte :
            #print('Block Consistency Check PASSED!')
            return True
        else:
            #print('Block Consistency Check FAILED!')
            return False
    
    def head_read(self,suppress=False):
        self.f.seek(0)
        start_byte = np.fromfile(self.f,dtype=np.uint32,count=1)[0]
        off_pos = self.f.tell()
        #print(start_byte)
        head_info = []
        npartThisFile = np.fromfile(self.f,dtype=np.uint32,count=6)
        head_info.append(npartThisFile)
        massTable     = np.fromfile(self.f,dtype=np.float64,count=6)
        head_info.append(massTable)
        time          = np.fromfile(self.f,dtype=np.float64,count=1)[0]
        head_info.append(time)
        redshift      = np.fromfile(self.f,dtype=np.float64,count=1)[0]
        head_info.append(redshift)
        flag_sfr      = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_sfr)
        flag_fb       = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_fb)
        npartTotal    = np.fromfile(self.f,dtype=np.uint32,count=6)
        head_info.append(npartTotal)
        flag_cool     = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_cool)
        nfiles        = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(nfiles)
        boxsize       = np.fromfile(self.f,dtype=np.float64,count=1)[0]
        head_info.append(boxsize)
        Omega0        = np.fromfile(self.f,dtype=np.float64,count=1)[0]
        head_info.append(Omega0)
        OmegaLambda   = np.fromfile(self.f,dtype=np.float64,count=1)[0]
        head_info.append(OmegaLambda)
        HubbleParam   = np.fromfile(self.f,dtype=np.float64,count=1)[0]
        head_info.append(HubbleParam)
        flag_age      = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_age)
        flag_metals   = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_metals)
        npartTotalHW  = np.fromfile(self.f,dtype=np.uint32,count=6)
        head_info.append(npartTotalHW)
        flag_entropy         = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_entropy)
        flag_doubleprecision = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_doubleprecision)
        flag_ic_info       = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_ic_info)
        lpt_scalingfactor             = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(lpt_scalingfactor)
        #print(self.f.tell())
        """
        flag_tmax            = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_tmax)
        flag_delaytime       = np.fromfile(self.f,dtype=np.int32,count=1)[0]
        head_info.append(flag_delaytime)
        H0  = HubbleParam * 100.
        Hz  = H0 * np.sqrt(OmegaLambda + Omega0 * (1.+redshift)**3)
        Hz /= 3.08567758e19  ## 1/s
        G   = 6.674e-8       ## cm^3/g/s^2
        rhocrit = 3. * Hz**2 / (8. * np.pi * G)
        """
        #Includes HighWord particles
        npart_HW = [npartTotalHW[i]<<32 for i in range(0,len(npartTotalHW))]        
        npart_inclHW = np.add(npartTotal, npart_HW)
        totnpart_inclHW = np.sum(npart_inclHW)
        
        self.header_vals = {'npartThisFile':npartThisFile,
                            'npartTotal':npartTotal,
                            'npartTotalHW':npartTotalHW,
                            'ngas':npartTotal[0],'ndm':npartTotal[1],
                            'ndisk':npartTotal[2],'nbulge':npartTotal[3],
                            'nstar':npartTotal[4],'nbh':npartTotal[5],
                            'massTable':massTable,
                            'time':time,
                            'nfiles':nfiles,
                            'redshift':redshift,
                            'boxsize':boxsize,
                            'Omega0':Omega0,
                            'Omegalambda':OmegaLambda,
                            'h':HubbleParam,
                            'flag_cooling':flag_cool,
                            'flag_sfr':flag_sfr,
                            'flag_fb':flag_fb,
                            #'flag_fh2':flag_fH2,
                            'flag_age':flag_age,
                            'flag_metals':flag_metals,
                            'flag_entropy':flag_entropy,
                            'flag_doubleprecision':flag_doubleprecision,
                            'flag_ic_info':flag_ic_info,
                            'lpt_scalingfactor':lpt_scalingfactor,
                            'npartTotal_inclHW':npart_inclHW,
                            'totalSimulationPart_inclHW':totnpart_inclHW}
        

        #self.pNames = {0:'GAS  ',
        #          1:'DM   ',
        #          2:'DISK ',
        #          3:'BULGE',
        #          4:'STAR ',
        #          5:'BNDRY'}
        if not(suppress):
            print('Header Reading Done!')
        check = self.consistent(start_byte,off_pos)
        self.f.close()
        return (self.header_vals,check)
"""
head = header('snapshot_068.0')
head.head_read(suppress=True)
import pprint
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(head.header_vals)
"""
