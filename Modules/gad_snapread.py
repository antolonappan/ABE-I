#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
MBII Blocks:

1.position (type 0, 1, 4, 5)
2.velocity (type 0, 1, 4, 5)
3.id       (type 0, 1, 4, 5)
4.masses   (type 0, 4, 5)
5.internal energy (type 0)
6.density (type 0)
7.electron abundance (type 0)
8.HI abundance (type 0)
9.smoothing length (type 0)
10.star formation rate (type 0)
11.stellar formation time (type 4)
12.metallicity (type 0, 4)
13.bh_mass (type 5)
14. bh_modot (type 5)
15.bh_nprogs (type 5)


1.position (type 0, 1, 4, 5)
2.velocity (type 0, 1, 4, 5)
3.id       (type 0, 1, 4, 5)
4.masses   (type 0, 4, 5)
5.internal energy (type 0)
6.electron abundance (type 0)
7.HI abundance (type 0)
8.smoothing length (type 0)
9.star formation rate (type 0)
10.stellar formation time (type 4)
11.metallicity (type 0, 4)
12.bh_mass (type 5)
13. bh_modot (type 5)
14.bh_nprogs (type 5)

so take bh_mass (as blackhole mass and it wont be discrete as you see in the mass field)
and bh_modot for accretion rate.

UnitMass_in_g = 1.989e+43            -->10^10 solar mass
UnitLength_in_cm = 3.085678e+21      -->1 kpc
UnitTime_in_s = 3.08568e16           -->this is derived from velocity and length
UnitVelocity_in_cm_per_s = 100000    --> km/sec

accretion rate is in mass/time unit in the above units. you will need to multiply by 
c^2 to get the bolometric luminosity.

"""
import numpy as np
import gad_head

class readsnap(object):
    def __init__(self,filename,suppress=False):
        self.filename = filename
        self.f = open(self.filename,'rb')
        self.h = gad_head.header(self.filename).head_read(suppress=1)[0]
        self.pNames = {0:'GAS',
                       1:'DM',
                       2:'DISK',
                       3:'BULGE',
                       4:'STAR',
                       5:'BH'}
        
        self.pBlocks = {0:'pos',
                        1:'vel',
                        2:'id',
                        3:'mass',
                        4:'U',
                        5:'density',
                        6:'ne',
                        7:'nHI',
                        8:'sml',
                        9:'sfr',
                        10:'sft',
                        11:'metallicity',
                        12:'M_BH',
                        13:'acr',
                        14:'nprogs_bh'}
        self.offsets = self.read_offsets(suppress=suppress)
    
    def read_offsets(self,suppress=True):
        self.f.seek(0)
        offsets = np.zeros((len(self.pBlocks)+1,),dtype=np.uint32)
        for i in range(0,len(offsets)):
            start_byte = np.fromfile(self.f,dtype=np.uint32,count=1)[0]
            #print(start_byte,self.f.tell())
            offsets[i] = self.f.tell()
            self.f.seek(start_byte,1)
            end_byte = np.fromfile(self.f,dtype=np.uint32,count=1)[0]
            #print(end_byte,self.f.tell())
            if int(start_byte)==int(end_byte)  and not(suppress):
                if i != 0:
                    print('Consistency check for Block ',self.pBlocks[i-1], ' PASSED!')
                else:
                    print('Consistency check for Block Header PASSED!')
            elif not(suppress):
                if i != 0:
                    print('Consistency check for Block ',self.pBlocks[i-1], ' FAILED!')
                else:
                    print('Consistency check for Block Header FAILED!')
            #else:
            #    print('Done!')
        self.f.seek(0)
        #self.f.seek(offsets[2])
        #print(np.fromfile(self.f,dtype=np.float64,count=3))
        return offsets

    def read_pos(self,ptype):
        self.f.seek(0)
        self.f.seek(self.offsets[1]-4)
        dataType, skipval = None, 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #off = self.f.tell()
        ntot = int(np.sum(self.h['npartThisFile']))
        if nbytes / (ntot * 4 * 3) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
            skipval = 4
        elif nbytes / (ntot * 8 * 3) == 1:
            dataType = np.float64
            skipval = 8
        else:
            print('Unrecognized Data Type!')
            return dataType
        pno = list(self.pNames.keys())[list(self.pNames.values()).index(ptype)]
        skipval = int(skipval* np.sum(self.h['npartThisFile'][:pno])*3)
        self.f.seek(skipval,1)  
        #off = self.f.tell()
        """ 
        if self.h['npartThisFile'][pno] == 0 : 
            print('Zero Particle of this type present in this File!')
            return np.array([])     
        """     
        pos = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][pno]*3)
        pos = np.reshape(pos,(self.h['npartThisFile'][pno],3))
        self.f.seek(0)
        if len(pos) == 0:
            print('No Particle of this type present!')
        return pos

    def read_vel(self,ptype):
        self.f.seek(0)
        self.f.seek(self.offsets[2]-4)
        dataType, skipval = None, 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #off = self.f.tell()
        ntot = int(np.sum(self.h['npartThisFile']))
        if nbytes / (ntot * 4 * 3) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
            skipval = 4
        elif nbytes / (ntot * 8 * 3) == 1:
            dataType = np.float64
            skipval = 8
        else:
            print('Unrecognized Data Type!')
            return dataType
        pno = list(self.pNames.keys())[list(self.pNames.values()).index(ptype)]
        skipval = int(skipval* np.sum(self.h['npartThisFile'][:pno])*3)
        self.f.seek(skipval,1)  
        #off = self.f.tell()          
        vel = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][pno]*3)
        vel = np.reshape(vel,(self.h['npartThisFile'][pno],3))
        self.f.seek(0)
        if len(vel) == 0:
            print('No Particle of this type present!')
        return vel
    
    def read_id(self,ptype):
        self.f.seek(0)
        self.f.seek(self.offsets[3]-4)
        dataType, skipval = None, 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #off = self.f.tell()
        ntot = int(np.sum(self.h['npartThisFile']))
        if nbytes / (ntot * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.uint32
            skipval = 4
        elif nbytes / (ntot * 8) == 1:
            dataType = np.uint64
            skipval = 8
        else:
            print('Unrecognized Data Type!')
            return dataType
        pno = list(self.pNames.keys())[list(self.pNames.values()).index(ptype)]
        skipval = int(skipval* np.sum(self.h['npartThisFile'][:pno]))
        self.f.seek(skipval,1)  
        #off = self.f.tell()          
        pid = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][pno])
        self.f.seek(0)
        if len(pid) == 0:
            print('No Particle of this type present!')
        return pid
    
    def read_mass(self,ptype):
        self.f.seek(self.offsets[4]-4)
        dataType, skipval, mass = None, 0, []
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #off = self.f.tell()
        #m_info = self.h['massTable']
        ntotm = 0
        for i in range(0,len(self.h['massTable'])):
            if self.h['massTable'][i] == 0:
                ntotm += self.h['npartThisFile'][i]
        if nbytes / (ntotm * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
            skipval = 4
        elif nbytes / (ntotm * 8) == 1:
            dataType = np.float64
            skipval = 8
        else:
            print('Unrecognized Data Type!')
            return dataType
        pno = list(self.pNames.keys())[list(self.pNames.values()).index(ptype)]
        if self.h['npartThisFile'][pno] == 0 : 
            print('No Particle of this type present!')
            self.f.seek(0)
            return mass
        else:
            if self.h['massTable'][pno] == 0:
                for i in range(0,pno):
                    if self.h['massTable'][i] == 0:
                        self.f.seek(int(skipval* self.h['npartThisFile'][i]),1)
            #mass = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][pno])
            #self.f.seek(0)
            #return mass
            else:
                mass = np.zeros((self.h['npartThisFile'][pno],),dtype=dataType)
                mass.fill(self.h['massTable'][pno])
                self.f.seek(0)
                print('All gas particles have same mass.')
                return mass
        #off = self.f.tell()          
        mass = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][pno])
        self.f.seek(0)
        return mass
    
    def read_U(self):
        self.f.seek(0)
        self.f.seek(self.offsets[5]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        ngas = self.h['npartThisFile'][0]
        if nbytes / (ngas * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (ngas * 8) == 1:
            dataType = np.float64
        else:
            print('Unrecognized Data Type!')
        U = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][0])
        self.f.seek(0)
        return U
   
    def read_density(self):
        self.f.seek(0)
        self.f.seek(self.offsets[6]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        ngas = self.h['npartThisFile'][0]
        if nbytes / (ngas * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (ngas * 8) == 1:
            dataType = np.float64
        else:
            print('Unrecognized Data Type!')
        rho = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][0])
        self.f.seek(0)
        return rho

    def read_ne(self):
        self.f.seek(0)
        self.f.seek(self.offsets[7]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        ngas = self.h['npartThisFile'][0]
        if nbytes / (ngas * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (ngas * 8) == 1:
            dataType = np.float64
        else:
            print('Unrecognized Data Type!')
        ne = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][0])
        self.f.seek(0)
        return ne
        
    def read_nHI(self):
        self.f.seek(0)
        self.f.seek(self.offsets[8]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        ngas = self.h['npartThisFile'][0]
        if nbytes / (ngas * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (ngas * 8) == 1:
            dataType = np.float64
        else:
            print('Unrecognized Data Type!')
        nHI = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][0])
        self.f.seek(0)
        return nHI
        
    def read_hsml(self):
        self.f.seek(0)
        self.f.seek(self.offsets[9]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        ngas = self.h['npartThisFile'][0]
        if nbytes / (ngas * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (ngas * 8) == 1:
            dataType = np.float64
        else:
            print('Unrecognized Data Type!')
        hsml = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][0])
        self.f.seek(0)
        return hsml
    
    def read_sfr(self):
        self.f.seek(self.offsets[10]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        ngas = self.h['npartThisFile'][0]
        if nbytes / (ngas * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (ngas * 8) == 1:
            dataType = np.float64
        else:
            print('Unrecognized Data Type!')
        sfr = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][0])
        self.f.seek(0)
        return sfr

    def read_sft(self):
        self.f.seek(self.offsets[11]-4)
        dataType = None
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #print(nbytes)
        nstar = self.h['npartThisFile'][4]
        if nbytes / (nstar * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (nstar * 8) == 1:
            dataType = np.float64
        else:
            #print(nbytes/(np.array(self.h['npartThisFile'])*8))
            print('Unrecognized Data Type!')
            return dataType
        sft = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][4])
        self.f.seek(0)
        return sft

    def read_metallicity(self,ptype):
        self.f.seek(0)
        self.f.seek(self.offsets[12]-4)
        dataType, skipval = None, 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #print(nbytes)
        npart = self.h['npartThisFile'][0] + self.h['npartThisFile'][4]
        if nbytes / (npart * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
            skipval = 4
        elif nbytes / (npart * 8) == 1:
            dataType = np.float64
            skipval = 8
        else:
            #print(nbytes/(np.array(self.h['npartThisFile'])*8))
            print('Unrecognized Data Type!')
            return dataType
        pno = list(self.pNames.keys())[list(self.pNames.values()).index(ptype)]
        if pno == self.h['npartThisFile'][0]:
            skipval = 0
        elif pno == self.h['npartThisFile'][4]:
            skipval = skipval*self.h['npartThisFile'][0]
        self.f.seek(skipval,1)
        met = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][pno])
        self.f.seek(0)
        return met
    
    def read_mBH(self):
        self.f.seek(0)
        self.f.seek(self.offsets[13]-4)
        dataType = None
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #print(nbytes)
        nBH = self.h['npartThisFile'][5]
        mBH = 0
        #print('4',nbytes/4,'8',nbytes/8)
        if nbytes / (nBH * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (nBH * 8) == 1:
            dataType = np.float64
        else:
            #print(nbytes/(np.array(self.h['npartThisFile'])*8))
            print('Unrecognized Data Type!')
            return dataType
            #return np.fromfile(self.f,dtype=np.float64,count=int(nbytes/8))
        mBH = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][5])
        self.f.seek(0)
        return mBH
    
    def read_mbhdot(self):
        self.f.seek(0)
        #print(self.offsets,self.offsets[12],len(self.offsets))
        self.f.seek(self.offsets[14]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #print(nbytes)
        nbh = self.h['npartThisFile'][5]
        #print('BH: ',nbh)
        if nbytes / (nbh * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.float32
        elif nbytes / (nbh * 8) == 1:
            dataType = np.float64
        else:
            #print(nbytes/(np.array(self.h['npartThisFile'])*8))
            print('Unrecognized Data Type!')
        mbhdot = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][5])
        mbhdot = mbhdot * (1.989e+43/3.08568e16) #g/s
        self.f.seek(0)
        return mbhdot

    def read_bhnprogs(self):
        self.f.seek(0)
        #print(self.offsets,self.offsets[12],len(self.offsets))
        self.f.seek(self.offsets[15]-4)
        dataType = 0
        nbytes = int(np.fromfile(self.f,dtype=np.uint32,count=1)[0])
        #print(nbytes)
        nbh = self.h['npartThisFile'][5]
        #print('BH: ',nbh)
        if nbytes / (nbh * 4) == 1:	#4 bytes = 32 bits times 3 for the 3 member array
            dataType = np.uint32
        elif nbytes / (nbh * 8) == 1:
            dataType = np.uint64
        else:
            #print(nbytes/(np.array(self.h['npartThisFile'])*8))
            print('Unrecognized Data Type!')
        bhnprogs = np.fromfile(self.f,dtype=dataType,count=self.h['npartThisFile'][5])
        self.f.seek(0)
        return bhnprogs


