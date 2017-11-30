#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 02:44:01 2017

@author: anto
"""
import matplotlib.pyplot as plt
import _pickle as pl
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

plt.rcParams['text.usetex'] = True


def pc(m, c):
        return pearsonr(m, c)


def fit_func(x, m, c):
    return m*x + c


with open('tab.p', 'rb') as f:
    tab = pl.load(f)
with open('stack.p', 'rb') as f:
    stack = pl.load(f)
acc_cov = 1.989e+43/3.08568e16
acc_cov = 1./acc_cov
acc_cov = acc_cov*10.2214

dat = np.loadtxt('bin.dat')
no = dat[:, 3]

lumst = []
for skey in stack.keys():
    p = np.log10(stack[skey]+.0000000000000000000000000001)
    lumst.append(np.average(p[0]+p[1]))

mass = []
acc = []
lum = []
Sigma = []
i=0
for key in tab.keys():
    dat = tab[key]
    lum.append(dat[0])
    mass.append(np.log10(dat[1]*1e10))
    acc.append(np.log10(dat[2]*acc_cov))
    Sigma.append(1./no[i])
    i += 1
# ##################### Accretion - Luminosity (Stack) ########################
params = curve_fit(fit_func, acc, lumst)
param = params[0]
perr = np.sqrt(np.diag(params[1]))
plt.errorbar(acc, lumst, yerr=Sigma, fmt='o')
bf = [((acc[i]*param[0])+param[1]) for i in range(len(acc))]
bf_l = [((acc[i]*param[0])+param[1]-perr[1]) for i in range(len(acc))]
bf_u = [((acc[i]*param[0])+param[1]+perr[1]) for i in range(len(acc))]
plt.plot(acc, bf, label='$Best\; fit$')
plt.fill_between(acc, bf_l, bf_u, color='grey', alpha='0.3')
plt.xlabel('$Log_{10}(Accretion\;M_\odot yr^{-1})$')
plt.ylabel('$Log_{10}(flux)$')
plt.text(-2.5, -33, '$PPMCC : %f$' % pc(acc, lumst)[0])
plt.text(-2.5, -33.25, '$Slope : %f$' % param[0])
plt.text(-2.5, -33.5, '$Intercept : %f$' % param[1])
plt.legend()
plt.title('$Accretion$-$Luminosity\;(Stack)$')
plt.savefig('acc_lum_stack.pdf', bbox_inches='tight')
plt.clf()

# ##################### Mass - Luminosity (Stack) ########################
params = curve_fit(fit_func, mass, lumst)
param = params[0]
perr = np.sqrt(np.diag(params[1]))
plt.errorbar(mass, lumst, yerr=Sigma, fmt='o')
bf = [((mass[i]*param[0])+param[1]) for i in range(len(acc))]
bf_l = [((mass[i]*param[0])+param[1]-perr[1]) for i in range(len(acc))]
bf_u = [((mass[i]*param[0])+param[1]+perr[1]) for i in range(len(acc))]
plt.plot(mass, bf, label='$Best\; fit$')
plt.fill_between(mass, bf_l, bf_u, color='grey', alpha='0.3')
plt.xlabel('$Log_{10}(Mass\;M_\odot)$')
plt.ylabel('$Log_{10}(flux)$')
plt.text(7.1, -33, '$PPMCC : %f$' % pc(mass, lumst)[0])
plt.text(7.1, -33.25, '$Slope : %f$' % param[0])
plt.text(7.1, -33.5, '$Intercept : %f$' % param[1])
plt.legend()
plt.title('$Mass$-$Luminosity\;(Stack)$')
plt.savefig('mass_lum_stack.pdf', bbox_inches='tight')
plt.clf()
# ##################### Accretion - Luminosity ########################
params = curve_fit(fit_func, acc, lum)
param = params[0]
perr = np.sqrt(np.diag(params[1]))
plt.errorbar(acc, lum, yerr=Sigma, fmt='o')
bf = [((acc[i]*param[0])+param[1]) for i in range(len(acc))]
bf_l = [((acc[i]*param[0])+param[1]-perr[1]) for i in range(len(acc))]
bf_u = [((acc[i]*param[0])+param[1]+perr[1]) for i in range(len(acc))]
plt.plot(acc, bf, label='$Best\; fit$')
plt.fill_between(acc, bf_l, bf_u, color='grey', alpha='0.3')
plt.xlabel('$Log_{10}(Accretion\;M_\odot yr^{-1})$')
plt.ylabel('$Log_{10}(flux)$')
plt.text(-2.5, -33, '$PPMCC : %f$' % pc(acc, lum)[0])
plt.text(-2.5, -33.5, '$Slope : %f$' % param[0])
plt.text(-2.5, -34, '$Intercept : %f$' % param[1])
plt.legend()
plt.title('$Accretion$-$Luminosity$')
plt.savefig('acc_lum.pdf', bbox_inches='tight')
plt.clf()
'''
'''
# ##################### Mass - Luminosity ########################
params = curve_fit(fit_func, mass, lum)
param = params[0]
perr = np.sqrt(np.diag(params[1]))
plt.errorbar(mass, lum, yerr=Sigma, fmt='o')
bf = [((mass[i]*param[0])+param[1]) for i in range(len(acc))]
bf_l = [((mass[i]*param[0])+param[1]-perr[1]) for i in range(len(acc))]
bf_u = [((mass[i]*param[0])+param[1]+perr[1]) for i in range(len(acc))]
plt.plot(mass, bf, label='$Best\; fit$')
plt.fill_between(mass, bf_l, bf_u, color='grey', alpha='0.3')
plt.xlabel('$Log_{10}(Mass\;M_\odot)$')
plt.ylabel('$Log_{10}(flux)$')
plt.text(7.1, -33, '$PPMCC : %f$' % pc(mass, lum)[0])
plt.text(7.1, -33.5, '$Slope : %f$' % param[0])
plt.text(7.1, -34, '$Intercept : %f$' % param[1])
plt.legend()
plt.title('$Mass$-$Luminosity$')
plt.savefig('mass_lum.pdf', bbox_inches='tight')
plt.clf()


# ##################### Mass - Accretion ########################
params = curve_fit(fit_func, mass, acc)
param = params[0]
perr = np.sqrt(np.diag(params[1]))
plt.errorbar(mass, acc, yerr=Sigma, fmt='o')
bf = [((mass[i]*param[0])+param[1]) for i in range(len(acc))]
bf_l = [((mass[i]*param[0])+param[1]-perr[1]) for i in range(len(acc))]
bf_u = [((mass[i]*param[0])+param[1]+perr[1]) for i in range(len(acc))]
plt.plot(mass, bf, label='$Best\; fit$')
plt.fill_between(mass, bf_l, bf_u, color='grey', alpha='0.3')
plt.xlabel('$Log_{10}(Mass\;M_\odot)$')
plt.ylabel('$Log_{10}(Accretion\;M_\odot yr^{-1})$')
plt.text(7.1, 1.5, '$PPMCC : %f$' % pc(mass, acc)[0])
plt.text(7.1, 1, '$Slope : %f$' % param[0])
plt.text(7.1, .5, '$Intercept : %f$' % param[1])
plt.legend()
plt.title('$Mass$-$Accretion$')
plt.savefig('mass_acc.pdf', bbox_inches='tight')
plt.clf()

