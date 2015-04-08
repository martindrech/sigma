# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:04:22 2015

@author: martin
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si


temp = np.linspace(0.01, 15, 10)
dispX, dispP = np.array([]), np.array([])
t = np.linspace(0,10, 100)
wc = 20
c = np.array([1])
g = 0.1

for T in temp:
    nu = np.sqrt(1-1/4 * g**2)
    
    sxx, spp, sxp = si.sigma_num(t, c, nu, g, wc, 0,T)
    
    sxx, spp = np.average([sxx]), np.average([spp]) 
    
    dispX = np.append(dispX, [sxx])
    dispP = np.append(dispP, [spp])
    print T


a = np.loadtxt('Temp.dat', skiprows=5)
au_T = a[:, 0]
au_dispX = a[:, 1]
au_dispP = a[:, 2]

plt.clf()
plt.plot(temp, dispX, 'ob', label = 'dispX', linewidth = 1)
plt.plot(temp, dispP, 'or',label = 'dispP', linewidth = 1)
plt.plot(au_T, au_dispX, 'ob', label = 'dispX') 
plt.plot(au_T, au_dispP, 'or', label = 'dispP')
plt.xlabel('$T$', size = 25)
plt.legend(loc=(0.1,.6))
plt.subplots_adjust(bottom=0.22)
#plt.savefig('/home/martin/Desktop/fig.jpg')
print 'done'