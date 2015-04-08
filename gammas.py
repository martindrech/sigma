# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:04:22 2015

@author: martin
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si

a = np.loadtxt('gamma2.dat', skiprows=0)
au_gammas = a[:, 0] * 2
au_dispX = a[:, 1]
au_dispP = a[:, 2]

gammas = au_gammas
dispX, dispP = np.array([]), np.array([])
t = np.linspace(0,10, 100)
wc = 20
c = np.array([1])
for g in gammas:
    nu = np.sqrt(1-.25 * g**2)
    
    sxx, spp, sxp = si.sigma_baja(t, c, nu, g, wc, 0)
    
    sxx, spp = np.average(sxx), np.average(spp) 
    
    dispX = np.append(dispX, [sxx])
    dispP = np.append(dispP, [spp])
    print g



plt.clf()
plt.plot(gammas , dispX, 'b', label = 'dispX', linewidth = 1)
plt.plot(gammas, dispP, 'r',label = 'dispP', linewidth = 1)
plt.plot(au_gammas, au_dispX, 'ob', label = 'dispX') 
plt.plot(au_gammas, au_dispP, 'or', label = 'dispP')
plt.xlabel('$\gamma$', size = 25)
plt.legend(loc=(0.1,.6))
plt.subplots_adjust(bottom=0.22)
#plt.savefig('/home/martin/Desktop/fig.jpg')
print 'done'
#dispX = dispX.real
#dispP = dispP.real

#A = np.matrix([gammas, dispX, dispP])
#A = np.transpose(A)
#np.savetxt('gammas.txt', A, fmt='%.3f')