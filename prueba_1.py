# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 18:49:50 2014

@author: martindrech
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si
#import fourier as f

c = np.loadtxt('0.75_0.24_20.txt', np.complex)
#c = np.loadtxt('.txt', np.complex)
nu = c[0]
c = c[1:]
i = 2
c = c[c.size//2-i:c.size//2+i+1]
e = -0.24
t = np.linspace(0,10, 30)
T1, T2, T3, T4 = .1, 1, 5, 10
wc = 50
g = 1
#
#d = np.array([1])
#s0e = si.sigma_baja(t, d, nu, g, wc, 0)
#s1e = si.sigma_num(t, d, nu, g, wc, 0, T1)
#s2e = si.sigma_num(t, d, nu, g, wc, 0, T2)
#s3e = si.sigma_num(t, d, nu, g, wc, 0, T3)
#s4e = si.sigma_num(t, d, nu, g, wc, 0, T4)
#
#print 'e'
#
#s0 = si.sigma_baja(t, c, nu, g, wc, e)
#s1 = si.sigma_num(t, c, nu, g, wc, e, T1)
#s2 = si.sigma_num(t, c, nu, g, wc, e, T2)
#s3 = si.sigma_num(t, c, nu, g, wc, e, T3)
#s4 = si.sigma_num(t, c, nu, g, wc, e, T4)

#plt.clf()
#plt.plot(t, s0[1], linewidth = 2, label = 'T=0')
#plt.plot(t, s1[1], linewidth = 2, label = 'T=0.1')
#plt.plot(t, s2[1], linewidth = 2, label = 'T=1')
#plt.plot(t, s3[1], linewidth = 2, label = 'T=5')
#plt.plot(t, s4[1], linewidth = 2, label = 'T=10')
#plt.xlabel('t', fontsize = 30), plt.ylabel('$\sigma_{pp}$', fontsize = 30)
#plt.legend(), plt.subplots_adjust(bottom=0.22)

plt.clf()
plt.plot(t, s2[0], 'b', linewidth = 2, label = '$\sigma_{xx}$')
plt.plot(t, s2e[0], 'b', linewidth = 2)
plt.plot(t, s2[1], 'r', linewidth = 2, label = '$\sigma_{pp}$')
plt.plot(t, s2e[1], 'r', linewidth = 2)
plt.plot(t, s2[2], 'g',linewidth = 2, label = '$\sigma_{xp}$')
plt.plot(t, s2e[2], 'g', linewidth = 2)
plt.legend(), plt.subplots_adjust(bottom=0.22)
plt.xlabel('t'), plt.title('T = 1')



print 'done'
