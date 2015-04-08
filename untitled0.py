# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:04:22 2015

@author: martin
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si


#c = np.loadtxt('0.75_-0.5.txt', np.complex)
c = np.loadtxt('0.85_0.1_10.txt', np.complex)

nu = c[0]
c = c[1:]
i = 2
c = c[c.size//2-i:c.size//2+i+1]

g = 1
e = 1
wc = 50
t = np.linspace(0,10, 80)

sxx, spp, sxp = si.sigma_baja(t, c, nu, g, wc, e)
sxxe, sppe, sxpe = si.sigma_baja_est(t, c, nu, g, wc)

plt.clf()
plt.plot(t, sxx, 'b', t, spp, 'r', t, sxp, 'g', linewidth = 3)
plt.plot(t, sxxe, 'b', t, sppe, 'r', t, sxpe, 'g', linewidth = 3)
