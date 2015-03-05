# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 18:49:50 2014

@author: martindrech
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si

def deriv(x, t):
    dt = t[1]-t[0]
    dx = np.diff(x)/dt
    dx = np.append([dx[0]], dx)
    return dx

c = np.loadtxt('0.85_0.1_10.txt', np.complex)
#c = np.loadtxt('.txt', np.complex)
nu = c[0]
c = c[1:]
i = 1
c = c[c.size//2-i:c.size//2+i+1]

t = np.linspace(0, 10, 100)

wc = 50
g = 2*np.sqrt(1-0.85)

sxx_baja = si.sigma_baja(t, c, nu, g, wc)
sxx_baja_est = si.sigma_baja(t, np.array([1]), nu, g, wc)

sxp = deriv(sxx_baja, t)
sxp_dt = deriv(sxp, t)



spp = (1+2*np.cos(2*t))*sxx_baja-si.S_pp_baja(t, c, nu, g, wc)+sxp_dt+g*sxp

spp = np.real(spp)
spp_est = sxx_baja_est-si.S_pp_baja(t, np.array([1]), nu, g, wc)

plt.clf()
plt.plot(t, np.real(sxx_baja), 'b', t, np.real(sxp), 'r',t, np.real(spp), 'g' )
plt.plot(t, np.real(sxx_baja_est), 'b', t, np.real(spp_est), 'g') 




plt.show()
print 'done'
