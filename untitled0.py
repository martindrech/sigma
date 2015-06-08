# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:04:22 2015

@author: martin
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si
import floquet as fl
from erres import R3_bajaT as R


def sigma_baja_nuevo(t, c, nu, g, wc):
    """
    Devuelve, para T=0 las dispersiones sxx, spp, sxp
    """
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    N,M,K,L, T = np.meshgrid(enes, enes, enes, enes, t)
    
    A1, A2, A3, A4, Ta = np.meshgrid(c, c, c,c, t)
    A = A1*A2*A3*A4
#    nu1 = np.conjugate(nu)
    nu1 = nu
    r =A*(R(wc, T, N, M, K, L, g, nu, nu1)-R(0, T, N, M, K, L, g, nu, nu1))
    
    sxx = np.sum(r, axis = (0,1,2,3)) *(si.B(c,nu)**2 * g / np.pi)
   
    return sxx.real

c_sci = np.loadtxt('0.75_-0.5.txt', np.complex)
c_sci = np.loadtxt('5_1.txt', np.complex)

nu_sci = c_sci[0]
c_sci = c_sci[1:]

a, q = .75, -5
nu = fl.mathieu_nu(a, q, 15)
c = fl.mathieu_coefs(a, q, nu, 11)

i = 3
c, c_sci = c[c.size//2-i:c.size//2+i+1], c_sci[c_sci.size//2-i:c_sci.size//2+i+1]

g = 1

wc = 50
t = np.linspace(0,10, 50)

plt.clf()
#sxx_sci = sigma_baja_nuevo(t, c_sci, nu_sci, g, wc)
sxx_new = sigma_baja_nuevo(t, c, nu, g, wc)
plt.plot(t, sxx_new, 'g')



