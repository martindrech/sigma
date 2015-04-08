# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 22:40:04 2015

@author: martin
"""

from __future__ import division
import numpy as np
import sigma as si
import pylab as plt

#c = np.loadtxt('hanggi.txt', np.complex)
c = np.loadtxt('0.85_0.1_10.txt', np.complex)

nu = c[0]
c = c[1:]
i = 0
#c = c[c.size//2-i:c.size//2+i+1]


def G_fou(s, t, nu, c):
    
     
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c)) 
    suma = np.array(np.zeros(len(s)))
    for n_i, n in enumerate(enes):
        for m_i, m in enumerate(enes):
            bn, bm = 2*n + nu, 2*m + nu 
            r = c[n_i]*c[m_i]*np.sin(bn*t-bm*s)
            suma = suma + r
    return si.B(c, nu, s) * suma

t = np.linspace(0,1)
a = si.phi1(t, c, nu)
b = si.phi2(t, c, nu)
n = 5
G = -a*b[n]+b*a[n]

plt.clf()
plt.plot(t, G, 'b', t, G_fou(t, t[n], nu, c), 'or')
