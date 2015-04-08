# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 10:52:57 2015

@author: martin
"""
import numpy as np
import sigma as si
import pylab as plt


c = np.loadtxt('hanggi.txt', np.complex)
nu = c[0]
c = c[1:]
i = 5
c = c[c.size//2-i:c.size//2+i+1]



def mi_G(t, s, c, nu):
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    
    suma = np.array(np.zeros(len(s))) + 0j
    for n_i, n in enumerate(enes):
        for m_i, m in enumerate(enes):
            bn, bm = 2*n + nu + 0j, 2*m + nu + 0j
            r = c[n_i]*c[m_i]*np.sin(bn*t-bm*s)
            suma = suma + r
    return si.B(c, nu) * suma
    
t = 2
s = np.linspace(0, 10, 1000)

plt.clf()
plt.plot(s, si.G(t, s, c, nu).real, 'b')
plt.plot(s, mi_G(t, s, c, nu).real, 'ro')
        
print 'done'