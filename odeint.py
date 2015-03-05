# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 18:50:56 2015

@author: martin
"""

from scipy.integrate import odeint
import numpy as np
import pylab as plt
import sigma as si

c = np.loadtxt('12_-5_100.txt', np.complex)
nu = c[0]
c = c[1:]
a = 12
q = -5

def V(t):
    return a-2*q*np.cos(2*t)
    
def harmonic(x_vec, t):
    x1, x2 = x_vec
    return [x2, -V(t)*x1]
    
ci = [1,0]
t = np.linspace(0, np.pi, 10000)
y_result = odeint(harmonic, ci, t)

plt.clf()
#plt.plot(t, y_result[:, 0])
plt.plot(t, y_result[:, 0],'b', t, si.phi2(t,c,nu), 'o')
#plt.plot(t, si.phi1(t,c,nu))



print 'done'