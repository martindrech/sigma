# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 18:50:56 2015

@author: martin
"""

from scipy.integrate import odeint
import numpy as np
import pylab as plt
import sigma as si
#
c = np.loadtxt('0.75_-0.5.txt', np.complex)
nu = c[0]
c = c[1:]
#i = 5
#c = c[c.size//2-i:c.size//2+i+1]
a = 0.75
q = -0.5
g = 1

def V(t):
    return a-2*q*np.cos(2*t)
    
def harmonic(x_vec, t):
    x1, x2 = x_vec
    return [x2, -V(t)*x1]
    
ci = [0,1]
t = np.linspace(0, 100, 1000)
y_result = odeint(harmonic, ci, t)
y = y_result[:, 0]
x = y*np.exp(-g*t/2)
plt.clf()
#plt.plot(t, y_result[:, 0])
plt.plot(t, x,'ob')
plt.plot(t, si.phi1(t,c,nu).real*np.exp(-t*g/2), 'b')



print 'done'