# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 19:00:54 2015

@author: martin
"""

from __future__ import division
import numpy as np
import floquet as fl
from matplotlib import pylab as plt
import sigma as si
import aes as aes
from numpy.linalg import det

w = 1
c_1 = .5
c_0 = 0
g = .005
T = 30

ca1, cq1 = w**2-g**2/4+c_0-2*g, -.5*c_1
ca2, cq2 = w**2-g**2/4-c_0-2*g, .5*c_1
nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
#print np.abs(nu1.imag/g), np.abs(nu2.imag/g)
i = 3
t = np.linspace(0,5, 30)
wc = 50

T1, T2 = T, T

def En(t, Mcov):
    var = np.array([])
    for i, ti in enumerate(t):
        
        A = Mcov[:2, :2, i]
        B = Mcov[2:, 2:, i]
        C = Mcov[2:, :2, i]
        
        delta = det(A)+det(B)-2*det(C)
        dete = det(Mcov[:, :, i])
        nmenos = np.sqrt((delta-np.sqrt(delta**2-4*dete))/2) 
#        print nmenos
        var = np.append(var, nmenos)
    var = -np.log2(2*var)
#    print var[0]
    
    var[var < 0] = 0

    return var
#plt.clf()
temps = np.arange(25, 35, 2)
negs = np.array([])
for T in temps:
    covarianzas = aes.cov(t, g, ca1, cq1, ca2, cq2, 25, T, wc, i)
    neg = En(t, covarianzas)
#    neg = np.average(neg)
#    negs = np.append(negs, neg)
    plt.plot(t, neg, 'o-', label = '$\Delta T=$'+str(T))    
    print T
plt.legend()
#plt.clf()
#plt.plot(temps, negs, 'o-', markersize=8, linewidth = 3)
plt.xlabel('t', size = 35)
plt.ylabel('$E_N$', size = 35)
#plt.text(30, 3,  str(np.abs(nu2.imag/g)//1) , fontsize = 30)
#print np.abs(nu2.imag/g)
#covarianzas = aes.cov(t, g, ca1, cq1, ca2, cq2, T1, T2, wc , i )
#neg = En(t, covarianzas)

#plt.clf()
#plt.plot(t, neg, 'o-')
#plt.axis([0,3, 5, 6])
#x1x2, x1p2, x2p1, x1x1 = covarianzas[0, 2, :], covarianzas[0, 3, :], covarianzas[1, 2, :], covarianzas[0, 0, :]
#plt.plot(t, x1x2)
#plt.plot(t, x1x2, 'o')
print 'done'
#phi1, dphi1, phi2, dphi2 = mathieu(ca1, cq1, t)
#phim1, dphim1, phim2, dphim2 = mathieu(ca2, cq2, t)
#plt.plot(t, phi1, t, phim1)
"""
temps = np.arange(0, 50, 3)
negs = np.array([])
for T in temps:
    covarianzas = aes.cov(t, g, ca1, cq1, ca2, cq2, T, T, wc, i)
    neg = En(t, covarianzas)
    neg = np.average(neg)
    negs = np.append(negs, neg)
    print T
plt.clf()
plt.plot(temps, negs, 'o-', linewidth = 2)
plt.text(30, .3, str(np.abs(nu2.imag/g)//1), fontsize = 30)
"""