# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 22:57:02 2015

@author: martin
"""

from __future__ import division
import floquet as fl
import numpy as np
import pylab as plt
from numpy.linalg import slogdet
###############################################################################
"""
La idea aca es resolver completamente el caso de 2 osciladores con acoplamiento c(t) = c0 + c1 cos(wd * t)
y principalmente, estudiar la evolucion temporal del entrelazamiento.

Lo parametros seran c0, c1, wd, w, y el estado inicial. 
"""
def En(t, Mcov):
    """
    Negatividad logaritmica
    """
    var = np.array([])
    def my_det(A):
        a1, a2 = slogdet(A)
        return a1*np.exp(a2)
        
    for i, ti in enumerate(t):
        
        A = Mcov[:2, :2, i]
        B = Mcov[2:, 2:, i]
        C = Mcov[:2, 2:, i]
        
        delta = my_det(A)+my_det(B)-2*my_det(C)
        dete = np.abs(my_det(Mcov[:, :, i]))
        nmenos2 = (delta-np.sqrt(delta**2-4*dete)) / 2
        nmenos = np.sqrt(nmenos2)
        var = np.append(var, nmenos)
        
    var = - np.log2(2*var)

    var[var < 0] = 0

    return var
    
def det_cov(t, Mcov):
    var = np.array([])
    def my_det(A):
        a1, a2 = slogdet(A)
        return a1*np.exp(a2)
    for i, ti in enumerate(t):
        dete = my_det(Mcov[:, :, i])
        var = np.append(var, dete)
    
    return var
        
            
            
def entrelazamiento(t, w, wd, c0, c1):
    
    wM0 = np.sqrt(w**2+c0+c1)
    wm0 = np.sqrt(w**2-c0-c1) 
    
    sxx0M = 1/2*wM0 
    spp0M = wM0/2
    sxx0m = 1/2*wm0 
    spp0m = wm0/2
    
    
    phi1M, dphi1M, phi2M, dphi2M = fl.mathieu(wd, w**2+c0, c1, t)
    phi1m, dphi1m, phi2m, dphi2m = fl.mathieu(wd, w**2-c0, -c1, t)
    
    sxxM = spp0M * phi1M**2 +sxx0M * phi2M**2
    sppM = dphi1M**2 * spp0M + dphi2M**2 * sxx0M
    sxpM = spp0M * phi1M * dphi1M + sxx0M * phi2M * dphi2M
    sxxm = spp0m * phi1m**2 +sxx0m * phi2m**2
    sppm = dphi1m**2 * spp0m + dphi2m**2 * sxx0m
    sxpm = spp0m * phi1m * dphi1m + sxx0m * phi2m * dphi2m
    
    x1x1 = 1/2 * (sxxM+sxxm)
    p1p1 = 1/2 * (sppM+sppm)
    x1p1 = 1/2 * (sxpM+sxpm)
    x2x2 = x1x1
    p2p2 = p1p1
    x2p2 = x1p1
    x1x2 = 1/2 * (sxxM-sxxm)
    p1p2 = 1/2 * (sppM-sppm)
    x1p2 = 1/2 * (sxpM-sxpm)
    x2p1 = x1p2
    
    Mcov= np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])
#    print Mcov[:, :, 500]
    neg = En(t, Mcov)
    detes = det_cov(t, Mcov) 
    
    return neg, sxxM, sxxm, sppM, sppm, detes
    
t = np.linspace(0, 100, 10000)
w = 1
wd = 2.05
c0 = 0.1
c1 = 0.11

nuM = fl.mathieu_nu(wd,  w**2+c0, c1)
num = fl.mathieu_nu(wd,  w**2-c0, -c1) 
print 'nu mas: ', nuM.imag
print 'nu menos: ', num.imag
neg, sxxM, sxxm, sppM, sppm, detes = entrelazamiento(t, w, wd, c0, c1)

plt.clf()
plt.subplot(1, 2, 1)
plt.plot(t, neg)
plt.grid()
plt.subplot(1, 2, 2)
plt.plot(t, sxxm, 'b',  t, sppm, 'r')
#plt.plot(t, detes)
#plt.axis([190, 210, 0, 100000000])
plt.grid()

print 'done'