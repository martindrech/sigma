# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 18:28:51 2015

@author: martin
"""

from __future__ import division
import numpy as np
import w_w as w_w
import floquet as fl
from matplotlib import pylab as plt
import sigma as si


ca1, cq1, g = 2, 0.5, 1
ca2, cq2 = 2, 10
nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
A1, A2 = fl.mathieu_coefs(ca1, cq1, nu1, 11), fl.mathieu_coefs(ca2, cq2, nu2, 11)
i = 4
A1, A2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
t = np.linspace(0,20, 100)
phi1, dphi1, phi2, dphi2 = fl.mathieu(ca1, cq1, t)    
phim1, dphim1, phim2, dphim2 = fl.mathieu(ca2, cq2, t)
f1 = phi1 * np.exp(-g*t/2)
df1 = (dphi1-g/2 * phi1) * np.exp(-g*t/2)
wc = 50

T1, T2 = 0,0

a11 = 1/2 * w_w.w1_w1(t, g, T1, nu1, A1, nu1 , A1, wc, phi1, phi1)
a12 = 1/2 * w_w.w1_w2(t, g, T1, nu1, A1, nu1 , A1, wc, phi1, phi1)
a22 = 1/2 * w_w.w2_w2(t, g, T1, nu1, A1, nu1 , A1, wc, phi1, phi1)

xx = 2*f1**2*a11
xp = 2*(f1*df1*a11+f1*a12)
pp = 2*(df1**2*a11+2*df1*a12+a22)


#plt.clf()
#plt.plot(t, pp, 'g', t, xx, 'b', t, xp, 'r')
#plt.axis([0, 20, 0, 10])
plt.plot(t, 1/2 * a22, 'b-o')
print 'done'