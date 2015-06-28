# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 22:40:04 2015

@author: martin
"""

from __future__ import division
import numpy as np
import sigma as si
import pylab as plt



def z(t, w, g, c, nu):
    def integr(s):
        return si.G(t, s, c, nu) * np.exp(g*s/2+1j*w*s)
    ret = si.complex_int(integr, 0, t)
    return ret

def zp(t, w, g, c, nu):
    def integr(s):
        return si.G(t, s, c, nu) * np.exp(g*s/2-1j*w*s)
    ret = si.complex_int(integr, 0, t)
    return ret
    
def w1w1(t, wc, g, c1, nu1, c2, nu2):
    ret = np.zeros(len(t)) + 0 * 1j
    for i, ti in enumerate(t):
        print ti
        def integr(w):
            return w*np.real(z(ti, w, g, c1, nu1) * z(ti, w, g, c2, nu2))
        int_ti = si.complex_int(integr, 0, wc)
        ret[i] = int_ti
    return ret

import floquet as fl
import w_w as wes
ca1, cq1 = .7, .5
nu1 = fl.mathieu_nu(ca1, cq1)
g = 1
ca2, cq2 = .8, .5
nu2 = fl.mathieu_nu(ca2, cq2)
A1, A2 = fl.mathieu_coefs(ca1, cq1, nu1, 11), fl.mathieu_coefs(ca2, cq2, nu2, 11)
i = 2
A1, A2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
wc = 50
t = np.linspace(0,10, 30)
phi1, dphi1, phi2, dphi2 = fl.mathieu(ca1, cq1, t)    
phim1, dphim1, phim2, dphim2 = fl.mathieu(ca2, cq2, t)


int_mia = wes.w1_w1(t, g, 0, nu1, A1, nu2, A2, wc, phi1, phim1)
int_num = w1w1(t, wc, g, A1, nu1, A2, nu2) 

plt.plot(t, int_mia, 'b')
plt.plot(t, int_num, 'g')
