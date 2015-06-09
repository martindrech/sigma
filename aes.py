# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 18:42:33 2015

@author: martin
"""

from __future__ import division
import numpy as np
import w_w as w_w
import floquet as fl
import time
from matplotlib import pylab as plt


def a1(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc):
    """
    Devuelve la matrix a1(t)    
    """
    w2w2t1 = w_w.w2_w2(t, g, temp1, nu1, c1, nu1 , c1, wc)
    w2mw2mt1 = w_w.w2_w2(t, g, temp1, nu2, c2, nu2 , c2, wc)
    w2w2mt1 = w_w.w2_w2(t, g, temp1, nu1, c1, nu2 , c2, wc)
    
    w2w2t2 = w_w.w2_w2(t, g, temp2, nu1, c1, nu1 , c1, wc)
    w2mw2mt2 =w_w.w2_w2(t, g, temp2, nu2, c2, nu2 , c2, wc)
    w2w2mt2 =w_w.w2_w2(t, g, temp2, nu1, c1, nu2 , c2, wc)
    
    a11 = w2w2t1+w2mw2mt1+2*w2w2mt1+w2w2t2+w2mw2mt2-2*w2w2mt2
    a12 = w2w2t1-w2mw2mt1+w2w2t2-w2mw2mt2
    a21 = a12
    a22 = w2w2t2+w2mw2mt2+2*w2w2mt2+w2w2t1+w2mw2mt1-2*w2w2mt1
    
    return 1/4 * np.array([[a11, a12], [a21, a22]])

def a2(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc):
    """
    Devuelve la matrix a2(t)    
    """
    w1w2t1 = w_w.w1_w2(t, g, temp1, nu1, c1, nu1, c1, wc)
    w1mw2mt1 = w_w.w1_w2(t, g, temp1, nu2, c2, nu2, c2, wc)   
    w1mw2t1 = w_w.w1_w2(t, g, temp1, nu2, c2, nu1, c1, wc)
    w1w2mt1 = w_w.w1_w2(t, g, temp1, nu1, c1, nu2, c2, wc)
    
    w1w2t2 = w_w.w1_w2(t, g, temp2, nu1, c1, nu1, c1, wc)
    w1mw2mt2 = w_w.w1_w2(t, g, temp2, nu2, c2, nu2, c2, wc)   
    w1mw2t2 = w_w.w1_w2(t, g, temp2, nu2, c2, nu1, c1, wc)
    w1w2mt2 = w_w.w1_w2(t, g, temp2, nu1, c1, nu2, c2, wc) 
    
    a11 = w1w2t1+w1mw2mt1+w1mw2t1+w1w2mt1 + w1w2t2+w1mw2mt2-w1mw2t2-w1w2mt2
    a12 = w1w2t1+w1mw2t1-w1w2mt1+w1mw2mt1 + w1w2t2+w1w2mt2-w1mw2t2-w1mw2mt2
    a21 = w1w2t2+w1mw2t2-w1w2mt2+w1mw2mt2 + w1w2t1+w1w2mt1-w1mw2t1-w1mw2mt1
    a22 = w1w2t2+w1mw2mt2+w1mw2t2+w1w2mt2 + w1w2t1+w1mw2mt1-w1mw2t1-w1w2mt1
    
    return 1/2 * np.array([[a11, a12], [a21, a22]])

def a3(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc):
    """
    Devuelve la matrix a3(t)    
    """
    w1w1t1 = w_w.w1_w1(t, g, temp1, nu1, c1, nu1 , c1, wc)
    w1mw1mt1 = w_w.w1_w1(t, g, temp1, nu2, c2, nu2 , c2, wc)
    w1w1mt1 = w_w.w1_w1(t, g, temp1, nu1, c1, nu2 , c2, wc)
    
    w1w1t2 = w_w.w1_w1(t, g, temp2, nu1, c1, nu1 , c1, wc)
    w1mw1mt2 =w_w.w1_w1(t, g, temp2, nu2, c2, nu2 , c2, wc)
    w1w1mt2 =w_w.w1_w1(t, g, temp2, nu1, c1, nu2 , c2, wc)
    
    a11 = w1w1t1+w1mw1mt1+2*w1w1mt1+w1w1t2+w1mw1mt2-2*w1w1mt2
    a12 = w1w1t1+w1w1mt1-2*w1w1mt1+w1w1t2+w1w1mt2-2*w1w1mt2
    a21 = a12
    a22 = w1w1t2+w1mw1mt2+2*w1w1mt2+w1w1t1+w1mw1mt1-2*w1w1mt1
    
    return 1/4 * np.array([[a11, a12], [a21, a22]])

def b(t, g, a1, q1, a2, q2):
    """
    Devuelve las matrix b1, b2, b3, b4
    """
    phi1, dphi1, phi2, dphi2 = fl.mathieu(a1, q1, t)
    phim1, dphim1, phim2, dphim2 = fl.mathieu(a2, q2, t)
    
    b1_11 = dphi1/phi1 + dphim1/phim1 - g 
    b1_12 = dphi1/phi1 - dphim1/phim1

    b2_11 = 1/phi1+1/phim1 
    b2_12 = 1/phi1-1/phim1

    b3_11 = 1/phi1+1/phim1 
    b3_12 = 1/phi1-1/phim1

    b4_11 =  g - phi2/phi1 - phim2/phim1
    b4_12 = phi2/phi1 + phim2/phim1
    
    b1 = 1/2 * np.array([[b1_11, b1_12], [b1_12, b1_11]])
    b2 = -1/2 * np.exp(-g*t/2) * np.array([[b2_11, b2_12], [b2_12, b2_11]])
    b3 = -1/2 * np.exp(g*t/2) * np.array([[b3_11, b3_12], [b3_12, b3_11]])    
    b4 = -1/2 * np.array([[b4_11, b4_12], [b4_12, b4_11]])
    
    return b1, b2, b3, b4
    
def cov(t, g, ca1, cq1, ca2, cq2, temp1, temp2, wc = 50, i = 3):
    """
    Devuelve <x1^2>, <x1^2>. 
    """    
    nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
    c1, c2 = fl.mathieu_coefs(ca1, cq1, nu1), fl.mathieu_coefs(ca2, cq2, nu2)
    c1, c2 = c1[c1.size//2-i:c1.size//2+i+1], c2[c2.size//2-i:c2.size//2+i+1]
    
    Ma1 = a1(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc)
    Ma2 = a2(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc)
    Ma3 = a3(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc)
    Mb1, Mb2, Mb3, Mb4 = b(t, g, ca1, cq1, ca2, cq2)
    
    det = Mb3[0][0]**2-Mb3[0][1]**2
    
    x1x1 = (1 / det**2) * ( Mb3[0][0]**2 * Ma3[0][0] -Mb3[0][0]*Mb3[0][1]*(Ma3[1][0]+Ma3[0][1]) +Mb3[1][0]**2 * Ma3[1][1] ) 
    x2x2 = (1 / det**2) * ( Mb3[0][1]**2 * Ma3[0][0] -Mb3[0][0]*Mb3[0][1]*(Ma3[1][0]+Ma3[0][1]) +Mb3[0][0]**2 * Ma3[1][1] ) 

    x1x2 = (1/ (2*det) ) * (1/2 * (Ma3[0][1]+Ma3[1][0]) - (Mb3[0][1]/Mb3[0][0]) * Ma3[1][1]) - (Mb3[0][1]/Mb3[0][0]) * x2x2
    
    
    x1p1 = (1/ (2*det) ) * (Mb3[1][0]*Ma2[0][1]-Mb3[0][0]*Ma2[0][0]) + Mb1[0][1] * x1x2 + Mb1[0][0]*x1x1   
    x2p2 = (1/ (2*det) ) * (Mb3[0][1]*Ma2[1][0]-Mb3[0][0]*Ma2[1][1]) + Mb1[0][1] * x1x2 + Mb1[0][0]*x2x2
    
    x1p2 = (1/ (2*det) ) * (Mb3[1][0]*Ma2[0][0]-Mb3[0][0]*Ma2[1][0]) + Mb1[0][1] * x1x1 + Mb1[0][0]*x1x2  
    x2p1 = (1/ (2*det) ) * (Mb3[0][1]*Ma2[0][0]-Mb3[0][0]*Ma2[0][1]) + Mb1[0][0] * x2x2 + Mb1[0][1]*x1x2
    
    p1p1 = Ma1[0][0] + Mb1[0][0]**2 * x1x1 + Mb1[0][1]**2 * x2x2 + 4*Mb1[0][0]*Mb1[0][1]*x1x2 + 2*Mb1[0][0] * x1p1 + 2*Mb1[0][1]*x2p1
    p2p2 = Ma1[1][1] + Mb1[1][0]**2 * x1x1 + Mb1[0][0]**2 * x2x2 + 4*Mb1[0][0]*Mb1[0][1]*x1x2 + 2*Mb1[0][0] * x2p2 + 2*Mb1[0][1]*x1p2
    
    p1p2 = 1/2 * (Ma1[0][1]+Ma1[1][0]) + 3*Mb1[0][0]*Mb1[0][1]* (x1x1+x2x2)+3*(Mb1[1][0]*Mb1[0][1]+Mb1[0][0]*Mb1[1][1])*x1x2
    
    plt.plot(t, Ma2[0][0])
    return x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2



    
start = time.time()  
ca1, cq1, g = .75, .5, 1
ca2, cq2 = .75, .5
nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
A1, A2 = fl.mathieu_coefs(ca1, cq1, nu1, 11), fl.mathieu_coefs(ca2, cq2, nu2, 11)
i = 3
A1, A2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
t = np.linspace(5,15, 100)

plt.clf()
T1, T2 = 10, 10
x1x1, x2x2,x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2 = cov(t, g, ca1, cq1, ca2, cq2, T1, T2, 50, i)
#plt.plot(t, x1x1, 'b-')
#plt.plot(t, x1x2, 'ro')
end = time.time()
print -start+end
#
print 'done'