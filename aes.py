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
    a12 = w1w1t1-w1mw1mt1+w1w1t2-w1mw1mt2
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
start = time.time()  
a, q, g = .75, .5, 1
nu = fl.mathieu_nu(a, q)
A = fl.mathieu_coefs(a, q, nu, 11)
i = 3
A = A[A.size//2-i:A.size//2+i+1]
t = np.linspace(0,10, 50)
#s1 = a1(t, 1, nu, A, 0, nu, A, 0, 50)
#s2 = a2(t, 1, nu, A, 0, nu, A, 0, 50)
s3 = a3(t, 1, nu, A, 0, nu, A, 0, 50)
b1, b2, b3, b4 = b(t, g, a, q, a, q)
plt.clf()
sxx =  ( 1/ (b3[0][0]**2) )**2 * b3[1][1]**2 * s3[0][0]
plt.plot(t, sxx)
end = time.time()
print -start+end

print 'done'