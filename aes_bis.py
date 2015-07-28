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
import sigma as si
#import fourier as fou
#
def impedir_peq(arr, eps):
    mascara = np.abs(arr) < eps
    arr[mascara] = eps*np.sign(arr[mascara])

def a1t1(t, g,  nu1, c1, temp, nu2, c2, wc, phi1, phim1):
    """
    Devuelve la matrix a1(t)    
    """

    
    w2w2t1 = w_w.w2_w2(t, g, temp, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w2mw2mt1 = w_w.w2_w2(t, g, temp, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w2w2mt1 = w_w.w2_w2(t, g, temp, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    
    
    a11 = w2w2t1+w2mw2mt1+2*w2w2mt1 
    a12 = w2w2t1-w2mw2mt1 
    a21 = a12
    a22 = w2w2t1+w2mw2mt1-2*w2w2mt1 

    
    return 1/4 * np.array([[a11, a12], [a21, a22]])
    
def a1t2(t, g,  nu1, c1, temp, nu2, c2, wc, phi1, phim1):
    """
    Devuelve la matrix a1(t)    
    """

    w2w2t2 = w_w.w2_w2(t, g, temp, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w2mw2mt2 =w_w.w2_w2(t, g, temp, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w2w2mt2 =w_w.w2_w2(t, g, temp, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    a11 =  w2w2t2+w2mw2mt2-2*w2w2mt2
    a12 =  w2w2t2-w2mw2mt2
    a21 = a12
    a22 =  w2w2t2+w2mw2mt2+2*w2w2mt2

    
    return 1/4 * np.array([[a11, a12], [a21, a22]])

def a2t1(t, g,  nu1, c1, temp, nu2, c2, wc, phi1, phim1):
    """
    Devuelve la matrix a2(t)    
    """
    w1w2t1 = w_w.w1_w2(t, g, temp, nu1, c1, nu1, c1, wc, phi1, phi1)
    w1mw2mt1 = w_w.w1_w2(t, g, temp, nu2, c2, nu2, c2, wc, phim1, phim1)   
    w1mw2t1 = w_w.w1_w2(t, g, temp, nu2, c2, nu1, c1, wc, phim1, phi1)
    w1w2mt1 = w_w.w1_w2(t, g, temp, nu1, c1, nu2, c2, wc, phi1, phim1)
    
     
    
    a11 = w1w2t1+w1w2mt1+w1mw2t1+w1mw2mt1 
    a12 = w1w2t1+w1mw2t1-w1w2mt1-w1mw2mt1 
    a21 = w1w2t1-w1mw2t1+w1w2mt1-w1mw2mt1 
    a22 = w1w2t1-w1w2mt1-w1mw2t1+w1mw2mt1
    
    return 1/2 * np.array([[a11, a12], [a21, a22]])

def a2t2(t, g,  nu1, c1, temp, nu2, c2, wc, phi1, phim1):
    """
    Devuelve la matrix a2(t)    
    """
    
    w1w2t2 = w_w.w1_w2(t, g, temp, nu1, c1, nu1, c1, wc, phi1, phi1)
    w1mw2mt2 = w_w.w1_w2(t, g, temp, nu2, c2, nu2, c2, wc, phim1, phim1)   
    w1mw2t2 = w_w.w1_w2(t, g, temp, nu2, c2, nu1, c1, wc, phim1, phi1)
    w1w2mt2 = w_w.w1_w2(t, g, temp, nu1, c1, nu2, c2, wc, phi1, phim1) 
    
    a11 =  w1w2t2-w1w2mt2-w1mw2t2+w1mw2mt2
    a12 =  w1w2t2-w1mw2t2+w1w2mt2-w1mw2mt2
    a21 =  w1w2t2+w1mw2t2-w1w2mt2-w1mw2mt2
    a22 =  w1w2t2+w1w2mt2+w1mw2t2+w1mw2mt2
    
    return 1/2 * np.array([[a11, a12], [a21, a22]])
    
def a3t1(t, g,  nu1, c1, temp, nu2, c2, wc, phi1, phim1):
    """
    Devuelve la matrix a3(t)    
    """
    w1w1t1 = w_w.w1_w1(t, g, temp, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w1mw1mt1 = w_w.w1_w1(t, g, temp, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w1w1mt1 = w_w.w1_w1(t, g, temp, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    
    
    a11 = w1w1t1+w1mw1mt1+2*w1w1mt1
    a12 = w1w1t1-w1mw1mt1
    a21 = a12
    a22 = w1w1t1+w1mw1mt1-2*w1w1mt1
    
    
    return 1/4 * np.array([[a11, a12], [a21, a22]])
    
def a3t2(t, g,  nu1, c1, temp, nu2, c2, wc, phi1, phim1):
    """
    Devuelve la matrix a3(t)    
    """
    
    
    w1w1t2 = w_w.w1_w1(t, g, temp, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w1mw1mt2 =w_w.w1_w1(t, g, temp, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w1w1mt2 =w_w.w1_w1(t, g, temp, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    a11 = w1w1t2+w1mw1mt2-2*w1w1mt2
    a12 = w1w1t2-w1mw1mt2
    a21 = a12
    a22 = w1w1t2+w1mw1mt2+2*w1w1mt2    
    
    return 1/4 * np.array([[a11, a12], [a21, a22]])


def b(t, g, phi1, dphi1, phi2, dphi2, phim1, dphim1, phim2, dphim2):
    """
    Devuelve las matrix b1, b2, b3, b4
    """

    f1 = phi1 * np.exp(-g*t/2)
    df1 = (dphi1- (g/2) * phi1) * np.exp(-g*t/2)
    g1 = phim1 * np.exp(-g*t/2)
    dg1 = (dphim1- (g/2) * phim1) * np.exp(-g*t/2)
    
    b1d = 1/2 * (df1/f1 + dg1/g1)
    b1c = 1/2 * (df1/f1 - dg1/g1)

    b2d = -1/2 * np.exp(-g*t) * (1/f1+1/g1)
    b2c =-1/2 * np.exp(-g*t) * (1/f1-1/g1)

    b3d = -1/2  * (1/f1+1/g1)
    b3c = -1/2  * (1/f1-1/g1)

    
    return b1d, b1c, b2d, b2c, b3d, b3c
    
def cov(t, g, ca1, cq1, ca2, cq2, temp1, temp2, wc = 50, i = 5):
    """
    Devuelve todos los valores medios para la matriz de covarianza. 
    """    
    nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
    c1, c2 = fl.mathieu_coefs(ca1, cq1, nu1), fl.mathieu_coefs(ca2, cq2, nu2)
    c1, c2 = c1[c1.size//2-i:c1.size//2+i+1], c2[c2.size//2-i:c2.size//2+i+1]
 
    phi1, dphi1, phi2, dphi2 = fl.mathieu(ca1, cq1, t)    
    phim1, dphim1, phim2, dphim2 = fl.mathieu(ca2, cq2, t)
    
        

    
    Ma1 = a1t1(t, g,  nu1, c1, temp1, nu2, c2, wc, phi1, phim1)+a1t2(t, g,  nu1, c1, temp2, nu2, c2, wc, phi1, phim1)
    Ma2 = a2t1(t, g,  nu1, c1, temp1, nu2, c2, wc, phi1, phim1)+a2t2(t, g,  nu1, c1, temp2, nu2, c2, wc, phi1, phim1)
    Ma3 = a3t1(t, g,  nu1, c1, temp1, nu2, c2, wc, phi1, phim1)+a3t2(t, g,  nu1, c1, temp2, nu2, c2, wc, phi1, phim1)
    
    Ma1, Ma2, Ma3 = 2*Ma1, 2*Ma2, 2*Ma3
    
    b1d, b1c, b2d, b2c, b3d, b3c = b(t, g, phi1, dphi1, phi2, dphi2, phim1, dphim1, phim2, dphim2)

    x1x1=(Ma3[1][1]*b3c**2-(Ma3[0][1]+Ma3[1][0])*b3c*b3d+Ma3[0][0]*b3d**2)/(b3c**2-b3d**2)**2
    x2x2=(Ma3[0][0]*b3c**2-(Ma3[0][1]+Ma3[1][0])*b3c*b3d+Ma3[1][1]*b3d**2)/(b3c**2-b3d**2)**2
    x1x2=((Ma3[0][1]+Ma3[1][0])*b3c**2-2*(Ma3[0][0]+Ma3[1][1])*b3c*b3d+(Ma3[0][1]+Ma3[1][0])*b3d**2)/(2*(b3c**2-b3d**2)**2)
    x1p1=(b3c**2*((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[1][1]*b1d-Ma2[0][1]*b3c)+b3c*(-2*(Ma3[0][0]+Ma3[1][1])*b1c-2*(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[0][0]*b3c)*b3d+((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[0][0]*b1d+Ma2[0][1]*b3c)*b3d**2-Ma2[0][0]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    x2p2=(b3c**2*((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[0][0]*b1d-Ma2[1][0]*b3c)+b3c*(-2*(Ma3[0][0]+Ma3[1][1])*b1c-2*(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[1][1]*b3c)*b3d+((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[1][1]*b1d+Ma2[1][0]*b3c)*b3d**2-Ma2[1][1]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    x1p2=(b3c**2*(2*Ma3[1][1]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d-Ma2[1][1]*b3c)+b3c*(-2*(Ma3[0][1]+Ma3[1][0])*b1c-2*(Ma3[0][0]+Ma3[1][1])*b1d+Ma2[1][0]*b3c)*b3d+(2*Ma3[0][0]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[1][1]*b3c)*b3d**2-Ma2[1][0]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    x2p1=(b3c**2*(2*Ma3[0][0]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d-Ma2[0][0]*b3c)+b3c*(-2*(Ma3[0][1]+Ma3[1][0])*b1c-2*(Ma3[0][0]+Ma3[1][1])*b1d+Ma2[0][1]*b3c)*b3d+(2*Ma3[1][1]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[0][0]*b3c)*b3d**2-Ma2[0][1]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    p1p1 = Ma1[0][0]+(b3c**2*(Ma3[0][0]*b1c**2+b1d*((Ma3[0][1]+Ma3[1][0])*b1c+Ma3[1][1]*b1d)-(Ma2[0][0]*b1c+Ma2[0][1]*b1d)*b3c)-b3c*(2*(Ma3[0][0]+Ma3[1][1])*b1c*b1d+Ma3[0][1]*(b1c**2+b1d**2)+Ma3[1][0]*(b1c**2+b1d**2)-(Ma2[0][1]*b1c+Ma2[0][0]*b1d)*b3c)*b3d+(Ma3[1][1]*b1c**2+Ma3[0][1]*b1c*b1d+Ma3[1][0]*b1c*b1d+Ma3[0][0]*b1d**2+Ma2[0][0]*b1c*b3c+Ma2[0][1]*b1d*b3c)*b3d**2-(Ma2[0][1]*b1c+Ma2[0][0]*b1d)*b3d**3)/(b3c**2-b3d**2)**2
    p2p2 = Ma1[1][1]+(b3c**2*(Ma3[1][1]*b1c**2+b1d*((Ma3[0][1]+Ma3[1][0])*b1c+Ma3[0][0]*b1d)-(Ma2[1][1]*b1c+Ma2[1][0]*b1d)*b3c)-b3c*(2*(Ma3[0][0]+Ma3[1][1])*b1c*b1d+Ma3[0][1]*(b1c**2+b1d**2)+Ma3[1][0]*(b1c**2+b1d**2)-(Ma2[1][0]*b1c+Ma2[1][1]*b1d)*b3c)*b3d+(Ma3[0][0]*b1c**2+Ma3[0][1]*b1c*b1d+Ma3[1][0]*b1c*b1d+Ma3[1][1]*b1d**2+Ma2[1][1]*b1c*b3c+Ma2[1][0]*b1d*b3c)*b3d**2-(Ma2[1][0]*b1c+Ma2[1][1]*b1d)*b3d**3)/(b3c**2-b3d**2)**2
    p1p2 = (Ma1[0][1]+Ma1[1][0]+(((Ma2[0][1]+Ma2[1][0])*b1c+(Ma2[0][0]+Ma2[1][1])*b1d)*b3c-((Ma2[0][0]+Ma2[1][1])*b1c+(Ma2[0][1]+Ma2[1][0])*b1d)*b3d)/(-b3c**2+b3d**2)+((b1c**2+b1d**2)*((Ma3[0][1]+Ma3[1][0])*b3c**2-2*(Ma3[0][0]+Ma3[1][1])*b3c*b3d+(Ma3[0][1]+Ma3[1][0])*b3d**2))/(b3c**2-b3d**2)**2+(2*b1c*b1d*((Ma3[0][0]+Ma3[1][1])*b3c**2-2*(Ma3[0][1]+Ma3[1][0])*b3c*b3d+(Ma3[0][0]+Ma3[1][1])*b3d**2))/(b3c**2-b3d**2)**2)/2
    
    for i in [x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2]:
        i[0] = i[1]
    return x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2
#    cov_matrix = np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])
#    return cov_matrix

    
#start = time.time()  
ca1, cq1, g = 3, 0.5, .01
ca2, cq2 = 2, 0.5
nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
#A1, A2 = fl.mathieu_coefs(ca1, cq1, nu1), fl.mathieu_coefs(ca2, cq2, nu2)
i = 3
#A1, A2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
t = np.linspace(0,5, 50)
#
#
wc = 50
#
T1, T2 = 25, 30
#
#

covarianzas = cov(t, g, ca1, cq1, ca2, cq2, T1, T2, wc, i)
x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2 = covarianzas
##t, x1x1, p1p1, x1p1 = t[1:], x1x1[1:], p1p1[1:], x1p1[1:]
##plt.figure(2)
plt.clf()
##plt.plot(t, x2p1, 'bo-', t, x1p2, 'go-')
#plt.plot(t, si.deriv(x1x1, t)/2, 'b-o')
#plt.plot(t, x1p1, 'go-')# -*- coding: utf-8 -*-

#plt.plot(t, x2x2, 'bo-')
plt.clf()
plt.plot(t, p1p2, 'go-')
plt.plot(t, p1p1, 'r-o')
plt.plot(t, p2p2, 'bo-')
#
##top = 10
#plt.axis([0, 30, 0, 1000])
##plt.plot(t, 1000*x1x1, 'g-o')
#phi1, dphi1, phi2, dphi2 = fl.mathieu(ca1, cq1, t)    
#phim1, dphim1, phim2, dphim2 = fl.mathieu(ca2, cq2, t)
#plt.plot(t, 1000*(phi1), 'b')
#plt.plot(t, 100*(phim1), 'r')
#plt.grid()
##
print 'done'