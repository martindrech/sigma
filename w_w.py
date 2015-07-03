# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 12:50:58 2015

@author: martin
"""
from __future__ import division
import numpy as np
import sigma as si
import erres as R
def impedir_peq(arr, eps):
    mascara = np.abs(arr) < eps
    arr[mascara] = eps*np.sign(arr[mascara])

def w2_w2(t, g, temp, nu1, c1, nu2 , c2, wc, phi1, phim1):
    """
    Devuelve la integral de w2 ruido w2
    """
    enes = np.linspace(-(len(c1)-1)/2, (len(c1)-1)/2, len(c1))
    N,K,T = np.meshgrid(enes, enes, t)
    A1, A2, Ta = np.meshgrid(c1, c2, t)
    A = si.C(c1, nu1)*si.C(c2, nu2)*A1*A2
    if temp == 0:
        r =A*(R.R1_bajaT(wc, T, N, K, g, nu1, nu2)-R.R1_bajaT(0, T, N, K, g, nu1, nu2))
    else:
        r =2*temp*A*(R.R1_altaT(wc, T, N, K, g, nu1, nu2)-R.R1_altaT(0, T, N, K, g, nu1, nu2))
    
    r = np.sum(r, axis = (0,1))
    phies = phi1*phim1
    
#    impedir_peq(phies, 0.1)
    return 1/phies * (g/np.pi)  * r.real
    
def w1_w2(t, g, temp, nu1, c1, nu2 , c2, wc, phi1, phim1):
    """
    Devuelve la integral de w2 ruido w2
    """
    enes = np.linspace(-(len(c1)-1)/2, (len(c1)-1)/2, len(c1))
    N,M,K,T = np.meshgrid(enes, enes, enes, t)
    A1, A2, A3, Ta = np.meshgrid(c1, c1, c2, t)
    A = si.B(c1, nu1)*si.C(c2, nu2)*A1*A2*A3
  
    if temp == 0:
        r =A*(R.R2_bajaT(wc,T,N,M,K,g,nu1,nu2)-R.R2_bajaT(0,T,N,M,K,g,nu1,nu2))
    else:
        r =2*temp*A*(R.R2_altaT(wc,T,N,M,K,g,nu1,nu2)-R.R2_altaT(0,T,N,M,K,g,nu1,nu2))
    
    
    r = np.sum(r, axis = (0,1,2))

    phies = phi1*phim1
#    phies = 1
#    impedir_peq(phies, 0.1)
    return 1/phies * np.exp(g*t/2)*(g/np.pi)  * r.real

def w1_w1(t, g, temp, nu1, c1, nu2 , c2, wc, phi1, phim1):
    """
    Devuelve la integral de w2 ruido w2
    """
    enes = np.linspace(-(len(c1)-1)/2, (len(c1)-1)/2, len(c1))
    N,M,K,L,T = np.meshgrid(enes, enes, enes, enes, t)
    A1, A2, A3, A4, Ta = np.meshgrid(c1, c1, c2, c2, t)
    A = si.B(c1, nu1)*si.B(c2, nu2)*A1*A2*A3*A4
  
    if temp == 0:
        r =A*(R.R3_bajaT(wc,T,N,M,K,L,g,nu1,nu2)-R.R3_bajaT(0,T,N,M,K,L,g,nu1,nu2))
    else:
        r =2*temp*A*(R.R3_altaT(wc,T,N,M,K,L,g,nu1,nu2)-R.R3_altaT(0,T,N,M,K,L,g,nu1,nu2))
    
    r = np.sum(r, axis = (0,1,2,3))
    phies = phi1*phim1
    
    return 1/phies * np.exp(g*t)*(g/np.pi) *  r.real
    
###############################################################################
#import pylab as plt
#import floquet as fl
##
#ca1, cq1, g = 2, .1, 1
#ca2, cq2 = 2, .1
#nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
#A1, A2 = fl.mathieu_coefs(ca1, cq1, nu1, 11), fl.mathieu_coefs(ca2, cq2, nu2, 11)
#i = 2
#c1, c2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
#t = np.linspace(0,20, 50)
#wc = 50
#plt.clf()
#temp1, temp2 = 0,0
#phi1, dphi1, phi2, dphi2 = fl.mathieu(ca1, cq1, t)
#phim1, dphim1, phim2, dphim2 = fl.mathieu(ca2, cq2, t)
#
#w1w2t1 = w1_w2(t, g, temp1, nu1, c1, nu1, c1, wc, phi1, phi1)
#w1mw2mt1 = w1_w2(t, g, temp1, nu2, c2, nu2, c2, wc, phim1, phim1)   
#w1mw2t1 = w1_w2(t, g, temp1, nu2, c2, nu1, c1, wc, phim1, phi1)
#w1w2mt1 = w1_w2(t, g, temp1, nu1, c1, nu2, c2, wc, phi1, phim1)
#    
#w1w2t2 = w1_w2(t, g, temp2, nu1, c1, nu1, c1, wc, phi1, phi1)
#w1mw2mt2 = w1_w2(t, g, temp2, nu2, c2, nu2, c2, wc, phim1, phim1)   
#w1mw2t2 = w1_w2(t, g, temp2, nu2, c2, nu1, c1, wc, phim1, phi1)
#w1w2mt2 = w1_w2(t, g, temp2, nu1, c1, nu2, c2, wc, phi1, phim1) 
#    
#a11 = w1w2t1+w1w2mt1+w1mw2t1+w1mw2mt1 + w1w2t2-w1w2mt2-w1mw2t2+w1mw2mt2
#a12 = w1w2t1+w1mw2t1-w1w2mt1-w1mw2mt1 + w1w2t2-w1mw2t2+w1w2mt2-w1mw2mt2
#a21 = w1w2t1-w1mw2t1+w1w2mt1-w1mw2mt1 + w1w2t2+w1mw2t2-w1w2mt2-w1mw2mt2
#a22 = w1w2t1-w1w2mt1-w1mw2t1+w1mw2mt1 + w1w2t2+w1w2mt2+w1mw2t2+w1mw2mt2
#
#plt.clf()
#plt.plot(t, a12, '-o')
##top = 100 
##plt.axis([0, 20, -top, top])
#print 'done'