# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 12:50:58 2015

@author: martin
"""
from __future__ import division
import numpy as np
import sigma as si
import floquet as fl
from matplotlib import pylab as plt
import erres as R
#para a1

def w2_w2(t, g, temp, nu1, c1, nu2 , c2, wc):
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
    return (g/np.pi) * (1/si.phi1(t, c1, nu1))*(1/si.phi1(t, c2, nu2)) * r.real
    
def w1_w2(t, g, temp, nu1, c1, nu2 , c2, wc):
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
    return np.exp(g*t/2)*(g/np.pi) * (1/si.phi1(t, c1, nu1))*(1/si.phi1(t, c2, nu2)) * r.real

def w1_w1(t, g, temp, nu1, c1, nu2 , c2, wc):
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
    return np.exp(g*t)*(g/np.pi) * (1/si.phi1(t, c1, nu1))*(1/si.phi1(t, c2, nu2)) * r.real
    
   
    
    
#a, q = .75, 5
#nu = fl.mathieu_nu(a, q)
#A = fl.mathieu_coefs(a, q, nu, 3)
#t = np.linspace(0,10, 1000)
#s = w1_w1(t, 1, 0, nu, A, nu, A, 50)
#plt.clf()
#plt.plot(t, s, '.-')
#
#print 'done'