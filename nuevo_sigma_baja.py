# -*- coding: utf-8 -*-
"""
Created on Tue May 26 16:36:15 2015

@author: martin
"""

from __future__ import division
import numpy as np
import floquet as fl
import sigma as si
import timeit
from matplotlib import pylab as plt
from time import time
from R_bajaT_nueva import R as int_baja
from R_bajaT_nueva import R_alta as int_alta 
from erres import *
from R_bajaT import R


def aij(t, i, j):
    return np.cos(2*t*(i-j))
def bij(t, i, j):
    return np.sin(2*t*(i-j))
def beta(nu, n):
    return 2*n+nu

def integrar_pdivq(a, b, P, Q, r):
    """
    Calcula la integral definida de P(x)/Q(x) entre a y b. Tiene que cumplirse gr P<gr Q y todas las r_Q diferentes. 
    Param:
    P, Q: lista con coeficientes.
    r: raices Q. 
    Observacion: funciona solo si las raices r tienen parte imaginaria no nula.
    """
    
    
#    Q_deriv = np.reshape(np.array([np.polyder(coefficients) for index, coefficients in np.ndenumerate(Q)]), np.shape(Q))
    Q_deriv = np.polyder(Q)
    fracc = ( np.polyval(P, r)/np.polyval(Q_deriv, r) )
    inte = ( 0.5 * np.log((b-r.real)**2+(r.imag)**2) + 1j * np.arctan((b-r.real)/(r.imag)) ) - ( 0.5 * np.log((a-r.real)**2+(r.imag)**2) + 1j * np.arctan((a-r.real)/(r.imag)) )    
      
    
    return np.sum(fracc*inte)






#def F_nmkl(t):
#    bm, bl, anm, akl, bnm, bkl = beta(nu, M), np.conjugate(beta(nu, L)), aij(t, N, M), aij(t, K, L),bij(t, N, M), bij(t, K, L) 
#
#    r = np.array([bm+1j*g/2, -bm+1j*g/2, -bl-1j*g/2, bl-1j*g/2, -bl-1j*g/2])
#    r = [r[0],r[1],r[2],r[3]]
#    p = np.array([bkl*bnm, 1j*(akl*bnm*bl-anm*bkl*bm), g**2/4*bkl*bnm+g/2*anm*bkl*bm+g/2*akl*bnm*bl+akl*anm*bm*bkl])
#    p = p[0]+p[1]+p[2]    
#    q = np.array([np.ones(np.shape(r)), np.zeros(np.shape(r)), g**2/2-bm**2-bl**2, g**2/16+g**2/4*(bm**2+bl**2)+1j*g*(bl**2-bm**2)+bm**2*bl**2])
#    q = q[0]+q[1]+q[2]+q[3]
#    print np.shape(q)
#    F = integrar_pdivq(0, wc, p, q, r)
#    return F
    
    
def F_mala(t, wc, n, m, k, l, nu, g):
    bm, bl, anm, akl, bnm, bkl = beta(nu, m), np.conjugate(beta(nu, l)), aij(t, n, m), aij(t, k, l),bij(t, n, m), bij(t, k, l)
    r = np.array([bm+1j*g/2, -bm+1j*g/2, -bl-1j*g/2, bl-1j*g/2, -bl-1j*g/2])
   
    p = np.array([bkl*bnm, 1j*(akl*bnm*bl-anm*bkl*bm), g**2/4*bkl*bnm+g/2*anm*bkl*bm+g/2*akl*bnm*bl+akl*anm*bm*bkl])
     
    q = np.array([np.ones(np.shape(r)), np.zeros(np.shape(r)), g**2/2-bm**2-bl**2, g**2/16+g**2/4*(bm**2+bl**2)+1j*g*(bl**2-bm**2)+bm**2*bl**2])
    
    return  integrar_pdivq(0, wc, p, q, r)
    

def sigma(t):
    enes = np.arange(-2, 3)
    suma = 0
    for n_i, n in enumerate(enes):
            for m_i, m in enumerate(enes):
                for k_i, k in enumerate(enes):
                    for l_i, l in enumerate(enes):
                        ter = c[n_i]*c[m_i]*np.conjugate(c[n_i]*c[m_i])*F_mala(t, wc, n, m, k, l, g, nu)
                        suma = suma + ter
    return suma
    
#t = np.linspace(0, 5, 100)
#s = np.array([])
#for tau in t: 
#    s = np.append(s, sigma(tau))
#    print tau


def sigma_baja_nuevo(t, c, nu, g, wc):
    """
    Devuelve, para T=0 la dispersion sxx
    """
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    N,M,K,L, T = np.meshgrid(enes, enes, enes, enes, t)
    
    A1, A2, A3, A4, Ta = np.meshgrid(c, c, np.conjugate(c), np.conjugate(c), t)
    A = A1*A2*A3*A4
    nu2 = np.conjugate(nu)
    r =A*(R3_bajaT(wc, T, N, M, K, L, g, nu, nu)-R3_bajaT(0, T, N, M, K, L, g, nu, nu))
    sxx = np.sum(r, axis = (0,1,2,3)) *(si.B(c,nu)**2  * g / np.pi)
   
    return np.real(sxx)

#def sigma_alta_nuevo(t, c, nu, g, wc):
#    """
#    Devuelve, para T=inf las dispersiones sxx
#    """
#    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
#    N,M,K,L, T = np.meshgrid(enes, enes, enes, enes, t)
#    
#    A1, A2, A3, A4, Ta = np.meshgrid(c, c, np.conjugate(c), np.conjugate(c), t)
#    A = A1*A2*A3*A4
#    r =A*(int_alta(wc, T, N, M, K, L, g, nu, nu)-int_alta(0, T, N, M, K, L, g, nu, nu))
#    sxx = np.sum(r, axis = (0,1,2,3)) * si.B(c,nu)**2 * (g / np.pi)
#   
#    return sxx.real    
#
    

    
    
    
#i = 2
#
#ci = c[5-i:5+1+i]
##cj = c[5-j:5+1+j]
#t = np.linspace(0, 5)
#beg1 = time()
#new = sigma_baja_nuevo(t, ci, nu, g, wc)
#end1 = time()
#beg2 = time()
#old, tr, ty = si.sigma_baja(t,  ci, nu, g, wc, 0)
#end2 = time()
#plt.clf()
#plt.plot(t, old )
#plt.plot(t, new, 'o')
#
#print 'new time: ', end1-beg1
#print 'old time: ', end2-beg2
#    
t = np.linspace(0, 10, 100)
plt.clf()
c_sci = np.loadtxt('0.75_-0.5.txt', np.complex)
nu_sci = c_sci[0]

c = c[1:]
i = 2
c = c[c.size//2-i:c.size//2+i+1]
a, q = .75, 1
g = 1
wc = 50
nu = fl.mathieu_nu(a, q, 15)
c = fl.mathieu_coefs(a, q, nu, 11)

c = c[c.size//2-i:c.size//2+i+1]
#old, cha, cha = si.sigma_baja(t, c, nu, g,  wc, 0)
new = sigma_baja_nuevo(t, c, nu, g, wc)
#plt.plot(t, old)
plt.plot(t, new, 'og')

#for wc in [50]:
#    s = sigma_alta_nuevo(t, ci, nu, g, wc)
#    plt.plot(t, s, label = str(wc))
#    plt.legend()
    
#for j in [0, .1, .3, .5, .6]:
#    i = 2    
#    a = 1
#    q = j
#    nu = fl.mathieu_nu(a, q, 15)
#    c = fl.mathieu_coefs(a, q, nu, 11)
#    ci =  c[5-i:5+1+i]
#    beg = time()
#    new = sigma_baja_nuevo(t, ci, nu, g, wc)
#    end = time()    
#    plt.plot(t, new.real, label = str(j))
#    plt.legend()
#    print i, end-beg
    
print 'done'

