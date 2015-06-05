# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:04:22 2015

@author: martin
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si
import floquet as fl

def R3_bajaT(w, t, n, m, k, l, g, nu1, nu2):
    """
    
    """
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    re = (1j/4)*(-(((np.cos(2*t*(n-m))+1j*np.sin(2*t*(n-m)))*(g-(2*1j)*(2*m+u+1j*v))*(np.sin(2*t*(k-l))*(g-(2*1j)*m-1j*u+v)+np.cos(2*t*(k-l))*(2*l+x+1j*y))*np.log((-1j/2)*g-2*m-u-1j*v+w))/(g**2+4*l**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*l*(x+1j*y)+(2*1j)*x*y-y**2))+((np.cos(2*t*(n-m))-1j*np.sin(2*t*(n-m)))*(g+(2*1j)*(2*m+u+1j*v))*(np.sin(2*t*(k-l))*(g+1j*(2*m+u+1j*v))+np.cos(2*t*(k-l))*(2*l+x+1j*y))*np.log((-1j/2)*g+2*m+u+1j*v+w))/(g**2+4*l**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*l*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos(2*t*(k-l))-1j*np.sin(2*t*(k-l)))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g+1j*(2*l+x+1j*y)))*(g+(2*1j)*(2*l+x+1j*y))*np.log(-2*l+w-x+(1j/2)*(g-2*y)))/(g**2-4*l**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*l*(x+1j*y)+(2*1j)*g*(2*l+x+1j*y)-(2*1j)*x*y+y**2)-((np.cos(2*t*(k-l))+1j*np.sin(2*t*(k-l)))*(g-(2*1j)*(2*l+x+1j*y))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g-(2*1j)*l-1j*x+y))*np.log(2*l+w+x+(1j/2)*(g+2*y)))/(g**2-4*l**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*l*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*l-1j*x+y)))
    
    return re 

def sigma_baja_nuevo(t, c, nu, g, wc):
    """
    Devuelve, para T=0 las dispersiones sxx, spp, sxp
    """
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    N,M,K,L, T = np.meshgrid(enes, enes, enes, enes, t)
    
    A1, A2, A3, A4, Ta = np.meshgrid(c, c, np.conjugate(c), np.conjugate(c), t)
    A = A1*A2*A3*A4
    nu1 = np.conjugate(nu)
    r =A*(R3_bajaT(wc, T, N, M, K, L, g, nu, nu1)-R3_bajaT(0, T, N, M, K, L, g, nu, nu1))
    
    sxx = np.sum(r, axis = (0,1,2,3)) *(np.abs(si.B(c,nu))**2 * g / np.pi)
   
    return sxx.real

#c = np.loadtxt('0.75_-0.5.txt', np.complex)
plt.clf()
#nu = c[0]
a, q = 0.75, 0
nu = fl.mathieu_nu(a, q, 15)
c = fl.mathieu_coefs(a, q, nu, 11)
#c = c[1:]
i = 3
c = c[c.size//2-i:c.size//2+i+1]

g = .5
e = 1
wc = 50
t = np.linspace(0,10, 50)
#
#sxx, spp, sxp = si.sigma_baja(t, c, nu, g, wc, e)
#sxxe, sppe, sxpe = si.sigma_baja_est(t, c, nu, g, wc)
#
plt.clf()
#plt.plot(t, sxx.real, 'b')
#plt.plot(t, sxxe.real, 'b')

sxx_new = sigma_baja_nuevo(t,c, nu, g, wc)
plt.plot(t, sxx_new, 'g')

print 'done'
