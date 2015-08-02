# -*- coding: utf-8 -*-
"""
Created on Wed May 13 20:07:08 2015

@author: martin
"""

from __future__ import division
from scipy.integrate import odeint
import numpy as np
import pylab as plt
from scipy.linalg import eig,logm

###############################################################################

"""
La ecuacion la tomo con V(t) = a + e *cos(wd * t)
"""

def mathieu(wd, a, e, t):
    """
    Devuelve las soluciones con ci [0,1] y [1, 0] de la ec. de Mathieu. 
    Donde la ecuacion es y'' + ( a + e * cos(wd * t) ) y = 0
    """
    def V(t):
        return a+e*np.cos(wd*t)
    
    def harmonic(x_vec, t):
        x1, x2 = x_vec
        return [x2, -V(t)*x1]
    
    
    
    phi1_result = odeint(harmonic, [0, 1], t)
    phi2_result = odeint(harmonic, [1, 0], t)
    phi1 = phi1_result[:, 0]   
    dphi1 = phi1_result[:, 1]   
    phi2 = phi2_result[:, 0]   
    dphi2= phi2_result[:, 1]
    
    return phi1, dphi1, phi2, dphi2



def mathieu_nu(wd, a, e):
    """
    Devuelve nu para los parametros a y e. 
    """
    
    na, ne = 4*a/wd**2, -2*e/wd**2
    T = np.pi
    n = 10000
    t = np.linspace(0, T, n)
    phi1, dphi1, phi2, dphi2 = mathieu(2, na, -2*ne, t)
    C = np.array([[phi2[n-1], phi1[n-1]], [dphi2[n-1], dphi1[n-1]]])
    B = -1j/T * logm(C)
    nu =  eig(B)[0][1] 
   
    
    return nu * wd/2
    
def H_nu(wd, e, nu, N):
    
    n = np.arange(-(N-1)/2, (N-1)/2+1)
    d = (wd*n+nu)**2 
    H = np.diagflat(d)
    H[0, 1] = -e/2
    H[N-1, N-2] = -e/2
    for i in range(1, N-1):
        H[i, i+1] = -e/2
        H[i, i-1] = -e/2
        
    return H

def mathieu_coefs(wd, a, e, nu, N=21):
    """
    Devuelve un vector c, con los coeficientes del desarrollo periodico de p(t).

    """
    
    A = H_nu(wd, e, nu, N) - a*np.diagflat(np.ones(N))
    autos = eig(A)

    autov_cero = np.argmin(np.abs(autos[0]))
    c = autos[1][:, autov_cero]
        
    return c
    
    
def p(t, c, wd):
    """
    Devuelve la funcion periodica P(t)
    """
    p = np.zeros(np.size(t)) + 0j
    c = c + 0j
    n = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    if np.size(t) > 1: 
        for i, tau in enumerate(t):
            p[i] = np.dot(c, np.exp(wd*tau*n*1j))
    else: 
        p = np.dot(c, np.exp(wd*t*n*1j))
    return p

def B(c, nu, wd):
    """
    Devuelve el coef B que aparece en las dispersiones
    """
    t = np.array([0, 0.00001])
    P = p(t, c, wd)
    dp = (P[1]-P[0])/(t[1]-t[0])
    p0 = P[0]
#    print p0
    return 1j/( p0 * (dp+nu*p0*1j))

def C(c, nu, wd):
    """
    Devuelve el coef C que aparece en las phi
    """
    t = np.array([0, 0.00001])
    P = p(t, c, wd)
    dp = (P[1]-P[0])/(t[1]-t[0])
    p0 = P[0]
#    print p0
    return 1j/( dp+nu*p0*1j )
    
def y1(t, c, nu, wd):
    """
    Solucion de eq. de mathieu sin ci
    """
    return np.exp(1j*nu*t)*p(t, c, wd)

def y2(t, c, nu, wd):
    """
    Solucion de eq. de mathieu sin ci
    """
    return np.exp(-1j*nu*t)*p(-t, c, wd)
        
def phi1_c(t, c, nu, wd):
    """
    Solucion de eq. de mathieu con ci = [0,1]
    """
    t_inf= np.array([0, 0.000001])
    dp = (p(t_inf[1], c, wd)-p(t_inf[0], c, wd))/(t_inf[1]-t_inf[0])
    cte = 2*(dp+1j*nu*p(t_inf[0], c, wd))
    
    return (y1(t,c,nu, wd)-y2(t,c,nu, wd))/cte
    
def phi2_c(t, c, nu, wd):
    """
    Solucion de eq. de mathieu con ci = [1,0]
    """
    
    cte = 2*p(0, c, wd)
    return (y1(t,c,nu, wd)+y2(t,c,nu, wd))/cte



    
###############################################################################
#plt.clf()
#
#a, e, wd = 1, .3, 2
#t = np.linspace(0, 100, 1000)
#phi1, dphi1, phi2, dphi2 = mathieu(wd, a, e, t)
#plt.plot(t, phi2, '-b')
#
#nu = mathieu_nu(wd, a, e)
#print nu
#c = mathieu_coefs(wd, a, e, nu)
#phi2c = phi2_c(t, c, nu, wd)
#plt.plot(t, phi2c, 'bo')

#plt.axis([90, 100, -2,2])




