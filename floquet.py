# -*- coding: utf-8 -*-
"""
Created on Wed May 13 20:07:08 2015

@author: martin
"""

from __future__ import division
from scipy.integrate import odeint
import numpy as np
import pylab as plt
import sigma as si
from scipy.linalg import eig, eigh, logm
import scipy
from scipy import linalg, matrix
import timeit
###############################################################################



def mathieu(a, q, t):
    """
    Devuelve las soluciones con ci [0,1] y [1, 0] de la ec. de Mathieu. 
    Donde la ecuacion es y'' + ( a-2q*cos(2t) ) y = 0
    """
    def V(t):
        return a-2*q*np.cos(2*t)
    
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

def mathieu_est(a_max, q_max, nro_a, nro_q):
    """
    Hace un plot medio cabeza de estabilidad.     
    """
    T = np.pi
    n = 1000
    t = np.linspace(0, T, n)
    A = np.zeros(5)
    for a in np.linspace(0, a_max, nro_a):
        print a
        for q in np.linspace(0, q_max, nro_q):
            phi1, dphi1, phi2, dphi2 = mathieu(a, q, t)
            C = np.array([[phi2[n-1], phi1[n-1]], [dphi2[n-1], dphi1[n-1]]])
            B = 1/T * logm(C)
            m1 = eig(B)[0][0]
            m2 = eig(B)[0][1]
            r1, r2 = np.abs(eig(C)[0][0]), np.abs(eig(C)[0][1])
            
            if np.abs(r1-1) < 10**-4 and np.abs(r2-1)< 10**-4:
                est = 'bo'
            else: 
                est = 'ro'
            A = np.vstack((A, np.array([a, q, m1, m2, est])))
            A = A[1:, :]
            for i in range(0, np.shape(A)[0]):
                plt.plot(A[i, 0], A[i, 1], A[i, 4], markersize = 6)
            plt.xlabel('a')
            plt.ylabel('q')
    return A

def mathieu_nu(a, q, N):
    """
    Devuelve nu para los parametros a y q. 
    """
    T = np.pi
    n = 1000
    t = np.linspace(0, T, n)
    phi1, dphi1, phi2, dphi2 = mathieu(a, q, t)
    C = np.array([[phi2[n-1], phi1[n-1]], [dphi2[n-1], dphi1[n-1]]])
    B = -1j/T * logm(C)
    nu =  eig(B)[0][0] 
   
    
    return nu
    
def H_nu(q, nu, N):
    
    n = np.arange(-(N-1)/2, (N-1)/2+1)
    d = (2*n+nu)**2 
    H = np.diagflat(d)
    H[0, 1] = q
    H[N-1, N-2] = q
    for i in range(1, N-1):
        H[i, i+1] = q
        H[i, i-1] = q
        
    return H

def mathieu_coefs(a, q, nu, N):
    """
    Devuelve un vector c, con los coeficientes del desarrollo periodico de p(t).

    """
    A = H_nu(q, nu, N) - a*np.diagflat(np.ones(N))
    autos = eig(A)

    autov_cero = np.argmin(np.abs(autos[0]))
    c = autos[1][:, autov_cero]
        
    return c



    
###############################################################################
#a_max, q_max, nro_a, nro_q = 10, 10, 10, 10
#mathieu_est(a_max, q_max, nro_a, nro_q)

#c_sci = np.loadtxt('3_-.5_10.txt', np.complex)
#nu_sci = c_sci[0]
#c_sci = c_sci[1:]
#
#a, q, N = 3, -.5, 37
#t = np.linspace(0, np.pi, 1000)
#nu= mathieu_nu(a, q, N)
#c = mathieu_coefs(a, q, nu, N)
#te = np.linspace(0, 20, 100)
#plt.clf()
#phi1_mio = si.phi1(te, c, nu)
#phi1_sci = si.phi1(te, c_sci, nu_sci)
#phi1_num = mathieu(a, q, te)[0]
#
#plt.plot(te, phi1_num.real, 'r-')
#plt.plot(te, phi1_mio.real, 'o')

"""
for N in np.arange(1, 100, 5):
    start = timeit.timeit()
    c = mathieu_coefs(a, q, nu, N)
    end = timeit.timeit()
    print "N: ", N , 'time: ', start-end
"""

#print 'done'
#a, q, N = .75, -.5, 11
#nu= mathieu_nu(a, q, N)
#
#for z, i in enumerate(np.arange(11, 80, 2)):
#    
#    c = mathieu_coefs(a, q, nu, i) 
#    c = c[:len(c)-2]
#    print i, np.argmax(np.abs(c))



#print 'done'