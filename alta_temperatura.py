# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 19:46:10 2015

@author: martin
"""


from __future__ import division
import numpy as np
import sigma as si
import floquet as fl
import aes as aes
from matplotlib import pylab as plt


def S_xx(t, g, temp, nu1, c1):
    
    enes = np.linspace(-(len(c1)-1)/2, (len(c1)-1)/2, len(c1))
    N,M,K,L,T = np.meshgrid(enes, enes, enes, enes, t)
    A1, A2, A3, A4, Ta = np.meshgrid(c1, c1, c1, c1, t)
    A = si.B(c1, nu1)**2*A1*A2*A3*A4
    
    def integ(t, g, nu1, n, m, k, l):
        return ((g*np.cos(((2*k+nu1)+(2*n+nu1))*t)-((2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)+(2*n+nu1))*t))/((2*l+nu1)**2+2*(2*l+nu1)*(2*m+nu1)+(2*m+nu1)**2+g**2)-(g*np.cos((2*k+nu1)*t-(2*n+nu1)*t)-((2*l+nu1)-(2*m+nu1))*np.sin((2*k+nu1)*t-(2*n+nu1)*t))/((2*l+nu1)**2-2*(2*l+nu1)*(2*m+nu1)+(2*m+nu1)**2+g**2))/2+(np.exp(g*t)*((g*np.cos((2*k+nu1)*t-(2*l+nu1)*t+(2*m+nu1)*t-(2*n+nu1)*t)-((2*l+nu1)-(2*m+nu1))*np.sin((2*k+nu1)*t-(2*l+nu1)*t+(2*m+nu1)*t-(2*n+nu1)*t))/((2*l+nu1)**2-2*(2*l+nu1)*(2*m+nu1)+(2*m+nu1)**2+g**2)-(g*np.cos(((2*l+nu1)+(2*m+nu1))*t-((2*k+nu1)+(2*n+nu1))*t)+((2*l+nu1)+(2*m+nu1))*np.sin(((2*l+nu1)+(2*m+nu1))*t-((2*k+nu1)+(2*n+nu1))*t))/((2*l+nu1)**2+2*(2*l+nu1)*(2*m+nu1)+(2*m+nu1)**2+g**2)))/2
        
    r =A*integ(T, g, nu1, N, M, K, L)
    r = np.sum(r, axis = (0,1,2,3))
    return 2*temp*g * r * np.exp(-g*t)
    
def sigma(t, g, temp, ca, cq, ci, V, i = 3):
    """
    V(t) = a-2q*cos(2*t)
    """
    nu = fl.mathieu_nu(ca, cq)
    C = fl.mathieu_coefs(ca, cq, nu)
    C = C[C.size//2-i:C.size//2+i+1]
    phi1, dphi1, phi2, dphi2 = fl.mathieu(ca, cq, t)
    sxx0, sxp0, spp0 = ci[0], ci[1], ci[2]
    
    
    sxx = np.exp(-g*t)*( (phi2-g/2*phi1)**2*sxx0 + 2*phi1*(phi2-g/2*phi1)*sxp0 + phi1**2 * spp0 )+S_xx(t, g, temp, nu, C)
    sxp = si.deriv(sxx, t)/2
    phi1, dphi1, phi2, dphi2 = fl.mathieu(ca1, cq1, t)    
#    f1 = phi1 * np.exp(-g*t/2)
#    df1 = (dphi1-g/2 * phi1) * np.exp(-g*t/2)

    spp = si.deriv(sxp, t)+g*sxp+(V-2*g)*sxx    
    
    return sxx, spp, sxp
    
w0 = 1
c = 0
delta = 0.1
g = .2
ca1, cq1 = w0**2-2*g-g**2/4+c, -.5 * delta
ca2, cq2 = w0**2-2*g-g**2/4-c, -.5 * delta
nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)

i = 2
#A1, A2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
t = np.linspace(20,30, 1000)
V = w0**2 + delta*np.cos(2*t)
wc = 50 
T = 100
T1, T2 = T, T


#x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2 = aes.cov(t, g, ca1, cq1, ca2, cq2, T1, T2, wc, i, unpacked=True)
#Mcov = np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])

ci = np.array([0, 0, 0])
sxx, spp, sxp = sigma(t, g, T, ca1, cq1, ci, V, i)

#spp2 = si.deriv(x1p1, t)+g*sxp+(V-2*g)*sxx

plt.clf()
plt.plot(t, spp, '-og')
#plt.plot(t, p1p1, 'r')
#plt.plot(t, spp2, 'b')
#plt.plot(t, -.18*np.cos(2*t), 'ro-')
#plt.axis([40, 50, 80, 120])

print 'Estable: ', nu1.imag <= g/2
print 'done'

