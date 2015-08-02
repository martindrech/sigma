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


def S_xx(t, g, temp, nu1, c1, deriv = 0):
    
    enes = np.linspace(-(len(c1)-1)/2, (len(c1)-1)/2, len(c1))
    N,M,K,L,T = np.meshgrid(enes, enes, enes, enes, t)
    A1, A2, A3, A4, Ta = np.meshgrid(c1, c1, c1, c1, t)
    A = si.B(c1, nu1)**2*A1*A2*A3*A4
    
    
    def integ(t, g, nu1, n, m, k, l):
        if deriv == 0:
            return (((-(g*np.cos(((2*k+nu1)-(2*n+nu1))*t))+((2*l+nu1)-(2*m+nu1))*np.sin(((2*k+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)+np.exp(g*t)*((g*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+(-(2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(g*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+((2*l+nu1)+(2*m+nu1))*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+(g*np.cos(((2*k+nu1)+(2*n+nu1))*t)-((2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)+(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))/(2*np.exp(g*t)))
        elif deriv ==1:
            return -(g*((-(g*np.cos(((2*k+nu1)-(2*n+nu1))*t))+((2*l+nu1)-(2*m+nu1))*np.sin(((2*k+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)+np.exp(g*t)*((g*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+(-(2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(g*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+((2*l+nu1)+(2*m+nu1))*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+(g*np.cos(((2*k+nu1)+(2*n+nu1))*t)-((2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)+(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2)))/(2*np.exp(g*t))+((((2*l+nu1)-(2*m+nu1))*((2*k+nu1)-(2*n+nu1))*np.cos(((2*k+nu1)-(2*n+nu1))*t)+((2*k+nu1)-(2*n+nu1))*g*np.sin(((2*k+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)+np.exp(g*t)*g*((g*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+(-(2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(g*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+((2*l+nu1)+(2*m+nu1))*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+np.exp(g*t)*(((-(2*l+nu1)+(2*m+nu1))*((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)-((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*g*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(((2*l+nu1)+(2*m+nu1))*(-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)-(-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*g*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+(-(((2*l+nu1)+(2*m+nu1))*((2*k+nu1)+(2*n+nu1))*np.cos(((2*k+nu1)+(2*n+nu1))*t))-((2*k+nu1)+(2*n+nu1))*g*np.sin(((2*k+nu1)+(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))/(2*np.exp(g*t))
        else:
            return (g**2*((-(g*np.cos(((2*k+nu1)-(2*n+nu1))*t))+((2*l+nu1)-(2*m+nu1))*np.sin(((2*k+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)+np.exp(g*t)*((g*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+(-(2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(g*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+((2*l+nu1)+(2*m+nu1))*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+(g*np.cos(((2*k+nu1)+(2*n+nu1))*t)-((2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)+(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2)))/(2*np.exp(g*t))+((((2*k+nu1)-(2*n+nu1))**2*g*np.cos(((2*k+nu1)-(2*n+nu1))*t)-((2*l+nu1)-(2*m+nu1))*((2*k+nu1)-(2*n+nu1))**2*np.sin(((2*k+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)+np.exp(g*t)*g**2*((g*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+(-(2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(g*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+((2*l+nu1)+(2*m+nu1))*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+np.exp(g*t)*((-(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))**2*g*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))-(-(2*l+nu1)+(2*m+nu1))*((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))**2*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(-((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))**2*g*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))-((2*l+nu1)+(2*m+nu1))*(-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))**2*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+2*np.exp(g*t)*g*(((-(2*l+nu1)+(2*m+nu1))*((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)-((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*g*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(((2*l+nu1)+(2*m+nu1))*(-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)-(-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*g*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+(-(((2*k+nu1)+(2*n+nu1))**2*g*np.cos(((2*k+nu1)+(2*n+nu1))*t))+((2*l+nu1)+(2*m+nu1))*((2*k+nu1)+(2*n+nu1))**2*np.sin(((2*k+nu1)+(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))/(2*np.exp(g*t))-(g*((((2*l+nu1)-(2*m+nu1))*((2*k+nu1)-(2*n+nu1))*np.cos(((2*k+nu1)-(2*n+nu1))*t)+((2*k+nu1)-(2*n+nu1))*g*np.sin(((2*k+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)+np.exp(g*t)*g*((g*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+(-(2*l+nu1)+(2*m+nu1))*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(g*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)+((2*l+nu1)+(2*m+nu1))*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+np.exp(g*t)*(((-(2*l+nu1)+(2*m+nu1))*((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*np.cos(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)-((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*g*np.sin(((2*k+nu1)-(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)-(2*m+nu1))**2+g**2)-(((2*l+nu1)+(2*m+nu1))*(-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*np.cos((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t)-(-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*g*np.sin((-(2*k+nu1)+(2*l+nu1)+(2*m+nu1)-(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2))+(-(((2*l+nu1)+(2*m+nu1))*((2*k+nu1)+(2*n+nu1))*np.cos(((2*k+nu1)+(2*n+nu1))*t))-((2*k+nu1)+(2*n+nu1))*g*np.sin(((2*k+nu1)+(2*n+nu1))*t))/(((2*l+nu1)+(2*m+nu1))**2+g**2)))/np.exp(g*t)
    r =A*integ(T, g, nu1, N, M, K, L)
    r = np.sum(r, axis = (0,1,2,3))
    
    return 2*temp*g * r * 0
    
def sigma(t, g, temp, a, e, wd, ci, V, i = 2):
    """
    V(t) = a + e * cos(wd*t)
    """
    nu = fl.mathieu_nu(wd, a, e)
    C = fl.mathieu_coefs(wd, a, e, nu)
    C = C[C.size//2-i:C.size//2+i+1]
    phi1, dphi1, phi2, dphi2 = fl.mathieu(wd, a, e, t)
    dt = t[1]-t[0]
    #d2phi1, d2phi2 = np.gradient(dphi1, dt), np.gradient(dphi2, dt)
    sxx0, sxp0, spp0 = ci[0], ci[1], ci[2]
#    dsxx0 = -((g*(spp0*phi1**2+2*sxp0*phi1*(-(g*phi1)/2+phi2)+sxx0*(-(g*phi1)/2+phi2)**2))/np.exp(g*t))+(2*spp0*phi1*dphi1+2*sxp0*(-(g*phi1)/2+phi2)*dphi1+2*sxp0*phi1*(-(g*dphi1)/2+dphi2)+2*sxx0*(-(g*phi1)/2+phi2)*(-(g*dphi1)/2+dphi2))/np.exp(g*t)
#    d2sxx0 = (g**2*(spp0*phi1**2+2*sxp0*phi1*(-(g*phi1)/2+phi2)+sxx0*(-(g*phi1)/2+phi2)**2))/np.exp(g*t)-(2*g*(2*spp0*phi1*dphi1+2*sxp0*(-(g*phi1)/2+phi2)*dphi1+2*sxp0*phi1*(-(g*dphi1)/2+dphi2)+2*sxx0*(-(g*phi1)/2+phi2)*(-(g*dphi1)/2+dphi2)))/np.exp(g*t)+(spp0*(2*dphi1**2+2*phi1*d2phi1)+2*sxp0*(2*dphi1*(-(g*dphi1)/2+dphi2)+(-(g*phi1)/2+phi2)*d2phi1+phi1*(-(g*d2phi1)/2+d2phi2))+sxx0*(2*(-(g*dphi1)/2+dphi2)**2+2*(-(g*phi1)/2+phi2)*(-(g*d2phi1)/2+d2phi2)))/np.exp(g*t) 
        
    
    sxx = np.exp(-g*t)*( (phi2-g/2*phi1)**2*sxx0 + 2*phi1*(phi2-g/2*phi1)*sxp0 + phi1**2 * spp0 )+S_xx(t, g, temp, nu, C)
    sxx = sxx.real
   
#    dsxx = dsxx0 + S_xx(t, g, temp, nu, C, deriv = 1)
#    d2sxx = d2sxx0 + S_xx(t, g, temp, nu, C, deriv = 2) 
    
    sxp = np.gradient(sxx, dt)/2
    spp = np.gradient(sxp, dt)+g*sxp+(V-2*g)*sxx    

#    sxp = 1/2 * dsxx
#    spp = 1/2 * d2sxx +g*sxp+(V-2*g)*sxx    

    
    return sxx.real, spp.real, sxp.real

def cov_galve(t, w0, g, c0, c1, T, ci1, ci2 ,wd, i = 2, unpacked = False):
    """
    2 osciladores, alta temperatura, ambos entornos a T
    """
    a1, e1 = w0**2-g**2/4+c0-2*g, c1
    a2, e2 = w0**2-g**2/4-c0-2*g, -c1
    nu1, nu2 = fl.mathieu_nu(wd, a1, e1), fl.mathieu_nu(wd, a2, e2)    
    
    
    print 'Limite: ', np.abs(nu1.imag/g)
    print 'Estable: ', nu1.imag - g/2 < 0, nu2.imag - g/2 < 0

    Vmas = w0**2 + c0 +  c1 * np.cos(wd*t)
    Vmenos = w0**2 - c0 - c1 * np.cos(wd*t) 

    sxxM, sppM, sxpM = sigma(t, g, T, a1, e1, wd, ci1, Vmas, i)
    sxxm, sppm, sxpm = sigma(t, g, T, a2, e2, wd, ci2, Vmenos, i)

    x1x1 = 1/2 * (sxxM+sxxm)
    p1p1 = 1/2 * (sppM+sppm)
    x1p1 = 1/2 * (sxpM+sxpm)
    x2x2 = x1x1
    p2p2 = p1p1
    x2p2 = x1p1
    x1x2 = 1/2 * (sxxM-sxxm)
    p1p2 = 1/2 * (sppM-sppm)
    x1p2 = 1/2 * (sxpM-sxpm)
    x2p1 = x1p2

    if unpacked:
        return x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2, nu1, nu2
    else:
        cov_matrix = np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])
        return cov_matrix, nu1, nu2
#w0 = 1
#c = 0
#delta = 0.1
#g = .2
#ca1, cq1 = w0**2-2*g-g**2/4+c, -.5 * delta
#ca2, cq2 = w0**2-2*g-g**2/4-c, -.5 * delta
#nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
#
#i = 2
##A1, A2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
#t = np.linspace(0,20, 1000)
#V = w0**2 + delta*np.cos(2*t)
#T = 100
#T1, T2 = T, T
#wc = 50

w0 = 1
wd = 1.94
c1 = 0.03
c0 = 0.1
g = 0.0005
T = 200
r = 0
wM, wm = np.sqrt(w0**2+c0+c1), np.sqrt(w0**2-c0-c1)
ci1 = np.array([1/2*wM, 0, wM/2])
ci2 = np.array([1/2*wm, 0, wm/2])
#ci1 = np.array([1,0,np.exp(-2*r)])
#ci2 = np.array([1,0,np.exp(2*r)])

t = np.linspace(0, 450, 1000)
x1x1g, x2x2g, x1x2g, x1p1g, x2p2g, x1p2g, x2p1g, p1p1g, p2p2g, p1p2g, nu1g, nu2g = cov_galve(t, w0, g, c0, c1, T, ci1, ci2 ,wd, i = 2, unpacked = True)
Mcov = np.array([[x1x1g, x1p1g, x1x2g, x1p2g], [x1p1g, p1p1g, x2p1g, p1p2g], [x1x2g, x2p1g, x2x2g, x2p2g], [x1p2g, p1p2g, x2p2g, p2p2g]])

neg = aes.En(t, Mcov)
plt.clf()
plt.plot(t, neg)

#ca1, cq1 = (4/wd**2)*(w0**2-g**2/4+c0-2*g), -(2/wd**2)*c1
#ca2, cq2 = (4/wd**2)*(w0**2-g**2/4-c0-2*g), (2/wd**2)*c1
#nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
#A1, 0A2 = fl.mathieu_coefs(ca1, cq1, nu1), fl.mathieu_coefs(ca2, cq2, nu2)
#i = 2
#t = np.linspace(0,100, 2000)
#t_mio = np.linspace(0, 30, 100)
#wc = 50
#T = 50
#T1, T2 = T, T

#r = 0
#ci1 = np.array([1,0,np.exp(-2*r)])
#ci2 = np.array([1,0,np.exp(2*r)])
#ci1, ci2 = np.array([0, 0, 0]), np.array([0, 0, 0])
#x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2 = aes.cov(t_mio, g, ca1, cq1, ca2, cq2, T1, T2, wc, i, unpacked=True)
#Mcov_mia = np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])
#x1x1g, x2x2g, x1x2g, x1p1g, x2p2g, x1p2g, x2p1g, p1p1g, p2p2g, p1p2g, nu1g, nu2g = cov_galve(t, w0, g, c0, c1, T, ci1, ci2, unpacked = True, i = 2)
#Mcov = np.array([[x1x1g, x1p1g, x1x2g, x1p2g], [x1p1g, p1p1g, x2p1g, p1p2g], [x1x2g, x2p1g, x2x2g, x2p2g], [x1p2g, p1p2g, x2p2g, p2p2g]])
#Mcov_galve = np.array([[x1x1g, x1x2g, x1p1g, x1p2g], [x1x2g, x2x2g, x2p1g, x2p2g], [x1p1g, x2p1g, p1p1g, p1p2g], [x1p2g, x2p2g, p1p2g, p2p2g]]) 
#t = 2/wd * t
#neg = aes.En(t, Mcov)
#dis = aes.discordia(t, Mcov)
#neg_mia = aes.En(t_mio, Mcov_mia)

#plt.clf()
#plt.subplot(1, 2, 1)
#plt.plot(t, neg, 'b')
#plt.plot(t_mio, neg_mia, 'g')
#plt.title('$Discordia$')
#plt.xlabel('t')
#text = '$\gamma$ = '+str(g) + '\n $T$ = ' + str(T) + '\n $c_1$ = ' + str(c1) + '\n $r$ = ' + str(r)
#plt.text(0.2, 2, text)
#plt.plot(t_mio, neg_mia, 'g')
#sxx, sxp, spp = Mcov[0, 0, :], Mcov[0, 1, :], Mcov[1, 1, :]
#xx, pp = Mcov[0, 2, :], Mcov[1, 3, :]
#plt.subplot(1, 3, 2)
#plt.plot(t, xx, 'b', t, pp, 'r')
#plt.subplot(1, 2, 2)
#plt.plot(t, sxx, label = '$\sigma_{xx}$') 
#plt.plot(t, spp, label = '$\sigma_{pp}$')
#plt.plot(t, sxp, label = '$\sigma_{xp}$')
#plt.legend(bbox_to_anchor=(0.5, 1))
#plt.xlabel('t')



#def En_galve(Mcov):
#    """
#    Negatividad logaritmica
#    """
#    S = np.array([[0,0,1,0], [0,0,0,1], [-1,0,0,0], [0,-1,0,0]])
#    M = -1j*np.dot(S, Mcov)
#    autov = np.linalg.eig(M)[0]
##    print 2*np.abs(autov)
#    autov[2*np.abs(autov) > 1] = 1
##    print autov
#    return -1/2 * np.sum(np.log2(autov))    
#    
##plt.clf()
#neg_galve = np.array([])
#for i in np.arange(len(t)):
#    neg_galve = np.append(neg_galve, En_galve(Mcov_galve[:, :, i]))
#plt.plot(t, neg_galve, 'g')


#sxx, spp, sxp = sigma(t, g, T, ca1, cq1, ci, V, i)

#spp2 = si.deriv(x1p1, t)+g*sxp+(V-2*g)*sxx
#spp = Mcov[1,1,:]
#p1p1 = Mcov_mia[1,1,:]

#plt.clf()
#plt.plot(t, spp.real, '-g')
#plt.plot(t, p1p1, 'r')
#from scipy.optimize import curve_fit
#def f(t, a, b, c):
#    return c+a*np.exp(b*t)
#popt, pcov = curve_fit(f, t, spp.real)
#a,b,c = popt[0], popt[1], popt[2]
#plt.plot(t, f(t, .03, 1, 0), '-')
#plt.axis([0, 100, 0, 500])

#plt.plot(t, spp2, 'b')
#plt.plot(t, -.18*np.cos(2*t), 'ro-')
#plt.axis([0, 45, 0, 120])
#print 'Limite: ', np.abs(nu1g.imag/g)
#print 'Estable: ', nu1.imag - g/2 < 0
#print 'Asymtotic: ', 1/np.abs(nu1.imag-g/2)
print 'done'

#plt.plot(t, f(t, .00001, 1, 0), '-')
