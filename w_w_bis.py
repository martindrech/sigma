# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 12:50:58 2015

@author: martin
"""
from __future__ import division
import numpy as np
import sigma as si
import erres as R
from scipy.special import expi
def impedir_peq(arr, eps):
    mascara = np.abs(arr) < eps
    arr[mascara] = eps*np.sign(arr[mascara])

def termino2(wc, t, n, m, k, g, nu1, nu2):
    
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    ret = (1j/4)*(2*k+x+1j*y)*(((np.cos(2*t*(n-m))+1j*np.sin(2*t*(n-m)))*(g-(2*1j)*(2*m+u+1j*v))*(expi((t*(g-(2*1j)*(2*m+u+1j*v)))/2)-expi((t*(g-(2*1j)*(2*m+u+1j*v-wc)))/2))*np.exp(-(t*(g-(2*1j)*(2*m+u+1j*v)))/2))/(g**2+4*k**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)-((np.cos(2*t*(n-m))-1j*np.sin(2*t*(n-m)))*(g+(2*1j)*(2*m+u+1j*v))*(expi((t*(g+(2*1j)*(2*m+u+1j*v)))/2)-expi((t*(g+(2*1j)*(2*m+u+1j*v+wc)))/2))*np.exp(-(t*(g+(2*1j)*(2*m+u+1j*v)))/2))/(g**2+4*k**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((g-(2*1j)*(2*k+x+1j*y))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g-(2*1j)*k-1j*x+y))*(expi(-(t*(g-(2*1j)*(2*k+x+1j*y)))/2)-expi(-(t*(g-(2*1j)*(2*k+wc+x+1j*y)))/2))*np.exp((t*(g-(2*1j)*(2*k+x+1j*y)))/2))/((2*k+x+1j*y)*(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)-(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2))-((np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g+1j*(2*k+x+1j*y)))*(g+(2*1j)*(2*k+x+1j*y))*(expi(-(t*(g+(2*1j)*(2*k+x+1j*y)))/2)-expi(-(t*(g+(2*1j)*(2*k-wc+x+1j*y)))/2))*np.exp((t*(g+(2*1j)*(2*k+x+1j*y)))/2))/((2*k+x+1j*y)*(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)+(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2)))
    return ret
def termino3(wc, t, n, m, k, g, nu1, nu2):
    
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    ret = (1j/4)*(((np.cos((2*n+u+1j*v)*t)+1j*np.sin((2*n+u+1j*v)*t))*(g-(2*1j)*(2*m+u+1j*v))*(np.sin((2*k+x+1j*y)*t)*(g-(2*1j)*m-1j*u+v)-np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*(expi(-(t*(g-(2*1j)*(2*m+u+1j*v)))/2)-expi(-(t*(g-(2*1j)*(2*m+u+1j*v-wc)))/2))*np.exp((t*(g-(2*1j)*(2*m+u+1j*v)))/2))/(g**2+4*k**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)-((np.cos((2*n+u+1j*v)*t)-1j*np.sin((2*n+u+1j*v)*t))*(g+(2*1j)*(2*m+u+1j*v))*(np.sin((2*k+x+1j*y)*t)*(g+1j*(2*m+u+1j*v))-np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*(expi(-(t*(g+(2*1j)*(2*m+u+1j*v)))/2)-expi(-(t*(g+(2*1j)*(2*m+u+1j*v+wc)))/2))*np.exp((t*(g+(2*1j)*(2*m+u+1j*v)))/2))/(g**2+4*k**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)-((np.cos((2*k+x+1j*y)*t)-1j*np.sin((2*k+x+1j*y)*t))*(g-(2*1j)*(2*k+x+1j*y))*(np.cos((2*n+u+1j*v)*t)*(2*m+u+1j*v)+np.sin((2*n+u+1j*v)*t)*(g-(2*1j)*k-1j*x+y))*(expi((t*(g-(2*1j)*(2*k+x+1j*y)))/2)-expi((t*(g-(2*1j)*(2*k+wc+x+1j*y)))/2))*np.exp(-(t*(g-(2*1j)*(2*k+x+1j*y)))/2))/(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*k-1j*x+y))+((np.cos((2*k+x+1j*y)*t)+1j*np.sin((2*k+x+1j*y)*t))*(np.cos((2*n+u+1j*v)*t)*(2*m+u+1j*v)+np.sin((2*n+u+1j*v)*t)*(g+1j*(2*k+x+1j*y)))*(g+(2*1j)*(2*k+x+1j*y))*(expi((t*(g+(2*1j)*(2*k+x+1j*y)))/2)-expi((t*(g+(2*1j)*(2*k-wc+x+1j*y)))/2))*np.exp(-(t*(g+(2*1j)*(2*k+x+1j*y)))/2))/(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)+(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2))
    return ret
def termino4(wc, t, n, m, k, g, nu1, nu2):
    
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    ret = (1j/8)*(2*k+x+1j*y)*((((-1j)*np.cos((2*n+u+1j*v)*t)+np.sin((2*n+u+1j*v)*t))*(g-(2*1j)*(2*m+u+1j*v))*(2*np.arctan((2*(2*m+u))/(g+2*v))-2*np.arctan((2*(2*m+u-wc))/(g+2*v))+1j*(np.log((4*(2*m+u)**2+(g+2*v)**2)/4)-np.log((g+2*v)**2/4+(2*m+u-wc)**2))))/(g**2+4*k**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos((2*n+u+1j*v)*t)-1j*np.sin((2*n+u+1j*v)*t))*((-1j)*g+4*m+2*u+(2*1j)*v)*(2*np.arctan((2*(2*m+u))/(g-2*v))-2*np.arctan((2*(2*m+u+wc))/(g-2*v))-1j*(np.log((4*(2*m+u)**2+(g-2*v)**2)/4)-np.log((g-2*v)**2/4+(2*m+u+wc)**2))))/(g**2+4*k**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)-((np.cos((2*n+u+1j*v)*t)*(2*m+u+1j*v)+np.sin((2*n+u+1j*v)*t)*(g+1j*(2*k+x+1j*y)))*(g+(2*1j)*(2*k+x+1j*y))*((2*1j)*np.arctan((2*(2*k+x))/(g-2*y))-(2*1j)*np.arctan((2*(2*k-wc+x))/(g-2*y))-np.log((2*k-wc+x)**2+(g-2*y)**2/4)+np.log((4*(2*k+x)**2+(g-2*y)**2)/4)))/((2*k+x+1j*y)*(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)+(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2))-((1j*np.cos((2*n+u+1j*v)*t)*(2*m+u+1j*v)+np.sin((2*n+u+1j*v)*t)*(1j*g+2*k+x+1j*y))*(1j*g+4*k+2*x+(2*1j)*y)*((-2*1j)*np.arctan((2*(2*k+x))/(g+2*y))+(2*1j)*np.arctan((2*(2*k+wc+x))/(g+2*y))-np.log((2*k+wc+x)**2+(g+2*y)**2/4)+np.log((4*(2*k+x)**2+(g+2*y)**2)/4)))/((2*k+x+1j*y)*(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)-(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2)))
    return ret
    
    
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
        r =2*temp*A*(R.R1_altaT(wc, T, N, K, g, nu1, nu2, temp)-R.R1_altaT(0, T, N, K, g, nu1, nu2, temp))
    
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
        integrales = ( np.exp(g*t)*(R.R2_bajaT(wc,T,N,M,K,g,nu1,nu2)-R.R2_bajaT(0,T,N,M,K,g,nu1,nu2)) ) +np.exp(-g*t/2) * ( termino2(wc,T,N,M,K,g,nu1,nu2)-termino3(wc,T,N,M,K,g,nu1,nu2) ) - termino4(wc,T,N,M,K,g,nu1,nu2)
        r =A*integrales
    else:
        r =2*temp*A*np.exp(g*t)*(R.R2_altaT(wc,T,N,M,K,g,nu1,nu2, temp)-R.R2_altaT(0,T,N,M,K,g,nu1,nu2, temp))
    
    
    r = np.sum(r, axis = (0,1,2))

    phies = phi1*phim1
     
#    impedir_peq(phies, 0.1)
    return 1/phies * np.exp(-g*t/2)*(g/np.pi)  * r.real

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
        r =2*temp*A*(R.R3_altaT(wc,T,N,M,K,L,g,nu1,nu2, temp)-R.R3_altaT(0,T,N,M,K,L,g,nu1,nu2, temp))
    
    r = np.sum(r, axis = (0,1,2,3))
    phies = phi1*phim1
    
    return 1/phies * np.exp(g*t)*(g/np.pi) *  r.real
    
