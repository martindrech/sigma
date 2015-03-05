# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 18:30:23 2014

@author: martindrech
"""


from __future__ import division
import numpy as np
from scipy.integrate import quad


'''
y''+(a-2q cos(2t))y = 0
nu es el exp de floquet, c el vector con la expansion de fourier de p(t)
'''
def p(t, c):
    """
    Devuelve la funcion periodica P(t)
    """
    p = np.zeros(np.size(t)) + 0j
    c = c + 0j
    n = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    if np.size(t) > 1: 
        for i, tau in enumerate(t):
            p[i] = np.dot(c, np.exp(2*tau*n*1j))
    else: 
        p = np.dot(c, np.exp(2*t*n*1j))
    return p

def B(c, nu, t):
    """
    Devuelve el coef B que aparece en las dispersiones
    """
    P = p(t, c)
    dp = (p(t[1], c)-p(t[0], c))/(t[1]-t[0])
    
    return 1j/( P[0]*(dp+nu*P[0]*1j))
    
def y1(t, c, nu):
    """
    Solucion de eq. de mathieu sin ci
    """
    return np.exp(1j*nu*t)*p(t, c)

def y2(t, c, nu):
    """
    Solucion de eq. de mathieu sin ci
    """
    return np.exp(-1j*nu*t)*p(-t, c)
        
def phi1(t, c, nu):
    """
    Solucion de eq. de mathieu con ci = [0,1]
    """
    dp = (p(t[1], c)-p(t[0], c))/(t[1]-t[0])
    cte = 2*(dp+1j*nu*p(t[0], c))
    
    return np.real((y1(t,c,nu)-y2(t,c,nu))/cte)
    
def phi2(t, c, nu):
    """
    Solucion de eq. de mathieu con ci = [1,0]
    """
    cte = 2*p(t[0], c)
    return np.real((y1(t,c,nu)+y2(t,c,nu))/cte)
    
### funciones auxiliares###
    

def complex_int(f, a, b, *args):
   """
   Integra numericamente funcion f entre a y b, siendo f compleja
   """
   
   def real_f(x, *args):
       return np.real(f(x, *args))       
   def imag_f(x, *args):
       return np.imag(f(x, *args))    
     
   real_int = quad(real_f, a, b, *args)[0]
   imag_int = quad(imag_f, a, b, *args)[0]
   
   return real_int + imag_int * 1j
   
def R(w, t, n, m, k, l, g, nu):
    """
    El R_nmkl que hay que integrar en frecuencia
    """
    bm, bl = 2*m+nu, 2*l+nu
    anm, bnm, akl, bkl = np.cos(2*t*(n-m)), np.sin(2*t*(n-m)), np.cos(2*t*(k-l)), np.sin(2*t*(k-l)) 
    g1, g2  = g/2+w*1j, g/2-w*1j
    
    return ((bm*anm+g1*bnm)/(g1**2+bm**2))*((bl*akl+g2*bkl)/(g2**2+bl**2))
    
    
def int_baja_num(wc, t, n, m, k, l, g, nu):
    """
    Integral numerica de alpha(w)*R_nmkl a bajas temperaturas
    """
    def integrando(w, t, n, m, k, l, g, nu):
        return w*R(w, t, n, m, k, l, g, nu)
        
    inte = np.array([])
    for ti in t:
        valor = complex_int(integrando, 0, wc, (ti, n, m, k, l, g, nu))
        inte = np.append(inte, [valor])
        print ti
    return inte
    
    

def int_baja(w, t, n, m, k, l, g, nu):
    """
    Integral de alpha(w)*R_nmkl a bajas temperaturas, la que me da wolfram por ahora
    """
    ret = ((-2*1j)*(np.cos(2*t*(k-l))*(2*l+nu)*(2*np.cos(2*t*(n-m))*(2*m+nu)*((2*l+nu)**2-(2*m+nu)**2)-np.sin(2*t*(n-m))*g*((2*l+nu)**2+3*(2*m+nu)**2+g**2))+np.sin(2*t*(k-l))*(np.cos(2*t*(n-m))*(2*m+nu)*g*(3*(2*l+nu)**2+(2*m+nu)**2+g**2)+np.sin(2*t*(n-m))*(2*(2*l+nu)**4-2*(2*l+nu)**2*(2*m+nu)**2+3*(2*l+nu)**2*g**2+(2*m+nu)**2*g**2+g**4)))*np.arctan((4*g*w)/(4*(2*l+nu)**2+g**2-4*w**2))+(2*1j)*(np.cos(2*t*(k-l))*(2*l+nu)*(np.cos(2*t*(n-m))*(-2*(2*l+nu)**2*(2*m+nu)+2*(2*m+nu)**3)+np.sin(2*t*(n-m))*g*((2*l+nu)**2+3*(2*m+nu)**2+g**2))+np.sin(2*t*(k-l))*(-(np.cos(2*t*(n-m))*(2*m+nu)*g*(3*(2*l+nu)**2+(2*m+nu)**2+g**2))+np.sin(2*t*(n-m))*(2*(2*m+nu)**4+3*(2*m+nu)**2*g**2+g**4+(2*l+nu)**2*(-2*(2*m+nu)**2+g**2))))*np.arctan((4*g*w)/(4*(2*m+nu)**2+g**2-4*w**2))+(4*1j)*(np.cos(2*t*(k-l))*(2*l+nu)*(2*np.sin(2*t*(n-m))*(2*m+nu)*((2*l+nu)**2-(2*m+nu)**2)+np.cos(2*t*(n-m))*g*((2*l+nu)**2+3*(2*m+nu)**2+g**2))+np.sin(2*t*(k-l))*(np.sin(2*t*(n-m))*(2*m+nu)*g*(3*(2*l+nu)**2+(2*m+nu)**2+g**2)+np.cos(2*t*(n-m))*(2*(2*m+nu)**4+3*(2*m+nu)**2*g**2+g**4+(2*l+nu)**2*(-2*(2*m+nu)**2+g**2))))*np.arctanh(((-1j)*g+2*w)/(2*(2*m+nu)))-(4*1j)*(np.sin(2*t*(k-l))*(2*l+nu)*(np.cos(2*t*(n-m))*(-2*(2*l+nu)**2*(2*m+nu)+2*(2*m+nu)**3)+np.sin(2*t*(n-m))*g*((2*l+nu)**2+3*(2*m+nu)**2+g**2))+np.cos(2*t*(k-l))*(np.cos(2*t*(n-m))*(2*m+nu)*g*(3*(2*l+nu)**2+(2*m+nu)**2+g**2)+np.sin(2*t*(n-m))*(2*(2*l+nu)**4-2*(2*l+nu)**2*(2*m+nu)**2+3*(2*l+nu)**2*g**2+(2*m+nu)**2*g**2+g**4)))*np.arctanh((1j*g+2*w)/(2*(2*l+nu)))+(np.cos(2*t*(k-l))*(2*l+nu)*(2*np.cos(2*t*(n-m))*(2*m+nu)*((2*l+nu)**2-(2*m+nu)**2)-np.sin(2*t*(n-m))*g*((2*l+nu)**2+3*(2*m+nu)**2+g**2))+np.sin(2*t*(k-l))*(np.cos(2*t*(n-m))*(2*m+nu)*g*(3*(2*l+nu)**2+(2*m+nu)**2+g**2)+np.sin(2*t*(n-m))*(2*(2*l+nu)**4-2*(2*l+nu)**2*(2*m+nu)**2+3*(2*l+nu)**2*g**2+(2*m+nu)**2*g**2+g**4)))*np.log(16*(2*l+nu)**4+8*(2*l+nu)**2*(g**2-4*w**2)+(g**2+4*w**2)**2)+(np.cos(2*t*(k-l))*(2*l+nu)*(np.cos(2*t*(n-m))*(-2*(2*l+nu)**2*(2*m+nu)+2*(2*m+nu)**3)+np.sin(2*t*(n-m))*g*((2*l+nu)**2+3*(2*m+nu)**2+g**2))+np.sin(2*t*(k-l))*(-(np.cos(2*t*(n-m))*(2*m+nu)*g*(3*(2*l+nu)**2+(2*m+nu)**2+g**2))+np.sin(2*t*(n-m))*(2*(2*m+nu)**4+3*(2*m+nu)**2*g**2+g**4+(2*l+nu)**2*(-2*(2*m+nu)**2+g**2))))*np.log(16*(2*m+nu)**4+8*(2*m+nu)**2*(g**2-4*w**2)+(g**2+4*w**2)**2))/(8*((2*l+nu)**4-2*(2*l+nu)**2*((2*m+nu)**2-g**2)+((2*m+nu)**2+g**2)**2))    
    return ret
    
def S_pp_baja_int(w, t, n, m, g, nu):
    """
    La integral que va en S_pp_baja
    """
    ret = (4*np.sin(2*t*(n-m))*(2*m+nu)*np.arctan((2*((2*m+nu)-w))/g)+2*np.cos(2*t*(n-m))*g*np.arctan((2*(-(2*m+nu)+w))/g)+4*np.sin(2*t*(n-m))*(2*m+nu)*np.arctan((2*((2*m+nu)+w))/g)-2*np.cos(2*t*(n-m))*g*np.arctan((2*((2*m+nu)+w))/g)-2*np.cos(2*t*(n-m))*(2*m+nu)*np.log(4*(2*m+nu)**2+g**2-8*(2*m+nu)*w+4*w**2)-np.sin(2*t*(n-m))*g*np.log(4*(2*m+nu)**2+g**2-8*(2*m+nu)*w+4*w**2)-2*np.cos(2*t*(n-m))*(2*m+nu)*np.log(4*(2*m+nu)**2+g**2+8*(2*m+nu)*w+4*w**2)-np.sin(2*t*(n-m))*g*np.log(4*(2*m+nu)**2+g**2+8*(2*m+nu)*w+4*w**2))/2
    return ret



def S_pp_baja(t,c,nu, g, wc):
    """
    El termino que hace falta para tener sigma_pp, para bajas temperaturas
    """
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    
    suma = np.array(np.zeros(len(t)))
    for n_i, n in enumerate(enes):
        for m_i, m in enumerate(enes):
            r = c[n_i]*c[m_i]*(S_pp_baja_int(wc, t, n, m, g, nu)-S_pp_baja_int(0, t, n, m, g, nu))
            suma = suma + r
    return suma* (B(c,nu,t) * g) / (2*np.pi)
   
    

    
def sigma_baja(t, c, nu, g, wc):
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    
    suma = np.array(np.zeros(len(t)))
    for n_i, n in enumerate(enes):
        for m_i, m in enumerate(enes):
            for k_i, k in enumerate(enes):
                for l_i, l in enumerate(enes):
                    print n, m, k, l
                    r = c[n_i]*c[m_i]*np.conj(c[k_i]*c[l_i])*(int_baja(wc, t, n, m, k, l, g, nu)-int_baja(0, t, n, m, k, l, g, nu))
                    suma = suma + r
    return suma* (np.abs(B(c,nu,t))**2 *g / np.pi)

def sigma_baja_num(t, c, nu, g, wc):
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    
    suma = np.array(np.zeros(len(t)))
    for n_i, n in enumerate(enes):
        for m_i, m in enumerate(enes):
            for k_i, k in enumerate(enes):
                for l_i, l in enumerate(enes):
                    print n, m, k, l
                    r = c[n_i]*c[m_i]*np.conj(c[k_i]*c[l_i])*int_baja_num(wc, t, n, m, k, l, g, nu)
                    suma = suma + r
    return suma* (np.abs(B(c,nu,t))**2 *g / np.pi)
    
    
    
def sigma_alta(t, c, nu, g, wc):
    enes = np.linspace(-(len(c)-1)/2, (len(c)-1)/2, len(c))
    
    suma = np.array(np.zeros(len(t)))
    for n_i, n in enumerate(enes):
        for m_i, m in enumerate(enes):
            for k_i, k in enumerate(enes):
                for l_i, l in enumerate(enes):
                    print n, m, k, l
                    r = c[n_i]*c[m_i]*c[k_i]*c[l_i]*int_alta(t, m, n, k, l, g, nu)
                    suma = suma + r
    return suma*np.abs(B(c,nu,t))**2
    
    
    


