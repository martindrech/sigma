# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 19:42:48 2015

@author: martin
"""

from __future__ import division
import numpy as np
from scipy.integrate import quad
import pylab as plt
#import sigma as si


def complex_int(f, a, b, *args):
   
   def real_f(x, *args):
     return np.real(f(x, *args))       
   def imag_f(x, *args):
     return np.imag(f(x, *args))    
     
   real_int = quad(real_f, a, b, *args)[0]
   imag_int = quad(imag_f, a, b, *args)[0]
   
   return real_int + imag_int * 1j
   

   
def a(n, m, t):
    return np.cos(2*t*(n-m))
def b(n, m, t):
    return np.sin(2*t*(n-m))
    
def r(w, n, m, t, g, nu):
    return ( (2*m+nu)*a(n,m,t)+(g/2+1j*w)*b(n,m,t) ) / ((g/2+w*1j)**2+(2*m+nu)**2)

def R(w, n, m, k, l, t, g, nu):
    return r(w, n, m, t, g, nu)*r(-w, k, l, t, g, nu)
    
def integrando(w, t):
    return np.sin(w*t)


def ArcCoth(x):
    
    return 1/2 * ( np.log(1+1/x)- np.log(1-1/x) )

#n, m, k, l = 1,0,1,0
#g = 1
#nu = 1.2
wc = 5
tiempo = np.linspace(1, 10, 100)

def una_int(tiempo, wc):
    inte = np.array([])
    for t in tiempo:
        print t
        valor = complex_int(integrando, 0, wc, (t))
        inte = np.append(inte, [valor])
    return inte

i = una_int(tiempo, wc)
plt.clf()
plt.plot(tiempo, i, 'b', tiempo, (1-np.cos(wc*tiempo))/tiempo, 'o')