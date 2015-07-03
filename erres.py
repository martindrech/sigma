# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 13:34:19 2015

@author: martin
"""

from __future__ import division
import numpy as np

def R1_bajaT(w, t, n, k, g, nu1, nu2):
    """
    R1_bajaT
    """
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    ret = (1j/8)*(((np.cos((2*n+u+1j*v)*t)-1j*np.sin((2*n+u+1j*v)*t))*(g-(2*1j)*(2*n+u+1j*v))*(np.sin((2*k+x+1j*y)*t)*(g-(2*1j)*n-1j*u+v)-np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*((-2*1j)*np.arctan((2*(2*n+u-w))/(g+2*v))+np.log((g+2*v)**2/4+(2*n+u-w)**2)))/(g**2+4*k**2-4*n**2-4*n*u-u**2-(4*1j)*n*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*n-1j*u+v)+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos((2*n+u+1j*v)*t)+1j*np.sin((2*n+u+1j*v)*t))*(g+(2*1j)*(2*n+u+1j*v))*(np.sin((2*k+x+1j*y)*t)*((-1j)*g+2*n+u+1j*v)+1j*np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*(2*np.arctan((2*(2*n+u+w))/(g-2*v))-1j*np.log((g-2*v)**2/4+(2*n+u+w)**2)))/(g**2+4*k**2-4*n**2-4*n*u-u**2+(2*1j)*g*(2*n+u+1j*v)-(4*1j)*n*v-(2*1j)*u*v+v**2+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos((2*k+x+1j*y)*t)+1j*np.sin((2*k+x+1j*y)*t))*(g+(2*1j)*(2*k+x+1j*y))*(1j*np.cos((2*n+u+1j*v)*t)*(2*n+u+1j*v)+np.sin((2*n+u+1j*v)*t)*((-1j)*g+2*k+x+1j*y))*(2*np.arctan((2*(2*k-w+x))/(g-2*y))-1j*np.log((2*k-w+x)**2+(g-2*y)**2/4)))/(g**2-4*k**2+4*n**2+4*n*u+u**2+(4*1j)*n*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)+(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2)+((np.cos((2*k+x+1j*y)*t)-1j*np.sin((2*k+x+1j*y)*t))*(g-(2*1j)*(2*k+x+1j*y))*(-(np.cos((2*n+u+1j*v)*t)*(2*n+u+1j*v))+np.sin((2*n+u+1j*v)*t)*(g-(2*1j)*k-1j*x+y))*((-2*1j)*np.arctan((2*(2*k+w+x))/(g+2*y))+np.log((2*k+w+x)**2+(g+2*y)**2/4)))/(g**2-4*k**2+4*n**2+4*n*u+u**2+(4*1j)*n*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*k-1j*x+y)))
    return ret
    
def R2_bajaT(w, t, n, m, k, g, nu1, nu2):
    """
    R2_bajaT
    """
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    ret = (1j/8)*(((np.cos(2*t*(n-m))+1j*np.sin(2*t*(n-m)))*(g-(2*1j)*(2*m+u+1j*v))*(np.sin((2*k+x+1j*y)*t)*(1j*g+2*m+u+1j*v)+np.cos((2*k+x+1j*y)*t)*((-2*1j)*k-1j*x+y))*(2*np.arctan((2*(2*m+u-w))/(g+2*v))+1j*np.log((g+2*v)**2/4+(2*m+u-w)**2)))/(g**2+4*k**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((1j*np.cos(2*t*(n-m))+np.sin(2*t*(n-m)))*(g+(2*1j)*(2*m+u+1j*v))*(np.sin((2*k+x+1j*y)*t)*(g+1j*(2*m+u+1j*v))-np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*(2*np.arctan((2*(2*m+u+w))/(g-2*v))-1j*np.log((g-2*v)**2/4+(2*m+u+w)**2)))/(g**2+4*k**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos((2*k+x+1j*y)*t)+1j*np.sin((2*k+x+1j*y)*t))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g+1j*(2*k+x+1j*y)))*((-1j)*g+4*k+2*x+(2*1j)*y)*(2*np.arctan((2*(2*k-w+x))/(g-2*y))-1j*np.log((2*k-w+x)**2+(g-2*y)**2/4)))/(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)+(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2)+((np.cos((2*k+x+1j*y)*t)-1j*np.sin((2*k+x+1j*y)*t))*(1j*g+4*k+2*x+(2*1j)*y)*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g-(2*1j)*k-1j*x+y))*(-2*np.arctan((2*(2*k+w+x))/(g+2*y))-1j*np.log((2*k+w+x)**2+(g+2*y)**2/4)))/(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*k-1j*x+y)))
    return ret


def R3_bajaT(w, t, n, m, k, l, g, nu1, nu2):
    """
    R3_bajaT
    """
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag

    re = (1j/8)*(((np.cos(2*t*(n-m))+1j*np.sin(2*t*(n-m)))*(1j*g+4*m+2*u+(2*1j)*v)*(np.sin(2*t*(k-l))*(g-(2*1j)*m-1j*u+v)+np.cos(2*t*(k-l))*(2*l+x+1j*y))*(2*np.arctan((2*(2*m+u-w))/(g+2*v))+1j*np.log((g+2*v)**2/4+(2*m+u-w)**2)))/(g**2+4*l**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*l*(x+1j*y)+(2*1j)*x*y-y**2)+((1j*np.cos(2*t*(n-m))+np.sin(2*t*(n-m)))*(g+(2*1j)*(2*m+u+1j*v))*(np.sin(2*t*(k-l))*(g+1j*(2*m+u+1j*v))+np.cos(2*t*(k-l))*(2*l+x+1j*y))*(2*np.arctan((2*(2*m+u+w))/(g-2*v))-1j*np.log((g-2*v)**2/4+(2*m+u+w)**2)))/(g**2+4*l**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*l*(x+1j*y)+(2*1j)*x*y-y**2)+((1j*np.cos(2*t*(k-l))+np.sin(2*t*(k-l)))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g+1j*(2*l+x+1j*y)))*(g+(2*1j)*(2*l+x+1j*y))*(2*np.arctan((2*(2*l-w+x))/(g-2*y))-1j*np.log((2*l-w+x)**2+(g-2*y)**2/4)))/(g**2-4*l**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*l*(x+1j*y)+(2*1j)*g*(2*l+x+1j*y)-(2*1j)*x*y+y**2)+(((-1j)*np.cos(2*t*(k-l))+np.sin(2*t*(k-l)))*(g-(2*1j)*(2*l+x+1j*y))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g-(2*1j)*l-1j*x+y))*(-2*np.arctan((2*(2*l+w+x))/(g+2*y))-1j*np.log((2*l+w+x)**2+(g+2*y)**2/4)))/(g**2-4*l**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*l*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*l-1j*x+y)))
    return re 
    
def R1_altaT(w, t, n, k, g, nu1, nu2):
    """
    R1_altaT
    """
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    ret = (((np.cos((2*n+u+1j*v)*t)-1j*np.sin((2*n+u+1j*v)*t))*(np.sin((2*k+x+1j*y)*t)*(g-(2*1j)*n-1j*u+v)-np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*((-2*1j)*np.arctan((2*(2*n+u-w))/(g+2*v))+np.log((g+2*v)**2/4+(2*n+u-w)**2)))/(g**2+4*k**2-4*n**2-4*n*u-u**2-(4*1j)*n*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*n-1j*u+v)+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos((2*n+u+1j*v)*t)+1j*np.sin((2*n+u+1j*v)*t))*(np.sin((2*k+x+1j*y)*t)*((-1j)*g+2*n+u+1j*v)+1j*np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*(2*np.arctan((2*(2*n+u+w))/(g-2*v))-1j*np.log((g-2*v)**2/4+(2*n+u+w)**2)))/(g**2+4*k**2-4*n**2-4*n*u-u**2+(2*1j)*g*(2*n+u+1j*v)-(4*1j)*n*v-(2*1j)*u*v+v**2+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+(((-1j)*np.cos((2*k+x+1j*y)*t)+np.sin((2*k+x+1j*y)*t))*(np.cos((2*n+u+1j*v)*t)*(2*n+u+1j*v)-np.sin((2*n+u+1j*v)*t)*(g+1j*(2*k+x+1j*y)))*(2*np.arctan((2*(2*k-w+x))/(g-2*y))-1j*np.log((2*k-w+x)**2+(g-2*y)**2/4)))/(g**2-4*k**2+4*n**2+4*n*u+u**2+(4*1j)*n*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)+(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2)+((np.cos((2*k+x+1j*y)*t)-1j*np.sin((2*k+x+1j*y)*t))*(np.cos((2*n+u+1j*v)*t)*((-2*1j)*n-1j*u+v)+np.sin((2*n+u+1j*v)*t)*(1j*g+2*k+x+1j*y))*(2*np.arctan((2*(2*k+w+x))/(g+2*y))+1j*np.log((2*k+w+x)**2+(g+2*y)**2/4)))/(g**2-4*k**2+4*n**2+4*n*u+u**2+(4*1j)*n*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*k-1j*x+y)))/4
    return ret
    
def R2_altaT(w, t, n, m, k, g, nu1, nu2):
    """
    R2_altaT
    """
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    ret = (((np.cos(2*t*(n-m))+1j*np.sin(2*t*(n-m)))*(np.sin((2*k+x+1j*y)*t)*(1j*g+2*m+u+1j*v)+np.cos((2*k+x+1j*y)*t)*((-2*1j)*k-1j*x+y))*(2*np.arctan((2*(2*m+u-w))/(g+2*v))+1j*np.log((g+2*v)**2/4+(2*m+u-w)**2)))/(g**2+4*k**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((1j*np.cos(2*t*(n-m))+np.sin(2*t*(n-m)))*(np.sin((2*k+x+1j*y)*t)*(g+1j*(2*m+u+1j*v))-np.cos((2*k+x+1j*y)*t)*(2*k+x+1j*y))*(2*np.arctan((2*(2*m+u+w))/(g-2*v))-1j*np.log((g-2*v)**2/4+(2*m+u+w)**2)))/(g**2+4*k**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*k*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos((2*k+x+1j*y)*t)+1j*np.sin((2*k+x+1j*y)*t))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g+1j*(2*k+x+1j*y)))*((2*1j)*np.arctan((2*(2*k-w+x))/(g-2*y))+np.log((2*k-w+x)**2+(g-2*y)**2/4)))/(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)+(2*1j)*g*(2*k+x+1j*y)-(2*1j)*x*y+y**2)+((1j*np.cos((2*k+x+1j*y)*t)+np.sin((2*k+x+1j*y)*t))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g-(2*1j)*k-1j*x+y))*(2*np.arctan((2*(2*k+w+x))/(g+2*y))+1j*np.log((2*k+w+x)**2+(g+2*y)**2/4)))/(g**2-4*k**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*k*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*k-1j*x+y)))/4
    return ret


def R3_altaT(w, t, n, m, k, l, g, nu1, nu2):
    """
    R3_altaT
    """
    u, v = nu1.real, nu1.imag
    x, y = nu2.real, nu2.imag
    re = ((((-1j)*np.cos(2*t*(n-m))+np.sin(2*t*(n-m)))*(np.sin(2*t*(k-l))*(g-(2*1j)*m-1j*u+v)+np.cos(2*t*(k-l))*(2*l+x+1j*y))*(-2*np.arctan((2*(2*m+u-w))/(g+2*v))-1j*np.log((g+2*v)**2/4+(2*m+u-w)**2)))/(g**2+4*l**2-4*m**2-4*m*u-u**2-(4*1j)*m*v-(2*1j)*u*v+v**2+2*g*((-2*1j)*m-1j*u+v)+x**2+4*l*(x+1j*y)+(2*1j)*x*y-y**2)+((1j*np.cos(2*t*(n-m))+np.sin(2*t*(n-m)))*(np.sin(2*t*(k-l))*(g+1j*(2*m+u+1j*v))+np.cos(2*t*(k-l))*(2*l+x+1j*y))*(2*np.arctan((2*(2*m+u+w))/(g-2*v))-1j*np.log((g-2*v)**2/4+(2*m+u+w)**2)))/(g**2+4*l**2-4*m**2-4*m*u-u**2+(2*1j)*g*(2*m+u+1j*v)-(4*1j)*m*v-(2*1j)*u*v+v**2+x**2+4*l*(x+1j*y)+(2*1j)*x*y-y**2)+((np.cos(2*t*(k-l))-1j*np.sin(2*t*(k-l)))*(np.cos(2*t*(n-m))*((-2*1j)*m-1j*u+v)+np.sin(2*t*(n-m))*((-1j)*g+2*l+x+1j*y))*(2*np.arctan((2*(2*l-w+x))/(g-2*y))-1j*np.log((2*l-w+x)**2+(g-2*y)**2/4)))/(g**2-4*l**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*l*(x+1j*y)+(2*1j)*g*(2*l+x+1j*y)-(2*1j)*x*y+y**2)+(((-1j)*np.cos(2*t*(k-l))+np.sin(2*t*(k-l)))*(np.cos(2*t*(n-m))*(2*m+u+1j*v)+np.sin(2*t*(n-m))*(g-(2*1j)*l-1j*x+y))*(2*np.arctan((2*(2*l+w+x))/(g+2*y))+1j*np.log((2*l+w+x)**2+(g+2*y)**2/4)))/(g**2-4*l**2+4*m**2+4*m*u+u**2+(4*1j)*m*v+(2*1j)*u*v-v**2-x**2-4*l*(x+1j*y)-(2*1j)*x*y+y**2+2*g*((-2*1j)*l-1j*x+y)))/4
    return re 
    


