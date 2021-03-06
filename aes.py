
"""
Created on Mon Jun  8 18:42:33 2015

@author: martin
"""

from __future__ import division
import numpy as np
import w_w as w_w
import floquet as fl
from sigma import deriv
import pylab as plt
from fourier import plot_ft as fft

def impedir_peq(arr, eps):
    mascara = np.abs(arr) < eps
    arr[mascara] = eps*np.sign(arr[mascara])

def a1(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc, phi1, phim1):
    """
    Devuelve la matrix a1(t)    
    """

    
    w2w2t1 = w_w.w2_w2(t, g, temp1, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w2mw2mt1 = w_w.w2_w2(t, g, temp1, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w2w2mt1 = w_w.w2_w2(t, g, temp1, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    w2w2t2 = w_w.w2_w2(t, g, temp2, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w2mw2mt2 =w_w.w2_w2(t, g, temp2, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w2w2mt2 =w_w.w2_w2(t, g, temp2, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    a11 = w2w2t1+w2mw2mt1+2*w2w2mt1 + w2w2t2+w2mw2mt2-2*w2w2mt2
    a12 = w2w2t1-w2mw2mt1 + w2w2t2-w2mw2mt2
    a21 = a12
    a22 = w2w2t1+w2mw2mt1-2*w2w2mt1 + w2w2t2+w2mw2mt2+2*w2w2mt2

    
    return 1/4 * np.array([[a11, a12], [a21, a22]])

def a2(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc, phi1, phim1):
    """
    Devuelve la matrix a2(t)    
    """
    w1w2t1 = w_w.w1_w2(t, g, temp1, nu1, c1, nu1, c1, wc, phi1, phi1)
    w1mw2mt1 = w_w.w1_w2(t, g, temp1, nu2, c2, nu2, c2, wc, phim1, phim1)   
    w1mw2t1 = w_w.w1_w2(t, g, temp1, nu2, c2, nu1, c1, wc, phim1, phi1)
    w1w2mt1 = w_w.w1_w2(t, g, temp1, nu1, c1, nu2, c2, wc, phi1, phim1)
    
    w1w2t2 = w_w.w1_w2(t, g, temp2, nu1, c1, nu1, c1, wc, phi1, phi1)
    w1mw2mt2 = w_w.w1_w2(t, g, temp2, nu2, c2, nu2, c2, wc, phim1, phim1)   
    w1mw2t2 = w_w.w1_w2(t, g, temp2, nu2, c2, nu1, c1, wc, phim1, phi1)
    w1w2mt2 = w_w.w1_w2(t, g, temp2, nu1, c1, nu2, c2, wc, phi1, phim1) 
    
    a11 = w1w2t1+w1w2mt1+w1mw2t1+w1mw2mt1 + w1w2t2-w1w2mt2-w1mw2t2+w1mw2mt2
    a12 = w1w2t1+w1mw2t1-w1w2mt1-w1mw2mt1 + w1w2t2-w1mw2t2+w1w2mt2-w1mw2mt2
    a21 = w1w2t1-w1mw2t1+w1w2mt1-w1mw2mt1 + w1w2t2+w1mw2t2-w1w2mt2-w1mw2mt2
    a22 = w1w2t1-w1w2mt1-w1mw2t1+w1mw2mt1 + w1w2t2+w1w2mt2+w1mw2t2+w1mw2mt2
    
    return  .5 * np.array([[a11, a12], [a21, a22]])

def a3(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc, phi1, phim1):
    """
    Devuelve la matrix a3(t)    
    """
    w1w1t1 = w_w.w1_w1(t, g, temp1, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w1mw1mt1 = w_w.w1_w1(t, g, temp1, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w1w1mt1 = w_w.w1_w1(t, g, temp1, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    w1w1t2 = w_w.w1_w1(t, g, temp2, nu1, c1, nu1 , c1, wc, phi1, phi1)
    w1mw1mt2 =w_w.w1_w1(t, g, temp2, nu2, c2, nu2 , c2, wc, phim1, phim1)
    w1w1mt2 =w_w.w1_w1(t, g, temp2, nu1, c1, nu2 , c2, wc, phi1, phim1)
    
    a11 = w1w1t1+w1mw1mt1+2*w1w1mt1+w1w1t2+w1mw1mt2-2*w1w1mt2
    a12 = w1w1t1-w1mw1mt1+w1w1t2-w1mw1mt2
    a21 = a12
    a22 = w1w1t2+w1mw1mt2+2*w1w1mt2+w1w1t1+w1mw1mt1-2*w1w1mt1
    
    
    return 1/4 * np.array([[a11, a12], [a21, a22]])


def b(t, g, phi1, dphi1, phi2, dphi2, phim1, dphim1, phim2, dphim2):
    """
    Devuelve las matrix b1, b2, b3, b4
    """

    f1 = phi1 * np.exp(-g*t/2)
    df1 = (dphi1- (g/2) * phi1) * np.exp(-g*t/2)
    g1 = phim1 * np.exp(-g*t/2)
    dg1 = (dphim1- (g/2) * phim1) * np.exp(-g*t/2)
    
    b1d = 1/2 * (df1/f1 + dg1/g1)
    b1c = 1/2 * (df1/f1 - dg1/g1)

    b2d = -1/2 * np.exp(-g*t) * (1/f1+1/g1)
    b2c =-1/2 * np.exp(-g*t) * (1/f1-1/g1)

    b3d = -1/2  * (1/f1+1/g1)
    b3c = -1/2  * (1/f1-1/g1)

    
    return b1d, b1c, b2d, b2c, b3d, b3c
    
def cov(t, g, ca1, cq1, ca2, cq2, temp1, temp2, wc = 50, i = 5, unpacked = False):
    """
    Devuelve todos los valores medios para la matriz de covarianza. Puede devolverlos como
    la matriz armada o unpacked.
    """    
    nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
    c1, c2 = fl.mathieu_coefs(ca1, cq1, nu1), fl.mathieu_coefs(ca2, cq2, nu2)
    c1, c2 = c1[c1.size//2-i:c1.size//2+i+1], c2[c2.size//2-i:c2.size//2+i+1]
 
    phi1, dphi1, phi2, dphi2 = fl.mathieu(ca1, cq1, t)    
    phim1, dphim1, phim2, dphim2 = fl.mathieu(ca2, cq2, t)
    
#    impedir_peq(phi1, .01)
#    impedir_peq(phim1, .01)

    
    Ma1 = a1(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc, phi1, phim1)
    Ma2 = a2(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc, phi1, phim1)
    Ma3 = a3(t, g,  nu1, c1, temp1, nu2, c2, temp2, wc, phi1, phim1)
    b1d, b1c, b2d, b2c, b3d, b3c = b(t, g, phi1, dphi1, phi2, dphi2, phim1, dphim1, phim2, dphim2)

    x1x1=(Ma3[1][1]*b3c**2-(Ma3[0][1]+Ma3[1][0])*b3c*b3d+Ma3[0][0]*b3d**2)/(b3c**2-b3d**2)**2
    x2x2=(Ma3[0][0]*b3c**2-(Ma3[0][1]+Ma3[1][0])*b3c*b3d+Ma3[1][1]*b3d**2)/(b3c**2-b3d**2)**2
    x1x2=((Ma3[0][1]+Ma3[1][0])*b3c**2-2*(Ma3[0][0]+Ma3[1][1])*b3c*b3d+(Ma3[0][1]+Ma3[1][0])*b3d**2)/(2*(b3c**2-b3d**2)**2)
    x1p1=(b3c**2*((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[1][1]*b1d-Ma2[0][1]*b3c)+b3c*(-2*(Ma3[0][0]+Ma3[1][1])*b1c-2*(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[0][0]*b3c)*b3d+((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[0][0]*b1d+Ma2[0][1]*b3c)*b3d**2-Ma2[0][0]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    x2p2=(b3c**2*((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[0][0]*b1d-Ma2[1][0]*b3c)+b3c*(-2*(Ma3[0][0]+Ma3[1][1])*b1c-2*(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[1][1]*b3c)*b3d+((Ma3[0][1]+Ma3[1][0])*b1c+2*Ma3[1][1]*b1d+Ma2[1][0]*b3c)*b3d**2-Ma2[1][1]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    x1p2=(b3c**2*(2*Ma3[1][1]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d-Ma2[1][1]*b3c)+b3c*(-2*(Ma3[0][1]+Ma3[1][0])*b1c-2*(Ma3[0][0]+Ma3[1][1])*b1d+Ma2[1][0]*b3c)*b3d+(2*Ma3[0][0]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[1][1]*b3c)*b3d**2-Ma2[1][0]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    x2p1=(b3c**2*(2*Ma3[0][0]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d-Ma2[0][0]*b3c)+b3c*(-2*(Ma3[0][1]+Ma3[1][0])*b1c-2*(Ma3[0][0]+Ma3[1][1])*b1d+Ma2[0][1]*b3c)*b3d+(2*Ma3[1][1]*b1c+(Ma3[0][1]+Ma3[1][0])*b1d+Ma2[0][0]*b3c)*b3d**2-Ma2[0][1]*b3d**3)/(2*(b3c**2-b3d**2)**2)
    p1p1 = Ma1[0][0]+(b3c**2*(Ma3[0][0]*b1c**2+b1d*((Ma3[0][1]+Ma3[1][0])*b1c+Ma3[1][1]*b1d)-(Ma2[0][0]*b1c+Ma2[0][1]*b1d)*b3c)-b3c*(2*(Ma3[0][0]+Ma3[1][1])*b1c*b1d+Ma3[0][1]*(b1c**2+b1d**2)+Ma3[1][0]*(b1c**2+b1d**2)-(Ma2[0][1]*b1c+Ma2[0][0]*b1d)*b3c)*b3d+(Ma3[1][1]*b1c**2+Ma3[0][1]*b1c*b1d+Ma3[1][0]*b1c*b1d+Ma3[0][0]*b1d**2+Ma2[0][0]*b1c*b3c+Ma2[0][1]*b1d*b3c)*b3d**2-(Ma2[0][1]*b1c+Ma2[0][0]*b1d)*b3d**3)/(b3c**2-b3d**2)**2
#    p1p1 = Ma1[0][0]+b1d**2*x1x1+b1c**2*x2x2+4*b1d*b1d*x1x2+2*b1d*x1p1+2*b1c*x2p1   
    p2p2 = Ma1[1][1]+(b3c**2*(Ma3[1][1]*b1c**2+b1d*((Ma3[0][1]+Ma3[1][0])*b1c+Ma3[0][0]*b1d)-(Ma2[1][1]*b1c+Ma2[1][0]*b1d)*b3c)-b3c*(2*(Ma3[0][0]+Ma3[1][1])*b1c*b1d+Ma3[0][1]*(b1c**2+b1d**2)+Ma3[1][0]*(b1c**2+b1d**2)-(Ma2[1][0]*b1c+Ma2[1][1]*b1d)*b3c)*b3d+(Ma3[0][0]*b1c**2+Ma3[0][1]*b1c*b1d+Ma3[1][0]*b1c*b1d+Ma3[1][1]*b1d**2+Ma2[1][1]*b1c*b3c+Ma2[1][0]*b1d*b3c)*b3d**2-(Ma2[1][0]*b1c+Ma2[1][1]*b1d)*b3d**3)/(b3c**2-b3d**2)**2
    p1p2 = (Ma1[0][1]+Ma1[1][0]+(((Ma2[0][1]+Ma2[1][0])*b1c+(Ma2[0][0]+Ma2[1][1])*b1d)*b3c-((Ma2[0][0]+Ma2[1][1])*b1c+(Ma2[0][1]+Ma2[1][0])*b1d)*b3d)/(-b3c**2+b3d**2)+((b1c**2+b1d**2)*((Ma3[0][1]+Ma3[1][0])*b3c**2-2*(Ma3[0][0]+Ma3[1][1])*b3c*b3d+(Ma3[0][1]+Ma3[1][0])*b3d**2))/(b3c**2-b3d**2)**2+(2*b1c*b1d*((Ma3[0][0]+Ma3[1][1])*b3c**2-2*(Ma3[0][1]+Ma3[1][0])*b3c*b3d+(Ma3[0][0]+Ma3[1][1])*b3d**2))/(b3c**2-b3d**2)**2)/2
    
    for i in [x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2]:
        i[0] = i[1]

    if unpacked:
        return x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2
    else:
        cov_matrix = np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])
        return cov_matrix

def En(t, Mcov):
    """
    Negatividad logaritmica
    """
    var = np.array([])
    for i, ti in enumerate(t):
        
        A = Mcov[:2, :2, i]
        B = Mcov[2:, 2:, i]
        C = Mcov[:2, 2:, i]
        
        delta = np.linalg.det(A)+np.linalg.det(B)-2*np.linalg.det(C)
        dete = np.linalg.det(Mcov[:, :, i])
        nmenos = np.sqrt((delta-np.sqrt(delta**2-4*dete))/2) 
        
        var = np.append(var, nmenos)
    var = -np.log2(2*var)

    var[var < 0] = 0

    return var
    

def n_menos(t, Mcov):
    """
    Menor autovalor simplectico
    """
    var = np.array([])
    for i, ti in enumerate(t):
        
        A = Mcov[:2, :2, i]
        B = Mcov[2:, 2:, i]
        C = Mcov[:2, 2:, i]
        
        delta = np.linalg.det(A)+np.linalg.det(B)+2*np.linalg.det(C)
        dete = np.linalg.det(Mcov[:, :, i])
        nmenos = np.sqrt((delta-np.sqrt(delta**2-4*dete))/2) 
        
        var = np.append(var, nmenos)
    
    return var


def heat(t, Mcov, c, V, mean = True):
    """
    Devuelve dQ1/dt, dQ2/dt
    """
#    Mcov = cov(t, g, ca1, cq1, ca2, cq2, temp1, temp2, wc, i)
    x1x1, x2x2 = Mcov[0, 0, :], Mcov[2, 2, :] 
    p1p1, p2p2 = Mcov[1, 1, :], Mcov[3, 3, :]
    x1p2, x2p1 = Mcov[0, 3, :], Mcov[2, 1, :]
    dQ1 = 1/2 * (deriv(p1p1, t)+V*deriv(x1x1, t))+c*x2p1 
    dQ2 = 1/2 * (deriv(p2p2, t)+V*deriv(x2x2, t))+c*x1p2
    if mean:
        return np.median(dQ1), np.median(dQ2)
    else:
        return dQ1, dQ2
    
def discordia(t, Mcov):
    """
    Gaussian discord
    """
    def f(x):
        return (x+1/2)*np.log2(x+1/2) - (x-1/2)*np.log2(x-1/2)
        
    discord = np.array([])
    for i, ti in enumerate(t):
        
        alpha = Mcov[:2, :2, i]
        beta = Mcov[2:, 2:, i]
        gamma = Mcov[:2, 2:, i] 
        
        A = np.linalg.det(alpha)
        B = np.linalg.det(beta)
        C = np.linalg.det(gamma)
        D = np.linalg.det(Mcov[:, :, i])
        
        delta = A+B+2*C
        
        nmenos = np.sqrt((delta-np.sqrt(delta**2-4*D))/2)
        nmas =   np.sqrt((delta+np.sqrt(delta**2-4*D))/2)
#        print nmenos
        ge = (D-A*B)**2-(1/4+B) * C**2 * (A+4*D)
#        print ge
        if ge <= 0:
            Emin =  (2*C**2+(-1/4+B)*(-A+4*D)+2*np.abs(C)*np.sqrt(C**2+(-1/4+B)*(-A+4*D))) / (4*(-1/4+B)**2)
        else:
            Emin = (A*B-C**2+D-np.sqrt(C**4+(-A*B+D)**2-2*C**2*(A*B+D))) / (2*B)
        
        dis_t = f(np.sqrt(B))-f(nmas)-f(nmenos)+f(np.sqrt(Emin))
        discord = np.append(discord, dis_t)
    
    return discord
"""

w0 = 1
wd = 2*w0
c_1 = 0.1*w0
c_0 = 0*w0
g = .0001*w0
ca1, cq1 = (4/wd**2)*(w0**2-g**2/4+c_0-2*g), -(2/wd**2)*c_1
ca2, cq2 = (4/wd**2)*(w0**2-g**2/4-c_0-2*g), (2/wd**2)*c_1
nu1, nu2 = fl.mathieu_nu(ca1, cq1), fl.mathieu_nu(ca2, cq2)
#A1, 0A2 = fl.mathieu_coefs(ca1, cq1, nu1), fl.mathieu_coefs(ca2, cq2, nu2)
i = 2
#A1, A2 = A1[A1.size//2-i:A1.size//2+i+1], A2[A2.size//2-i:A2.size//2+i+1]
t = np.linspace(0, 3, 50)
#
#
wc = 50
T = 100
T1, T2 = 40, 200
#
c = c_0+c_1*np.cos(wd*t)
V = w0**2-2*g
#covarianzas = cov(t, g, ca1, cq1, ca2, cq2, T1, T2, wc, i)
#x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2 = cov(t, g, ca1, cq1, ca2, cq2, T1, T2, wc, i, unpacked=True)
#Mcov = np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])
#heat1, heat2 = heat(t, Mcov, c, V, False)

plt.clf()
temps = np.arange(0, 20, 5)
#Q1, Q2 = np.array([]), np.array([]) 
for T in temps:
    T1, T2 = 10, T
    x1x1, x2x2, x1x2, x1p1, x2p2, x1p2, x2p1, p1p1, p2p2, p1p2 = cov(t, g, ca1, cq1, ca2, cq2, T1, T2, wc, i, unpacked=True)
    Mcov = np.array([[x1x1, x1p1, x1x2, x1p2], [x1p1, p1p1, x2p1, p1p2], [x1x2, x2p1, x2x2, x2p2], [x1p2, p1p2, x2p2, p2p2]])
    plt.plot(t, discordia(t, Mcov), label = str(T))
#    heat1, heat2 = heat(t, Mcov, c, V)
#    Q1 = np.append(Q1, heat1)
#    Q2 = np.append(Q2, heat2)
    
    print T
plt.legend()
#plt.clf()
#plt.plot(temps, Q2+Q1, 'b')
#plt.plot(temps, Q2, 'r')    

#plt.clf()
#Q1, Q2 = heat(t, Mcov, c, V, False)

#plt.axis([0, 5, -100, 100])
#t = 2/wd * t
#neg = En(t, Mcov)
#dis = discordia(t, Mcov)


#plt.clf()
#plt.plot(t, dis, 'o-')
#plt.plot(t, np.median(dis)*np.ones(len(t)))








#plt.clf()
#plt.subplot(1, 3, 1)
#plt.plot(t, x1x1, 'bo-', t, x2x2, 'go-', t, x1x2, 'ro-')
#plt.subplot(1, 2, 1)
#plt.plot(t, p1p1, 'bo-', t, p2p2, 'go-', t, p1p2, 'ro-')
#plt.subplot(1, 2, 2)
#plt.plot(t, x1p2, 'ro-', t, x2p1, 'mo-')
#plt.plot(t, np.average(heat1)*np.ones(len(t)), 'o-b', t, np.average(heat2)*np.ones(len(t)), 'o-r')

#
#dis = discordia(t, covarianzas)
#neg = En(t, covarianzas)
#plt.clf()
#plt.plot(t, dis,'-xb', label = 'discord', linewidth = 3)
#plt.plot(t, neg,'-g', label = 'En', linewidth = 3)
#plt.axis([0, np.max(t), -1, .5])
#plt.legend()

#print 'Estable: ', nu1.imag <= g/2, nu2.imag <= g/2

#plt.legend()

print 'Estable: ', nu1.imag <= g/2, nu2.imag <= g/2
print 'Limite clasico: ', np.abs(nu1.imag/g)

print 'done'

"""