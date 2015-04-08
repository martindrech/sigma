# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 10:36:49 2015

@author: martin
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si
import time

#c = np.loadtxt('hanggi.txt', np.complex)
c = np.loadtxt('0.75_-0.5.txt', np.complex)

nu = c[0]
c = c[1:]

i = 2
c = c[c.size//2-i:c.size//2+i+1]

e = 1
t = np.linspace(0,5, 30)
wc = 50
g = 1

#tic = time.clock()
s1 = si.sigma_baja(t, c, nu, g, wc, e)
s2 = si.sigma_num(t, c, nu, g, wc, e, .1)
s3 = si.sigma_num(t, c, nu, g, wc, e, 1)
s4 = si.sigma_num(t, c, nu, g, wc, e, 5)
s5 = si.sigma_num(t, c, nu, g, wc, e, 10)


#sn = si.sigma_baja_num(t, c, nu, g, wc, e)
#se = si.sigma_baja_est(t, np.sqrt(.75), g, wc)[0]

#tac = time.clock()
#ti = tac-tic  
#print 'time: ', ti

#A = np.matrix([s1, s2, s3, s4, s5])

#sxxn = sn[0]
#sxxe, sxxe = se[0], se[1]

#n = 1
#sxx = s1
plt.clf()
plt.plot(t, s4[0])

#plt.plot(t, ex.imag, 'b',t, num.imag, 'ro')
#plt.plot(t, ex, 'b', t, np.imag(num), 'ro-')


print 'done'