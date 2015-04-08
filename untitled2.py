# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 10:52:57 2015

@author: martin
"""
import numpy as np
from R_bajaT import R as Rint
import pylab as plt
import sigma as si



#c = np.loadtxt('hanggi.txt', np.complex)
c = np.loadtxt('0.85_0.1_10.txt', np.complex)

nu = c[0]
c = c[1:]
i = 5
c = c[c.size//2-i:c.size//2+i+1]

e = 1
t = np.linspace(0,5, 50)
wc = 50
g = 1

n, m, k, l = 1,2,3,4
ex = Rint(wc, t, n, m, k, l, g, nu)-Rint(0, t, n, m, k, l, g, nu)
num = si.int_baja_num(wc, t, n, m, k, l, g, nu)

plt.clf()
#v = [0, 5, -5, 5]
plt.subplot(1, 2, 1)
plt.plot(t, ex.imag, t, num.imag, 'o')
plt.title('imag')
#plt.axis(v)
plt.subplot(1, 2, 2)
plt.plot(t, ex.real, t, num.real, 'o')
plt.title('real')
#plt.axis(v)
print 'done'