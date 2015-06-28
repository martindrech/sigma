# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 13:04:22 2015

@author: martin
"""

from __future__ import division
import numpy as np
import pylab as plt
import sigma as si
import floquet as fl
import erres as R

t = np.linspace(0, 30, 1000)
a1, q1 = 1, .5
a2, q2 = 1.1, .5
nu1, nu2 = fl.mathieu_nu(a1, q1), fl.mathieu_nu(a2, q2)
c1, c2 = fl.mathieu_coefs(a1, q1, nu1, 3), fl.mathieu_coefs(a1, q1, nu1, 3)
phi1, dphi1, phi2, dphi2 = fl.mathieu(a1, q1, t)
phim1, dphim1, phim2, dphim2 = fl.mathieu(a2, q2, t)
g = 1
plt.clf()
plt.plot(t, phi1, 'g-', t, phi2, 'b-')



