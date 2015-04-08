# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 19:30:30 2015

@author: martin
"""
import numpy as np
import pylab as plt



def plot_ft(t, x):
    """
    asd
    """
    dt = t[1]-t[0]
    f = np.fft.fft(x)
    freqs = np.fft.fftfreq(x.size)*2*np.pi/dt
    
    f = f[:t.size//2]
    freqs = freqs[:t.size//2]
    plt.plot(freqs, np.abs(f))
    plt.show()
    
    return 0

#t = np.linspace(0,10,10000)
#x = np.sinc(t)
#
#plt.figure(1)
#plt.clf()
#plt.plot(t, x)
#plt.figure(2)
#plt.clf()
#plot_ft(t, x)

#plt.plot(t, x)
#
#dt = t[1]-t[0]
#f = np.fft.fft(x)
#freqs = np.fft.fftfreq(x.size)*2*np.pi/dt
#
#plt.figure()
#plt.plot(freqs[:50], f[:50], 'o')