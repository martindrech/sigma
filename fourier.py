# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 19:30:30 2015

@author: martin
"""
import numpy as np
import pylab as plt



def plot_ft(t, x):
    """
    Plotea fft de x(t)
    """
    dt = t[1]-t[0]
    f = np.fft.fft(x)
    freqs = np.fft.fftfreq(x.size)*2*np.pi/dt
    
    f = f[:t.size//2]
    freqs = freqs[:t.size//2]
    plt.plot(freqs, np.abs(f), 'o')
    plt.show()
    
    return freqs[np.argmax(np.abs(f[1:]))+1]
    
def fourier_series_coeff_numpy(f, T, N, return_complex=False):
    """Calculates the first 2*N+1 Fourier series coeff. of a periodic function.
    
    Given a periodic, function f(t) with period T, this function returns the
    coefficients a0, {a1,a2,...},{b1,b2,...} such that:

    f(t) ~= a0/2+ sum_{k=1}^{N} ( a_k*cos(2*pi*k*t/T) + b_k*sin(2*pi*k*t/T) )

    If return_complex is set to True, it returns instead the coefficients
    {c0,c1,c2,...}
    such that:

    f(t) ~= sum_{k=-N}^{N} c_k * exp(i*2*pi*k*t/T)

    where we define c_{-n} = complex_conjugate(c_{n})

    Refer to wikipedia for the relation between the real-valued and complex
    valued coeffs at http://en.wikipedia.org/wiki/Fourier_series.

    Parameters
    ----------
    f : the periodic function, a callable like f(t)
    T : the period of the function f, so that f(0)==f(T)
    N_max : the function will return the first N_max + 1 Fourier coeff.

    Returns
    -------
    if return_complex == False, the function returns:

    a0 : float
    a,b : numpy float arrays describing respectively the cosine and sine coeff.

    if return_complex == True, the function returns:

    c : numpy 1-dimensional complex-valued array of size N+1

    """
    # From Shanon theoreom we must use a sampling freq. larger than the maximum
    # frequency you want to catch in the signal. In this case, we know in
    # advance that it is 50 Hz because we are manufacturing the signal
    f_sample = 2 * N
    # we also need to use an integer sampling frequency, or the
    # points will not be equispaced between 0 and 1. We then add +2 to f_sample
    t, dt = np.linspace(0, T, f_sample + 2, endpoint=False, retstep=True)

    y = np.fft.rfft(f(t)) / t.size

    if return_complex:
        return y
    else:
        y *= 2
        return y[0].real, y[1:-1].real, -y[1:-1].imag

def coeficientes(f, T, N):
    """
    Devuelve los coeficientes para la reconstruccion: desde el -(N+1) hasta el N+1.
    """
    c = fourier_series_coeff_numpy(f, T, N, True)
    d = c[1:][::-1]
    c = np.append(d, c)
    return c

def rebuild(f, T, N, t):
    """
    Reconstruye la funcion f, de período T, usando 2(N+1)+1 terminos de fourier.
    Es decir, devuelve y = f(t) ~= sum_{k=-(N+1)}^{N+1} c_k * exp(i*2*pi*k*t/T)
    """
    k = np.arange(-N-1, N+2)
    y = np.zeros(np.size(t))
    c = coeficientes(f, T, N)
    for i, tau in enumerate(t):
            y[i] = np.dot(c, np.exp(k*1j*2*np.pi*tau/T))
    return y

###############################################################################


