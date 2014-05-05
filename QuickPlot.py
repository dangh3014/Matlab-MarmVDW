# -*- coding: utf-8 -*-
"""
Some plotting commands, calculation of radiated pressure, and fft
"""

# File setup
import numpy as np
# import scipy as sp # (Not needed)
import matplotlib as mpl
mpl.rc('font', family='serif', size=16)
mpl.rc('xtick',labelsize='small')
mpl.rc('ytick',labelsize='small')

import matplotlib.pyplot as plt


# Some test data for plotting
t = np.linspace(0,1e-6,num=1000)      # Time vector
frequency = 5e6
P = 1e6*np.sin(2*np.pi*frequency*t)      # Pressure pulse
R = 1e-6*np.sin(2*np.pi*frequency*t)      # Radius
V = 1e-6*(2*np.pi*frequency)*np.cos(2*np.pi*frequency*t)     # Velocity
A = -1e-6*(2*np.pi*frequency)**2*np.sin(2*np.pi*frequency*t)     # Acceleration
P_rad = (1000/(2e-6))*(A*R**2+2*R*V**2)     # Radiated pressure


# Compute Fourier transform / power spectrum
data = R

xincr = (t[-1]-t[0])/len(t) # Sample time
fs = 1/xincr     # Sampling frequency
NFFT = 2**14    # Number of points in FFT
vectfs = fs*np.arange(NFFT/2)/NFFT # Frequency axis

# Zero padded fft / ps
ft_data = np.fft.fft(data,NFFT)/NFFT      # Normalized fft   
ps_data = (abs(ft_data[0:(NFFT/2)])**2)

"""
# Exact fft
vectfs2 = np.fft.fftfreq(data.size,xincr)
ft_data2 = np.fft.fft(data)
"""


# 4 Sample Plots
# Acoustic Pulse
plt.figure(1)
plt.plot(t*10**6,P*10**-6)
plt.xlabel(r'Time [$10^{-6}$ s]')  # ,fontsize=16
plt.ylabel(r'Pressure [MPa]')

"""
# Bubble Radius-Time
plt.figure(2)
plt.plot(t*10**6,R*10**6)
plt.xlabel(r'Time [$10^{-6}$ s]')  # ,fontsize=16
plt.ylabel(r'Radius [$10^{-6}$ m]')

# Bubble Radiated Pressure
plt.figure(3)
plt.plot(t*10**6,P_rad*10**-6)
plt.xlabel(r'Time [$10^{-6}$ s]')  # ,fontsize=16
plt.ylabel(r'Pressure [MPa]')
"""

# Bubble Power Spectrum
plt.figure(4)
plt.plot(vectfs*10**-6,ps_data)
plt.xlim(0,10)
plt.xlabel(r'Frequency [MHz]')
plt.ylabel(r'Power Spectrum [A.U.]')

