# -*- coding: utf-8 -*-
"""
Mass-Damper-Spring ODE solver
Test run for bubble solver
"""
import numpy as np
import scipy as sp
import numpy.lib.scimath as sm
import matplotlib as mpl
mpl.rc('font', family='serif', size=16)
mpl.rc('xtick',labelsize='small')
mpl.rc('ytick',labelsize='small')
mpl.rc('legend',fontsize='small')
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Defining the ODE to be integrated
def ForcedMassDamperSpring(ICs,t,Constants,F):
    # Unpack the state vector
    x = ICs[0]
    v = ICs[1]
    
    # Unpack the constants
    m = Constants[0]
    b = Constants[1]
    k = Constants[2]
    
    F_int = F #sp.interpolate.interp1d(t,F(t))?
      
    # Here are the equations being solved 
    xd = v
    vd = (-b/m)*v+(-k/m)*x+F_int
    return [xd,vd]


# Set initial conditions
x0 = 1       # Initial position
v0 = 0       # Initial velocity
ICs = np.array([x0,v0])


# Here the constants are defined 
m = 1
b = 0.1
k = 1
Constants = np.array([m,b,k])


# Time frame for integration
t = np.arange(0.0, 20.0, 0.1)


# Here, the forcing term is defined
F = 0 # 0.1*np.sin(2*np.pi*0.5*t)


# Integration; yinfo turned on for error checking; note the possible need for stricter tolerances than default
RelTol = 1e-6
AbsTol = RelTol/100
numsolution, yinfo = odeint(ForcedMassDamperSpring, ICs, t, args = (Constants, F,), rtol = RelTol, atol = AbsTol,full_output = 1)
xnum = numsolution[:,0]
vnum = numsolution[:,1]


# Exact solution, provided F = 0
omega_param = sm.sqrt(k/m)
gamma_param = sm.sqrt((b**2-4*k*m))/(2*m)
xexact = np.real((x0*np.exp(gamma_param*t)+v0*np.exp(-gamma_param*t))*np.exp(-b*t/(2*m)))

"""
# Compare errors
mse = np.sqrt((xexact-xnum)**2)
"""


# Plot of results
plt.figure(1)
plt.plot(t,xnum,'bx')
plt.plot(t,xexact,'g')
plt.legend(('Numerical','Exact'))
plt.xlabel(r'Time [s]')
plt.ylabel(r'Position [m]')

