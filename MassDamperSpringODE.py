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
from scipy.interpolate import interp1d


# Defining the ODE to be integrated
def ForcedMassDamperSpring(ICs,t,Constants,F,T):
    # Unpack the state vector
    x = ICs[0]
    v = ICs[1]
    
    # Unpack the constants
    m = Constants[0]
    b = Constants[1]
    k = Constants[2]
    
    # Intepolate the force
    set_interp = sp.interpolate.interp1d(T,F)
    F_int = set_interp(t)      
     
    # Here are the equations being solved 
    xd = v
    vd = (-b/m)*v+(-k/m)*x+F_int
    return [xd,vd]


# Set initial conditions
x0 = 1      # Initial position
v0 = 0      # Initial velocity
ICs = np.array([x0,v0])

# Define the constants 
m = 0.1
b = 0.1
k = 1
Constants = np.array([m,b,k])


#####
# This next section is primarily for error checking
# Calculate some parameters for the exact solutions
freq = 0.25 # Frequency of sinusoidal forcing in Hz
ang_freq = 2*np.pi*freq # Angular frequency
F0 = 0.5 # Amplitude of sinusoidal forcing

omega0_param = sm.sqrt(k/m) # Undamped angular frequency
zeta_param = b/(2*sm.sqrt(k*m)) # Damping ratio
if zeta_param >= 1:
    exit('Please choose a smaller zeta parameter (<1)...exiting')
else:
    pass
gamma_param1 = 0.5*sm.sqrt(4*omega0_param**2-b**2)
gamma_param = sm.sqrt((b**2-4*k*m))/(2*m)
Zm_param = sm.sqrt((2*omega0_param*zeta_param)**2+((omega0_param**2-ang_freq**2)**2)/(ang_freq**2))
psi_param = np.arctan((b*ang_freq)/(k-m*ang_freq**2))-np.pi

A_param = (F0)/(np.sqrt((omega0_param**2-ang_freq**2)**2)+4*(b/(2*m))**2*ang_freq**2)
psih_param = np.arctan(np.sqrt(omega0_param**2-(b/(2*m))**2)*(x0-A_param*np.cos(psi_param))/(v0+(b/(2*m))*(x0-A_param*np.cos(psi_param))-A_param*ang_freq*np.sin(psi_param)))
Ah_param = (x0-A_param*np.cos(psi_param))/np.sin(psih_param)
#####



# Time frame for integration (t) and interpolation (T)
t_beg = 0.0
t_end = 40.0
t_step = 0.1

t = np.arange(t_beg, t_end, t_step)
T = np.arange(t_beg, 2*t_end, t_step) # Should be ~2x as long for odeint


# Here, the forcing term is defined
Ftype = 'Zero' # Options: 'Zero', 'Step', 'Sine'

if Ftype=='Zero': # Zero forcing
    F = 0*T 
elif Ftype=='Step': # Step forcing
    F = omega0_param**2*np.ones(len(T)) # Sinusoidal forcing 
elif Ftype=='Sine':
    F = F0*np.cos(ang_freq*T+np.pi) 
else:
    exit('Undefined forcing type...exiting')


#####
# Again, for error checking
# Here are the exact solutions for several cases
if Ftype=='Zero': 
    # xexact = np.real((x0*np.exp(gamma_param*t)+v0*np.exp(-gamma_param*t))*np.exp(-b*t/(2*m)))
    xexact = np.exp(-b/(2*m)*t)*(x0*np.cos(np.sqrt(k/m-(b/(2*m))**2)*t)+(b*x0/(2*gamma_param1)+v0/gamma_param1)*np.sin(np.sqrt(k/m-(b/(2*m))**2)*t))
elif Ftype=='Step': # If zeta < 1, step input, the system response
    xexact = 1-np.exp(-zeta_param*omega0_param*t)*np.sin(sm.sqrt(1-zeta_param**2)*omega0_param*t+np.arccos(zeta_param))/np.sin(np.arccos(zeta_param))
elif Ftype=='Sine': # If zeta < 1, sinusoidal input, the steady state response
    # xexact = F0/(m*Zm_param*ang_freq)*np.sin(ang_freq*t+psi_param)
    xexact = Ah_param*np.exp(-(b*t/(2*m)))*np.sin(np.sqrt(omega0_param**2-(b/(2*m))**2)*t+psih_param)+A_param*np.cos(ang_freq*t-psi_param)
else:
    pass
#####

# Integration; yinfo turned on for error checking; note the possible need for stricter tolerances than default
RelTol = 1e-5
AbsTol = RelTol/100
numsolution, yinfo = odeint(ForcedMassDamperSpring, ICs, t, args = (Constants, F, T), rtol = RelTol, atol = AbsTol,full_output = 1)
xnum = numsolution[:,0]
vnum = numsolution[:,1]


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

