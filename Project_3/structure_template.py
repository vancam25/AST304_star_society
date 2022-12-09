""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

<Star Society, Abram Anderson, Steven Vancamp, Sanskriti Verma, Hannah Sullivan>
"""

import numpy as np
from eos import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi
from reactions import pp_rate

def central_thermal(m,r,mu):
    """ 
    Computes the central pressure, density, and temperature from the polytropic
    relations for n = 3/2.

    Arguments
        m
            mass in solar units
        r
            radius is solar units
        mu
            mean molecular weight
    Returns
        Pc, rhoc, Tc
            central pressure, density, and temperature in solar units
    """
    # fill this in
    G = 6.6743*10**-11
    Pc = 0.77*G*(m**2)/(r**4)
    rhoc = 5.99*3*m/(4*np.pi*r**3)
    Tc = 0.54*(mu*1.66*10**-27/(1.380649*10**-23))*(G*m/r)
    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
def stellar_derivatives():
    pass

def central_values():
    pass

def lengthscales():
    pass

def integrate():
    pass
