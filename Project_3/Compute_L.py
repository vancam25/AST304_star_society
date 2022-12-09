"""
####Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
#### AST 304, Fall 2022
#### Michigan State University
"""

## INSTRUCTIONS: 

## Write a function that computes Lnuc(R) − 4πR**2σSBTeff**4. Here Lnuc(R) is
## computed by integrating over the star’s structure, and Teff is found by
## interpolating from the table of masses and effective temperatures. The
## root of this equation is the star’s main-sequence radius. 
## Solve for this radius using a rootfind routine, either scipy.optimize.bisect or its more
## efficient cousin scipy.optimize.brentq.


from zams import Teff 
from astro_const import fourpi, sigmaSB, G, Msun, Rsun, Lsun, kB, m_u
from structure import central_thermal, stellar_derivatives, central_values, lengthscales, integrate
from eos_template import get_rho_and_T, mean_molecular_weight
from ode import rk4
from reactions_template import pp_rate
from scipy import optimize 



def Compute_L(Mwant, Pc, delta_m, delta_r, eta, xi, mu, XH, pp_factor):
    
    """
    
    Computes Lnuc(R) - 4πR**2σSBTeff**4, via Teff found by interpolation in
    function Teff. Root of this function is the star's main-sequence radius, 
    found via bisectioning. 
    
    Parameters
        Teff (float, scalar or array)
            Interpolated effective temperature of the star in Kelvin [K].
            
    Returns
        Lnuc(R) - 4πR**2σSBTeff**4 (float, scalar or array)
            (FINISH THIS DESCRIPTION) 
    
    """
    
    
    Teff = Teff(Mwant)
    
    get_integrate = integrate(Pc, delta_m, delta_r, eta, xi, mu, XH, pp_factor)
    Lnuc = get_integrate[3] 
    r_test = get_integrate[1]
    
    subtraction = Lnuc - fourpi*(R**2)*sigmaSB*(Teff**4)
    function = (Lnuc/(fourpi*sigmaSB*(Teff**4)))**(1/2)
    
    main_sequence_radius = scipy.optimize.brentq(function, 0, r_test)
    
    return main_sequence_radius 