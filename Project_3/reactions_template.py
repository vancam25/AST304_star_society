# 
# <Star Society, Abram Anderson, Steven Vancamp, Sanskriti Verma, Hannah Sullivan>
# 

import numpy as np

def pp_rate(T,rho,XH,pp_factor=1.0):
    """
    Specific heating rate from pp chain hydrogen burning. Approximate rate 
    taken from Hansen, Kawaler, & Trimble.
    
    Arguments
        T, rho
            temperature [K] and density [kg/m**3]
        XH
            mass fraction hydrogen
        pp_factor
            multiplicative factor for rate
    Returns
        heating rate from the pp-reaction chain [W/kg]
    """
    
    # fill this out
    rate = pp_factor*((0.0024*rho*XH**2)/((T/10**9)**(2/3))*np.exp(-3.38/(T/10**9)**(1/3)))
    return rate
