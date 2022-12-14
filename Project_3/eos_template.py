"""
Routines to compute an adiabatic equation of state.

<Star Society, Abram Anderson, Steven Vancamp, Sanskriti Verma, Hannah Sullivan>
"""

import numpy as np

def mean_molecular_weight(Z,A,X):
    """Computes the mean molecular weight for a fully ionized plasma with an 
    arbitrary mixture of species
    
    Arguments
        Z, A, X (either scaler or array)
            charge numbers, atomic numbers, and mass fractions
            The mass fractions must sum to 1
    """
    Zs = np.array(Z)
    As = np.array(A)
    Xs = np.array(X)
    assert np.sum(Xs) == 1.0
    
    
    # compute value of mean molecular weight
    mu = np.sum(Xs*(1+Zs)/As)  #formula of mean molecular weight for a fully ionized plasma
    mu =mu**(-1) #since the formula above is the formula for the inverse of mu, we inverse it over here to return the correct computation
    return mu
    
def get_rho_and_T(P,P_c,rho_c,T_c):
    """
    Compute density and temperature along an adiabat of index gamma given a 
    pressure and a reference point (central pressure, density, and temperature).
    
    Arguments
        P (either scalar or array-like)
            value of pressure
    
        P_c, rho_c, T_c
            reference values; these values should be consistent with an ideal 
            gas EOS
    
    Returns
        density, temperature
    """

    # replace with computed values
    rho = rho_c*(P/P_c)**(3/5) 
    T = T_c*((P/P_c)**(1-1/(5/3)))

    return rho, T
