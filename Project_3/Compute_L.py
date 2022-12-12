"""
####Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
#### AST 304, Fall 2022
#### Michigan State University
"""

# import scipy
import zams as zms
import structure as stc
import astro_const as ac

# from astro_const import Rsun
from scipy.optimize import brentq

def Compute_L(Mwant, Rwant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, mx):
    
    """
    
    Computes Lnuc(R) - 4πR**2σSBTeff**4, via Teff found by interpolation in
    function Teff. Root of this function is the star's main-sequence radius, 
    found via bisectioning. 
    
    Parameters
        Teff (float, scalar or array)
            Interpolated effective temperature of the star in Kelvin [K].
            
    Returns
        main_sequence_radius (float, scalar) #???
            Main sequence radius of the star found via bisectioning of Lnuc(R)-4πR**2σSBTeff**4. 
            In [] ##what are the units dear god help 
    
    """
    
    
    teff = zms.Teff(Mwant)
    
    def bisect_input_func(r_test, teff):        
        L_nuc = stc.integrate(Mwant, r_test, def_delta_m, def_delta_r, 
                              def_eta, def_xi, comp_array, def_pp_factor, 
                              mx, err_max_step=True)[3][-1]
        
        surf_lum = zms.surface_luminosity(teff,r_test*ac.Rsun)
        
        subtraction = L_nuc - surf_lum
        
        return subtraction
    
    a = Rwant*1E-8
    b = 0
    b_res = 0
    b_offset = 0
    
    while b_res >= 0:
        b_offset += 1
        b = Rwant * b_offset
        b_res = bisect_input_func(b, teff)
    
    main_sequence_radius = brentq(bisect_input_func, a, b, args = (teff))
    
    return main_sequence_radius