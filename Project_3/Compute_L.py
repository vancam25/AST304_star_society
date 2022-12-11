"""
####Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
#### AST 304, Fall 2022
#### Michigan State University
"""


import scipy
import zams as zms
import structure as stc

from astro_const import Rsun


def Compute_L(Mwant, RWant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, mx):
    
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
    
    get_integrate = stc.integrate(Mwant, RWant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, mx)
    Lnuc = get_integrate[3] 
    r_test = get_integrate[1]
    
    #surf_lum = zms.surface_luminosity(teff,r_test)
    #subtraction = Lnuc - surf_lum
    
    def bisect_input(Lnuc, r_test, teff):
        
        surf_lum = zms.surface_luminosity(teff,r_test)
        subtraction = Lnuc - surf_lum
        
        return subtraction
    
    
    main_sequence_radius = scipy.optimize.brentq(bisect_input, a = 0, b = def_delta_r, args = ((Lnuc,r_test),teff))
    
    return main_sequence_radius 