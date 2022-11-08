########################################################################
# Team <your team name>: <names>
# AST304, Fall 2020
# Michigan State University
########################################################################

"""
<Description of this module goes here: what it does, how it's used.>
"""

import astro_const as ac

def pressure(rho, mue):
    """
    Arguments
        rho
            mass density (kg/m**3)
        mue
            baryon/electron ratio
    
    Returns
        electron degeneracy pressure (Pascal)
    """
    
    p = ac.K_e * (rho/mue)**(5/3) # Eq.(1) of "instructions-1.pdf"
    
    return(p)

def density(p, mue):
    """
    Arguments
        p
            electron degeneracy pressure (Pascal)
        mue
            baryon/electron ratio
        
    Returns
        mass density (kg/m**3)
    """
    
    rho = (p/ac.K_e)**(3/5) * mue # Eq.(1) of "instructions-1.pdf"
    
    return(rho)
