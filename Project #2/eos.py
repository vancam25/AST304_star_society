########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################

"""
Description: 

pressure: A function that calculates the electron degeneracy pressure based on 
    inputs rho (mass density) and mue (baryon/electron ratio). Follows from 
    our first equation in "instructions-1.pdf".

density: A function that computes the mass density of our white dwarf star. 
    Takes inputs of electron degeneracy pressure rho and baryon/electron ratio mue. 
    Follows from a restructuring of our first equation in "instructions-1.pdf".

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
