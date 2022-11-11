########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################


#####FILL OUT DESCRIPTION OF FUNCTIONS ######
"""
<Description of this module goes here: what it does, how it's used.>
"""

import numpy as np
import eos
import astro_const as ac

from ode import rk4

def stellar_derivatives(m, z, mue):
    """
    Description:
        RHS of Lagrangian differential equations for radius and pressure
    
    Arguments:
        m (float): 
            Current mass value; units [kg]
            
        z (array): 
            Current radius & pressure values; like [radius, pressure]; units [m],[Pa]
            
        mue (float):
            Nucleon/electron ratio
        
    Returns:
        dzdm (array):
            Lagrangian derivatives dr/dm, dP/dm; like [dr/dm, dP/dm]
    """
    
    rho = eos.density(z[1], mue) # current density
    
    drdm = 1 / (4 * ac.pi * z[0]**2 * rho) # Eq.(4) of "instructions-1.pdf"
    dPdm = -(ac.G * m) / (4 * ac.pi * z[0]**4) # Eq.(5) of "instructions-1.pdf"
    dzdm = np.array([drdm, dPdm])
    
    return(dzdm)

def central_values(Pc, delta_m, mue):
    """
    Description:
        Constructs the boundary conditions at the edge of a small, constant density 
        core of mass delta_m with central pressure P_c
    
    Arguments:
        Pc (float):
            Central pressure; units [Pa]
            
        delta_m (float):
            "Core" mass; units [kg]
            
        mue (float):
            Nucleon/electron ratio
    
    Returns:
        z (array):
            Current radius & pressure values; like [radius, pressure]; units [m],[Pa]
    """
    z = np.zeros(2)
    
    rho_i = eos.density(Pc,mue) # calculate the density of the "core"
    r_i = ((3 * delta_m) / (4 * ac.pi * rho_i))**(1/3) # Eq.(9) of "instructions-1.pdf"
    
    z[0] = r_i
    z[1] = Pc
    
    return(z)
    
def lengthscales(m, z, mue):
    """
    Description:
        Computes the radial length scale H_r and the pressure length H_P
    
    Arguments:
        m (float): 
            Current mass value; units [kg]
            
        z (array):
            Current radius & pressure values; like [radius, pressure]; units [m],[Pa]
           
        mue (float)
            Mean electron weight ???
    
    Returns:
        z/|dzdm| (array):
            Current lengthscale for radius and pressure; like [Hr, Hp]
    """
    
    rho = eos.density(z[1], mue)
    
    Hr = 4 * ac.pi * z[0]**3 * rho # Eq.(10) of "instructions-1.pdf"
    Hp = (4 * ac.pi * z[0]**4 * z[1]) / (ac.G * m) # Eq.(11) of "instructions-1.pdf"
    
    Hr_Hp = np.array([Hr, Hp])
    
    return(Hr_Hp)
    
def integrate(Pc, delta_m, eta, xi, mue, max_steps=10000, err_max_step = False):
    """
    Description:
        Integrates the scaled stellar structure equations

    Arguments:
        Pc (float):
            Central pressure; units [Pa]
            
        delta_m (float):
            Core mass; units [kg]
            
        eta (float):
            The integration stops when P < eta * Pc
            
        xi (float):
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
            
        mue (float):
            mean electron mass
            
        max_steps (int):
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
            
        err_max_step (bool):
            boolean, when true dont through exception if the number of steps in the 
            integration loop 'Nsteps' has surpased 'max_steps'
            
    Returns:
        m_step, r_step, p_step 
            arrays containing mass coordinates, radii and pressures during 
            integration; units [kg], [m], [Pa]
    """
    
    # initialize parameter arrays with a length equal to the set max number of steps
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    
    # set starting conditions using central values 
    z = central_values(Pc, delta_m, mue)
    
    Nsteps = 0
    max_step_reached = 0
    
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        
        # are we at the surface?
        if (pressure < eta*Pc):
            break
           
        # store the step, aka current values of m, r, p
        m_step[step] = delta_m
        r_step[step] = radius
        p_step[step] = pressure
        
        # set the stepsize
        stepsize = xi * min(lengthscales(m_step[step], z, mue))
        
        # take a step
        # USE RK4
        
        z = rk4(stellar_derivatives, delta_m, z, stepsize, mue) 
        
        delta_m = m_step[step] + stepsize
        
        # increment the counter
        Nsteps += 1
        
    # if the loop runs to max_steps, then signal an error
    else:
        if err_max_step:
            max_step_reached = 1
        else:
            max_step_reached = 1
            raise Exception('too many iterations')
        
    return m_step[0:Nsteps], r_step[0:Nsteps], p_step[0:Nsteps], max_step_reached

def pressure_guess(m, mue):
    """
    Description:
        Computes a guess for the central pressure based on virial theorem and mass-
        radius relation. 
    
    Arguments:
        m (float): 
            Total mass value; units [kg]
        
        mue (float):
            Mean electron mass 
    
    Returns:
        P_guess (float):
            Guess for pressure; units [Pa]
    """
    
    P_guess = (ac.G**5 / ac.K_e**4) * (m * mue**2)**(10/3) # Eq.(16) of "instructions-1.pdf"
    
    return P_guess
