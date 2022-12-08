""" 
Routines for computing structure of fully convective star of a given mass and 
radius.

<Team name, members go here>
"""

import numpy as np
import astro_const as ac

from eos_template import get_rho_and_T, mean_molecular_weight
from ode import rk4
from astro_const import G, Msun, Rsun, Lsun, kB, m_u, fourpi
from reactions_template import pp_rate

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
    Pc = 0.77*ac.G*m**2/r**4
    rhoc = 5.99*3*m/(4*ac.pi*r**3)
    Tc = 0.54*(mu*1.66*(10**-27)/(1.380649*10**-23))*(G*m/r)
    
    return Pc, rhoc, Tc

# The following should be modified versions of the routines you wrote for the 
# white dwarf project
def stellar_derivatives(m, z, rho, ep_nuc):
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
    
    drdm = 1 / (4 * ac.pi * z[0]**2 * rho) # Eq.(4) of "instructions-1.pdf"
    dPdm = -(ac.G * m) / (4 * ac.pi * z[0]**4) # Eq.(5) of "instructions-1.pdf"
    dLdm = ep_nuc
    
    dzdm = np.array([drdm, dPdm, dLdm])
    
    return(dzdm)

def central_values(P_c, Rho_c, T_c, delta_m, XH, pp_factor):
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
    z = np.zeros(3)
    
    r_i = ((3 * delta_m) / (4 * ac.pi * Rho_c))**(1/3) # Eq.(9) of "instructions-1.pdf"
    
    z[0] = r_i
    z[1] = P_c
    z[2] = pp_rate(T_c, Rho_c, XH, pp_factor) * delta_m
    
    return(z)
    
def lengthscales(m, z, rho, ep_nuc):
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
    
    Hr = (4/3) * ac.pi * z[0]**3 * rho # Eq.(10) of "instructions-1.pdf"
    Hp = (4 * ac.pi * z[0]**4 * z[1]) / (ac.G * m) # Eq.(11) of "instructions-1.pdf"
    Hl = z[2] / ep_nuc # Eq.(14) of "instructions-3.pdf", the nuclear heating rate; units [W kg^-1]
    
    Hr_Hp_Hl = np.array([Hr, Hp, Hl])
    
    return(Hr_Hp_Hl)
    
def integrate(Pc, delta_m, delta_r, eta, xi, mu, XH, pp_factor, max_steps=10000, err_max_step = False):
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
    l_step = np.zeros(max_steps)
    
    # set starting conditions using central values 
    P_c, Rho_c, T_c = central_thermal(delta_m, delta_r, mu)
    z = central_values(P_c, Rho_c, T_c, delta_m, XH, pp_factor)
    
    Nsteps = 0
    max_step_reached = 0
    
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        luminosity = z[2]
        
        rho, T = get_rho_and_T(pressure, P_c, Rho_c, T_c)
        
        ep_nuc = pp_rate(T, rho, XH, pp_factor)
        
        # are we at the surface?
        if (pressure < eta*Pc):
            break
           
        # store the step, aka current values of m, r, p
        m_step[step] = delta_m
        r_step[step] = radius
        p_step[step] = pressure
        l_step[step] = luminosity
        
        # set the stepsize
        stepsize = xi * min(lengthscales(m_step[step], z, rho, ep_nuc))
        
        # take a step
        # USE RK4     
        z = rk4(stellar_derivatives, delta_m, z, stepsize, rho, ep_nuc) 
        
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
        
    return m_step[0:Nsteps], r_step[0:Nsteps], p_step[0:Nsteps], l_step[0:Nsteps], max_step_reached