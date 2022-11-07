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
from eos import pressure, density
from ode import fEuler, rk2, rk4
import astro_const as ac 

def stellar_derivatives(m,z,mue):
    """
    RHS of Lagrangian differential equations for radius and pressure
    
    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure)
        mue
            ratio, nucleons to electrons.  For a carbon-oxygen white dwarf, 
            mue = 2.
        
    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm
    """
    r = z[0:2]
    row = z[2:4]
    
    dzdm = np.zeros_like(z)
    
    # evaluate dzdm
    
    drdm = 1/(4*np.pi*r**2*row)
    dPdm = (ac.G*m/4*np.pi*r**4)
    dzdm = np.concatenate((drdm,dPdm))
    
    return dzdm

def central_values(Pc,delta_m,mue):
    """
    Constructs the boundary conditions at the edge of a small, constant density 
    core of mass delta_m with central pressure P_c
    
    Arguments
        Pc
            central pressure (units = ?)
        delta_m
            core mass (units = ?)
        mue
            nucleon/electron ratio
    
    Returns
        z = array([ r, p ])
            central values of radius and pressure (units = ?)
    """
    z = np.zeros(2)
    
    # compute initial values of z = [ r, p ]
    
    #From the document: 
    #m = delta_m
    #P = Pc(m)
    #p(m) = p(Pc)
    #r(m) = ((3*m)/(4*np.pi*p))**(1/3)
    
    #FROM ABRAM:
    #z[0] = Pc
    #z[1] = ((3*delta_m)/(4*np.pi*Pc))**(1/3)
    
    
    rho_i = eos.density(Pc,mue)
    r_i = ((3*m)/(4*ac.pi*rho_i))**(1/3)
    
    z[0] = r_i
    z[1] = Pc
    
    return z
    
def lengthscales(m,z,mue):
    """
    Computes the radial length scale H_r and the pressure length H_P
    
    Arguments
        m
            current mass coordinate (units = ?)
        z (array)
           [ r, p ] (units = ?)
        mue
            mean electron weight
    
    Returns
        z/|dzdm| (units = ?)
    """

    ##### COME BACK TO; NEED TO USE MUE? ###########
    
    H_r = 4*ac.pi*m
    H_p = (4*ac.pi*(z[0]**4)*z[1])/(ac.G*m)
    
    #Returns z/|dzdm| (units = ?)
    return [(H_r, H_p)]
    
def integrate(Pc,delta_m,eta,xi,mue,max_steps=10000):
    """
    Integrates the scaled stellar structure equations

    Arguments
        Pc
            central pressure (units = ?)
        delta_m
            initial offset from center (units = ?)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
        max_steps
            solver will quit and throw error if this more than max_steps are 
            required (default is 10000)
                        
    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during 
            integration (what are the units?)
    """
        
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)
    
    # set starting conditions using central values 
    z = central_values(Pc, delta_m, mue)
    
    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        # are we at the surface?
        if (pressure < eta*Pc):
            break
            
        # store the step, aka current values of m, r, p
        m_step[Nstep] = delta_m
        r_step[Nstep] = radius
        p_step[Nstep] = pressure
        
        # set the stepsize
        # remember xi: The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|
        #might change to xi[step]
        
        l = lengthscale(m_step[step], z, mue)
        stepsize = xi*min(pressure/l[1], radius/l[0])
        
        # take a step
        # USE RK4
        
        z[step] = rk4(stellar_derivatives, delta_m, z[step-1], stepsize, mue) 
        
        delta_m = m_step[Nstep] + stepsize
        
        radius = z[0]
        pressure = z[1]
        
        # increment the counter
        Nsteps += 1
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')
        
    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]

def pressure_guess(m,mue):
    """
    Computes a guess for the central pressure based on virial theorem and mass-
    radius relation. 
    
    Arguments
        m
            mass of white dwarf (units are ?)
        mue
            mean electron mass 
    
    Returns
        P
            guess for pressure
    """
    # fill this in
    #take stuff from integrate? Use pressure eq (16???) from instructions document
    #do we need all values of M from m_step??? no idea <3 
    
    Pguess = (ac.G**5/ac.Ke**4)*(M*mue**2)**10/3
    
    return Pguess
