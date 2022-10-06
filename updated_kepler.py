########################################################################
# Team <your team name>: <names>
# AST 304, Fall 2022
# Michigan State University
# This header (minus this line) should go at the top of all code files.
########################################################################

"""
<Description of this module goes here: what it does, how it's used.>
"""

import numpy as np
import math
# ode.py is your routine containing integrators
from ode import fEuler, rk2, rk4

# use this to identify methods
integration_methods = {
    'Euler':fEuler,
    'RK2':rk2,
    'RK4':rk4
    }
    
# energies
def kinetic_energy(v):
    '''
    description:
        Calculate the kinetic energy of the system
        
    arguments:
        v: velocity vector of the orbiting body; np.array([v_x,v_y]), in units []
            
    returns:
        KE: kinetic energy per unit mass; float, in units []
    '''
    
    return(0.5*np.dot(v,v))

def potential_energy(r,m):
    '''
    description:
        Calculate the potential energy of the system
        
    arguments:
        r: position vector of the orbiting body; np.array([x,y]), in units [au]
        
        m: total mass of the system; float, in units [solar mass']
        
    returns:
        PE: potential energy per unit mass; float, in units []
    '''
    
    return(-m/np.linalg.norm(r))

def total_energy(z,m):
    '''
    description:
        Calculate the total energy of the system
        
    arguments:
        z: compound position and velocity vector; np.array([x,y,v_x,v_y]), in respective units
            
        m: total mass of the system; float, in units [solar mass'']
                
    returns:
        total_energy: total energy per unit mass; float, in units []
    '''
    
    r = z[0:2]  # get position vector
    v = z[2:4]  # get velocity vector 
    
    return(kinetic_energy(v) + potential_energy(r, m))


def derivs(t,z,m):
    """
    Computes derivatives of position and velocity for Kepler's problem 
    
    Arguments
        <t> The time
        
        <z> This is an array containing the positions of x and y, along with the derivatives of x and y
        
        <m> This is the reduced mass of the Earth - Sun system.
    
    Returns
        numpy array dzdt with components [ dx/dt, dy/dt, dv_x/dt, dv_y/dt ]
    """
    # Fill in the following steps
    # 1. split z into position vector and velocity vector (see total_energy for example)
    r = z[0:2]
    v = z[2:4]
    # 2. compute the norm of position vector, and use it to compute the force
    norm_r = np.linalg.norm(r)
    
    # 3. compute drdt (array [dx/dt, dy/dt])
    drdt = v
    
    # 4. compute dvdt (array [dvx/dt, dvy/dt])
    dvdt = (-m*r) / (norm_r**3)
    
    # join the arrays
    dzdt = np.concatenate((drdt,dvdt))
    
    return(dzdt)

def integrate_orbit(z0, m, tend, h, method='RK4'):
    """
    Integrates orbit starting from an initial position and velocity from t = 0 
    to t = tend.
    
    Arguments:
        z0
            < The initial positions of x and y and the corresponding derivatives >

        m
            < This is the reduced mass of the Earth - Sun system. >
    
        tend
            < The ending time that the model will run to. >
    
        h
            < The stepsize being applied to the positions and derivatives. >
    
        method ('Euler', 'RK2', or 'RK4')
            identifies which stepper routine to use (default: 'RK4')

    Returns
        ts, Xs, Ys, KEs, PEs, TEs := arrays of time, x postions, y positions, 
        and energies (kin., pot., total) 
    """
    
    # set the initial time and phase space array
    t = 0.0

    # expected number of steps
    Nsteps = int(tend/h)+1

    # arrays holding t, x, y, kinetic energy, potential energy, and total energy
    z = np.zeros((Nsteps,4))
    ts = np.zeros(Nsteps)
    Xs = np.zeros(Nsteps)
    Ys = np.zeros(Nsteps)
    KEs = np.zeros(Nsteps)
    PEs = np.zeros(Nsteps)
    TEs = np.zeros(Nsteps)
    
    # store the initial point
    z[0] = z0
    ts[0] = t
    Xs[0] = z[0][0]
    Ys[0] = z[0][1]
    
    # now extend this with KEs[0], PEs[0], TEs[0]
    KEs[0] = kinetic_energy(z0[2:4])
    PEs[0] = potential_energy(z0[0:2], m)
    TEs[0] = KEs[0] + PEs[0]

    # select the stepping method
    advance_one_step = integration_methods[method]
    
    # run through the steps
    for step in range(1,Nsteps):
        z[step] = advance_one_step(derivs,t,z[step-1],h,args=m)
        
        # insert statement here to increment t by the stepsize h
        t += h
        
        # store values
        ts[step] = t
        
        # fill in with assignments for Xs, Ys, KEs, PEs, TEs
        Xs[step] = z[step][0]
        Ys[step] = z[step][1]
        KEs[step] = kinetic_energy(z[step][2:4])
        PEs[step] = potential_energy(z[step][0:2],m)
        TEs[step] = KEs[step] + PEs[step]
        
    return ts, Xs, Ys, KEs, PEs, TEs
    
def set_initial_conditions(a, m, e):
    """
    set the initial conditions for the orbit.  The orientation of the orbit is 
    chosen so that y0 = 0.0 and vx0 = 0.0.
    
    Arguments
        a
            semi-major axis in AU
        m
            total mass in Msun
        e
            eccentricity ( x0 = (1+e)*a )
    
    Returns:
    x0, y0, vx0, vy0, eps0, Tperiod := initial position and velocity, energy, 
        and period
    """
        
    # total energy per unit mass
    eps0 = -m/(2*a)
    
    # period of motion
    Tperiod = (math.pi/2**(1/2))*m*(np.abs(eps0)**(-3/2))

    # initial position
    x0 = (1 + e) * a
    y0 = 0.0

    # initial velocity is in y-direction; we compute it from the energy
    vx0 = 0.0
    vy0 = (2 * eps0 + (2 * m) / x0)**(1/2)
    
    return np.array([x0,y0,vx0,vy0]), eps0, Tperiod
