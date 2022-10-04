########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
Description
This module is the creation of three integration methods- Forward Euler (fEuler), Runge-Kutta 2nd Order (rk2), 
and Runge-Kutta 4th order (rk4). These functions will be used to estimate the slope or path at a time of t+h, 
with h being the step-size from our initial value of t. 

All three of these integration methods are a Taylor expansion of z(t+h). As an example, we know for fEuler that
z(t_h) = z(t) + h*(dz/dt) evaluated at t (omitting higher terms). dz/dt = f(t,z) is a known function,
so that means that z(t+h) approximates to z+h*f(t,z) at time t. The same goes for the other two methods, except 
we evaluate them at higher orders of the Taylor Expansion for more accuracy. 

The module can be used to call the functions defined below. Each function works by entering a function f, 
the time we're looking to integrate over t, the initial conditions z (which may change as it is dependant on t), and
the step size away from our value of t is h. There are additional arguements you may pass if your function needs them, 
and they can be passed last as a tuple of 'args = (arguments here)'. 

Each function will return znew with different accuracies, which reflect the different order of integration that we are using. 

To use these functions, make sure this file is in the same location as your notebook. Import the functions from this file.
To actually use the functions, you can call the function as:
x = fEuler(function, time, initial condition, time step h, args =(additional arguements)). 
rk2 and rk4 can be called and stored as variables in much of the same way. 

"""

def fEuler(f,t,z,h,args=()):
    """
    Description
        fEuler is a function that caluclates the Forward Euler method of taking a step h from initial time t of
        a differential equation dz/dt = f(t+h,z). This equation has been expanded into a Taylor expansion.
        This is a first order method, though will have more errors. This is due to the step size
        being reduced by 2, but the integration error only goes down by 2 as well. 
    
        The function takes a function f that relies on time and conditions z, (the RHS) the time length, 
        our initial conditions z, and any additional arugments. 
    
        You can call the function as x = fEuler(function, time, initial condition, time step h, args =(additional arguements)). 
    
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        t
            time length that we are evaluating over. Typically linspace. 
            
        z
            initial condition of z at time t. 
        
        h
            time step h that we will be looking for the value of. 
    
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
    
    # The following trick allows us to pass additional parameters to f
    # first we make sure that args is of type tuple; if not, we make it into
    # that form
    if not isinstance(args,tuple):
        args = (args,)
    
    # when we call f, we use *args to pass it as a list of parameters.
    # for example, if elsewhere we define f like
    # def f(t,z,x,y):
    #    ...
    # then we would call this routine as
    # znew = fEuler(f,t,z,h,args=(x,y))
    #
    znew = z + h*f(t,z,*args)
    
    return znew


def rk2(f,t,z,h,args=()):
    """
    Description
    rk2 is the Runge-Kutta method of integration for a differential equation at time t+h. 
    This method is more accurate as we take many small steps and are computing to a higher degree
    of h in the Taylor Expansion. 
    
    Using Forward Euler, we take a step to the midpoint between t and t+h and calculate zp, an estimate of the solution
    at t + h/2. By taking this midpoint, we can gleam a better look at the solution for znew at z(t+h) in full. 
    We've reduced the stepsize by a factor or two, so the integration error is expected to go down by a factor of four. 
    
    The function takes a function f that relies on time and conditions z, (the RHS) the time length, 
    our initial conditions z, and any additional arugments. 
    
    You can call the function as x = fEuler(function, time, initial condition, time step h, args =(additional arguements)).
    
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        t
            time length that we are evaluating over. Typically linspace. 
            
        z
            initial condition of z at time t. 
        
        h
            time step h that we will be looking for the value of. 
    
    
        args (tuple, optional)
            additional arguments to pass to f
  
    Returns
        znew = z(t+h)
    """
    # The following trick allows us to pass additional parameters to f
    # first we make sure that args is of type tuple; if not, we make it into
    # that form
    if not isinstance(args,tuple):
        args = (args,)
    
    #Calculate zp, or the slope of the midpoint between t and t+h of z. Calculates the same way as Forward Euler. 
    zp =  z + (h/2)*f(t,z,*args)
    
    #Returning z of full step h calculated from the midpoint. 
    znew = z + h*f(t+(h/2),zp,*args)
    
    return znew

def rk4(f,t,z,h,args=()):
    """
    Description
    rk4 is the integration method of Runge-Kutta 4, which is a fourth order integration. This method 
    created a weighted approximation from different averages and peeks at the slope at the midpoints 
    in order to create the solution. Once again, we are computing at a higher degree of the Taylor expansion. 
    
    The function works as taking four estimates; k1 is the estimate of simply returning the slope at time t of the initial conditions.
    k2 is estimating the slope of the midpoint and constructing the solution based on k1. 
    k3 estimates the slope of the midpoint based on the construction of the solution with k2. 
    And k4 estimates the slope in full time stepp t+h based on our function plus the timestep multiplied by our construction of k3. 
    
    We have reduced our step size by 2**4, so our integration error is expected to go down by a factor of 16. This makes this 
    function the most accurate out of our three functions. 
    
    The function takes a function f that relies on time and conditions z, (the RHS) the time length, 
    our initial conditions z, and any additional arugments. 
    
    You can call the function as x = fEuler(function, time, initial condition, time step h, args =(additional arguements)).
    
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
    
        t
            time length that we are evaluating over. Typically linspace. 
            
        z
            initial condition of z at time t. 
        
        h
            time step h that we will be looking for the value of.
    
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
   
    # The following trick allows us to pass additional parameters to f
    # first we make sure that args is of type tuple; if not, we make it into
    # that form
    if not isinstance(args,tuple):
        args = (args,)
    
    #Calculating k1, or our RHS at the initial conditions. In otherwords, the slope at t. 
    k1 = f(t,z, *args)
    #Calculating k2, or the midpoint from k1. Or, calculating slope from k1 and updating solution.
    k2 = f(t+(h/2),z+((h/2)*k1), *args)
    #Calculating k3; updating solution from created slope at k2. 
    k3 = f(t+(h/2),z+((h/2)*k2), *args)
    #Calulating k4; in otherwords, updating our solution with a slope of k3 created from the two prior midpoints. 
    #Taking the timestep in full now. 
    k4 = f(t+h,z + h*k3, *args)
    
    #Putting everything together. Will return znew. 
    znew = z + (h/6)*(k1+2*k2+2*k3+k4)
    
    return znew
