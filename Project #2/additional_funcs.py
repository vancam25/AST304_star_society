########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################

import random as rnd
import numpy as np
import structure as stc
import astro_const as ac

from scipy.optimize import bisect

import configparser

confile = configparser.ConfigParser()
confile.read('proj_2_config.ini')

def_delta_m = float(confile['Parameters']['delta_m'])
def_eta = float(confile['Parameters']['eta'])
def_xi = float(confile['Parameters']['xi'])
def_max_step = int(confile['Parameters']['max_steps'])

def sample_param_normal(sample_size, mu, sigma):
    '''
    Description:
        This function generates an array of of values randomly sampled from a normal
        distribution centered at 'mu' with a standard deviation
        of 'sigma'
        
        Uses the Box-Muller transform method
            https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
        
    Arguments:
        sample_size (int):
            number of samples to generate
            
        mu (float):
            mean value of the normal distribution
            
        sigma (float):
            standard deviation of the normal distribution
            
    Returns:
        z_0 (array):
            array of values randomly sampled from a normal 
            distribution centered at 'mu' with a standard 
            deviation of 'sigma'
            
        z_1 (array):
            array of values randomly sampled from a normal 
            distribution centered at 'mu' with a standard 
            deviation of 'sigma'
    '''
    # initialize arrays of size 'sample_size'
    z_0 = np.zeros(sample_size)
    z_1 = np.zeros(sample_size)
    
    # generate random values from a normal distribution
    for i in range(sample_size):
        u1 = rnd.random()
        u2 = rnd.random()
        
        mag = sigma * np.sqrt(-2.0 * np.log(u1)) # calculate the magnitude of the randomly sampled values
        
        z_0[i] = mag * np.cos(2 * np.pi * u2) + mu
        z_1[i] = mag * np.sin(2 * np.pi * u2) + mu
    
    return(z_0, z_1)

def calc_structure_routine(Pc, delta_m = def_delta_m, eta = def_eta,
                           xi = def_xi, max_steps = def_max_step, 
                           err_max_step = True):
    '''
    Description:
        runs the integration loop with the parameters 'delta_m', 'eta', 'xi' and a 'Pc' as set 
        by 'M_test' and returns the final values of the integration results for mass and radius
        
    Arguments:
        Pc (float):
            central pressure; units [Pa]
            
        delta_m (float):
            "Core" mass; units [kg]
        
        eta (float):
            the integration stops when P < eta * Pc
            
        xi (float):
            the stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
            
    Returns:
        m_r_array (array):
            array containing the final values for mass, and radius from 
            the integration loop; like [mass, radius]; units [kg, m]
            
    '''
    m_res, r_res, p_res, max_step_reached = stc.integrate(Pc, delta_m, eta, xi, ac.mue, max_steps=def_max_step, err_max_step=err_max_step) # run integration loop
    
    m_r_array = np.array([m_res[-1], r_res[-1]]) # save the final value  from the results of the integration loop
    
    return(m_r_array)

def test_param_routine(Pc, delta_m = def_delta_m, eta = def_eta, 
                       xi = def_xi, max_steps = def_max_step, 
                       err_max_step = True):
    '''
    Description:
        runs the integration loop with the parameters 'delta_m', 'eta', 'xi' and a 'Pc' as set 
        by 'M_test' and returns the final values of the integration results
        
    Arguments:
        Pc (float):
            central pressure; units [Pa]
            
        delta_m (float):
            "Core" mass; units [kg]
        
        eta (float):
            the integration stops when P < eta * Pc
            
        xi (float):
            the stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
            
    Returns:
        m_r_p_array (array):
            array containing the final values for mass, radius, and pressure from 
            the integration loop as well as a value defining if the integration 
            went over 'max_steps'; ike [mass, radius, pressure, bool]; units [kg, m, Pa, none]
            
    '''
    m_res, r_res, p_res, max_step_reached = stc.integrate(Pc, delta_m, eta, xi, ac.mue, max_steps=def_max_step, err_max_step=err_max_step) # run integration loop
    
    m_r_p_array = np.array([m_res[-1], r_res[-1], p_res[-1], max_step_reached]) # save the final value  from the results of the integration loop
    
    return(m_r_p_array)

def shooting_func(Pc, m_want, delta_m = def_delta_m, eta = def_eta, 
                       xi = def_xi, max_steps = def_max_step, 
                       err_max_step = True):
    '''
    Description:
        function used as apart of the shooting routine to choose
        a 'better' Pc
        
    Arguments:
        Pc (float):
            central pressure; units [Pa]
            
        m_want (float):
            the desired final mass value: units [kg]
        
        delta_m (float):
            "Core" mass; units [kg]
        
        eta (float):
            the integration stops when P < eta * Pc
            
        xi (float):
            the stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
            
    Returns:
        m_res[-1] - m_want (float):
            the difference in the desired mass 'm_want' and the mass calculated with 'Pc'
    '''
    
    m_res, r_res, p_res, max_step_reached = stc.integrate(Pc, delta_m, eta, xi, ac.mue, max_steps=10000, err_max_step=err_max_step) # run integration loop
    
    return(m_res[-1] - m_want)

def shooting_routine(m_want):
    '''
    Description:
        function that uses the shooting method (as described in the 'instructions-1.pdf') 
        to improve the central pressure guess for a star of mass 'm_want'. As written this 
        assumes 'm_want' is relatively close to the mass of the sun
        
    Arguments:
        m_want (float):
            the desired final mass value; units [kg]
            
    Returns:
        bisect_res (float):
            the guess for 'Pc' found using the shooting method; units [Pa]
    '''
    Pc_m_want = stc.pressure_guess(m_want, ac.mue) # make initial guess
    
    Pc_low_high = [0,0]
    
    shoot_result_1 = shooting_func(Pc_m_want, m_want)
    
    shoot_result_2 = 0
    
    if shoot_result_1 >= 0:
        Pc_low_high[1] = Pc_m_want
        
        div_iter = 1
        
        while shoot_result_2 >= 0:
            shoot_result_2 = shooting_func(Pc_m_want/div_iter, m_want)
            
            div_iter += 1
            
            if div_iter >= 1000:
                raise Exception('Too Many Iterations')
                
        Pc_low_high[0] = Pc_m_want/div_iter
    
    else:
        Pc_low_high[0] = Pc_m_want 
        
        mult_iter = 1
        
        while shoot_result_2 >= 0:
            shoot_result_2 = shooting_func(Pc_m_want*mult_iter, m_want)
            
            mult_iter += 1
            
            if mult_iter >= 1000:
                raise Exception('Too Many Iterations')
        
        Pc_low_high[1] = Pc_m_want*mult_iter
    
    bisect_res = bisect(shooting_func, Pc_low_high[1], Pc_low_high[0], args=(m_want))
        
    return(bisect_res)
