########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################

import structure as stc
import astro_const as ac
import numpy as np
import pandas as pd
import random as rnd

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

def test_routine(M_test, delta_m, eta, xi):
    '''
    Description:
        runs the integration loop with the parameters 'delta_m', 'eta', 'xi' and a 'Pc' as set 
        by 'M_test' and returns the final values of the integration results
        
    Arguments:
        M_test (float):
            total mass of the white dwarf; units [kg]
        
        delta_m (float):
            "Core" mass; units [kg]
        
        eta (float):
            the integration stops when P < eta * Pc
            
        xi (float):
            the stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
            
    Returns:
        m_r_p_array (array):
            array containing the final values for mass, radius, and pressure from the integration loop 
            as well as a value defining if the integration went over 'max_steps';
            like [mass, radius, pressure, bool]; units [kg, m, Pa, none]
    '''
    
    Pc = stc.pressure_guess(M_test, ac.mue) # make an initial guess of central pressure 'Pc'
    
    m_res, r_res, p_res, max_step_reached = stc.integrate(Pc, delta_m, eta, xi, ac.mue, max_steps=10000, err_max_step=True) # run integration loop
    
    m_r_p_array = np.array([m_res[-1], r_res[-1], p_res[-1], max_step_reached]) # save the final value  from the results of the integration loop
    
    return(m_r_p_array)

num_samples = 1000 # number of randomly sampled parameter values to generate

delta_m = ac.Msun # test mass

eta = 10**(-10) # initial guess for 'eta'
std_eta = np.sqrt(eta) # assumed standard deviation when sampling parameter values from normal dist

xi = 0.1 # initial guess for 'xi'
std_xi = np.sqrt(xi) # assumed standard deviation when sampling parameter values from normal dist

# get randomly sampled parameter values with a mean value of 'xi' and 'eta' with 
    # standard deviations set as the square root of the mean value
# take the absolute value of the returned arrays to ensure parameter space is within
    # the correct bounds
xi_sample_1, xi_sample_2 = np.abs(sample_param_normal(num_samples, xi, std_xi))
eta_sample_1, eta_sample_2 = np.abs(sample_param_normal(num_samples, eta, std_eta))

# initialize arrays of shape ['num_samples', 4]
xi_sample_1_m_r_p = np.zeros((num_samples,4))
xi_sample_2_m_r_p = np.zeros((num_samples,4))
eta_sample_1_m_r_p = np.zeros((num_samples,4))
eta_sample_2_m_r_p = np.zeros((num_samples,4))

# run integrations to get integration results for all randomly generated parameters
for i in range(num_samples):
    xi_sample_1_m_r_p[i] = test_routine(ac.Msun, delta_m, eta, xi_sample_1[i])
    xi_sample_2_m_r_p[i] = test_routine(ac.Msun, delta_m, eta, xi_sample_2[i])
    eta_sample_1_m_r_p[i] = test_routine(ac.Msun, delta_m, eta_sample_1[i], xi)
    eta_sample_2_m_r_p[i] = test_routine(ac.Msun, delta_m, eta_sample_2[i], xi)

# save integration results to csv files
param_analysis_dir = 'ParamAnalysis/'

pd.DataFrame(xi_sample_1).to_csv(param_analysis_dir+'xi_sample_1_indep.csv')
pd.DataFrame(xi_sample_2).to_csv(param_analysis_dir+'xi_sample_2_indep.csv')
pd.DataFrame(eta_sample_1).to_csv(param_analysis_dir+'eta_sample_1_indep.csv')
pd.DataFrame(eta_sample_2).to_csv(param_analysis_dir+'eta_sample_2_indep.csv')

pd.DataFrame(xi_sample_1_m_r_p).to_csv(param_analysis_dir+'xi_sample_1.csv')
pd.DataFrame(xi_sample_2_m_r_p).to_csv(param_analysis_dir+'xi_sample_2.csv')
pd.DataFrame(eta_sample_1_m_r_p).to_csv(param_analysis_dir+'eta_sample_1.csv')
pd.DataFrame(eta_sample_2_m_r_p).to_csv(param_analysis_dir+'eta_sample_2.csv')

#pdf document with final answers and plots showing this
#where to find it in the code