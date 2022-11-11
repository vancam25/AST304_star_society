########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################

import structure as stc
import astro_const as ac
import additional_funcs as af
import numpy as np
import pandas as pd

num_samples = 1000 # number of randomly sampled parameter values to generate

delta_m = 1E-8 # test mass

eta = 10**(-20) # initial guess for 'eta'

xi = 0.1 # initial guess for 'xi'
std_xi = np.sqrt(xi) # assumed standard deviation when sampling parameter values from normal dist

# central pressure guess based on expexted final mass
Pc_guess = stc.pressure_guess(ac.Msun, ac.mue)

# get randomly sampled parameter values with a mean value of 'xi' and 'eta' with 
    # standard deviations set as the square root of the mean value
# take the absolute value of the returned arrays to ensure parameter space is within
    # the correct bounds
xi_sample_1, xi_sample_2 = np.abs(af.sample_param_normal(num_samples, xi, std_xi))

# initialize arrays of shape ['num_samples', 4]
xi_sample_1_m_r_p = np.zeros((num_samples,4))
xi_sample_2_m_r_p = np.zeros((num_samples,4))

# run integrations to get integration results for all randomly generated parameters
for i in range(num_samples):
    xi_sample_1_m_r_p[i] = af.test_param_routine(Pc_guess, delta_m, eta, xi_sample_1[i])
    xi_sample_2_m_r_p[i] = af.test_param_routine(Pc_guess, delta_m, eta, xi_sample_2[i])

# save integration results to csv files
param_analysis_dir = ''

pd.DataFrame(xi_sample_1).to_csv(param_analysis_dir+'xi_sample_1_indep.csv')
pd.DataFrame(xi_sample_2).to_csv(param_analysis_dir+'xi_sample_2_indep.csv')

pd.DataFrame(xi_sample_1_m_r_p).to_csv(param_analysis_dir+'xi_sample_1.csv')
pd.DataFrame(xi_sample_2_m_r_p).to_csv(param_analysis_dir+'xi_sample_2.csv')