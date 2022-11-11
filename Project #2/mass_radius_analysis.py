########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################

import configparser
import eos
import astro_const as ac
import additional_funcs as af
import numpy as np
import pandas as pds

confile = configparser.ConfigParser()
confile.read('proj_2_config.ini')

def_delta_m = float(confile['Parameters']['delta_m'])
def_eta = float(confile['Parameters']['eta'])
def_xi = float(confile['Parameters']['xi'])
def_max_step = int(confile['Parameters']['max_steps'])

results_array = np.zeros((10,6))
mass_frac = np.arange(0.1,1.1,0.1)

for i in range(len(mass_frac)):
    
    results_array[i,2] = af.shooting_routine(ac.Msun*mass_frac[i])
    
for i in range(len(results_array)):
    integration_res = af.calc_structure_routine(results_array[i,2])
    
    results_array[i,0] = integration_res[0]
    results_array[i,1] = integration_res[1]

results_array[:,3] = results_array[:,2] / (ac.G * results_array[:,0]**2 * 
                                           results_array[:,1]**(-4))

results_array[:,4] = eos.density(results_array[:,2], ac.mue)
results_array[:,5] = results_array[:,4] / (3 * results_array[:,0] / 
                                           (4 * np.pi * results_array[:,1]**3))

results_array[:,0] /= ac.Msun
results_array[:,1] /= ac.Rsun

results_col_names = np.array(['M/M_solar', 'R/R_solar', 'P_c [MKS]',
                              'P_c/(G*M^{2}*R^{-4})', 'rho_c [MKS]',
                              'rho_c/[3*M/(3*pi*R^{3})]'])

results_df = pds.DataFrame(data = results_array, 
                           columns = results_col_names)

results_df.to_csv('mass_radius_relation_table.csv')