import configparser
import numpy as np
import pandas as pds
import matplotlib.pyplot as plt

import structure as stc
import astro_const as ac

confile = configparser.ConfigParser()
confile.read('proj_3_config.ini')

def_delta_m = float(confile['Parameters']['delta_m'])
def_delta_r = float(confile['Parameters']['delta_r'])
def_eta = float(confile['Parameters']['eta'])
def_xi = float(confile['Parameters']['xi'])
def_max_step = int(confile['Parameters']['max_steps'])

Z_array = np.array([1,2,7])
A_array = np.array([1,2,7])
XH_array = np.array([0.706, 0.275, 0.019])

comp_array = np.vstack((Z_array,A_array,XH_array))
pp_factor = 1

Mwant = 0.3
RWant = 0.33
def_delta_m = 1E-12
def_delta_r = 1E-12

m, r, p, l, rho, T, mx = stc.integrate(Mwant, RWant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, pp_factor, 1000)

print(m[-1] / ac.Msun)
print(r[-1] / ac.Rsun)
print(p[-1])
print(l[-1] / ac.Lsun)
print(rho[-1])
print(T[-1])
print(mx)

# plt.scatter(np.arange(mx),m/ac.Msun)
# plt.scatter(np.arange(mx),r/ac.Msun)
