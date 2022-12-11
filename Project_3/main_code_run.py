import configparser
import numpy as np
import pandas as pds
import matplotlib.pyplot as plt

import structure as stc
import astro_const as ac
import Compute_L as cl
import zams as zms

confile = configparser.ConfigParser()
confile.read('proj_3_config.ini')

def_delta_m = float(confile['Parameters']['delta_m'])
def_delta_r = float(confile['Parameters']['delta_r'])
def_eta = float(confile['Parameters']['eta'])
def_xi = float(confile['Parameters']['xi'])
def_max_step = int(confile['Parameters']['max_steps'])
def_pp_factor = float(confile['Parameters']['pp_factor'])

Z_array = np.array([int(i) for i in confile['Comparray']['Zs'].split(',')])
A_array = np.array([int(i) for i in confile['Comparray']['As'].split(',')])
XH_array = np.array([float(i) for i in confile['Comparray']['Xs'].split(',')])

comp_array = np.vstack((Z_array,A_array,XH_array))

Mwant = 0.3
RWant = 0.33

m, r, p, l, rho, T, mx = stc.integrate(Mwant, RWant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, 1000)

#print('m', m[-1] / ac.Msun)
#print('r', r[-1] / ac.Rsun)
#print('p', p[-1])
#print('l', l[-1] / ac.Lsun)
#print('rho', rho[-1])
#print('T', T[-1])
#print('mx', mx)

#test1 = cl.Compute_L(Mwant, RWant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, 1000)
#print(test1)

teff = zms.Teff(Mwant)

surf_lum = zms.surface_luminosity(teff,r)

print(l-surf_lum)



