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

# Mwant = 0.4
# Rwant = 0.33

# m, r, p, l, rho, T, mx = stc.integrate(Mwant, Rwant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, 1000)

#print('m', m[-1] / ac.Msun)
#print('r', r[-1] / ac.Rsun)
#print('p', p[-1])
#print('l', l[-1] / ac.Lsun)
#print('rho', rho[-1])
#print('T', T[-1])
#print('mx', mx)

# test1 = cl.Compute_L(Mwant, Rwant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, 10000)
# print(test1)

test_masses = np.linspace(0.1, 0.3, 20)
Rwant = 0.33

main_sequence_radii = np.zeros_like(test_masses)
main_sequence_Teff = np.zeros_like(test_masses)
main_sequence_Lsurf = np.zeros_like(test_masses)

for i in range(len(test_masses)):
    main_sequence_radii[i] = cl.Compute_L(test_masses[i], Rwant, def_delta_m, def_delta_r, def_eta, def_xi, comp_array, def_pp_factor, 10000)[0]
    main_sequence_Teff[i] = zms.Teff(test_masses[i])
    main_sequence_Lsurf[i] = zms.surface_luminosity(main_sequence_Teff[i], main_sequence_radii[i])

print('R', main_sequence_radii)
print('Teff', main_sequence_Teff)
print('Lsurf', main_sequence_Lsurf)

fig, ax1 = plt.subplots()

ax1.scatter(main_sequence_Teff, main_sequence_Lsurf)

ax1.set_xlim(np.max(main_sequence_Teff)*1.2, np.min(main_sequence_Teff)*0.8)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_xlabel('Teff')
ax1.set_ylabel('Lsurf')

ax1.grid()

# ax2.scatter(main_sequence_Teff, main_sequence_Lsurf)

# ax2.set_xlim(np.max(main_sequence_Teff)*1.2, np.min(main_sequence_Teff)*0.8)

# ax2.set_xscale('log')
# ax2.set_yscale('log')

# ax2.set_xlabel('Teff')
# ax2.set_ylabel('Lsurf')

# ax2.grid()

# fig.savefig('FinalPlots/Mass_Radius_Plot.svg')
# plt.clf()