import configparser
import numpy as np
import pandas as pds
import matplotlib.pyplot as plt
import matplotlib

import structure as stc
import astro_const as ac
import Compute_L as cl
import zams as zms
import eos

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

mu = eos.mean_molecular_weight(Z_array, A_array, XH_array)

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
main_sequence_PC = np.zeros_like(test_masses)
main_sequence_RhoC = np.zeros_like(test_masses)
main_sequence_TC = np.zeros_like(test_masses)

for i in range(len(test_masses)):
    main_sequence_radii[i] = cl.Compute_L(test_masses[i], Rwant, def_delta_m, def_delta_r, 
                                          def_eta, def_xi, comp_array, def_pp_factor, 10000)
    main_sequence_Teff[i] = zms.Teff(test_masses[i])
    main_sequence_Lsurf[i] = zms.surface_luminosity(main_sequence_Teff[i], main_sequence_radii[i]*ac.Rsun)/ac.Lsun
    main_sequence_PC[i], main_sequence_RhoC[i], main_sequence_TC[i] = stc.central_thermal(test_masses[i], main_sequence_radii[i], mu)

'''
Plotting
'''
fig, (ax1, ax2) = plt.subplots(2,1)
fig.set_size_inches(6,10)

log_10_Teff = np.log10(main_sequence_Teff)
log_10_Lsurf = np.log10(main_sequence_Lsurf)

ax1.scatter(log_10_Teff, log_10_Lsurf, c=test_masses)

mass_norm = matplotlib.colors.Normalize(vmin=np.min(test_masses), vmax=np.max(test_masses))
mass_sm = matplotlib.cm.ScalarMappable(mass_norm, cmap='viridis')

plt.colorbar(mass_sm, ax=ax1, label='Mass [$M_{\odot}$]')

ax1.set_title("HR Diagram of Low Mass Main Sequence Stars")
ax1.set_xlabel("$log_{10}(T_{eff}/K)$")
ax1.set_ylabel("$log_{10}(L/L_{\odot})$")

ax1.grid()

ax1.invert_xaxis()


log_10_RhoC = np.log10(main_sequence_RhoC * 0.001) #unit conversion
log_10_TC = np.log10(main_sequence_TC)

ax2.scatter(log_10_RhoC, log_10_TC, c=test_masses)

mass_norm = matplotlib.colors.Normalize(vmin=np.min(test_masses), vmax=np.max(test_masses))
mass_sm = matplotlib.cm.ScalarMappable(mass_norm, cmap='viridis')

plt.colorbar(mass_sm, ax=ax2, label='Mass [$M_{\odot}$]')

ax2.set_title("LogLog of Central Temperature vs Central Density")
ax2.set_xlabel("$log_{10}(Rho/g cm^{-3}))$")
ax2.set_ylabel("$log_{10}(T_{c}/K)$")

ax2.grid()

plt.tight_layout()

fig.savefig('FinalPlots/Part_10_plots.svg')
plt.clf()