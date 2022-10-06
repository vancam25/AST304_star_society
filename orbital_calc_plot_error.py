# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 17:24:52 2022

@author: crazy
"""

import matplotlib.pyplot as plt
import numpy as np

import updated_kepler as kpr

#TO RUN PROGRAM ON YOUR COMPUTER CHANGE THE FOLLOWING VARIABLES
save_plot_directory = 'C:/Users/crazy/MSU_Work/Fall_2022/AST_304/Project1_Plots/'
#TO RUN PROGRAM ON YOUR COMPUTER CHANGE THE ABOVE VARIABLES

def calc_plot_orbit(n, ode_method = 'RK4', save_TF = False, save_directory = '', show_fig = True):
    '''
    description:
        calculate and plot the orbits defined by the initial parameters bellow; 
        'mass_reduced', 'semi_major_axis', 'eccentricity'. Calculate the orbits 
        for 'n' number of orbits.
    
    arguments:
        n: number of orbital periods to integrate over; int, in units [Hz]
        
        ode_method: ode solving method to use; str, {'RK2', 'RK4', 'Euler'}
        
        save_TF: conditional for saving figures to a directory; bool, {True: save} {False: don't save}
        
        save_directory: directory of where to save the resulting figures; str
    
        show_fig: conditional for showing figure inline; bool, {True: show} {False: don't show}
    
    returns:
        energy_percent_diff: % change in energy from the first (analytical) solution to the final time; int
    
    '''
    #Calculating orbital mechanics
    mass_reduced = 1 # reduced mass of the sun-earth system; float, in units [earth mass']
    semi_major_axis = 1 # semi major axis of the orbiting system; float, in units [au] 
    eccentricity = 0.5 # eccentricity of the orbiting system; float, found with eq -> (x0 = (1+e)*a), unitless

    z0, eps0, Tperiod = kpr.set_initial_conditions(semi_major_axis, mass_reduced, eccentricity)

    time_integration = Tperiod*3 #calculate the total integration time
    time_step = (Tperiod*0.1)/(2**n) #calculate the time step for integration
    
    ts, Xs, Ys, KEs, PEs, TEs = kpr.integrate_orbit(z0, mass_reduced, time_integration, time_step, ode_method)

    time_day = ts*(1/6.28)*(365) #convert time in periods to time in days 

    #Calculating the percent difference in the total energy value at the begining of the first
        #period and the end of the third period
    energy_abs_diff = np.abs(TEs[-1] - TEs[0])
    energy_average = np.abs((TEs[-1] + TEs[0]) / 2)
    
    energy_percent_diff = 100 * (energy_abs_diff/energy_average)
    
    #Plotting
    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(8,16))
    
    plt.rcParams['font.size'] = '15'
    
    fig.suptitle('Orbital Calculations when $h_{{0}}={:.6f}$ $method:'.format(time_step) + ode_method + '$')
    
    ax1.plot(Xs,Ys)
    ax1.scatter(Xs[0], Ys[0], s=100, c='Blue', label='Start')
    ax1.scatter(Xs[-1], Ys[-1], s=100, c='Red', label='End')
    
    ax1.set_xlabel('X [au]')
    ax1.set_ylabel('Y [au]')
    ax1.set_title('Orbit of Earth')
    ax1.set_xlim(-1,2)
    ax1.set_ylim(-1.5,1.5)
    ax1.legend()
    
    ax2.plot(time_day,TEs, label = 'Total Energy')
    ax2.plot(time_day,KEs, label = 'Kinetic Energy')
    ax2.plot(time_day,PEs, label = 'Potential Energy')
    for i in range(0,4):
        ax2.axvline(365*i, c='black')
    
    ax2.set_title('Energy of Earth')
    ax2.set_xlabel('Time [day]')
    ax2.set_ylabel('Energy/Mass [$GM_{\odot}/1au$]')
    ax2.legend()

    fig.tight_layout()
    
    file_name = 'Orbit_Calc_M_' + ode_method + '_H0_{:.6f}'.format(time_step) + '.pdf'
    
    if save_TF:
        plt.savefig(save_directory + file_name)
    
    if show_fig:
        plt.show(fig)
        
    plt.close(fig)
    
    return(energy_percent_diff, time_step)


#run orbital calculations for various step sizes, and methods
#save total energy errors & step size to arrays for use in error analysis plotting
method_list = ['RK2', 'RK4', 'Euler']
n_array = np.arange(0,11)

TEs_err = np.zeros((len(method_list), len(n_array)))
h0_array = np.zeros(TEs_err.shape)

for i in range(len(method_list)):
    for j in range(len(n_array)):
        TEs_err[i,j], h0_array[i,j] = calc_plot_orbit(n_array[j], method_list[i], save_TF=True, save_directory=save_plot_directory, show_fig=False)
        
        
#plot error analysis results as a loglog plot
fig, ax1 = plt.subplots(figsize=(8,8))
    
plt.rcParams['font.size'] = '15'

ax1.plot(h0_array[0], TEs_err[0], label='RK2')
ax1.plot(h0_array[1], TEs_err[1], label='RK4')
ax1.plot(h0_array[2], TEs_err[2], label='Euler')

ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_xlabel('$h_0$')
ax1.set_ylabel('$\sigma_{E}$')

ax1.set_title('% Error in Total Energy')
ax1.legend()

plt.savefig(save_plot_directory + 'error_analysis.pdf')
plt.show(fig)
plt.close(fig)