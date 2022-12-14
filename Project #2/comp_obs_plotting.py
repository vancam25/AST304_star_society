########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################

"""
Description: 

This file creates our mass_radius plots. We first import in our MassRadiusObservations and
     our make_observation_plot from our observations file. We then take our values and plot them
     as both a line plot and a scatter plot (our values and our errors are plotted here). 
     We then save the plot for our writeup.

"""

import pandas as pds
import matplotlib.pyplot as plt

from IPython.display import display

from observations import MassRadiusObservations
from observations import make_observation_plot

mass_radius_table = pds.read_csv('mass_radius_relation_table.csv', index_col=(0))

display(mass_radius_table)

obs = MassRadiusObservations()

fig, ax1 = plt.subplots()

ax1 = make_observation_plot(ax1, obs)

ax1.plot(mass_radius_table.iloc[:,0], mass_radius_table.iloc[:,1]*100, c='g', linestyle='--')
ax1.scatter(mass_radius_table.iloc[:,0], mass_radius_table.iloc[:,1]*100, c='g', marker='*',
            s=80, label='Our Results')

ax1.legend()
ax1.grid()

fig.savefig('FinalPlots/Mass_Radius_Plot.svg')
plt.clf()
