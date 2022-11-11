########################################################################
# Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
# AST304, Fall 2022
# Michigan State University
########################################################################

import astro_const as ac
import matplotlib.pyplot as plt
import pandas as pd

param_analysis_dir = ''
final_plot_dir = '../FinalPlots/'

xi_sample_1_indep = pd.read_csv(param_analysis_dir+'xi_sample_1_indep.csv')
xi_sample_2_indep = pd.read_csv(param_analysis_dir+'xi_sample_2_indep.csv')
eta_sample_1_indep = pd.read_csv(param_analysis_dir+'eta_sample_1_indep.csv')
eta_sample_2_indep = pd.read_csv(param_analysis_dir+'eta_sample_2_indep.csv')
delta_m_sample_1_indep = pd.read_csv(param_analysis_dir+'delta_m_sample_1_indep.csv')
delta_m_sample_2_indep = pd.read_csv(param_analysis_dir+'delta_m_sample_2_indep.csv')

xi_sample_1 = pd.read_csv(param_analysis_dir+'xi_sample_1.csv')
xi_sample_2 = pd.read_csv(param_analysis_dir+'xi_sample_2.csv')
eta_sample_1 = pd.read_csv(param_analysis_dir+'eta_sample_1.csv')
eta_sample_2 = pd.read_csv(param_analysis_dir+'eta_sample_2.csv')
delta_m_sample_1 = pd.read_csv(param_analysis_dir+'delta_m_sample_1.csv')
delta_m_sample_2 = pd.read_csv(param_analysis_dir+'delta_m_sample_2.csv')


fig, [[ax1, ax2, ax3], [ax4, ax5, ax6]] = plt.subplots(2,3)
fig.set_size_inches(15,10)

plt.rcParams['lines.markersize'] = 5

# xi parameter analysis
ax1.scatter(xi_sample_1_indep.iloc[:,1], xi_sample_1.iloc[:,1]/ac.Msun, c='b', label='Run 1')
ax1.scatter(xi_sample_2_indep.iloc[:,1], xi_sample_2.iloc[:,1]/ac.Msun, c='r', label='Run 2')

ax1.set_ylim(2.975, 3.1)

ax1.set_xlabel('xi')
ax1.set_ylabel('Final mass [Solar Masses]')
ax1.set_title('Final mass as a function of xi', pad=20)
ax1.legend()

ax4.scatter(xi_sample_1_indep.iloc[:,1], xi_sample_1.iloc[:,2], c='b', label='Run 1')
ax4.scatter(xi_sample_2_indep.iloc[:,1], xi_sample_2.iloc[:,2], c='r', label='Run 2')

ax4.set_ylim(5E6, 6.5E6)

ax4.set_xlabel('xi')
ax4.set_ylabel('Final radius [m]')
ax4.set_title('Final radius as a function of xi', pad=20)
ax4.legend()


# eta parameter analysis
ax2.scatter(eta_sample_1_indep.iloc[:,1], eta_sample_1.iloc[:,1]/ac.Msun, c='b', label='Run 1')
ax2.scatter(eta_sample_2_indep.iloc[:,1], eta_sample_2.iloc[:,1]/ac.Msun, c='r', label='Run 2')

ax2.set_xlabel('eta')
ax2.set_ylabel('Final mass [Solar Masses]')
ax2.set_title('Final mass as a function of eta', pad=20)
ax2.legend()

ax5.scatter(eta_sample_1_indep.iloc[:,1], eta_sample_1.iloc[:,2], c='b', label='Run 1')
ax5.scatter(eta_sample_2_indep.iloc[:,1], eta_sample_2.iloc[:,2], c='r', label='Run 2')

ax5.set_xlabel('eta')
ax5.set_ylabel('Final radius [m]')
ax5.set_title('Final radius as a function of eta', pad=20)
ax5.legend()


# delta_m parameter analysis
ax3.scatter(delta_m_sample_1_indep.iloc[:,1], delta_m_sample_1.iloc[:,1]/ac.Msun, c='b', label='Run 1')
ax3.scatter(delta_m_sample_2_indep.iloc[:,1], delta_m_sample_2.iloc[:,1]/ac.Msun, c='r', label='Run 2')

ax3.set_xlabel('delta_m')
ax3.set_ylabel('Final mass [Solar Masses]')
ax3.set_title('Final mass as a function of delta_m', pad=20)
ax3.legend()

ax6.scatter(delta_m_sample_1_indep.iloc[:,1], delta_m_sample_1.iloc[:,2], c='b', label='Run 1')
ax6.scatter(delta_m_sample_2_indep.iloc[:,1], delta_m_sample_2.iloc[:,2], c='r', label='Run 2')

ax6.set_xlabel('delta_m')
ax6.set_ylabel('Final radius [m]')
ax6.set_title('Final radius as a function of delta_m', pad=20)
ax6.legend()


fig.tight_layout()

fig.savefig(final_plot_dir+'ParameterAnalysis.svg')