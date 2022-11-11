########################################################################
# Team <your team name>: <names>
# AST304, Fall 2020
# Michigan State University
########################################################################

"""
Reads in table 4 of Joyce et al. (2018) containing masses and radii
of white dwarf stars. You should not need to alter this file. The routines
`gen_latex_table` and particularly `make_observation_plot` may be of use.
"""

from numpy import genfromtxt,array

class WhiteDwarf:
    def __init__(self,source,instrument,R,R_err,M,M_err):
        self.source = source
        self.instrument = instrument
        self.radius = R
        self.radius_error = R_err
        self.mass = M
        self.mass_error = M_err      

from observations import MassRadiusObservations
obs = MassRadiusObservations()
print(obs.masses)
print(obs.radii)
print(obs.radius_errors)

for source, info in obs.sources.items():
    print('{0:20} M =  {1:5.3f} +/- {2:5.3f} Msun'.format(source, info.mass, info.mass_error))


    def __init__(self):
        self._data = genfromtxt(\
            'Joyce.txt',
            delimiter=[16,12,7,7,12,8],
            names=['source','instrument','R','R_err','M','M_err'],
            dtype=['S16','S12','f8','f8','f8','f8'],
            autostrip=True)
        self._observations = {}
        sources = array([str(s,'utf-8')\
            for s in self._data['source']])
        self._instruments = array([str(s,'utf-8')\
            for s in self._data['instrument']])
        for s, i, m, r, m_e, r_e in zip(sources,
            self._instruments,
            self._data['M'], self._data['R'],
            self._data['M_err'], 
            self._data['R_err']):
            name = '-'.join([s,i])
            self._observations[name] = WhiteDwarf(s,i,r,r_e,m,m_e)

    #def sources(self):
        #return self._observations

   
    #def masses(self):
        #return self._data['M']

    
    #def radii(self):
        #return self._data['R']

   
    #def mass_errors(self):
        #return self._data['M_err']
        
    
    #def radius_errors(self):
        #return self._data['R_err']
    
    
    #def instruments(self):
        #return self._instruments

def gen_latex_table(observations):
    """
    Makes a nice latex table of the observations
    
    Arguments
        observations
    Returns
        lines
            an array containing the formatted table
    """
    lines=[ r'\begin{tabular}{llrrrr}' ]
    lines.append(\
    r"source     & instrument & $M$   & $\Delta M$ & $R$ & $\Delta R$ \\")
    lines.append(\
    r"           &     & \multicolumn{2}{c}{$(\Msun)$} & \multicolumn{2}{c}{$(0.01\,\Rsun)$}\\")
    lines.append(r"\hline")

    for obj,wd in observations.sources.items():
        lines.append(r'{0:16s} & {1:4s} & {2:5.3f} & {3:5.3f} & {4:5.3f} & {5:5.3f}\\'.format(\
            wd.source,
            wd.instrument,
            wd.mass,
            wd.mass_error,
            wd.radius,
            wd.radius_error))
    lines.append(r'\end{tabular}')
    return lines

def make_observation_plot(ax,observations):
    """
    Adds the observed data to an axis instance along a legend
    
    Arguments
        ax
            An instance of a matplotlib axes
        observations
            An instance of MassRadiusObservations
    Returns
        ax
            The axes instance so that further items can be added to the plot
    """
    
    ax.set_xlabel(r'$M/M_\odot$')
    ax.set_ylabel(r'$R/(0.01\,R_\odot)$')
    ax.minorticks_on()
    
    # split observations into FUSE and HST observations
    fuse = observations.instruments=='FUSE'
    hst = observations.instruments=='HST'
    ax.errorbar(observations.masses[fuse],observations.radii[fuse],
        yerr=observations.radius_errors[fuse],
        xerr=observations.mass_errors[fuse],
        fmt='ro',markersize=6)
    ax.errorbar(observations.masses[hst],observations.radii[hst],
        yerr=observations.radius_errors[hst],
        xerr=observations.mass_errors[hst],
        fmt='b^',markersize=6)
    # make legend with observation type: seems we need to make two fake points,
    # one for each category
    ax.plot([],[],color='blue',marker='^',markersize=8,
        linestyle='none',label='HST')
    ax.plot([],[],color='red',marker='o',markersize=8,
        linestyle='none',label='FUSE')
    ax.legend(frameon=False,labelcolor='mfc')
    return ax

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    obs = MassRadiusObservations()
    plotfile = 'MR_Joyce.png'
    tabfile = 'Joyce-table4.tex'
    print('making plot {}'.format(plotfile))

    width, aspect = 8, 0.62
    fig = plt.figure(figsize=(width,width*aspect))
    ax = fig.add_subplot(111)
    ax = make_observation_plot(ax,obs)
    fig.savefig(plotfile,bbox_inches='tight')

    print('writing {}'.format(tabfile))
    with open(tabfile,'w') as fancytab:
        for line in gen_latex_table(obs):
            fancytab.write(line+'\n')

        
    
