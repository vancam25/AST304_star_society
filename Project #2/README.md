White Dwarf
===========

Project directory for a model of white dwarf. 

Using multiple functions and routines, we model the pressure, mass, and radius of white dwarf stars. Additionally, we compare this data to that of Joyce et al. (2018). To get acceptable results we create and test functions/routines that model the equation of state for the white dwarf, and finds converging values for parameters δ, η, and ξ. Our results are stored in a few directories for organization and ease of access, and these results are used/discussed in our WriteUp.

Contents
--------

0. `README.md`: This file.
1. `astro_const.py`: A module containing physical constants, uses `astropy`.
2. `eos.py`: The starter code for the equation of state.
3. `test_eos.py`: Unit test for the equation of state; do not modify this file.
4. `eos_table.txt`: The comparison data for the equation of state, used by `test_eos.py`.
5. `ode.py`: A module of integration funtions, including fEuler, rk2, and rk4. 
6. `structure.py`: Starter code to integrate stellar structure equations.
7. `observations.py`: A module containing class for reading and storing tabulated
    white dwarf masses, radii, and associated uncertainties. 
8. `Joyce.txt`: Table 4, Joyce et al. (2018). Data provided for `observations.py`.
9. `Provencal.txt`: A table of masses and radii and their error for multiple white dwarf stars, which was provided to us for comparison.
10. `additional_funcs.py`: A module of additional functions that creates the shooting function and routine for use in the `mass_radius_analysis.py` file.
11. `mass_radius_analysis.py`: A module for calculating our mass_radius table for analysis. 
12. `mass_radius_relation_table.csv`: A csv of the mass_radius relation saved as a table for observation/comparison. 
13. `comp_obs_plotting.py`: Creates and saves the figure FinalPlots/Mass_Radius_Plot.svg.
14. `proj_2_config.ini`: Final parameters of xi, eta, delta_m, and max_steps to be used in `mass_radius_analysis.py`. 
15. `Project_2_WriteUp.md`: A markdown document of our final write up and results from the project.
16. `FinalPlots`: A directory that holds our final plots for the project.
17. `ParamAnalysis`: A directory filled with csv files as a result of testing for our δ, η, and ξ convergences. 
18. `rubric.md`: The rubric provided to us. 
19. `GroupReadMe.txt`: A text file with the breakdown of project duties, meant for us to designate work. 

Report/closeout
---------------

1. As mentioned in the contents, file `Project_2_WriteUp.md` is where you can find our report and writeup of the project. 

2. How to run our code: 

    Our code is first run through the file `mass_radius_analysis.py`. This generates the mass_radius table and where the shooting function (which is a function routine that to choose a more adept Pc) occurs. These functions are called from our file `additional_funcs.py` which describes and creates the shooting function and routine to compute our mass values. Additionally, our final graphical results in the file `comp_obs_plotting.py` are generated from the table created in `mass_radius_analysis.py`. 
    
    Secondly, the file `comps_obs_plotting.py` is where our comparison happens. There we compare our values to the observations plot. 
    
    Third, the folder `ParamAnalysis` holds all our parameter analysis, as the name suggests. 

Resources
---------

S. R. G. Joyce, M. A. Barstow, S. L. Casewell, M. R. Burleigh, J. B. Holberg,
and H. E. Bond. Testing the white dwarf mass-radius relation and compar-
ing optical and far-UV spectroscopic results with Gaia DR2, HST, and FUSE.
MNRAS, 479:1612–1626, September 2018. doi: 10.1093/mnras/sty1425.
