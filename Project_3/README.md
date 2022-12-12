### Project Three: Radius and Luminosity of Low Mass Stars
#### Team Flat Star Society: Hannah Sullivan, Steven Vancamp, Abram Anderson, Sanskriti Verma
#### AST 304, Fall 2022
#### Michigan State University
===========

In this project, we calculate radius and luminosity values for low mass stars, as well as experimenting with the p+p reaction rate. Changing this rate has an impact on the structure of the star, which in change influences the luminosity. In order to explore these tasks, we were able to take the code from Project 2 on white dwarfs and add in needed equations for luminosity and reaction rate. From this, we were able to generate the needed values and plots to answer the questions set out for us. 

Contents
===========

0. `README.md`: This document. 
1. `rubric.md`: The rubric for which we followed along the duration of this project. 
2. `astro_const.py`: A file with astrological constants used in this project. 
3. `zams.py`: A file for finding Teff and Surface_Luminosity. 
4. `tests`: A folder containing all the testing files for ensuring accurate functions. 
5. `Compute_L.py`: A file for computing the luminosity of the stars and finding radius. 
6. `ode.py`: A module of integration funtions, including fEuler, rk2, and rk4.
7. `main_code_run.py`: A file where we ran our code to complete the objectives set out by the instructions. 
8. `structure.py`: A module for the structure of the star. 
9. `reactions_template.py`: A module for the stellar reactions happening inside the star. 
10. `proj_3_config.ini`: A file of the initial parameters used. 
11. `additional_funcs.py`: Additional functions used to help make effiencient code for the project. 
12. `Project3_WriteUp.md`: The final write up of the project. 
13. `FinalPlots`: A folder holding the final plots, .csv files and excel files of comparison data. 


Report/closeout
---------------
1. You can find our write up in `Project3_WriteUp.md`. 
2. How to run our code: 

        Run `main_code_run.py` to get the run for parts 2.10, 2.11, and should also be for part 3. Folder `FinalPlots` holds the graph outputs and also the comparisons for part 3. 


Resources
---------

G. Chabrier, I. Baraffe, F. Allard, and P. Hauschildt. Evolutionary Models for
  Very Low-Mass Stars and Brown Dwarfs with Dusty Atmospheres. ApJ, 542:
  464–472, October 2000. doi: 10.1086/309513.
  
Bill Paxton, Lars Bildsten, Aaron Dotter, Falk Herwig, Pierre Lesaffre, and
  Frank Timmes. Modules for experiments in stellar astrophysics (MESA).
  ApJS, 192:3, January 2011.
