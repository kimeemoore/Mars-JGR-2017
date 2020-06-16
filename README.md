# Mars-JGR-2017
This Matlab code is to facilitate the use and interpretation of the Mars crustal magnetic field models included in supplementary information of:

Moore, K.M. & Bloxham, J. (2017). The construction of sparse models of Mars’ crustal magnetic field. JGR Planets. DOI: 10.1002/2016JE005238

Use of these files should cite the above reference. 

Note that these Mars models are given as a physical grid of Br values (using an equal-area tiling of a sphere, included in the paper's supplementary information), as opposed to the standard method of a spherical harmonic model. The field can be upward continued to points of the user's choice using a Green's function. 

File descriptions:
1. test_plotting1.m  --> this file quickly plots the Mars magnetic field models
2. test_plotting_new_locations.m --> this fileshows how to use the coordinates, mesh, and models included in the suppplementary information. It provides an in-depth description of all variables, and shows how to calculate Mars' magnetic field at arbitrary points in space using my models, and the Greens Function file listed below (3).
3. calcGreensFunction_weighted_correct_areas.m --> this function calculates a Green's function to upward continue my Mars magnetic field models from the planet's surface (defined here as 3393.5km) to arbitrary points of the user's choice. 
