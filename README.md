# LatticeKMC
This is the code developed by me, Joshua Molloy, and Luke Pimlott for our Msci final year project at the University of Bath.
The project was to create a lattice KMC simulation that could be used to simulate the water gas shift reaction (WGSR) on a Palladium surface, Pd(100).
The simulation allows the user to vary the temperature, pressure, partial pressure and lattice dimensions for simulating WGSR.
Includes stiffness scaling methods to throttle fast events. 

Included are the python scripts used for the simulation as well as for the analysis of the data from the simulation. 
Lattie_KMC.py is the code containing the simulation itself 
  NOTE - Lattice_KMC.py can take a long time to run, so decrease lattice dimension from 50
Data Analysis 1.1 plots the reaction frequencies and makes text files of the H2 and CO2 gas amounts, to be used in the Further Data Analysis code
