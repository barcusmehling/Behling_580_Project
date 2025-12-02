# Behling_580_Project
Run "main.m" to produce the three figures (also saved as .mats in the folder)

Loads in Env and FRF files in folder. 

A flight environment was simulated by applying 14 forces to a finite element model and calculating the 
resulting accelerations as a function of frequency at various locations. Many other forces were used to
calculate a frequency response function matrix relating the force to response. Various locations of
forces and "training data" accelerometers were used to calculate the model forces and estimate the error 
at "test data" accelerometers. Leave-one-out cross validation was also compared to this approach.
