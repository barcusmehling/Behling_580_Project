# Behling_580_Project
Download zip file and extract. Set current directory in MATLAB to the resulting folder. Run "main.m", which 
loads in Env and FRF files in folder. This script produces the three figures on the poster (also included as .fig files here).

A flight vibration environment was simulated by applying 14 broadband forces to a finite element model and calculating the 
resulting accelerations as a function of frequency at various locations. This environment represents the true
environment that a component / vehicle experienced in flight. In practice, one might like to find a set of modeled
forces that accurately represents this environment. To this end, a different set of forces was used to
calculate a frequency response function matrix relating the model forces to the response. The goal of
this project is to understand how the number and placement of accelerometers and forces affects the accuracy of the
resulting model (the model being the calculated forces). Various locations of forces and "training data" 
accelerometers were used to calculate the model forces and estimate the error at "test data" accelerometers. 
Leave-one-out cross validation was also compared to this approach.

See the poster for more background and interpretation of results.
