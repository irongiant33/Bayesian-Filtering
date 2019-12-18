# Linear Probabilistic Data Association (PDA) Filter

* getClutteredObservations: the purpose of this file is to inject clutter (observations not originating from the true target) into the set of observations. Do not modify this program.
* getError.m: the purpose of this file is to calculate the error between the true object position and the estimated object position. Do not overwrite this program.
* getObservations.m: the purpose of this file is to inject noise on the true target position to simulate the noisy observations of a sensor. Note that these observations are linearly related to the true target position. Do not modify this program.
* getTrueTrack.m: the purpose of this file is to generate a random object track using the parameters defined in the main file. Do not modify this program.
* pdaf.m: the purpose of this file is to take the cluttered observations as an input and output the estimated state of the target. Assume linear measurement model.
* main.m: the purpose of this file is to define the parameters necessary for simulation, call the necessary functions to simulate, plot the results and display the error. You may partially modify this program.

The only programs that should be modified are the main.m and pdaf.m files. The main.m file should only be modified if the user wants to experiment
with what happens to the error as certain parameters are changed. 
