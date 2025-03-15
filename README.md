# Ising-Model-2D

The code in this repository is a simulation of a 2D Ising model, together with related figures and LaTeX code.

## This repository contains the following code: 

### C++:
- Ising.cpp: Is an object that stores four attributes associated with Ising model. The constructor takes in lattice length and initial state. There are some associated functions to for instance, change temperature or solve the Ising model using Marco chain Monte Carlo. This code uses armadillo library and has to be built with the flag -larmadillo.

- utils.cpp: This code contains two functions, both of which are only generalized, such that code is not repeated in main file. This code, too includes armadillo and should be built with an appropriate flag. There is a commented out section of this code that attempts to use OpenMP parallelization, which it seems dosn't work as desired.

- main.cpp: This is a main function running the simulations. It consists of multiple section, each of which has a headline to explain what it does. This code requires to be linked with both Ising.cpp and utils.cpp to function properly. Also, this code uses armadillo and needs an appropriate flag while building. The use of -O3 flag while building this code, has resulted in a speeding up of the code.


### Python:
- plotting.py: This code interprets data collected by the above-mentioned code. Its main purpose is to plot the results, but it is also used to find variance in solutions, determine uncertainty and use linear regression.

### LaTeX:
- main.tex: This is the source code for the report I handed in. Running this code properly requires figures that are stored in the directory LaTeX/figures.
