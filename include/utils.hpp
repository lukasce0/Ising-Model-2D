
#ifndef __utils_hpp__  
#define __utils_hpp__



#include <iostream>
#include <random>
#include <string>
#include <map>
#include <cmath>
#include <math.h>
#include <armadillo>
#include <omp.h>
#include <chrono>







void temperatureEvolution(int mc_cycles, std::string filename, int L, double temperature, std::string orientation);


void findExpectationValues(std::string filename, int L, std::string state, int burn_in, int save_every, int temp_points, int steps_per_T, double min_temperature, double max_temperature, int equolibrium_steps);








#endif