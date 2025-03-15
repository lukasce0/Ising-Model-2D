


#include <iostream>
#include <random>
#include <string>
#include <map>
#include <cmath>
#include <math.h>
#include <armadillo>
#include <omp.h>
#include <chrono>

#include "Ising.hpp"
#include "utils.hpp"

using namespace std::chrono;





int main() {
	
	


	/*
	Now we will run the simulation for 2x2 lattice and differant temperatures
	we will save the expectation values, to see how they grow with temperature.
	These will be plotted using Python against the theoretical values.
	*/

	int L = 2; //lattice length and widht
	int poins_per_T = (int)1e6; //number of steps per given temperature
	int temp_points = 20; //number of differant temperatures
	int equo_points = 10; //number of steps to reach equolibrium from initial state
	int T_change_equo = 0; //number of steps to reach equolbrium after temperature change
	int save_every = 4; //save values (eps, m, c_V and chi) at every 5th step (corresponds to saving every MC cycle)
	double min_temperature = 1; //minimal temperature in the simulation [J/kb]
	double max_temperature = 5; //maximal temperature in the simulation [J/kb]

	findExpectationValues("expe_val_4b.txt", L, "random", T_change_equo, save_every, temp_points, poins_per_T, min_temperature, max_temperature, equo_points);




	/*
	Now we will study how the expectation values converges to their analytical solutions. We will
	run simulatins at 2 differant temperatures, and save the expectation values after every Monte
	Carlo cycle (N=4 steps). We keep the min and max temperatures. Here it is important that all
	temperatures have same starting point.
	*/
	

	double T; //temperature

	std::string state = "random";
	
	
	L = 2;


	T = 1.;
	temperatureEvolution(1000, "expe_val_4c_1.txt", L, T, state); //1000 MC cycles

	T = 2.4;
	temperatureEvolution(1000, "expe_val_4c_24.txt", L, T, state);








	
	/*
	Now we will look at a simular thing, but with L=20. Here we are interessted in the burn-in
	numbr of Monte Carlo cycles, based on two differant temperatures and initial spin orientation.
	*/
	
	
	L = 20;

	state = "ordered";
	T = 1;
	temperatureEvolution((int) 2e3, "expe_val_5a_1_o.txt", L, T, state); //2e3 MC cycles
		
	T = 2.4;
	temperatureEvolution((int) 2e3, "expe_val_5a_24_o.txt", L, T, state);

	state = "random";
	T = 1;
	temperatureEvolution((int) 2e3, "expe_val_5a_1_r.txt", L, T, state);
	
	T = 2.4;
	temperatureEvolution((int) 2e3, "expe_val_5a_24_r.txt", L, T, state);





	/*
	Here we will be solving the model for L=20 lattice, two differant temperatures, and we want
	to save the individual epsilon measurements and temperature. The purpose is to later (in Python)
	estimate the pdf p(eps;T).
	*/


	state = "random";


	T = 1;
	Ising instance_2020(L, "random", T); //creates an instance
    instance_2020.solve((int) 1e4, (int) 2e3*L*L, L*L);
	std::vector<double> energy_1_vec = instance_2020.energy_pp; //energy measurements

	T = 2.4;
	instance_2020.changeTemperature(T); //changing temperature, no need to create a new instance
	instance_2020.randomConfig();
	instance_2020.clearResults(); //clears the result vectors
	instance_2020.solve((int) 1e4, (int)2e3*L*L, L*L);
	std::vector<double> energy_24_vec = instance_2020.energy_pp; //energy measurements


	//Writing the values to files:


	std::string energy_6_1 = "energy_val_6_1.txt";
	std::string energy_6_24 = "energy_val_6_24.txt";
	std::ofstream oenergy_6_1;
	std::ofstream oenergy_6_24;
	oenergy_6_1.open(energy_6_1);
	oenergy_6_24.open(energy_6_24);


	for (int i = 0; i < energy_24_vec.size(); i++) {
		oenergy_6_1 << energy_1_vec.at(i) << std::endl;
		oenergy_6_24 << energy_24_vec.at(i) << std::endl;
	}

	oenergy_6_1.close();
	oenergy_6_24.close();





	/*
	Here we will run a simulation for multiple temperatures and multiple lattice sizes, to map the
	region around the critical temperature.
	*/

	L = 40; //lattice length and widht
	poins_per_T = (int)1e9; //number of steps per given temperature
	temp_points = 100; //number of differant temperatures
	equo_points = (int)2000*pow(L,2); //number of steps to reach equolibrium from initial state (2000 MC cycles)
	T_change_equo = 10*pow(L, 2); //number of steps to reach equolbrium after temperature change
	save_every = L*L; //save values (eps, m, c_V and chi) at every 5th step (corresponds to saving every MC cycle)
	min_temperature = 2.1; //minimal temperature in the simulation [J/kb]
	max_temperature = 2.4; //maximal temperature in the simulation [J/kb]

	findExpectationValues("expe_val_40_8.txt", L, "random", T_change_equo, save_every, temp_points, poins_per_T, min_temperature, max_temperature, equo_points);

	L = 60;
	equo_points = (int)5 * pow(L, 3);
	T_change_equo = 1e-2 * pow(L, 3);
	save_every = L * L;

	findExpectationValues("expe_val_60_8.txt", L, "random", T_change_equo, save_every, temp_points, poins_per_T, min_temperature, max_temperature, equo_points);


	L = 80;
	equo_points = (int)5 * pow(L, 3);
	T_change_equo = 1e-2 * pow(L, 3);
	save_every = L * L;

	findExpectationValues("expe_val_80_8.txt", L, "random", T_change_equo, save_every, temp_points, poins_per_T, min_temperature, max_temperature, equo_points);

	L = 100;
	equo_points = (int)5 * pow(L, 3);
	T_change_equo = 1e-2 * pow(L, 3);
	save_every = L * L;

	findExpectationValues("expe_val_100_8.txt", L, "random", T_change_equo, save_every, temp_points, poins_per_T, min_temperature, max_temperature, equo_points);


	return 0;
}
