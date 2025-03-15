

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








void temperatureEvolution(int mc_cycles, std::string filename, int L, double temperature, std::string orientation) {

	std::string file = filename; //"expe_val_4c_1.txt"
	std::ofstream ofile;
	ofile.open(file);

	Ising expe_evolution_instance(L, orientation, temperature);
	std::vector<double> expe_val;

	//This will slowly extend the results vectors:
	for (int j = 0; j < mc_cycles; j++) {
		expe_evolution_instance.solve(L * L, 0, 1); //mente carlo cycle, save every step,


		expe_val = expe_evolution_instance.findExpectation(); //findes expectation values (unefficent as we recalculate the first exp vals, also we compute 6, save 4 although we only need 1)


		for (int value = 0; value < 4; value++) {
			ofile << expe_val.at(value) << "   "; //writes expectation values into file
		}
		ofile << std::endl;

	}
	ofile.close();
}








void findExpectationValues(std::string filename, int L, std::string state, int burn_in, int save_every, int temp_points, int steps_per_T, double min_temperature, double max_temperature, int equolibrium_steps) {

	std::vector<std::vector<double>> expe_val(temp_points);

	//Here we parallelize the code for differant temperatures.
	//Uses the matrix from the end of one simulation as a beginning og the next

	Ising instance(L, state, min_temperature); //creates a lattice
	instance.solve(0, equolibrium_steps, 1); //get a equolibrium for fisrt temperature

	auto start = std::chrono::high_resolution_clock::now();

	//#pragma opm parallel for num_threads(8)
	for (int i = 0; i < temp_points; i++) {


		std::vector<double> vector_element;
		double temperature_i = i * (max_temperature - min_temperature) / (double)(temp_points - 1) + min_temperature; //provides temperature T_i

		instance.changeTemperature(temperature_i);
		instance.clearResults(); //clears the result vectors
		instance.solve(steps_per_T, burn_in, save_every); //steps_per_T steps, skip saving first burn_in steps (to get equolibrium at new temoerature), afterwards save values every save_every step,


		vector_element = instance.findExpectation(); //findes expectation values
		vector_element.push_back(temperature_i);

		expe_val[i] = vector_element;

	}

	auto stop = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	std::cout << duration.count() << std::endl;


		std::string file = filename;
	std::ofstream ofile;
	ofile.open(file);

	for (int i = 0; i < expe_val.size(); i++) {
		ofile << expe_val[i][4] << "   "; //writes temperature to file

		for (int value = 0; value < 4; value++) {
			ofile << expe_val[i][value] << "   "; //writes expectation values into file
		}
		ofile << std::endl;
	}
	ofile.close();
	

}






