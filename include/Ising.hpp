#ifndef __Ising_hpp__
#define __Ising_hpp__


#include <iostream>
#include <random>
#include <string>
#include <map>
#include <armadillo>




class Ising
{
public:
	int L; //lattice length (size along one axis)
	double N; //lattice size N=L^2
	double T; //temperature
	arma::mat M; //lattice (matrix with elements +-1)
	std::vector<int> indices; //due to periodic boundary conditions
	std::map<int, double> transitions; //precalculating Boltzmann's factors for possible energy transitions
	std::vector<double> magnetization_pp; //absolute value of magnetization per particle
	std::vector<double> energy_pp; //energy per particle

	std::mt19937 generator;
	std::uniform_real_distribution<double> double_dis;
	std::uniform_int_distribution<int> int_dis; //an alternative would be to multiply and round the double distribution up





	Ising(int lengdth, std::string init, double temperature);

		

	void changeTemperature(double temperature);



	void sameConfig();



	void randomConfig();



	void clearResults();


	//steps: number of flip proposals (not including burn-in steps), burn_in: number of initial steps (to avoid burn-in)
	//save_point: how many steps between each reading of expectation values 
	void solve(int steps, int burn_in, int save_point);


	void metropolis(int steps);


	void saveMag();


	void saveEnergy();


	int findDE(int i, int j);


	std::vector<double> findExpectation();

};










#endif