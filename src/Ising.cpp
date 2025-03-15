#include <iostream>
#include <random>
#include <string>
#include <map>
#include <cmath>
#include <math.h>
#include <armadillo>


#include "Ising.hpp"








	Ising::Ising(int lengdth, std::string init, double temperature) {

		L = lengdth;
		N = pow(L, 2);
		T = temperature;

		//As we are dealing with square matrix, we may create one index vector for both x and y
		indices.push_back(L - 1);
		for (int i = 0; i < L; i++) {
			indices.push_back(i);
		}
		indices.push_back(0);


		if (init == "random") { //makes it easy to implement alternative initial configurations
			randomConfig(); //uniform distribution
		}
		else if (init == "ordered") {
			sameConfig(); //creates LxL matrix with homogenous spin.
		}


		transitions[4] = exp(-4. / T);
		transitions[8] = exp(-8. / T);


		std::mt19937 generator(std::clock()); //RNG mt19937 using clock time as seed (which will vary)
		std::uniform_real_distribution<double> double_dis(0.0, 1.0); //for rejecting proposed spin flips
		std::uniform_int_distribution<int> int_dis(0, L - 1); //for choosing random positions in the lattice

	}


	void Ising::changeTemperature(double temperature) {
		T = temperature;
		transitions[4] = exp(-4. / T);
		transitions[8] = exp(-8. / T);
	}



	void Ising::sameConfig() {
		M = arma::mat(L, L).fill(1.); //creates LxL matrix with spin up (as we don't have external B field, there is no differance in spin down).
	}



	void Ising::randomConfig() {
		arma::arma_rng::set_seed(123); //chooses seed for randu() (may be set to random using set_seed_random(), this uses current time)
		M = arma::mat(L, L).randu(); //creates LxL matrix with elements random elements in [0,1] (unifrom distribution)
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {

				if (M(i, j) < 0.5) { //changes the value to spin down
					M(i, j) = -1.;
				}
				else {
					M(i, j) = 1.; //changes the value to spin up
				}
			}
		}
	}



	void Ising::clearResults() {
		energy_pp.clear(); //clears the energy measurements
		magnetization_pp.clear(); //clears the magnetization measurements
	}




	//steps: number of flip proposals (not including burn-in steps), burn_in: number of initial steps (to avoid burn-in)
	//save_point: how many steps between each reading of expectation values 
	void Ising::solve(int steps, int burn_in, int save_point) {

		metropolis(burn_in); //Do burn in steps

		int d_step = steps / save_point; //Steps between each expectation value calculation (integer division, so this may not be exactly what user asked for)

		for (int i = 0; i < d_step; i++) {
			metropolis(save_point);
			saveMag();
			saveEnergy();
		}
	}




	void Ising::metropolis(int steps) {

		int i; //random x value on the lattice
		int j; //random y value on the lattice
		int DE; //Delta energy

		for (int step = 0; step < steps; step++) {
			i = int_dis(generator) % L;
			j = int_dis(generator) % L;
			DE = findDE(i, j);

			if (DE <= 0) {
				M(i, j) = -M(i, j);
			}
			else if (transitions.at(DE) >= double_dis(generator)) {
				M(i, j) = -M(i, j);
			}
		}
	}






	void Ising::saveMag() {
		//Magnetization:
		double mag_tot = 0.; //total magnetization
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {
				mag_tot += M(i, j);
			}
		}
		magnetization_pp.push_back(fabs(mag_tot / N));
	}


	void Ising::saveEnergy() {
		//This will be a bit more difficoult, as we need to avoid double counting:
		double E_tot = 0.; //total energy
		int spin_ij;

		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {
				//Here we consider all conections within the lattice togheter with right and lower boundary.
				//By the periodic conditions we don't need to consider upper and left boundary, as they are 
				//identical to lower and its only half of total boundary energy that belongs to our lattice.
				spin_ij = M(i, j); //Here we dont need the index vector
				E_tot -= spin_ij * M(indices[i + 2], indices[j + 1]);
				E_tot -= spin_ij * M(indices[i + 1], indices[j + 2]);
			}
		}
		energy_pp.push_back(E_tot / N);
	}







	int Ising::findDE(int i, int j) { //finds the change in energy in units [J]
		int DE = 0;
		int spin_ij = M(i, j); //index vector is shifted as first element is L-1

		DE += spin_ij * M(indices[i], j);
		DE += spin_ij * M(indices[i + 2], j);
		DE += spin_ij * M(i, indices[j]);
		DE += spin_ij * M(i, indices[j + 2]);

		DE *= 2; //Delta E is minus twice the original energy, but the energy has a prefactor -

		return DE;
	}



	std::vector<double> Ising::findExpectation() { //Returns expectation value of epsilon, magnetization and estimates chi and thermal capacity
		double E_eps = 0; //<eps>
		double E_eps2 = 0; //<eps^2> which we need for c_V
		double E_m = 0; //<m>
		double E_m2 = 0; //<m^2> which we need for chi
		double chi; //suseptibility per spin
		double c_V; //Heat capacity per spin
		int lengdth = energy_pp.size();
		std::vector<double> return_vec;

		for (int i = 0; i < lengdth; i++) {
			E_eps += energy_pp.at(i);
			E_eps2 += pow(energy_pp.at(i), 2);
			E_m += magnetization_pp.at(i);
			E_m2 += pow(magnetization_pp.at(i), 2);
		}
		E_eps /= (double)lengdth;
		E_eps2 /= (double)lengdth;
		E_m /= (double)lengdth;
		E_m2 /= (double)lengdth;

		return_vec.push_back(E_eps);
		return_vec.push_back(E_m);

		c_V = N * 1 / pow(T, 2) * (E_eps2 - pow(E_eps, 2));
		chi = N * 1 / T * (E_m2 - pow(E_m, 2));

		return_vec.push_back(c_V);
		return_vec.push_back(chi);

		return return_vec;
	}


