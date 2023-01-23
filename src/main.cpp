#include <iostream>
#include <random>
#include <omp.h>
#include <fstream>
#include <string>
#include "settings.h"
#include "func_declaration.h"
#include "random_generators.h"

int main()
{
	double spin[n + 2][n + 2];
	double temperatures[nT + 1];
	double magnetizations[nT + 1];
	double susceptibilities[nT + 1];
	double internalEnergies[nT + 1];
	double specificHeat[nT + 1];
	double acceptRatio[nT + 1];
	double magnetizationValue;
	double internalEnergy;
	double normalConst;

	for (int i = 0; i < nT + 1; i++) {
		temperatures[i] = Tfirst + i * Tstep;
		magnetizations[i] = 0;
		susceptibilities[i] = 0;
		internalEnergies[i] = 0;
		specificHeat[i] = 0;
		acceptRatio[i] = 0;
	}

	omp_set_num_threads(N_THREAD);
	#pragma omp parallel for private(spin, magnetizationValue, internalEnergy) shared(temperatures, magnetizations, susceptibilities, internalEnergies, specificHeat, acceptRatio) 
	for (int iT = 0; iT < nT + 1; iT++) {
		double trash = 0;

		for (int k = 0; k < N_INDEP; k++) { //thermalozation
			initialize(spin, 1);
			for (int i = 0; i < N_TSWEEP; i++) {
				trash = sweep(spin, temperatures[iT], internalEnergies[iT]);
			}

			internalEnergy = total_hamiltonian(spin);

			for (int i = 0; i < N_SWEEP; i++) { //simulation
				normalConst = double(double(k) + 1) * N_SWEEP + double(double(i) + 1);
				acceptRatio[iT] += sweep(spin, temperatures[iT], internalEnergy);
				internalEnergies[iT] += internalEnergy;
				magnetizationValue = magnetization(spin);
				magnetizations[iT] += magnetizationValue;
				susceptibilities[iT] += susceptibility(spin, magnetizations[iT] / normalConst, temperatures[iT]);
				specificHeat[iT] += specific_heat(spin, internalEnergies[iT] / normalConst, temperatures[iT]);
			}
		}
		acceptRatio[iT] /= normalization;
		magnetizations[iT] /= normalization;
		internalEnergies[iT] /= normalization;
		susceptibilities[iT] /= normalization;
		specificHeat[iT] /= normalization;

		susceptibilities[iT] = (susceptibilities[iT] - magnetizations[iT] * magnetizations[iT]) / (temperatures[iT] * temperatures[iT]);
		specificHeat[iT] = (specificHeat[iT] - internalEnergies[iT] * internalEnergies[iT])/(temperatures[iT] * temperatures[iT]);
	}

	std::ofstream myfile;
	std::string fileName = "output/data_" + std::to_string(n) /*+ "_finer" */ + ".dat";
	myfile.open(fileName);

	myfile << "#temperature \t magnetization \t susceptibility \t specific heat \t internal energy \t accept ratio";
	for (int i = 0; i < nT + 1; i++)
		myfile << std::to_string(temperatures[i]) + "\t" + std::to_string(magnetizations[i]) + "\t" + std::to_string(susceptibilities[i]) + "\t" + std::to_string(specificHeat[i]) + "\t" + std::to_string(internalEnergies[i]) + "\t" + std::to_string(acceptRatio[i]) + "\n";
	
	myfile.close();
}

double local_hamiltonian(double spin[n + 2][n + 2], int i, int j)
{
	return -double(J) * spin[i][j] * (spin[i - 1][j - 1] + spin[i - 1][j + 1] + spin[i + 1][j - 1] + spin[i + 1][j + 1]);
}

double total_hamiltonian(double spin[n + 2][n + 2])
{
	double totalEnergy = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			totalEnergy += -double(J) * spin[i][j] * (spin[i + 1][j] + spin[i][j + 1]);
		}
	}
	totalEnergy -= -double(J) * spin[0][0] * (spin[1][0] + spin[0][1]);
	return totalEnergy;
}


void change_spin(double spin[n + 2][n + 2], int i, int j)
{
	spin[i][j] = -spin[i][j];
	//rest of this funcion ensures periodic boundary conditions
	if (i == 1)
		spin[n + 1][j] = spin[i][j];
	else if (i == n)
		spin[0][j] = spin[i][j];
	if (j == 1)
		spin[i][n + 1] = spin[i][j];
	else if (j == n)
		spin[i][0] = spin[i][j];
}

double magnetization(double spin[n + 2][n + 2])
{
	double result = 0;
	for (int i = 1; i < n + 1; i++) {
		for (int j = 1; j < n + 1; j++) {
			result += spin[i][j];
		}
	}
	return result;
}

double susceptibility(double spin[n + 2][n + 2], double magnetizationValue, double T)
{
	double suscept = 0;
	for (int i = 1; i < n + 1; i++) {
		for (int j = 1; j < n + 1; j++) {
			suscept += spin[i][j] * spin[i][j];
		}
	}
	return suscept; // we will compute the real susceptibility at the end of simulation
	//return (suscept - magnetizationValue * magnetizationValue) / (T * T);
}

double specific_heat(double spin[n + 2][n + 2], double energy, double T)
{
	double heat = 0;
	double subheat;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			subheat = -double(J) * spin[i][j] * (spin[i + 1][j] + spin[i][j + 1]);
			heat += (subheat * subheat);
		}
	}
	subheat = -double(J) * spin[0][0] * (spin[1][0] + spin[0][1]);
	heat -= (subheat * subheat);

	//return (heat - energy * energy) / (T * T);
	return heat; // we will compute the real specific heat at the end of simulation
}

void initialize(double spin[n + 2][n + 2], double state)
{
	if (state != 0)
	{
		for (int i = 0; i < n + 2; i++) {
			for (int j = 0; j < n + 2; j++) {
				spin[i][j] = state;
			}
		}
	}
	else
	{
		int upDown[2] = { 1, -1 };
		int idx;
		for (int i = 0; i < n + 2; i++) {
			for (int j = 0; j < n + 2; j++) {
				idx = distInt01(rng);
				spin[i][j] = upDown[idx];
			}
		}
	}
}

double sweep(double spin[n + 2][n + 2], double temperature, double& internalEnergy) //makes one sweep through all spins and return number of succesful configuration changes(necessary for acceptance ratio)
{
	double accepted = 0;
	double prevEnergy;
	double newEnergy;
	for (int i = 1; i < n + 1; i++) {
		for (int j = 1; j < n + 1; j++) {
			prevEnergy = local_hamiltonian(spin, i, j);
			change_spin(spin, i, j);
			newEnergy = local_hamiltonian(spin, i, j);
			if (metropolis_condition(prevEnergy, newEnergy, temperature)) {
				internalEnergy += newEnergy - prevEnergy;
				accepted++;
			}
			else
				change_spin(spin, i, j); // change the spin back if we dont accept new configuration
		}
	}
	return accepted;
}


bool metropolis_condition(double prevEnergy, double newEnergy, double temp)
{
	return dist(rng) < std::exp((prevEnergy - newEnergy) / temp);
}

void show_spin(double spin[n + 2][n + 2])
{
	for (int i = 1; i < n + 1; i++) {
		for (int j = 1; j < n + 1; j++) {
			std::cout << spin[i][j];
		}
		std::cout << std::endl;
	}
}