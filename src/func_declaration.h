#pragma once

double local_hamiltonian(double spin[n + 2][n + 2], int i, int j);
double total_hamiltonian(double spin[n + 2][n + 2]);
void change_spin(double spin[n + 2][n + 2], int i, int j);
double magnetization(double spin[n + 2][n + 2]);
double susceptibility(double spin[n + 2][n + 2], double magnetizationValue, double T);
double specific_heat(double spin[n + 2][n + 2], double energy, double T);
void initialize(double spin[n + 2][n + 2], double state);
double sweep(double spin[n + 2][n + 2], double temperature, double& internalEnergy);
bool metropolis_condition(double prevEnergy, double newEnergy, double temp);
void show_spin(double spin[n + 2][n + 2]);
