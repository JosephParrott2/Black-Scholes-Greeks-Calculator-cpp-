#include "black_scholes.h"
#include <iostream>

using namespace std;

void monte_carlo_call_price(const int num_sims, const double S, const double K, const double r, const double v, const double T, const double delta_S, double& price_Sp, double& price_S, double& price_Sm) {

	// We wish to use the same Gaussian random draws for each path, so need three seperate adjusted stock paths
	double Sp_adjust = (S + delta_S) * exp(T * (r - 0.5 * v * v));
	double S_adjust = S * exp(T * (r - 0.5 * v * v));
	double Sm_adjust = (S - delta_S) * exp(T * (r - 0.5 * v * v));

	// To store all three current prices during simulation
	double Sp_cur = 0.0;
	double S_cur = 0.0;
	double Sm_cur = 0.0;

	// The three seperate payoff sums for final price
	double payoff_sum_p = 0.0;
	double payoff_sum = 0.0;
	double payoff_sum_m = 0.0;

	for (int i = 0; i < num_sims; i++) {

		double gauss_bm = gaussian_box_muller();	// Random gauss draw
		
		// Adjust three stock paths
		double expgauss = exp(sqrt(v * v * T) * gauss_bm);
		Sp_cur = Sp_adjust * expgauss;
		S_cur = S_adjust * expgauss;
		Sm_cur = Sm_adjust * expgauss;

		// Calculate payoff sums
		payoff_sum_p += max(Sp_cur - K, 0.0);
		payoff_sum += max(S_cur - K, 0.0);
		payoff_sum_m += max(Sm_cur - K, 0.0);
	}

	price_Sp = (payoff_sum_p / static_cast<double>(num_sims)) * exp(-r * T);
	price_S = (payoff_sum / static_cast<double>(num_sims)) * exp(-r * T);
	price_Sm = (payoff_sum_m / static_cast<double>(num_sims)) * exp(-r * T);
}

// Delta
double call_delta_mc(const int num_sims, const double S, const double K, const double r, const double v, const double T, const double delta_S) {

	double price_Sp = 0.0;
	double price_S = 0.0;
	double price_Sm = 0.0;

	monte_carlo_call_price(num_sims, S, K, r, v, T, delta_S, price_Sp, price_S, price_Sm);
	return (price_Sp - price_S) / delta_S;
}

// Gamma
double call_gamma_mc(const int num_sims, const double S, const double K, const double r, const double v, const double T, const double delta_S) {

	double price_Sp = 0.0;
	double price_S = 0.0;
	double price_Sm = 0.0;

	monte_carlo_call_price(num_sims, S, K, r, v, T, delta_S, price_Sp, price_S, price_Sm);
	return (price_Sp - 2 * price_S + price_Sm) / (delta_S * delta_S);
}

int main(int argc, char** argv) {

	// Parameter list
	double S = 100.0;				// Asset price
	double delta_S = 0.001;			// Asset price increment
	double K = 100.0;				// Strike price
	double r = 0.05;				// Risk-free rate
	double v = 0.2;					// Volatility
	double T = 1.0;					// Time to expiry in years
	int num_sims = 10000000;		// Number of simulations for Monte Carlo
	
	// Calculating Delta and Gamma
	double call_delta_m = call_delta_mc(num_sims, S, K, r, v, T, delta_S);
	double call_gamma_m = call_gamma_mc(num_sims, S, K, r, v, T, delta_S);

	// Ouputting answers
	cout << "Number of simulations: " << num_sims << endl;
	cout << "Underlying asset price: " << S << endl;
	cout << "Increment size: " << delta_S << endl;
	cout << "Strike price " << K << endl;
	cout << "Risk-free rate: " << r * 100 << "%" << endl;
	cout << "Volatility: " << v * 100 << "%" << endl;
	cout << "Time to expiry " << T << " years" << endl;

	cout << "Call Delta: " << call_delta_m << endl;
	cout << "Call Gamma: " << call_gamma_m << endl;


}




