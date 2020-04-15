#include "black_scholes.h"

using namespace std;

int main(int argc, char** argv) {

	//Parameter list
	double S = 100.0;	// Asset price
	double K = 100.0;	// Strike price
	double r = 0.05;	// Risk-free rate
	double v = 0.2;		// Volatility 
	double T = 1.0;		// Time to maturity in years

	// Calculating call/put values and the Greeks
	double call = call_price(S, K, r, v, T);
	double call_delta_v = call_delta(S, K, r, v, T);
	double call_gamma_v = call_gamma(S, K, r, v, T);
	double call_vega_v = call_vega(S, K, r, v, T);
	double call_theta_v = call_theta(S, K, r, v, T);
	double call_rho_v = call_rho(S, K, r, v, T);

	double put = put_price(S, K, r, v, T);
	double put_delta_v = put_delta(S, K, r, v, T);
	double put_gamma_v = put_gamma(S, K, r, v, T);
	double put_vega_v = put_vega(S, K, r, v, T);
	double put_theta_v = put_theta(S, K, r, v, T);
	double put_rho_v = put_rho(S, K, r, v, T);

	// Outputting parameters and prices
	cout << "Asset price " << S << endl;
	cout << "Strike price " << K << endl;
	cout << "Risk-free rate: " << 100 * r << "%" << endl;
	cout << "Volatility: " << 100 * v << "%" << endl;
	cout << "Time to maturity: " << T << endl;

	cout << "Call price: " << call << endl;
	cout << "Call delta: " << call_delta_v << endl;
	cout << "Call gamma: " << call_gamma_v << endl;
	cout << "Call vega: " << call_vega_v << endl;
	cout << "Call theta: " << call_theta_v << endl;
	cout << "Call rho " << call_rho_v << endl;

	cout << "Put price: " << put << endl;
	cout << "Put delta: " << put_delta_v << endl;
	cout << "Put gamma: " << put_gamma_v << endl;
	cout << "Put vega: " << put_vega_v << endl;
	cout << "Put theta: " << put_theta_v << endl;
	cout << "Put rho " << put_rho_v << endl;

	return 0;
}