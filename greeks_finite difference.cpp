#include "black_scholes.h"

using namespace std;

int main(int argc, char** argv) {

	// Parameters
	double S = 100.0;			// Asset price
	double delta_S = 0.001;		// Asset price increment
	double K = 100.0;			// Strike price
	double r = 0.05;			// Risk-free rate
	double delta_r = 0.00001;	// Risk-free rate increment
	double v = 0.2;				// Volatility
	double delta_v = 0.00001;	// Volatility increment
	double T = 1.0;				// Time to maturity in years
	double delta_T = 0.001;		// Maturity increment

	// Calculating Greeks
	double call_delta_f = call_delta_fdm(S, K, r, v, T, delta_S);
	double call_gamma_f = call_gamma_fdm(S, K, r, v, T, delta_S);
	double call_vega_f = call_vega_fdm(S, K, r, v, T, delta_v);
	double call_theta_f = call_theta_fdm(S, K, r, v, T, delta_T);
	double call_rho_f = call_rho_fdm(S, K, r, v, T, delta_r);

	double put_delta_f = put_delta_fdm(S, K, r, v, T, delta_S);
	double put_gamma_f = put_gamma_fdm(S, K, r, v, T, delta_S);
	double put_vega_f = put_vega_fdm(S, K, r, v, T, delta_v);
	double put_theta_f = put_theta_fdm(S, K, r, v, T, delta_T);
	double put_rho_f = put_rho_fdm(S, K, r, v, T, delta_r);

	// Outputting parameters and prices
	cout << "Asset price " << S << endl;
	cout << "Strike price " << K << endl;
	cout << "Risk-free rate: " << 100 * r << "%" << endl;
	cout << "Volatility: " << 100 * v << "&" << endl;
	cout << "Time to maturity: " << T << endl;

	cout << "Call delta: " << call_delta_f << endl;
	cout << "Call gamma: " << call_gamma_f << endl;
	cout << "Call vega: " << call_vega_f << endl;
	cout << "Call theta: " << call_theta_f << endl;
	cout << "Call rho: " << call_rho_f << endl;

	cout << "Put delta: " << put_delta_f << endl;
	cout << "Put gamma: " << put_gamma_f << endl;
	cout << "Put vega: " << put_vega_f << endl;
	cout << "Put theta: " << put_theta_f << endl;
	cout << "Put rho: " << put_rho_f << endl;

	return 0;
}
