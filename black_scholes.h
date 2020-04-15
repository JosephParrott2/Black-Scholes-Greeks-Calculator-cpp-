#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <algorithm>


// ===================
// ANALYTICAL SOLUTION
// ===================

// Standard normal probability density function
double norm_pdf(const double& x) {
	return (1.0 / pow(2 * M_PI, 0.5)) * exp(-0.5 * x * x);
}

// Approximation to the cumulative distribution function
double norm_cdf(const double& x) {
	double k = 1.0 / (1.0 + 0.2316419 * x);
	double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

	if (x >= 0.0) {
		return (1.0 - (1.0 / (pow(2 * M_PI, 0.5)) * exp(-0.5 * x * x) * k_sum));
	}
	else {
		return 1.0 - norm_cdf(-x);
	}
}

// Calculating d_j for analytical solution
double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
	return (log(S / K) + (r + (pow(-1, j - 1)) * 0.5 * v * v) * T) / (v * (pow(T, 0.5)));
}

// Calculating Greeks for call option

// Call Price
double call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {

	return S * norm_cdf(d_j(1, S, K, r, v, T)) - K * exp(-r * T) * norm_cdf(d_j(2, S, K, r, v, T));
}

// Delta
double call_delta(const double S, const double K, const double r, const double v, const double T) {
	return norm_cdf(d_j(1, S, K, r, v, T));
}

// Gamma
double call_gamma(const double S, const double K, const double r, const double v, const double T) {
	return norm_pdf(d_j(1, S, K, r, v, T)) / (S * v * sqrt(T));
}

// Vega
double call_vega(const double S, const double K, const double r, const double v, const double T) {
	return S * norm_pdf(d_j(1, S, K, r, v, T)) * sqrt(T);
}

// Theta
double call_theta(const double S, const double K, const double r, const double v, const double T) {
	return -S * norm_pdf(d_j(1, S, K, r, v, T)) * v / (2 * sqrt(T)) - r * K * exp(-r * T) * norm_cdf(d_j(2, S, K, r, v, T));
}

// Rho
double call_rho(const double S, const double K, const double r, const double v, const double T) {
	return K * T * exp(-r * T) * norm_cdf(d_j(2, S, K, r, v, T));
}

// Calculating Greeks for put option

// Put Price
double put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {

	return -S * norm_cdf(-d_j(1, S, K, r, v, T)) + K * exp(-r * T) * norm_cdf(-d_j(2, S, K, r, v, T));
}

// Delta
double put_delta(const double S, const double K, const double r, const double v, const double T) {
	return norm_cdf(d_j(1, S, K, r, v, T)) - 1;
}

// Gamma
double put_gamma(const double S, const double K, const double r, const double v, const double T) {
	return norm_pdf(d_j(1, S, K, r, v, T)) / (S * v * sqrt(T));
}

// Vega
double put_vega(const double S, const double K, const double r, const double v, const double T) {
	return S * norm_pdf(d_j(1, S, K, r, v, T)) * sqrt(T);
}

// Theta
double put_theta(const double S, const double K, const double r, const double v, const double T) {
	return -S * norm_pdf(d_j(1, S, K, r, v, T)) * v / (2 * sqrt(T)) + r * K * exp(-r * T) * norm_cdf(-d_j(2, S, K, r, v, T));
}

// Rho
double put_rho(const double S, const double K, const double r, const double v, const double T) {
	return -K * T * exp(-r * T) * norm_cdf(-d_j(2, S, K, r, v, T));
}

// ========================
// FINITE DIFFERENCE METHOD
// ========================

// Call Options

// Delta
double call_delta_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_S) {
	return (call_price(S + delta_S, K, r, v, T) - call_price(S, K, r, v, T)) / delta_S;
}

// Gamma
double call_gamma_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_S) {
	return (call_price(S + delta_S, K, r, v, T) + call_price(S - delta_S, K, r, v, T) - 2 * call_price(S, K, r, v, T)) / (delta_S * delta_S);
}

// Vega
double call_vega_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_v) {
	return (call_price(S, K, r, v + delta_v, T) - call_price(S, K, r, v, T)) / delta_v;
}

// Theta
double call_theta_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_T) {
	return (call_price(S, K, r, v, T + delta_T) - call_price(S, K, r, v, T)) / delta_T;
}

// Rho
double call_rho_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_r) {
	return (call_price(S, K, r + delta_r, v, T) - call_price(S, K, r, v, T)) / delta_r;
}

// Put Options

// Delta
double put_delta_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_S) {
	return (put_price(S + delta_S, K, r, v, T) - put_price(S, K, r, v, T)) / delta_S;
}

// Gamma
double put_gamma_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_S) {
	return (put_price(S + delta_S, K, r, v, T) + put_price(S - delta_S, K, r, v, T) - 2 * put_price(S, K, r, v, T)) / (delta_S * delta_S);
}

// Vega
double put_vega_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_v) {
	return (put_price(S, K, r, v + delta_v, T) - put_price(S, K, r, v, T)) / delta_v;
}

// Theta
double put_theta_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_T) {
	return (put_price(S, K, r, v, T + delta_T) - put_price(S, K, r, v, T)) / delta_T;
}

// Rho
double put_rho_fdm(const double S, const double K, const double r, const double v, const double T, const double delta_r) {
	return (put_price(S, K, r + delta_r, v, T) - put_price(S, K, r, v, T)) / delta_r;
}


// ============
// MONTE CARLO
// ============

// Box-muller algorithm to generate guassian random numbers
double gaussian_box_muller() {
	double x = 0.0;
	double y = 0.0;
	double euclid_sq = 0.0;

	// Continue generating two uniform random variables until the square of their euclidean distance is less than one
	do {
		x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
		euclid_sq = x * x + y * y;
	} while (euclid_sq >= 1.0);


	return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}



