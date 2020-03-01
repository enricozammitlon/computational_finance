#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <functional>

using namespace std;

void path_dependant_option()
{

	mt19937 rng;
	normal_distribution<> ND(0, 1.);

	double S0 = 64000, sigma = 0.42, r = 0.03, T = 2, X = 64000, D = 0.02;
	int K = 20;

	int N = 100000;
	double sum = 0.;
	for (int i = 0; i < N; i++)
	{

		double dt = T / K;
		vector<double> stockPath(K + 1);
		// initialise first value
		stockPath[0] = S0;
		for (int k = 1; k <= K; k++)
		{
			double phi = ND(rng);
			stockPath[k] = stockPath[k - 1] * exp((r - D - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * phi);
		}

		// now calculate the value of A
		double A = 0.;
		for (int k = 1; k <= K; k++)
			A += fabs(stockPath[k]);
		A /= K;
		sum += max(A - X, 0.);
	}
	cout << " V = " << sum / N * exp(-r * T) << endl;
}

double european_options_monteCarlo(const double stock0, const double strikePrice, const double D,
								   const double interestRate, const double sigma, const double maturity,
								   const function<double(double, double)> payoff,
								   const int N)
{
	// declare the random number generator
	static mt19937 rng;
	normal_distribution<> ND(0, 1);

	double sum = 0.;
	for (int i = 0; i < N; i++)
	{
		double phi = ND(rng);
		double ST = stock0 * exp((interestRate - D - 0.5 * sigma * sigma) * maturity + phi * sqrt(maturity) * sigma);
		sum += payoff(ST, strikePrice);
	}
	return sum / N * exp(-interestRate * maturity);
}

int main()
{
	//path_dependant_option();
	double S0 = 100., T = 0.5, sigma = 0.41, r = 0.03, D = 0.04, X1 = 70., X2 = 100., N = 100000;

	auto long_call_payoff = [](const double S, const double X) { return max(S - X, 0.); };
	auto short_call_payoff = [](const double S, const double X) { return -max(S - X, 0.); };
	auto binary_call_payoff = [](const double S, const double X) { if(S<=X){return 0.;}else{return 1.;} };

	double pi = 0.;
	pi += european_options_monteCarlo(S0, X1, D, r, sigma, T, short_call_payoff, N);
	pi += european_options_monteCarlo(S0, X2, D, r, sigma, T, long_call_payoff, N);
	pi += 2 * X2 * european_options_monteCarlo(S0, X2, D, r, sigma, T, binary_call_payoff, N);
	pi += european_options_monteCarlo(S0, 0., D, r, sigma, T, long_call_payoff, N);

	cout << pi;
}