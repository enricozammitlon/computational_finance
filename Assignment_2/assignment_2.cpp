#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

double normalDistribution(double x)
{
	return 0.5*erfc(-x / sqrt(2.));
}
double qFunc(double t, double T, double kappa, double theta, double sigma)
{
	return (0.25*sigma*sigma/(kappa*kappa))*pow(1-exp(-kappa*(T-t)),3.);
}
double kSquaredFunc(double t, double T, double kappa, double theta, double sigma)
{
	return 0.5*sigma*sigma/pow(kappa,3.)
		*( 2.*exp(-kappa*(T-t)) - 3.* exp(-2.*kappa*(T-t))
		+ 4.*kappa*(T-t) + 1. );
}
double nFunc(double r, double t, double T, double kappa, double theta, double sigma)
{
	return r*(T-t) + (theta-r)*(1.-exp(-kappa*(T-t)))/(2. * kappa);
}
double mFunc(double r, double t, double T, double kappa, double theta, double sigma)
{
	return r*exp(-kappa*(T-t)) + (1-exp(-kappa*(T-t)))*theta;
}
double vSquaredFunc(double t, double T, double kappa, double theta, double sigma)
{
	return (sigma*sigma/(3.*kappa))*(1-exp(-3.*kappa*(T-t)));
}
double fFunc(double r, double t, double T, double kappa, double theta, double sigma)
{
	return mFunc(r,t,T,kappa,theta,sigma) - 0.5*qFunc(t,T,kappa,theta,sigma);
}
//H function used by normalising z=xr-f/v where xr is N(f,v^2)
double hFunc(double r, double t, double T, double kappa, double theta, double sigma,double Xr)
{
	return (Xr-fFunc(r,t,T,kappa,theta,sigma))/sqrt(vSquaredFunc(t,T,kappa,theta,sigma));
}
double PFunc(double r, double t, double T, double kappa, double theta, double sigma)
{
	return exp((2/3)*kSquaredFunc(t,T,kappa,theta,sigma) - 0.25*nFunc(r,t,T,kappa,theta,sigma));
}
double VFunc(double r, double t, double T, double kappa, double theta, double sigma,double Xr)
{
	return PFunc(r, t, T, kappa, theta, sigma)*normalDistribution(hFunc(r, t, T, kappa, theta, sigma,Xr));
}

int main()
{
	double r0 = 0.052, t = 0., T = 2., kappa = 0.0944, theta = 0.0616, sigma = 0.0317, Xr=0.06;
	
    double rMin = 0., rMax = 0.2;
	double n = 100.;
	double dr = (rMax - rMin) / n;
    cout << VFunc(r0, t, T, kappa, theta, sigma,Xr) << endl;
	ofstream output("./Assignment_2/test.csv");
	for (int i = 0; i <= 100; i++)
	{
		double r = rMin + i*dr;
		output << r << "," << PFunc(r, t, T, kappa, theta, sigma);
		output << "," << VFunc(r, t, T, kappa, theta, sigma,Xr) << endl;
	}
}