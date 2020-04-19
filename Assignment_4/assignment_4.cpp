#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;

/*
 * ON INPUT:
 * a, b and c -- are the diagonals of the matrix
 * rhs        -- is the right hand side
 * x          -- is the initial guess
 * iterMax    -- is maximum iterations
 * tol        -- is the tolerance level
 * omega      -- is the relaxation parameter
 * sor        -- not used
 * ON OUTPUT:
 * a, b, c, rhs        -- unchanged
 * x                   -- solution to Ax=b
 * iterMax, tol, omega -- unchanged
 * sor                 -- number of iterations to converge
 */
void sorSolve(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c, const std::vector<double> &rhs,
              std::vector<double> &x, int iterMax, double tol, double omega, int &sorCount)
{
  // assumes vectors a,b,c,d,rhs and x are same size (doesn't check)
  int n = a.size() - 1;
  // sor loop
  for (sorCount = 0; sorCount < iterMax; sorCount++)
  {
    double error = 0.;
    // implement sor in here
    {
      double y = (rhs[0] - c[0] * x[1]) / b[0];
      x[0] = x[0] + omega * (y - x[0]);
    }
    for (int j = 1; j < n; j++)
    {
      double y = (rhs[j] - a[j] * x[j - 1] - c[j] * x[j + 1]) / b[j];
      x[j] = x[j] + omega * (y - x[j]);
    }
    {
      double y = (rhs[n] - a[n] * x[n - 1]) / b[n];
      x[n] = x[n] + omega * (y - x[n]);
    }
    // calculate residual norm ||r|| as sum of absolute values
    error += std::fabs(rhs[0] - b[0] * x[0] - c[0] * x[1]);
    for (int j = 1; j < n; j++)
      error += std::fabs(rhs[j] - a[j] * x[j - 1] - b[j] * x[j] - c[j] * x[j + 1]);
    error += std::fabs(rhs[n] - a[n] * x[n - 1] - b[n] * x[n]);
    // make an exit condition when solution found
    if (error < tol)
      break;
  }
}
std::vector<double> thomasSolve(const std::vector<double> &a, const std::vector<double> &b_, const std::vector<double> &c, std::vector<double> &d)
{
  int n = a.size();
  std::vector<double> b(n), temp(n);
  // initial first value of b
  b[0] = b_[0];
  for (int j = 1; j < n; j++)
  {
    b[j] = b_[j] - c[j - 1] * a[j] / b[j - 1];
    d[j] = d[j] - d[j - 1] * a[j] / b[j - 1];
  }
  // calculate solution
  temp[n - 1] = d[n - 1] / b[n - 1];
  for (int j = n - 2; j >= 0; j--)
    temp[j] = (d[j] - c[j] * temp[j + 1]) / b[j];
  return temp;
}
/* Template code for the Crank Nicolson Finite Difference
 */
double crank_nicolson(double S0, double X, double F, double T, double r, double sigma,
                      double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, int &sorCount)
{
  // declare and initialise local variables (ds,dt)
  double dS = S_max / jMax;
  double dt = T / iMax;
  // create storage for the stock price and option price (old and new)
  vector<double> S(jMax + 1), vOld(jMax + 1), vNew(jMax + 1);
  // setup and initialise the stock price
  for (int j = 0; j <= jMax; j++)
  {
    S[j] = j * dS;
  }
  // setup and initialise the final conditions on the option price
  for (int j = 0; j <= jMax; j++)
  {
    vOld[j] = max(F, R * S[j]);
    vNew[j] = max(F, R * S[j]);
  }
  // start looping through time levels
  for (int i = iMax - 1; i >= 0; i--)
  {
    // declare vectors for matrix equations
    vector<double> a(jMax + 1), b(jMax + 1), c(jMax + 1), d(jMax + 1);
    // set up matrix equations a[j]=
    double theta = (1 + mu) * X * exp(mu * i * dt);
    a[0] = 0;
    b[0] = (-1 / dt) - (r / 2) - (kappa * theta / dS);
    c[0] = (kappa * theta / dS);
    d[0] = (-C * exp(-alpha * i * dt)) + (vOld[0] * (-(1 / dt) + (r / 2)));
    for (int j = 1; j <= jMax - 1; j++)
    {
      //
      a[j] = (pow(sigma, 2) * pow(j * dS, 2 * beta) / (4 * pow(dS, 2))) - (kappa * (theta - j * dS) / (4 * dS));
      b[j] = (-1 / dt) - ((pow(sigma, 2.) * pow(j * dS, 2. * beta)) / (2. * pow(dS, 2))) - (r / 2.);
      c[j] = ((pow(sigma, 2.) * pow(j * dS, 2. * beta)) / (4. * pow(dS, 2.))) + ((kappa * (theta - j * dS)) / (4. * dS));
      d[j] = (-vOld[j] / dt) - ((pow(sigma, 2.) * pow(j * dS, 2. * beta) / (4. * pow(dS, 2.))) * (vOld[j + 1] - 2. * vOld[j] + vOld[j - 1])) - (((kappa * (theta - j * dS)) / (4. * dS)) * (vOld[j + 1] - vOld[j - 1])) + ((r / 2.) * vOld[j]) - (C * exp(-alpha * dt * i));
    }
    double A = R * exp((kappa + r) * (i * dt - T));
    double B = -X * A + C * exp(-alpha * i * dt) / (alpha + r) + X * R * exp(r * (i * dt - T)) - C * exp(-(alpha + r) * T + r * i * dt) / (alpha + r);
    a[jMax] = 0;
    b[jMax] = 1;
    c[jMax] = 0;
    d[jMax] = jMax * dS * A + B;
    // solve matrix equations with SOR
    sorSolve(a, b, c, d, vNew, iterMax, tol, omega, sorCount);
    //vNew = thomasSolve(a, b, c, d);

    if (sorCount == iterMax)
      return -1;

    // set old=new
    vOld = vNew;
  }
  // finish looping through time levels

  // output the estimated option price
  double optionValue;

  int jStar = S0 / dS;
  double sum = 0.;
  sum += (S0 - S[jStar]) / (dS)*vNew[jStar + 1];
  sum += (S[jStar + 1] - S0) / (dS)*vNew[jStar];
  optionValue = sum;

  return optionValue;
}

int main()
{
  //
  double T = 2., F = 95., R = 2., r = 0.0229, kappa = 0.125, altSigma = 0.416,
         mu = 0.0213, X = 47.66, C = 1.09, alpha = 0.02, beta = 0.486, sigma = 3.03, tol = 1.e-7, omega = 1., S_max = 10 * X;
  //
  /*
  double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.083333333, altSigma = 0.369,
         mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 1., sigma = 3.73, S_max = 10 * X, tol = 1.e-7, omega = 1.;
  */
  int iterMax = 10000, iMax = 100, jMax = 100;
  //Create graph of varying S0 and beta and bond
  int length = 300;
  double S_range = 3 * X;
  int sor;
  /*
  std::ofstream outFile7("./data/varying_s_beta.csv");

  for (double beta = 0; beta < 1.3; beta += 0.1)
  {
    for (int j = 1; j <= length - 1; j++)
    {
      outFile7 << beta << " , " << altSigma << " , " << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, altSigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor) << "\n";
    }
  }
  outFile7.close();
  */
  /*
  std::ofstream outFile8("./data/varying_s_sigma.csv");
  beta = 1;
  for (double altSigma = 0; altSigma < 3.5; altSigma += 0.1)
  {
    for (int j = 1; j <= length - 1; j++)
    {
      outFile8 << beta << " , " << altSigma << " , " << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, altSigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor) << "\n";
    }
  }
  outFile8.close();
  */
  /*
  std::ofstream outFile9("./data/varying_s_sigma_beta.csv");
  for (double altSigma = 0; altSigma < 3.5; altSigma += 0.1)
  {
    for (double beta = 0; beta < 1.3; beta += 0.1)
    {
      double S0 = X;
      outFile9 << beta << " , " << altSigma << " , " << S0 << " , " << crank_nicolson(S0, X, F, T, r, altSigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor) << "\n";
    }
  }
  outFile9.close();
  */
  /*

  std::ofstream outFile1("./data/varying_s_beta_1.csv");
  std::ofstream outFile2("./data/varying_s_beta_0_4.csv");
  for (int j = 1; j <= length - 1; j++)
  {
    vector<double> gamma(jMax + 1);
    outFile1 << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, altSigma, R, kappa, mu, C, alpha, 1, iMax, jMax, S_max, tol, omega, iterMax, sor, gamma) << "\n";
    outFile2 << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor, gamma) << "\n";
  }
  outFile1.close();
  outFile2.close();
  */
  /*
  std::ofstream outFile3("./data/varying_imax.csv");
  jMax = 100;
  for (iMax = 1; iMax <= 500; iMax += 1)
  {
    double S = X;
    auto t1 = std::chrono::high_resolution_clock::now();
    double result = crank_nicolson(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_taken =
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
            .count();
    outFile3 << S_max << "," << iMax << "," << jMax << "," << S << " , " << std::fixed << result << "," << time_taken << "\n";
  }
  outFile3.close();
  */
  /*
  for (int s_Mult = 10; s_Mult <= 10; s_Mult += 1)
  {
    double S = X;
    S_max = s_Mult * X;
    string title = "./data/smax_jmax/" + to_string(s_Mult) + "_varying_jmax.csv";
    std::ofstream outFile4(title);
    iMax = 40;
    for (jMax = 1; jMax <= 500; jMax += 1)
    {
      auto t1 = std::chrono::high_resolution_clock::now();
      double result = crank_nicolson(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor);
      auto t2 = std::chrono::high_resolution_clock::now();
      auto time_taken =
          std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
              .count();
      outFile4 << S_max << "," << iMax << "," << jMax << "," << S << " , " << std::fixed << result << "," << time_taken << "\n";
    }
    outFile4.close();
  }
  */
  /*
  std::ofstream outFile5("./data/varying_smax.csv");
  iMax = 25;
  tol = 1.e-7;
  for (int s_Mult = 10; s_Mult <= 50; s_Mult += 1)
  {
    jMax = s_Mult * 10;
    double S = X;
    S_max = s_Mult * X;
    int sorCount;
    auto t1 = std::chrono::high_resolution_clock::now();
    double result = crank_nicolson(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sorCount);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_taken =
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
            .count();
    outFile5 << S_max << "," << iMax << "," << jMax << "," << S << " , " << std::fixed << result << "," << time_taken << "\n";
  }
  outFile5.close();
  */
  /*
  S_max = 10 * X;
  std::ofstream outFile6("./data/analytic.csv");
  iMax = 25, jMax = 100;

  for (int j = 1; j <= length - 1; j++)
  {
    vector<double> gamma(jMax + 1);

    outFile6 << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, altSigma, R, 0, mu, C, alpha, 1, iMax, jMax, S_max, tol, omega, iterMax, sor, gamma) << "\n";
  }
  outFile6.close();
  */
  S_max = 10 * X;
  iMax = 50, jMax = 400;
  double S0 = X;
  std::cout << std::fixed << crank_nicolson(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor);
}