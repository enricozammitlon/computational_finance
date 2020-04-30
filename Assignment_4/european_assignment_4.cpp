#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
using namespace std;

/* Code for the Crank Nicolson Finite Difference
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
    int sor;
    for (sor = 0; sor < iterMax; sor++)
    {
      double error = 0.;
      // implement sor in here
      {
        double y = (d[0] - c[0] * vNew[1]) / b[0];
        y = vNew[0] + omega * (y - vNew[0]);
        error += (y - vNew[0]) * (y - vNew[0]);
        vNew[0] = y;
      }
      for (int j = 1; j < jMax; j++)
      {
        double y = (d[j] - a[j] * vNew[j - 1] - c[j] * vNew[j + 1]) / b[j];
        y = vNew[j] + omega * (y - vNew[j]);
        error += (y - vNew[j]) * (y - vNew[j]);
        vNew[j] = y;
      }
      {
        double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
        y = vNew[jMax] + omega * (y - vNew[jMax]);
        error += (y - vNew[jMax]) * (y - vNew[jMax]);
        vNew[jMax] = y;
      }
      // make an exit condition when solution found
      if (error < tol * tol)
      {
        sorCount += sor;
        break;
      }
    }
    if (sor >= iterMax)
    {
      std::cout << " Error NOT converging within required iterations\n";
      std::cout.flush();
      throw;
    }
    vOld = vNew;
  }
  // finish looping through time levels

  // output the estimated option price
  double optionValue;
  //linear interp
  int jStar = S0 / dS;
  double sum = 0.;
  sum += (S0 - S[jStar]) / (dS)*vNew[jStar + 1];
  sum += (S[jStar + 1] - S0) / (dS)*vNew[jStar];
  optionValue = sum;
  // alternatively
  //optionValue = lagrangeInterpolation(vNew, S, S0, 4);
  return optionValue;
}

int main()
{
  // Initial condition
  double T = 2., F = 95., R = 2., r = 0.0229, kappa = 0.125, altSigma = 0.416,
         mu = 0.0213, X = 47.66, C = 1.09, alpha = 0.02, beta = 0.486, sigma = 3.03, tol = 1.e-8, omega = 1.4;

  int iterMax = 100000, iMax = 200, jMax = 200, S_max = 6 * X;
  int length = 300;
  double S_range = 3 * X;
  int sor;

  // Run to obtain 3d graph
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

  // Run to obtain varying configurations of beta, sigma graph
  std::ofstream outFile1("./data/varying_s_beta_1.csv");
  std::ofstream outFile2("./data/varying_s_beta_0_4.csv");
  for (int j = 1; j <= length - 1; j++)
  {
    vector<double> gamma(jMax + 1);
    outFile1 << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, altSigma, R, kappa, mu, C, alpha, 1, iMax, jMax, S_max, tol, omega, iterMax, sor) << "\n";
    outFile2 << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor) << "\n";
  }
  outFile1.close();
  outFile2.close();

  // Run to obtain varying imax graph
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

  // Run to obtain varying smax per varying jmax graph
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

  //Run to obtain graph of varying smax
  std::ofstream outFile5("./data/varying_smax.csv");
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

  double oldResult = 0;
  double oldDiff = 0;
  for (int N = 100; N < 3200; N *= 2)
  {
    jMax = N;
    iMax = N;
    double S = X;
    double S_max = int(N / 30) * X;
    int sorCount{0};
    auto t1 = std::chrono::high_resolution_clock::now();
    double result = crank_nicolson(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sorCount);
    double diff = result - oldResult;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_taken =
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
            .count();
    cout << S_max << "," << iMax << "," << jMax << "," << S << " , " << setprecision(10) << result << "," << time_taken << "," << setprecision(3) << oldDiff / diff << "," << sorCount << "\n";
    oldDiff = diff;
    oldResult = result;
  }

  // Run to obtain graph to check with analytic value
  std::ofstream outFile6("./data/analytic.csv");
  S_max = 6 * X;
  iMax = 200, jMax = 200;
  for (int j = 1; j <= length - 1; j++)
  {
    outFile6 << j * S_range / length << " , " << crank_nicolson(j * S_range / length, X, F, T, r, altSigma, R, 0, mu, C, alpha, 1, iMax, jMax, S_max, tol, omega, iterMax, sor) << "\n";
  }
  outFile6.close();

  // Run to obtain final accurate value
  S_max = 58 * X;
  iMax = 800, jMax = 800;
  double S0 = X;
  auto t1 = std::chrono::high_resolution_clock::now();
  double result = crank_nicolson(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor);
  auto t2 = std::chrono::high_resolution_clock::now();
  auto time_taken =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
          .count();
  std::cout << setprecision(10) << result << "," << time_taken << endl;
}