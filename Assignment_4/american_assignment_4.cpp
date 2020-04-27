#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>

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
              std::vector<double> &x, int iterMax, double tol, double omega, int &sorCount, double rs)
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
      x[0] = std::max(x[0], rs);
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

double lagrangeInterpolation(const vector<double> &y, const vector<double> &x, double x0, unsigned int n)
{
  if (x.size() < n)
    return lagrangeInterpolation(y, x, x0, x.size());
  if (n == 0)
    throw;
  int nHalf = n / 2;
  int jStar;
  double dx = x[1] - x[0];
  if (n % 2 == 0)
    jStar = int((x0 - x[0]) / dx) - (nHalf - 1);
  else
    jStar = int((x0 - x[0]) / dx + 0.5) - (nHalf);
  jStar = std::max(0, jStar);
  jStar = std::min(int(x.size() - n), jStar);
  if (n == 1)
    return y[jStar];
  double temp = 0.;
  for (unsigned int i = jStar; i < jStar + n; i++)
  {
    double int_temp;
    int_temp = y[i];
    for (unsigned int j = jStar; j < jStar + n; j++)
    {
      if (j == i)
      {
        continue;
      }
      int_temp *= (x0 - x[j]) / (x[i] - x[j]);
    }
    temp += int_temp;
  }
  // end of interpolate
  return temp;
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

void psorSolve(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c, const std::vector<double> &rhs,
               std::vector<double> &x, int iterMax, double tol, double omega, int &sorCount)
{
  double P = 100.;
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
      error += std::fabs(rhs[0] - b[0] * x[0] - c[0] * x[1]);
      x[0] = std::max(P, x[0]);
    }
    for (int j = 1; j < n; j++)
    {
      double y = (rhs[j] - a[j] * x[j - 1] - c[j] * x[j + 1]) / b[j];
      x[j] = x[j] + omega * (y - x[j]);
      error += std::fabs(rhs[j] - a[j] * x[j - 1] - b[j] * x[j] - c[j] * x[j + 1]);
      x[j] = std::max(P, x[j]);
    }
    {
      double y = (rhs[n] - a[n] * x[n - 1]) / b[n];
      x[n] = x[n] + omega * (y - x[n]);
      error += std::fabs(rhs[n] - a[n] * x[n - 1] - b[n] * x[n]);
      x[n] = std::max(P, x[n]);
    }

    // make an exit condition when solution found
    if (error < tol)
      break;
  }
}

/* Template code for the Crank Nicolson Finite Difference
 */
double crank_nicolson1(double S0, double X, double F, double T, double r, double sigma,
                       double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, int &sorCount, double t0)
{
  // declare and initialise local variables (ds,dt)
  double P = 100.;
  double dS = S_max / jMax;
  double f = (T - t0) / T;
  double dt = (T - t0) / (iMax * f);
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
  for (int i = iMax; i >= 0; i--)
  {

    if (i * dt < t0)
    {
      dt = t0 / (iMax * (1 - f));
    }

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
    double penalty = 1.e8;
    int q;
    for (q = 0; q < 100000; q++)
    {
      vector<double> bHat(b), dHat(d);
      for (int j = 1; j < jMax; j++)
      {
        if (i * dt < t0)
        {
          if (vNew[j] < max(R * S[j], P))
          {
            bHat[j] = b[j] - penalty;
            dHat[j] = d[j] - penalty * max(R * S[j], P);
          }
        }
        else
        {
          // turn on penalty if V < RS
          if (vNew[j] < R * S[j])
          {
            bHat[j] = b[j] - penalty;
            dHat[j] = d[j] - penalty * R * S[j];
          }
        }
      }
      // solve matrix equations with SOR
      vector<double> y = thomasSolve(a, bHat, c, dHat);
      // calculate difference from last time
      double error = 0.;
      for (int j = 0; j <= jMax; j++)
        error += fabs(vNew[j] - y[j]);
      vNew = y;
      if (error < 1.e-8)
      {
        sorCount += q;
        break;
      }
    }
    if (q == 100000)
    {
      std::cout << " Error NOT converging within required iterations\n";
      std::cout.flush();
      throw;
    }

    // set old=new
    vOld = vNew;
  }
  // finish looping through time levels

  // output the estimated option price
  double optionValue;
  /*
  int jStar = S0 / dS;
  double sum = 0.;
  sum += (S0 - S[jStar]) / (dS)*vNew[jStar + 1];
  sum += (S[jStar + 1] - S0) / dS * vNew[jStar];
  optionValue = sum;
  */
  optionValue = lagrangeInterpolation(vNew, S, S0, 4);

  return optionValue;
}

double crank_nicolson2(double S0, double X, double F, double T, double r, double sigma,
                       double R, double kappa, double mu, double C, double alpha, double beta, int iMax, int jMax, int S_max, double tol, double omega, int iterMax, int &sorCount, double t0)
{
  // declare and initialise local variables (ds,dt)
  double P = 100.;
  double dS = S_max / jMax;
  double f = (T - t0) / T;
  double dt = (T - t0) / (iMax * f);
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
  for (int i = iMax; i >= 0; i--)
  {
    if (i * dt < t0)
    {
      dt = t0 / (iMax * (1 - f));
    }
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
        if (i * dt < t0)
        {
          y = std::max(std::max(y, R * S[0]), P);
        }
        else
        {
          y = std::max(y, R * S[0]);
        }
        error += (y - vNew[0]) * (y - vNew[0]);
        vNew[0] = y;
      }
      for (int j = 1; j < jMax; j++)
      {
        double y = (d[j] - a[j] * vNew[j - 1] - c[j] * vNew[j + 1]) / b[j];
        y = vNew[j] + omega * (y - vNew[j]);
        if (i * dt < t0)
        {
          y = std::max(std::max(y, R * j * dS), P);
        }
        else
        {
          y = std::max(y, R * j * dS);
        }
        error += (y - vNew[j]) * (y - vNew[j]);
        vNew[j] = y;
      }
      {
        double y = (d[jMax] - a[jMax] * vNew[jMax - 1]) / b[jMax];
        y = vNew[jMax] + omega * (y - vNew[jMax]);
        if (i * dt < t0)
        {
          y = std::max(std::max(y, R * jMax * dS), P);
        }
        else
        {
          y = std::max(y, R * jMax * dS);
        }
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

    if (sorCount == iterMax)
      return -1;

    // set old=new
    vOld = vNew;
  }
  // finish looping through time levels

  // output the estimated option price
  double optionValue;
  /*
  int jStar = S0 / dS;
  double sum = 0.;
  sum += (S0 - S[jStar]) / (dS)*vNew[jStar + 1];
  sum += (S[jStar + 1] - S0) / dS * vNew[jStar];
  optionValue = sum;
  */
  optionValue = lagrangeInterpolation(vNew, S, S0, 4);

  return optionValue;
}

int main()
{
  double T = 2., F = 95., R = 2., r = 0.0229, kappa = 0.125, altSigma = 0.416,
         mu = 0.0213, X = 47.66, C = 1.09, alpha = 0.02, beta = 0.486, sigma = 3.03, tol = 1.e-8, omega = 1., S_max = 20 * X;
  int iMax = 600;
  int jMax = 600;
  double t0 = 0.57245;
  //
  /*
  double T = 3., F = 56., R = 1., r = 0.0038, kappa = 0.083333333, altSigma = 0.369,
         mu = 0.0073, X = 56.47, C = 0.106, alpha = 0.01, beta = 0.425, sigma = 3.73, S_max = 10 * X, tol = 1.e-7, omega = 1.;
  */
  int iterMax = 100000;
  //Create graph of varying S0 and beta and bond
  int length = 300;
  double S_range = 3 * X;
  int sor;
  /*
  std::ofstream outFile1("./data/no_put_american_varying_s_beta_0_4.csv");
  std::ofstream outFile2("./data/american_varying_s_beta_0_4.csv");
  for (int j = 1; j <= length - 1; j++)
  {
    std::cout << j << std::endl;
    outFile1 << j * S_range / length << " , " << crank_nicolson1(j * S_range / length, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor, 0.) << "\n";
    outFile2 << j * S_range / length << " , " << crank_nicolson1(j * S_range / length, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor, t0) << "\n";
    outFile1.flush();
    outFile2.flush();
  }
  outFile1.close();
  outFile2.close();
  */
  /*
  std::ofstream outFile4("./data/american_varying_s_kappa_625.csv");
  std::ofstream outFile5("./data/american_varying_s_kappa_125.csv");
  std::ofstream outFile6("./data/american_varying_s_kappa_187.csv");

  for (int j = 1; j <= length - 1; j++)
  {
    std::cout << j << std::endl;
    double result1 = crank_nicolson1(j * S_range / length, X, F, T, r, sigma, R, 0.0625, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor, t0);
    double result2 = crank_nicolson1(j * S_range / length, X, F, T, r, sigma, R, 0.125, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor, t0);
    double result3 = crank_nicolson1(j * S_range / length, X, F, T, r, sigma, R, 0.1875, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor, t0);
    outFile4 << j * S_range / length << " , " << result1 << "\n";
    outFile5 << j * S_range / length << " , " << result2 << "\n";
    outFile6 << j * S_range / length << " , " << result3 << "\n";
    outFile4.flush();
    outFile5.flush();
    outFile6.flush();
  }
  outFile4.close();
  outFile5.close();
  outFile6.flush();
  */
  double S0 = X;
  iMax = 1300;
  jMax = 1300;
  S_max = jMax * X / 20;
  auto t1 = std::chrono::high_resolution_clock::now();
  double result = crank_nicolson1(S0, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sor, t0);
  auto t2 = std::chrono::high_resolution_clock::now();
  auto time_taken =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
          .count();
  cout << fixed << result << "," << time_taken << endl;
  /*
  std::ofstream outFile7("./data/american_varying_smax_penalty.csv");
  double oldResult = 0, oldDiff = 0;
  double S = X;
  iMax = 100;
  jMax = 100;
  for (int n = 100; n <= 10000; n *= 2)
  {
    iMax = n;
    jMax = n;
    S_max = n / 20 * X;
    int sorCount{0};
    //t0 = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    double result = crank_nicolson1(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sorCount, t0);
    double diff = result - oldResult;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_taken =
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
            .count();
    double extrap = (4 * result - oldResult) / 3.;
    outFile7 << S_max << "," << iMax << "," << jMax << "," << S << "," << setprecision(10) << result << "," << time_taken << "," << extrap << "," << setprecision(3) << oldDiff / diff << "," << sorCount << "\n";
    oldDiff = diff;
    oldResult = result;
  }
  outFile7.close();

  std::ofstream outFile8("./data/american_varying_smax_sor.csv");
  oldResult = 0;
  oldDiff = 0;
  S = X;
  iMax = 100;
  jMax = 100;
  for (int n = 100; n <= 10000; n *= 2)
  {
    iMax = n;
    jMax = n;
    S_max = n / 20 * X;
    int sorCount{0};
    //t0 = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    double result = crank_nicolson2(S, X, F, T, r, sigma, R, kappa, mu, C, alpha, beta, iMax, jMax, S_max, tol, omega, iterMax, sorCount, t0);
    double diff = result - oldResult;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time_taken =
        std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
            .count();
    double extrap = (4 * result - oldResult) / 3.;
    outFile8 << S_max << "," << iMax << "," << jMax << "," << S << "," << setprecision(10) << result << "," << time_taken << "," << extrap << "," << setprecision(3) << oldDiff / diff << "," << sorCount << "\n";
    oldDiff = diff;
    oldResult = result;
  }
  outFile8.close();
  */
}