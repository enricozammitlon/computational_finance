#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

double generateGaussianNoise(double mu, double sigma, int i) {
  auto halton_seq = [](int index, int base = 2) {
    double f = 1, r = 0;
    while (index > 0) {
      f = f / base;
      r = r + f * (index % base);
      index = index / base;
    }
    return r;
  };
  static const double two_pi = 2.0 * 3.14159265358979323846;

  double u1, u2;
  u1 = halton_seq(i, 2);
  u2 = halton_seq(i, 3);

  double z1;
  thread_local bool generate;
  generate = !generate;

  if (!generate)
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
  return z1 * sigma + mu;

  double z0;
  z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
  return z0 * sigma + mu;
}

double normalDistribution(double x) { return 0.5 * erfc(-x / sqrt(2.)); }

double european_options_monteCarlo(
    const double stock0, const double strikePrice, const double D,
    const double interestRate, const double sigma, const double maturity,
    const std::function<double(double, double)> payoff, const int N) {
  // declare the random number generator
  static std::mt19937 rng;
  std::normal_distribution<> ND(0, 1);

  double sum = 0.;
  for (int i = 0; i < N; i++) {
    double phi = ND(rng);
    double ST =
        stock0 * exp((interestRate - D - 0.5 * sigma * sigma) * maturity +
                     phi * sqrt(maturity) * sigma);
    sum += payoff(ST, strikePrice);
  }
  return sum / N * exp(-interestRate * maturity);
}

double european_options_halton(
    const double stock0, const double strikePrice, const double D,
    const double interestRate, const double sigma, const double maturity,
    const std::function<double(double, double)> payoff, const int N) {
  double sum = 0.;
  for (int i = 1; i <= N; i++) {
    double phi = generateGaussianNoise(0, 1, i);
    double ST =
        stock0 * exp((interestRate - D - 0.5 * sigma * sigma) * maturity +
                     phi * sqrt(maturity) * sigma);
    sum += payoff(ST, strikePrice);
  }
  return sum / N * exp(-interestRate * maturity);
}

double european_options_monteCarlo_antithetic(
    const double stock0, const double strikePrice, const double D,
    const double interestRate, const double sigma, const double maturity,
    const std::function<double(double, double)> payoff, const int N) {
  // declare the random number generator
  static std::mt19937 rng;
  std::normal_distribution<> ND(0, 1);

  double sum = 0.;
  for (int i = 0; i < int(N / 2); i++) {
    double phi = ND(rng);
    double ST1 =
        stock0 * std::exp((interestRate - D - 0.5 * sigma * sigma) * maturity +
                          phi * std::sqrt(maturity) * sigma);
    double ST2 =
        stock0 * std::exp((interestRate - D - 0.5 * sigma * sigma) * maturity +
                          (-phi) * std::sqrt(maturity) * sigma);
    sum += payoff(ST1, strikePrice);
    sum += payoff(ST2, strikePrice);
  }
  return sum / N * std::exp(-interestRate * maturity);
}

double european_options_monteCarlo_moment_match(
    const double stock0, const double strikePrice, const double D,
    const double interestRate, const double sigma, const double maturity,
    const std::function<double(double, double)> payoff, const int N) {
  // declare the random number generator
  static std::mt19937 rng;
  std::normal_distribution<> ND(0, 1);
  std::vector<double> PI_N(N);
  double sum = 0.0;
  for (int i = 0; i < int(N / 2); i += 1) {
    double phi = ND(rng);
    PI_N[i] = phi;
    sum += 2 * std::pow(phi, 2);
  }
  double sample_variance = std::sqrt(sum / (N - 1));
  sum = 0.0;
  for (int i = 0; i < int(N / 2); i++) {
    double phi = PI_N[i] / sample_variance;
    double ST1 =
        stock0 * exp((interestRate - D - 0.5 * sigma * sigma) * maturity +
                     phi * sqrt(maturity) * sigma);
    double ST2 =
        stock0 * exp((interestRate - D - 0.5 * sigma * sigma) * maturity +
                     (-phi) * sqrt(maturity) * sigma);
    sum += payoff(ST1, strikePrice);
    sum += payoff(ST2, strikePrice);
  }
  return sum / N * exp(-interestRate * maturity);
}

double european_options_analytic(
    const double stock0, const double strikePrice, const double D,
    const double interestRate, const double sigma, const double maturity,
    const std::function<double(double, double, double, double, double, double,
                               double, double)>
        payoff,
    const double time) {

  const double d1 =
      (std::log(stock0 / strikePrice) +
       ((interestRate - D + (pow(sigma, 2) / 2)) * (maturity - time))) /
      (std::sqrt(maturity - time) * sigma);
  const double d2 = d1 - sigma * std::sqrt(maturity - time);
  return payoff(strikePrice, interestRate, maturity, time, d1, d2, D, stock0);
}

double path_dependent_options_monteCarlo_momentmatching(double T, double K,
                                                        double N, double S0,
                                                        double r, double D,
                                                        double sigma,
                                                        double X) {
  static std::mt19937 rng;
  std::normal_distribution<> ND(0, 1);
  std::vector<std::vector<double>> phi_vector(int(N / 2),
                                              std::vector<double>(int(K)));
  double sum = 0.0;
  for (int i = 0; i < int(N / 2); i += 1) {
    for (int k = 0; k < K; k++) {
      double phi = ND(rng);
      phi_vector[i][k] = phi;
      sum += 2 * std::pow(phi, 2);
    }
  }
  double sample_variance = std::sqrt(sum / ((N)*2 * (K)));
  sum = 0.;
  for (int i = 0; i < int(N / 2); i++) {
    double dt = T / K;
    std::vector<double> stockPath1(K + 1);
    std::vector<double> stockPath2(K + 1);
    // initialise first value
    stockPath1[0] = S0;
    stockPath2[0] = S0;
    double A1 = 0.;
    double A2 = 0.;
    for (int k = 1; k <= K; k++) {
      double phi = phi_vector[i][(k - 1)];
      stockPath1[k] =
          stockPath1[k - 1] *
          exp((r - D - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * phi);
      stockPath2[k] =
          stockPath2[k - 1] *
          exp((r - D - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * (-phi));
    }
    for (int k = 1; k <= K; k++) {
      A1 += stockPath1[k];
      A2 += stockPath2[k];
    }
    A1 /= K;
    A2 /= K;
    sum += std::max(A1 - X, 0.);
    sum += std::max(A2 - X, 0.);
  }
  double V_MC = (sum / (N)) * exp(-r * T);
  return V_MC;
}

double path_dependent_options_monteCarlo(double T, double K, double N,
                                         double S0, double r, double D,
                                         double sigma, double X) {
  static std::mt19937 rng;
  std::normal_distribution<> ND(0., 1.);
  double sum = 0.;
  for (int n = 0; n < N; n++) {
    // now create a path
    double dt = T / K;
    std::vector<double> stockPath(K + 1);
    stockPath[0] = S0;
    for (int i = 1; i <= K; i++) {
      double phi = ND(rng);
      stockPath[i] = stockPath[i - 1] * exp((r - D - 0.5 * sigma * sigma) * dt +
                                            phi * sigma * sqrt(dt));
    }
    // and calculate A
    double A = 0.;
    for (int i = 1; i <= K; i++) {
      A += (stockPath[i]);
    }
    A /= K;
    sum += std::max(A - X, 0.);
  }
  return sum / N * exp(-r * T);
}

void path_independent_options() {
  // path_dependant_option();
  const double T = 0.5, sigma = 0.41, r = 0.03, D = 0.04, X1 = 70., X2 = 100.;

  auto long_call_payoff = [](const double S, const double X) {
    return std::max(S - X, 0.);
  };
  auto short_call_payoff = [](const double S, const double X) {
    return -std::max(S - X, 0.);
  };
  auto binary_call_payoff = [](const double S, const double X) {
    if (S <= X) {
      return 0.;
    } else {
      return 1.;
    }
  };
  auto analytic_long_call_payoff =
      [](const double X, const double r, const double T, const double t,
         const double d1, const double d2, const double D, const double S) {
        return -X * std::exp(-r * (T - t)) * normalDistribution(d2) +
               S * std::exp(-D * (T - t)) * normalDistribution(d1);
      };
  auto analytic_short_call_payoff =
      [](const double X, const double r, const double T, const double t,
         const double d1, const double d2, const double D, const double S) {
        return +X * std::exp(-r * (T - t)) * normalDistribution(d2) -
               S * std::exp(-D * (T - t)) * normalDistribution(d1);
      };
  auto analytic_long_binary_call_payoff =
      [](const double X, const double r, const double T, const double t,
         const double d1, const double d2, const double D, const double S) {
        return std::exp(-r * (T - t)) * normalDistribution(d2);
      };

  const double min_N = 1000;
  const double max_N = 100000;
  const size_t N_data_points = 100;
  const double N_interval = (max_N - min_N) / (N_data_points - 1);
  const double min_S0 = 70.;
  const double max_S0 = 100.;
  const size_t S0_data_points = 2;
  const double S0_interval = (max_S0 - min_S0) / (S0_data_points - 1);
  const double M = 100;
  double S0 = min_S0;

  std::ofstream output1("./Assignment_3/outputs/analytic.task.2.1.csv");
  for (size_t S{1}; S < 170; S += 1) {
    if (S == X1 || S == X2) {
      continue;
    }
    double analytic_pi_0 = 0;
    analytic_pi_0 += european_options_analytic(S, X1, D, r, sigma, T,
                                               analytic_short_call_payoff, 0);
    analytic_pi_0 += european_options_analytic(S, X2, D, r, sigma, T,
                                               analytic_long_call_payoff, 0);
    analytic_pi_0 +=
        2 * X2 *
        european_options_analytic(S, X2, D, r, sigma, T,
                                  analytic_long_binary_call_payoff, 0);
    analytic_pi_0 += european_options_analytic(S, 0., D, r, sigma, T,
                                               analytic_long_call_payoff, 0);
    double analytic_pi_T = 0;
    analytic_pi_T += european_options_analytic(S, X1, D, r, sigma, T,
                                               analytic_short_call_payoff, T);
    analytic_pi_T += european_options_analytic(S, X2, D, r, sigma, T,
                                               analytic_long_call_payoff, T);
    analytic_pi_T +=
        2 * X2 *
        european_options_analytic(S, X2, D, r, sigma, T,
                                  analytic_long_binary_call_payoff, T);
    analytic_pi_T += european_options_analytic(S, 0., D, r, sigma, T,
                                               analytic_long_call_payoff, T);
    output1 << S << "," << analytic_pi_0 << "," << analytic_pi_T << std::endl;
  }

  std::ofstream output2("./Assignment_3/outputs/momentmatch.task.2.1.csv");
  for (size_t i{}; i < S0_data_points; i += 1) {
    double N = min_N;
    for (size_t j{}; j < N_data_points; j += 1) {
      // Carry out M calculations of montecarlo simulations for the portfolio
      // with N random generations

      std::vector<double> PI_N(M);
      double sum = 0;
      double variance = 0;
      double sample_variance = 0;
      double lower_confidence_limit = 0;
      double upper_confidence_limit = 0;
      double mean = 0;
      auto t1 = std::chrono::high_resolution_clock::now();
      for (size_t k{}; k < M; k += 1) {
        double montecarlo_pi = 0;
        montecarlo_pi += european_options_monteCarlo_moment_match(
            S0, X1, D, r, sigma, T, short_call_payoff, N);
        montecarlo_pi += european_options_monteCarlo_moment_match(
            S0, X2, D, r, sigma, T, long_call_payoff, N);
        montecarlo_pi += 2 * X2 *
                         european_options_monteCarlo_moment_match(
                             S0, X2, D, r, sigma, T, binary_call_payoff, N);
        montecarlo_pi += european_options_monteCarlo_moment_match(
            S0, 0., D, r, sigma, T, long_call_payoff, N);
        sum += montecarlo_pi;
        PI_N[k] = montecarlo_pi;
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      auto time_taken =
          std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
              .count();
      mean = sum / M;
      sum = 0.;
      for (int l{0}; l < M; l += 1) {
        sum += std::pow((PI_N[l] - mean), 2);
      }
      sample_variance = sum / M * (M - 1);
      lower_confidence_limit = mean - 2 * std::sqrt(sample_variance);
      upper_confidence_limit = mean + 2 * std::sqrt(sample_variance);
      // Carry out analytic value of portfolio
      double analytic_pi = 0;
      analytic_pi += european_options_analytic(S0, X1, D, r, sigma, T,
                                               analytic_short_call_payoff, 0);
      analytic_pi += european_options_analytic(S0, X2, D, r, sigma, T,
                                               analytic_long_call_payoff, 0);
      analytic_pi +=
          2 * X2 *
          european_options_analytic(S0, X2, D, r, sigma, T,
                                    analytic_long_binary_call_payoff, 0);
      analytic_pi += european_options_analytic(S0, 0., D, r, sigma, T,
                                               analytic_long_call_payoff, 0);
      // Output to file
      output2 << time_taken << "," << N << "," << S0 << "," << mean << ","
              << std::sqrt(sample_variance) << "," << lower_confidence_limit
              << "," << upper_confidence_limit << "," << analytic_pi
              << std::endl;
      N += N_interval;
    }

    S0 += S0_interval;
  }
}

void path_dependent_options() {

  std::mt19937 rng;
  std::normal_distribution<> ND(0, 1.);

  double S0 = 64000, sigma = 0.42, r = 0.03, T = 2, X = 64000, D = 0.02;
  int K = 20;
  const double min_N = 1000;
  const double max_N = 50000;
  const size_t N_data_points = 3;
  double N_Values[3] = {9000, 10000, 8000};
  double M_Values[3] = {300, 500, 600};
  const double N_interval = (max_N - min_N) / (N_data_points - 1);
  // const double N_interval = 0;
  std::ofstream output1("./Assignment_3/outputs/final.task.2.2.csv");
  double K_Values[4] = {70, 150};
  size_t max_paths = 70;
  double N = min_N;

  for (size_t j{}; j < N_data_points; j += 1) {
    N = N_Values[j];
    // for (size_t current_k{1}; current_k < 100; current_k += 1) {
    for (size_t current_path{0}; current_path < 3; current_path += 1) {
      // double K = K_Values[current_k];
      K = 20.0;
      double sum = 0;
      double variance = 0;
      double sample_variance = 0;
      double lower_confidence_limit = 0;
      double upper_confidence_limit = 0;
      auto t1 = std::chrono::high_resolution_clock::now();
      size_t paths = M_Values[current_path];
      std::vector<double> PI_N(paths);
      for (size_t p{0}; p < paths; p += 1) {
        double montecarlo = 0;
        montecarlo = path_dependent_options_monteCarlo_momentmatching(
            T, K, N, S0, r, D, sigma, X);
        sum += montecarlo;
        PI_N[p] = montecarlo;
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      auto time_taken =
          std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
              .count();
      double mean = sum / paths;
      sum = 0.;
      for (int l{0}; l < paths; l += 1) {
        sum += std::pow((PI_N[l] - mean), 2);
      }
      sample_variance = sum / (paths * (paths - 1));
      lower_confidence_limit = mean - 2 * std::sqrt(sample_variance);
      upper_confidence_limit = mean + 2 * std::sqrt(sample_variance);
      // Output to file
      output1 << time_taken << "," << N << "," << K << "," << paths << ","
              << mean << "," << std::sqrt(sample_variance) << ","
              << lower_confidence_limit << "," << upper_confidence_limit
              << std::endl;
      // N += N_interval;
    }
    //}
  }
}

int main() {
  path_dependent_options();
  // path_independent_options();
}