//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "TVS_header.h"
#include <random>
#include <chrono>
#include <numeric>
#include <thread>
#include <functional> // for std::hash
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
//
// [[Rcpp::export]]
double beta_neg_lk_cpp(
    arma::vec beta_lk, double beta0_lk, double nu_lk, double ga_lk, arma::vec betaPRE, double t0, double t1,
    arma::vec Y_lk, arma::mat X_lk,
    double theta_lk) {

  int n = Y_lk.n_elem;
  int p = X_lk.n_cols;  // Get number of predictors

  // Check if theta_lk was provided
  double theta_val = (theta_lk < 0) ? 1.0 / p : theta_lk;

  // Ensure beta_lk size matches X_lk columns
  if (static_cast<int>(beta_lk.n_elem) != p) {
    Rcpp::stop("Size of beta coefficients does not match predictor dimensions");
  }
  // Compute residuals
  arma::vec residuals = Y_lk - beta0_lk - X_lk * beta_lk;

  arma::vec ga_pow(n);
  for (int i = 0; i < n; ++i) {
    double sign_res = (residuals[i] >= 0) ? -1.0 : 1.0;
    ga_pow[i] = std::pow(ga_lk, 2 * sign_res);
  }
  /* arma::vec ratio = pow(residuals, 2) / nu_lk * ga_pow; */
  arma::vec ratio = (pow(residuals, 2) / nu_lk) % ga_pow;

  arma::vec inverseprob = 1 + (t0 * (1 - theta_val) / (t1 * theta_val)) * exp(-abs(betaPRE) * (t0 - t1));

  double bneglog = arma::dot(abs(beta_lk), (t0 - (t0 - t1) / inverseprob)) + (nu_lk/2 + 0.5) * arma::sum(arma::log1p(ratio));

  return  bneglog;
}

// [[Rcpp::export]]
arma::vec beta_neg_gradient_cpp(
    arma::vec beta_lk, double beta0_lk, double ga_lk, double nu_lk,
    arma::vec betaPRE, double t0, double t1, arma::vec Y_lk, arma::mat X_lk,
    double theta_lk) {

  int n = Y_lk.n_elem;
  int p = X_lk.n_cols;  // Get number of predictors

  // Check if theta_lk was provided
  double theta_val = (theta_lk < 0) ? 1.0 / p : theta_lk;

  // Ensure beta_lk size matches X_lk columns
  if (static_cast<int>(beta_lk.n_elem) != p) {
    Rcpp::stop("Size of beta coefficients does not match predictor dimensions");
  }
  // Compute residuals
  arma::vec residuals = Y_lk - beta0_lk - X_lk * beta_lk;

  arma::vec ga_pow(n);
  for (int i = 0; i < n; ++i) {
    double sign_res = (residuals[i] >= 0) ? -1.0 : 1.0;
    ga_pow[i] = std::pow(ga_lk, 2 * sign_res);
  }
  arma::vec ratio = (pow(residuals, 2) / nu_lk) % ga_pow;

  arma::vec inverseprob = 1 + (t0 * (1 - theta_val) / (t1 * theta_val)) * exp(-abs(betaPRE) * (t0 - t1));

  arma::vec grad_denom = (2 * residuals % ga_pow / nu_lk) / (1 + ratio);

  arma::vec signbeta_lk = arma::sign(beta_lk);

  arma::vec gradient = signbeta_lk % (t0 - (t0 - t1) / inverseprob) -(nu_lk/2 + 0.5) * X_lk.t() * grad_denom;

  return Rcpp::NumericVector(gradient.begin(), gradient.end());
  //return gradient;
}

// [[Rcpp::export]]
double beta0_neg_lk_cpp(double beta0_lk, arma::vec beta_lk, double nu_lk, double gamma_lk,
                        double hyper_mu_beta0, double hyper_sigma_beta0,
                        arma::vec Y_lk, arma::mat X_lk) {

  int n = Y_lk.n_elem;
  int p = X_lk.n_cols;

  // Ensure beta_lk size matches X_lk columns
  if (static_cast<int>(beta_lk.n_elem) != p) {
    Rcpp::stop("Size of beta coefficients does not match predictor dimensions");
  }

  // Compute residuals
  arma::vec residuals = Y_lk - beta0_lk - X_lk * beta_lk;

  // Create index vector
  arma::vec index(n);
  for(int i = 0; i < n; i++) {
    index[i] = (residuals[i] >= 0) ? -1.0 : 1.0;
  }

  // Calculate sum in the formula
  double sum_term = 0.0;
  for(int i = 0; i < n; i++) {
    double ga_term = std::pow(gamma_lk, 2.0 * index[i]);
    sum_term += std::log(1.0 + residuals[i] * residuals[i] / nu_lk * ga_term);
  }

  return  0.5 * std::pow(beta0_lk - hyper_mu_beta0, 2) / hyper_sigma_beta0 +
    (nu_lk / 2.0 + 0.5) * sum_term;
}


// [[Rcpp::export]]
double nu_neg_lk_cpp(double nu_lk, double ga_lk, Rcpp::NumericVector error_lk,
                     double hyper_mu, double hyper_sigma) {
  int n = error_lk.size();
  Rcpp::NumericVector index(n);

  // Create index vector (ifelse(error_lk>=0,-1,1))
  for(int i = 0; i < n; i++) {
    index[i] = (error_lk[i] >= 0) ? -1.0 : 1.0;
  }

  // Calculate sum in the formula
  double sum_term = 0.0;
  for(int i = 0; i < n; i++) {
    double ga_term = std::pow(ga_lk, 2.0 * index[i]);
    sum_term += std::log(1.0 + error_lk[i] * error_lk[i] / nu_lk * ga_term);
  }

  return (std::pow(std::log(nu_lk) - hyper_mu, 2) / (2.0 * std::pow(hyper_sigma, 2))) +
    (n / 2.0 + 1.0) * std::log(nu_lk) -
    n * (R::lgammafn(nu_lk / 2.0 + 0.5) - R::lgammafn(nu_lk / 2.0)) +
    (nu_lk / 2.0 + 0.5) * sum_term;
}


// [[Rcpp::export]]
double gamma_neg_lk_cpp(double ga_lk, double nu_lk, Rcpp::NumericVector error_lk,
                        double hyper_c, double hyper_d) {
  int n = error_lk.size();
  Rcpp::NumericVector index(n);

  // Create index vector (ifelse(error_lk>=0,-1,1))
  for(int i = 0; i < n; i++) {
    index[i] = (error_lk[i] >= 0) ? -1.0 : 1.0;
  }

  // Calculate sum in the formula
  double sum_term = 0.0;
  for(int i = 0; i < n; i++) {
    double ga_term = std::pow(ga_lk, 2.0 * index[i]);
    sum_term += std::log(1.0 + error_lk[i] * error_lk[i] / nu_lk * ga_term);
  }

  return -(n + hyper_c - 1.0) * std::log(ga_lk) + hyper_d * ga_lk +
    n * std::log(std::pow(ga_lk, 2) + 1.0) + (nu_lk / 2.0 + 0.5) * sum_term;
}
