//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "TVS_header.h"
#include <random>
#include <chrono>
#include <numeric>
#include <thread>
#include <functional> // for std::hash
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
//
// [[Rcpp::export]]
arma::vec TVS_cpp(
    Rcpp::Function wrapper_beta,
    Rcpp::Function wrapper_beta0,
    Rcpp::Function wrapper_nu,
    Rcpp::Function wrapper_gamma,
    Rcpp::List dataXY,
    arma::vec init_beta_TVS,
    int B,
    double init_beta0_TVS,
    double init_nu_TVS,
    double init_gamma_TVS,
    double init_theta_TVS,
    double SS_t0_TVS,
    double SS_t1_TVS,
    double hyper_mu_beta0_TVS,
    double hyper_sigma_beta0_TVS,
    double hyper_mu_nu_TVS,
    double hyper_sigma_nu_TVS,
    double hyper_c_gamma_TVS,
    double hyper_d_gamma_TVS,
    double hyper_a_theta_TVS,
    double hyper_b_theta_TVS,
    int max_iter_TVS,
    double tol_TVS,
    double add_correc_CiS) {

  arma::vec dat_Y = dataXY["Y"];
  arma::mat dat_X = dataXY["X"];
  int p = dat_X.n_cols;

  // Check if hyper_b_theta was provided
  double hyper_b_theta_val = (hyper_b_theta_TVS < 0) ? p : hyper_b_theta_TVS;

  // Original CiS
  Rcpp::List em_orig = TVS_EM_cpp(wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma,
                                  dataXY, init_beta_TVS, init_beta0_TVS, init_nu_TVS, init_gamma_TVS, init_theta_TVS,
                                  SS_t0_TVS, SS_t1_TVS, hyper_mu_beta0_TVS, hyper_sigma_beta0_TVS, hyper_mu_nu_TVS, hyper_sigma_nu_TVS,
                                  hyper_c_gamma_TVS, hyper_d_gamma_TVS, hyper_a_theta_TVS, hyper_b_theta_val, max_iter_TVS, tol_TVS);

  arma::vec pvals(p, arma::fill::zeros);

  for (int j = 0; j < p; ++j) {
    arma::uvec test_j = {static_cast<unsigned int>(j)};
    int j_index_cpp = j + 1;

    //observed CiSj
    double CiS_orig = CiS_group_fun_cpp(test_j+1,
                                        Rcpp::as<arma::vec>(em_orig["beta"]),
                                        Rcpp::as<double>(em_orig["beta0"]),
                                        Rcpp::as<double>(em_orig["nu"]),
                                        Rcpp::as<double>(em_orig["gamma"]),
                                        dataXY,
                                        add_correc_CiS);

    int count = 0;

    for (int b = 0; b < B; ++b) {

      Rcpp::Rcout << "Running perm #" << b+1 << ", predictor " << j_index_cpp << std::endl;

      double CiS_perm = per_fun_cpp(j_index_cpp, wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma,
                                    dataXY, init_beta_TVS, init_beta0_TVS, init_nu_TVS, init_gamma_TVS, init_theta_TVS,
                                    SS_t0_TVS, SS_t1_TVS, hyper_mu_beta0_TVS, hyper_sigma_beta0_TVS, hyper_mu_nu_TVS, hyper_sigma_nu_TVS,
                                    hyper_c_gamma_TVS, hyper_d_gamma_TVS, hyper_a_theta_TVS, hyper_b_theta_val, max_iter_TVS, tol_TVS, add_correc_CiS
      );

      if (CiS_perm >= CiS_orig) count++;
    }

    pvals(j) = static_cast<double>(count) / (B); // add 1 for smoothing
  }

  return pvals;
}


// [[Rcpp::export]]
Rcpp::List TVS_multi_stage_cpp(
    Rcpp::Function wrapper_beta,
    Rcpp::Function wrapper_beta0,
    Rcpp::Function wrapper_nu,
    Rcpp::Function wrapper_gamma,
    Rcpp::List dataXY,
    arma::vec init_beta_TVS,
    int group_B,
    int indiv_B,
    int B_final,
    double group_cutoff,
    double indiv_cutoff,
    int group_size,
    double init_beta0_TVS,
    double init_nu_TVS,
    double init_gamma_TVS,
    double init_theta_TVS,
    double SS_t0_TVS,
    double SS_t1_TVS,
    double hyper_mu_beta0_TVS,
    double hyper_sigma_beta0_TVS,
    double hyper_mu_nu_TVS,
    double hyper_sigma_nu_TVS,
    double hyper_c_gamma_TVS,
    double hyper_d_gamma_TVS,
    double hyper_a_theta_TVS,
    double hyper_b_theta_TVS,
    int max_iter_TVS,
    double tol_TVS,
    double add_correc_CiS) {

  arma::mat dat_X = dataXY["X"];
  int p = dat_X.n_cols;
  int G = std::ceil((double)p / group_size);
  arma::vec final_pvals(p, arma::fill::value(-1.0)); // initialize with -1
  std::vector<int> final_selected;

  // Check if hyper_b_theta was provided
  double hyper_b_theta_val = (hyper_b_theta_TVS < 0) ? p : hyper_b_theta_TVS;

  // Original model fit
  Rcpp::List em_orig = TVS_EM_cpp(wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma,
                                  dataXY, init_beta_TVS, init_beta0_TVS, init_nu_TVS, init_gamma_TVS, init_theta_TVS,
                                  SS_t0_TVS, SS_t1_TVS, hyper_mu_beta0_TVS, hyper_sigma_beta0_TVS, hyper_mu_nu_TVS, hyper_sigma_nu_TVS,
                                  hyper_c_gamma_TVS, hyper_d_gamma_TVS, hyper_a_theta_TVS, hyper_b_theta_val, max_iter_TVS, tol_TVS);
  arma::vec beta_orig = em_orig["beta"];
  double beta0_orig = em_orig["beta0"];
  double nu_orig = em_orig["nu"];
  double gamma_orig = em_orig["gamma"];

  // ==== STEP 1: Group-Level Screening ====
  std::vector<arma::uvec> selected_groups;
  for (int g = 0; g < G; ++g) {
    int start = g * group_size;
    int end = std::min((g + 1) * group_size, p);
    arma::uvec group_indices = arma::regspace<arma::uvec>(start + 1, end); // 1-based

    double CiS_orig = CiS_group_fun_cpp(group_indices, beta_orig, beta0_orig, nu_orig, gamma_orig, dataXY, add_correc_CiS);

    int count = 0;
    for (int b = 0; b < group_B; ++b) {
      double perm_CiS = per_group_fun_cpp(group_indices, wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma,
                                          dataXY, init_beta_TVS, init_beta0_TVS, init_nu_TVS, init_gamma_TVS, init_theta_TVS,
                                          SS_t0_TVS, SS_t1_TVS, hyper_mu_beta0_TVS, hyper_sigma_beta0_TVS, hyper_mu_nu_TVS, hyper_sigma_nu_TVS,
                                          hyper_c_gamma_TVS, hyper_d_gamma_TVS, hyper_a_theta_TVS, hyper_b_theta_val, max_iter_TVS, tol_TVS, add_correc_CiS
      );
      if (perm_CiS >= CiS_orig) count++;
    }

    double pval_group = static_cast<double>(count + 1) / (group_B + 1);
    if (pval_group <= group_cutoff) {
      selected_groups.push_back(group_indices);
    }
  }
  Rcpp::Rcout << "[Info] Group screening done.\n";

  if (selected_groups.empty()) {
    Rcpp::Rcout << "[Info] No groups passed the screening step.\n";
    return Rcpp::List::create(
      Rcpp::Named("selected_indices") = Rcpp::IntegerVector(0),
      Rcpp::Named("p_values") = final_pvals
    );
  }

  // ==== STEP 2: Individual-Level Screening ====
  arma::uvec predictors_to_screen;
  for (const auto& group : selected_groups) {
    predictors_to_screen = arma::join_vert(predictors_to_screen, group);
  }

  predictors_to_screen = arma::unique(predictors_to_screen);
  predictors_to_screen = arma::sort(predictors_to_screen);

  std::vector<int> predictors_pass;

  for (unsigned int i = 0; i < predictors_to_screen.n_elem; ++i) {
    int j_idx = predictors_to_screen[i]; // 1-based

    arma::uvec j_vec = {static_cast<unsigned int>(j_idx)};
    double CiS_orig = CiS_group_fun_cpp(j_vec, beta_orig, beta0_orig, nu_orig, gamma_orig, dataXY, add_correc_CiS);

    int count = 0;
    for (int b = 0; b < indiv_B; ++b) {
      double perm_CiS = per_fun_cpp(j_idx, wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma,
                                    dataXY, init_beta_TVS, init_beta0_TVS, init_nu_TVS, init_gamma_TVS, init_theta_TVS,
                                    SS_t0_TVS, SS_t1_TVS, hyper_mu_beta0_TVS, hyper_sigma_beta0_TVS, hyper_mu_nu_TVS, hyper_sigma_nu_TVS,
                                    hyper_c_gamma_TVS, hyper_d_gamma_TVS, hyper_a_theta_TVS, hyper_b_theta_val, max_iter_TVS, tol_TVS, add_correc_CiS
      );
      if (perm_CiS >= CiS_orig) count++;
    }

    double pval = static_cast<double>(count + 1) / (indiv_B + 1);
    if (pval <= indiv_cutoff) {
      predictors_pass.push_back(j_idx - 1); // save as 0-based index
    }
  }
  Rcpp::Rcout << "[Info] Individual screening done.\n";

  if (predictors_pass.empty()) {
    Rcpp::Rcout << "[Info] No individual predictors passed the screening step.\n";
    return Rcpp::List::create(
      Rcpp::Named("selected_indices") = Rcpp::IntegerVector(0),
      Rcpp::Named("p_values") = final_pvals
    );
  }

  // ==== STEP 3: Final Full Permutation ====
  for (int j : predictors_pass) {
    int j_idx = j + 1;

    arma::uvec j_vec = {static_cast<unsigned int>(j_idx)};
    double CiS_orig = CiS_group_fun_cpp(j_vec, beta_orig, beta0_orig, nu_orig, gamma_orig, dataXY, add_correc_CiS);

    int count = 0;
    for (int b = 0; b < B_final; ++b) {
      double perm_CiS = per_fun_cpp(j_idx, wrapper_beta, wrapper_beta0, wrapper_nu, wrapper_gamma,
                                    dataXY, init_beta_TVS, init_beta0_TVS, init_nu_TVS, init_gamma_TVS, init_theta_TVS,
                                    SS_t0_TVS, SS_t1_TVS, hyper_mu_beta0_TVS, hyper_sigma_beta0_TVS, hyper_mu_nu_TVS, hyper_sigma_nu_TVS,
                                    hyper_c_gamma_TVS, hyper_d_gamma_TVS, hyper_a_theta_TVS, hyper_b_theta_val, max_iter_TVS, tol_TVS, add_correc_CiS
      );
      if (perm_CiS >= CiS_orig) count++;
    }

    double final_pval = static_cast<double>(count + 1) / (B_final + 1);
    final_pvals(j) = final_pval;
    final_selected.push_back(j);
  }
  Rcpp::Rcout << "[Info] Final step done.\n";

  arma::uvec arma_vec_final_selected = arma::conv_to<arma::uvec>::from(final_selected);
  arma::uvec final_selected_1based = arma_vec_final_selected + 1;

  return Rcpp::List::create(
    Rcpp::Named("selected_indices") = final_selected_1based,
    Rcpp::Named("p_values") = final_pvals
  );
}
