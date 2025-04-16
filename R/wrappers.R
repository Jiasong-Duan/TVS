#' Wrapper function for beta optimization
#'
#' @param beta Initial beta parameter
#' @param beta0 Beta0 parameter
#' @param nu Nu parameter
#' @param gamma Gamma parameter
#' @param betaPRE Previous beta value
#' @param t0 Initial temperature for simulated annealing (default: 10)
#' @param t1 Final temperature for simulated annealing (default: 1)
#' @param Y_lk Response vector
#' @param X_lk Predictor matrix
#' @param theta_lk Theta parameter (default: -2)
#' @return Optimized beta parameter
#' @export
wrapper_beta <- function(beta, beta0, nu, gamma, betaPRE, t0=10, t1=1, Y_lk, X_lk, theta_lk=-2) {
  func_lk <- function(x) beta_neg_lk_cpp(beta_lk=x, beta0_lk=beta0, nu_lk=nu, ga_lk=gamma, betaPRE=betaPRE,
                                         t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk)
  func_gradient <- function(x) beta_neg_gradient_cpp(beta_lk=x, beta0_lk=beta0, nu_lk=nu, ga_lk=gamma, betaPRE=betaPRE,
                                                     t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk)
  result <- tryCatch(
    optim(par = beta, fn = func_lk, gr= func_gradient, method = "L-BFGS-B", lower = min(Y_lk), upper = max(Y_lk)),
    error = function(e) {
      warning("beta optimization failed: ", conditionMessage(e))
      list(par = NA)
    }
  )
  return(result$par)
}

#' Wrapper function for beta0 optimization
#'
#' @param beta0.lk Initial beta0 parameter
#' @param beta.lk Beta parameter
#' @param nu.lk Nu parameter
#' @param gamma.lk Gamma parameter
#' @param hyper.mu.beta0 Hyperparameter for beta0 mean (default: 0)
#' @param hyper.sigma.beta0 Hyperparameter for beta0 variance (default: 10^6)
#' @param Y.lk Response vector
#' @param X.lk Predictor matrix
#' @return Optimized beta0 parameter
#' @export
wrapper_beta0 <- function(beta0.lk, beta.lk, nu.lk, gamma.lk, hyper.mu.beta0=0, hyper.sigma.beta0=10^6, Y.lk, X.lk){
  func <- function(x) beta0_neg_lk_cpp(beta0_lk=x, beta_lk=beta.lk, nu_lk=nu.lk, gamma_lk=gamma.lk,
                                       hyper_mu_beta0= hyper.mu.beta0, hyper_sigma_beta0= hyper.sigma.beta0, Y_lk=Y.lk, X_lk=X.lk)
  result <- tryCatch(
    optim(par = beta0.lk, fn = func, method = "L-BFGS-B", lower = -100, upper = 100),
    error = function(e) {
      warning("beta0 optimization failed: ", conditionMessage(e))
      list(par = NA)
    }
  )
  return(result$par)
}

#' Wrapper function for nu optimization
#'
#' @param nu Initial nu parameter
#' @param gamma Gamma parameter
#' @param error_lk Error vector
#' @param hyper_mu Hyperparameter for nu mean
#' @param hyper_sigma Hyperparameter for nu variance
#' @return Optimized nu parameter
#' @export
wrapper_nu <- function(nu, gamma, error_lk, hyper_mu, hyper_sigma) {
  func <- function(x) nu_neg_lk_cpp(nu_lk=x, ga_lk=gamma, error_lk, hyper_mu, hyper_sigma)
  result <- tryCatch(
    optim(par = nu, fn = func, method = "L-BFGS-B", lower = 0.01, upper = 100),
    error = function(e) {
      warning("nu optimization failed: ", conditionMessage(e))
      list(par = NA)
    }
  )
  return(result$par)
}

#' Wrapper function for gamma optimization
#'
#' @param gamma Initial gamma parameter
#' @param nu Nu parameter
#' @param error_lk Error vector
#' @param hyper_c Hyperparameter for gamma shape
#' @param hyper_d Hyperparameter for gamma rate
#' @return Optimized gamma parameter
#' @export
wrapper_gamma <- function(gamma, nu, error_lk, hyper_c, hyper_d) {
  func <- function(x) gamma_neg_lk_cpp(ga_lk=x, nu_lk=nu, error_lk, hyper_c, hyper_d)
  result <- tryCatch(
    optim(par = gamma, fn = func, method = "L-BFGS-B", lower = 0.01, upper = 100),
    error = function(e) {
      warning("gamma optimization failed: ", conditionMessage(e))
      list(par = NA)
    }
  )
  return(result$par)
}

#' EM algorithm for TVS
#'
#' @param wrapper_beta Function wrapper for beta
#' @param wrapper_beta0 Function wrapper for beta0
#' @param wrapper_nu Function wrapper for nu
#' @param wrapper_gamma Function wrapper for gamma
#' @param dataXY List containing data X and Y
#' @param init_beta Initial vector for beta
#' @param init_beta0 Initial value for beta0 (default: 1)
#' @param init_nu Initial value for nu (default: 1)
#' @param init_gamma Initial value for gamma (default: 1)
#' @param init_theta Initial value for theta (default: 0.5)
#' @param SS_t0 Initial temperature for simulated annealing (default: 10.0)
#' @param SS_t1 Final temperature for simulated annealing (default: 1.0)
#' @param hyper_mu_beta0 Hyperparameter for beta0 mean (default: 0)
#' @param hyper_sigma_beta0 Hyperparameter for beta0 variance (default: 1e6)
#' @param hyper_mu_nu Hyperparameter for nu mean (default: 1)
#' @param hyper_sigma_nu Hyperparameter for nu variance (default: 1)
#' @param hyper_c_gamma Hyperparameter for gamma shape (default: 0.0001)
#' @param hyper_d_gamma Hyperparameter for gamma rate (default: 0.0001)
#' @param hyper_a_theta Hyperparameter for theta shape (default: 1)
#' @param hyper_b_theta Hyperparameter for theta shape (default: -1.0)
#' @param max_iter Maximum number of iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param conv_type Convergence type (default: "param")
#' @return List of results from the EM algorithm
#' @export
TVS_EM <- function(
    wrapper_beta,
    wrapper_beta0,
    wrapper_nu,
    wrapper_gamma,
    dataXY,
    init_beta = rep(1,8),
    init_beta0 = 1,
    init_nu = 1,
    init_gamma = 1,
    init_theta = 0.5,
    SS_t0 = 10.0,
    SS_t1 = 1.0,
    hyper_mu_beta0 = 0,
    hyper_sigma_beta0 = 1e6,
    hyper_mu_nu = 1,
    hyper_sigma_nu = 1,
    hyper_c_gamma = 0.0001,
    hyper_d_gamma = 0.0001,
    hyper_a_theta = 1,
    hyper_b_theta = -1.0,
    max_iter = 100,
    tol = 1e-6,
    conv_type = "param") {
  TVS_EM_cpp(
    wrapper_beta = wrapper_beta,
    wrapper_beta0 = wrapper_beta0,
    wrapper_nu = wrapper_nu,
    wrapper_gamma = wrapper_gamma,
    dataXY = dataXY,
    init_beta = init_beta,
    init_beta0 = init_beta0,
    init_nu = init_nu,
    init_gamma = init_gamma,
    init_theta = init_theta,
    SS_t0 = SS_t0,
    SS_t1 = SS_t1,
    hyper_mu_beta0 = hyper_mu_beta0,
    hyper_sigma_beta0 = hyper_sigma_beta0,
    hyper_mu_nu = hyper_mu_nu,
    hyper_sigma_nu = hyper_sigma_nu,
    hyper_c_gamma = hyper_c_gamma,
    hyper_d_gamma = hyper_d_gamma,
    hyper_a_theta = hyper_a_theta,
    hyper_b_theta = hyper_b_theta,
    max_iter = max_iter,
    tol = tol,
    conv_type = conv_type
  )
}

#' Calculate slope error
#'
#' @param err Error vector
#' @param nu_slo Nu parameter for slope
#' @param ga_slo Gamma parameter for slope
#' @return Vector of slope errors
#' @export
slope_error <- function(err, nu_slo, ga_slo) {
  slope_error_cpp(err = err, nu_slo = nu_slo, ga_slo = ga_slo)
}

#' Calculate curvature error
#'
#' @param err Error vector
#' @param nu_cur Nu parameter for curvature
#' @param ga_cur Gamma parameter for curvature
#' @return Vector of curvature errors
#' @export
curvature_error <- function(err, nu_cur, ga_cur) {
  curvature_error_cpp(err = err, nu_cur = nu_cur, ga_cur = ga_cur)
}

#' Calculate CiS for a single index
#'
#' @param test_index Index to test
#' @param beta_opt Optimized beta vector
#' @param beta0_opt Optimized beta0 value
#' @param nu_opt Optimized nu value
#' @param ga_opt Optimized gamma value
#' @param dataXY List containing data X and Y
#' @param add_correc_CiS Correction factor for CiS (default: 0.001)
#' @return CiS value
#' @export
CiS_j_fun <- function(test_index,
                      beta_opt,
                      beta0_opt,
                      nu_opt,
                      ga_opt,
                      dataXY,
                      add_correc_CiS = 0.001) {

  CiS_j_fun_cpp(
    test_index = test_index,
    beta_opt = beta_opt,
    beta0_opt = beta0_opt,
    nu_opt = nu_opt,
    ga_opt = ga_opt,
    dataXY = dataXY,
    add_correc_CiS = add_correc_CiS
  )
}

#' Calculate permutation value for a single index
#'
#' @param j_index Index for permutation
#' @param wrapper_beta Function wrapper for beta
#' @param wrapper_beta0 Function wrapper for beta0
#' @param wrapper_nu Function wrapper for nu
#' @param wrapper_gamma Function wrapper for gamma
#' @param dataXY List containing data X and Y
#' @param init_beta_per Initial vector for beta
#' @param init_beta0_per Initial value for beta0 (default: 1)
#' @param init_nu_per Initial value for nu (default: 1)
#' @param init_gamma_per Initial value for gamma (default: 1)
#' @param init_theta_per Initial value for theta (default: 0.5)
#' @param SS_t0_per Initial temperature for simulated annealing (default: 10.0)
#' @param SS_t1_per Final temperature for simulated annealing (default: 1.0)
#' @param hyper_mu_beta0_per Hyperparameter for beta0 mean (default: 0)
#' @param hyper_sigma_beta0_per Hyperparameter for beta0 variance (default: 1e6)
#' @param hyper_mu_nu_per Hyperparameter for nu mean (default: 1)
#' @param hyper_sigma_nu_per Hyperparameter for nu variance (default: 1)
#' @param hyper_c_gamma_per Hyperparameter for gamma shape (default: 0.0001)
#' @param hyper_d_gamma_per Hyperparameter for gamma rate (default: 0.0001)
#' @param hyper_a_theta_per Hyperparameter for theta shape (default: 1)
#' @param hyper_b_theta_per Hyperparameter for theta shape (default: -1.0)
#' @param max_iter_per Maximum number of iterations (default: 100)
#' @param tol_per Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 0.001)
#' @return Permutation function value
#' @export
per_fun <- function(j_index,
                    wrapper_beta,
                    wrapper_beta0,
                    wrapper_nu,
                    wrapper_gamma,
                    dataXY,
                    init_beta_per = rep(1,8),
                    init_beta0_per = 1,
                    init_nu_per = 1,
                    init_gamma_per = 1,
                    init_theta_per = 0.5,
                    SS_t0_per = 10.0,
                    SS_t1_per = 1.0,
                    hyper_mu_beta0_per = 0,
                    hyper_sigma_beta0_per = 1e6,
                    hyper_mu_nu_per = 1,
                    hyper_sigma_nu_per = 1,
                    hyper_c_gamma_per = 0.0001,
                    hyper_d_gamma_per = 0.0001,
                    hyper_a_theta_per = 1,
                    hyper_b_theta_per = -1.0,
                    max_iter_per = 100,
                    tol_per = 1e-6,
                    add_correc_CiS = 0.001) {

  per_fun_cpp(
    j_index = j_index,
    wrapper_beta = wrapper_beta,
    wrapper_beta0 = wrapper_beta0,
    wrapper_nu = wrapper_nu,
    wrapper_gamma = wrapper_gamma,
    dataXY = dataXY,
    init_beta_per = init_beta_per,
    init_beta0_per = init_beta0_per,
    init_nu_per = init_nu_per,
    init_gamma_per = init_gamma_per,
    init_theta_per = init_theta_per,
    SS_t0_per = SS_t0_per,
    SS_t1_per = SS_t1_per,
    hyper_mu_beta0_per = hyper_mu_beta0_per,
    hyper_sigma_beta0_per = hyper_sigma_beta0_per,
    hyper_mu_nu_per = hyper_mu_nu_per,
    hyper_sigma_nu_per = hyper_sigma_nu_per,
    hyper_c_gamma_per = hyper_c_gamma_per,
    hyper_d_gamma_per = hyper_d_gamma_per,
    hyper_a_theta_per = hyper_a_theta_per,
    hyper_b_theta_per = hyper_b_theta_per,
    max_iter_per = max_iter_per,
    tol_per = tol_per,
    add_correc_CiS = add_correc_CiS
  )
}

#' Calculate CiS for a group of indices
#'
#' @param test_indices Vector of indices to test
#' @param beta_opt Optimized beta vector
#' @param beta0_opt Optimized beta0 value
#' @param nu_opt Optimized nu value
#' @param ga_opt Optimized gamma value
#' @param dataXY List containing data X and Y
#' @param add_correc_CiS Correction factor for CiS (default: 0.001)
#' @return Group CiS value
#' @export
CiS_group_fun <- function(test_indices,
                          beta_opt,
                          beta0_opt,
                          nu_opt,
                          ga_opt,
                          dataXY,
                          add_correc_CiS = 0.001) {

  CiS_group_fun_cpp(
    test_indices = test_indices,
    beta_opt = beta_opt,
    beta0_opt = beta0_opt,
    nu_opt = nu_opt,
    ga_opt = ga_opt,
    dataXY = dataXY,
    add_correc_CiS = add_correc_CiS
  )
}

#' Calculate permutation value for a group of indices
#'
#' @param j_indices Vector of indices for permutation
#' @param wrapper_beta Function wrapper for beta
#' @param wrapper_beta0 Function wrapper for beta0
#' @param wrapper_nu Function wrapper for nu
#' @param wrapper_gamma Function wrapper for gamma
#' @param dataXY List containing data X and Y
#' @param init_beta_per Initial vector for beta
#' @param init_beta0_per Initial value for beta0 (default: 1)
#' @param init_nu_per Initial value for nu (default: 1)
#' @param init_gamma_per Initial value for gamma (default: 1)
#' @param init_theta_per Initial value for theta (default: 0.5)
#' @param SS_t0_per Initial temperature for simulated annealing (default: 10.0)
#' @param SS_t1_per Final temperature for simulated annealing (default: 1.0)
#' @param hyper_mu_beta0_per Hyperparameter for beta0 mean (default: 0)
#' @param hyper_sigma_beta0_per Hyperparameter for beta0 variance (default: 1e6)
#' @param hyper_mu_nu_per Hyperparameter for nu mean (default: 1)
#' @param hyper_sigma_nu_per Hyperparameter for nu variance (default: 1)
#' @param hyper_c_gamma_per Hyperparameter for gamma shape (default: 0.0001)
#' @param hyper_d_gamma_per Hyperparameter for gamma rate (default: 0.0001)
#' @param hyper_a_theta_per Hyperparameter for theta shape (default: 1)
#' @param hyper_b_theta_per Hyperparameter for theta shape (default: -1.0)
#' @param max_iter_per Maximum number of iterations (default: 100)
#' @param tol_per Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 0.001)
#' @return Group permutation function value
#' @export
per_group_fun <- function(j_indices,
                          wrapper_beta,
                          wrapper_beta0,
                          wrapper_nu,
                          wrapper_gamma,
                          dataXY,
                          init_beta_per = rep(1,8),
                          init_beta0_per = 1,
                          init_nu_per = 1,
                          init_gamma_per = 1,
                          init_theta_per = 0.5,
                          SS_t0_per = 10.0,
                          SS_t1_per = 1.0,
                          hyper_mu_beta0_per = 0,
                          hyper_sigma_beta0_per = 1e6,
                          hyper_mu_nu_per = 1,
                          hyper_sigma_nu_per = 1,
                          hyper_c_gamma_per = 0.0001,
                          hyper_d_gamma_per = 0.0001,
                          hyper_a_theta_per = 1,
                          hyper_b_theta_per = -1.0,
                          max_iter_per = 100,
                          tol_per = 1e-6,
                          add_correc_CiS = 0.001) {

  per_group_fun_cpp(
    j_indices = j_indices,
    wrapper_beta = wrapper_beta,
    wrapper_beta0 = wrapper_beta0,
    wrapper_nu = wrapper_nu,
    wrapper_gamma = wrapper_gamma,
    dataXY = dataXY,
    init_beta_per = init_beta_per,
    init_beta0_per = init_beta0_per,
    init_nu_per = init_nu_per,
    init_gamma_per = init_gamma_per,
    init_theta_per = init_theta_per,
    SS_t0_per = SS_t0_per,
    SS_t1_per = SS_t1_per,
    hyper_mu_beta0_per = hyper_mu_beta0_per,
    hyper_sigma_beta0_per = hyper_sigma_beta0_per,
    hyper_mu_nu_per = hyper_mu_nu_per,
    hyper_sigma_nu_per = hyper_sigma_nu_per,
    hyper_c_gamma_per = hyper_c_gamma_per,
    hyper_d_gamma_per = hyper_d_gamma_per,
    hyper_a_theta_per = hyper_a_theta_per,
    hyper_b_theta_per = hyper_b_theta_per,
    max_iter_per = max_iter_per,
    tol_per = tol_per,
    add_correc_CiS = add_correc_CiS
  )
}

#' Calculate TVS values
#'
#' @param wrapper_beta Function wrapper for beta
#' @param wrapper_beta0 Function wrapper for beta0
#' @param wrapper_nu Function wrapper for nu
#' @param wrapper_gamma Function wrapper for gamma
#' @param dataXY List containing data X and Y
#' @param init_beta_TVS Initial vector for beta
#' @param B Number of iterations (default: 300)
#' @param init_beta0_TVS Initial value for beta0 (default: 1)
#' @param init_nu_TVS Initial value for nu (default: 1)
#' @param init_gamma_TVS Initial value for gamma (default: 1)
#' @param init_theta_TVS Initial value for theta (default: 0.5)
#' @param SS_t0_TVS Initial temperature for simulated annealing (default: 10.0)
#' @param SS_t1_TVS Final temperature for simulated annealing (default: 1.0)
#' @param hyper_mu_beta0_TVS Hyperparameter for beta0 mean (default: 0)
#' @param hyper_sigma_beta0_TVS Hyperparameter for beta0 variance (default: 1e6)
#' @param hyper_mu_nu_TVS Hyperparameter for nu mean (default: 1)
#' @param hyper_sigma_nu_TVS Hyperparameter for nu variance (default: 1)
#' @param hyper_c_gamma_TVS Hyperparameter for gamma shape (default: 0.0001)
#' @param hyper_d_gamma_TVS Hyperparameter for gamma rate (default: 0.0001)
#' @param hyper_a_theta_TVS Hyperparameter for theta shape (default: 1)
#' @param hyper_b_theta_TVS Hyperparameter for theta shape (default: -1.0)
#' @param max_iter_TVS Maximum number of iterations (default: 100)
#' @param tol_TVS Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 0.001)
#' @return Vector of TVS values
#' @export
TVS <- function(wrapper_beta,
                wrapper_beta0,
                wrapper_nu,
                wrapper_gamma,
                dataXY,
                init_beta_TVS = rep(1,8),
                B = 300,
                init_beta0_TVS = 1,
                init_nu_TVS = 1,
                init_gamma_TVS = 1,
                init_theta_TVS = 0.5,
                SS_t0_TVS = 10.0,
                SS_t1_TVS = 1.0,
                hyper_mu_beta0_TVS = 0,
                hyper_sigma_beta0_TVS = 1e6,
                hyper_mu_nu_TVS = 1,
                hyper_sigma_nu_TVS = 1,
                hyper_c_gamma_TVS = 0.0001,
                hyper_d_gamma_TVS = 0.0001,
                hyper_a_theta_TVS = 1,
                hyper_b_theta_TVS = -1.0,
                max_iter_TVS = 100,
                tol_TVS = 1e-6,
                add_correc_CiS = 0.001) {

  TVS_cpp(
    wrapper_beta = wrapper_beta,
    wrapper_beta0 = wrapper_beta0,
    wrapper_nu = wrapper_nu,
    wrapper_gamma = wrapper_gamma,
    dataXY = dataXY,
    init_beta_TVS = init_beta_TVS,
    B = B,
    init_beta0_TVS = init_beta0_TVS,
    init_nu_TVS = init_nu_TVS,
    init_gamma_TVS = init_gamma_TVS,
    init_theta_TVS = init_theta_TVS,
    SS_t0_TVS = SS_t0_TVS,
    SS_t1_TVS = SS_t1_TVS,
    hyper_mu_beta0_TVS = hyper_mu_beta0_TVS,
    hyper_sigma_beta0_TVS = hyper_sigma_beta0_TVS,
    hyper_mu_nu_TVS = hyper_mu_nu_TVS,
    hyper_sigma_nu_TVS = hyper_sigma_nu_TVS,
    hyper_c_gamma_TVS = hyper_c_gamma_TVS,
    hyper_d_gamma_TVS = hyper_d_gamma_TVS,
    hyper_a_theta_TVS = hyper_a_theta_TVS,
    hyper_b_theta_TVS = hyper_b_theta_TVS,
    max_iter_TVS = max_iter_TVS,
    tol_TVS = tol_TVS,
    add_correc_CiS = add_correc_CiS
  )
}

#' Multi-stage TVS algorithm
#'
#' @param wrapper_beta Function wrapper for beta
#' @param wrapper_beta0 Function wrapper for beta0
#' @param wrapper_nu Function wrapper for nu
#' @param wrapper_gamma Function wrapper for gamma
#' @param dataXY List containing data X and Y
#' @param init_beta_TVS Initial vector for beta
#' @param group_B Number of iterations for group stage (default: 20)
#' @param indiv_B Number of iterations for individual stage (default: 20)
#' @param B_final Number of iterations for final stage (default: 300)
#' @param group_cutoff Cutoff value for group stage (default: 3.0/20.0)
#' @param indiv_cutoff Cutoff value for individual stage (default: 3.0/20.0)
#' @param group_size Size of groups (default: 4)
#' @param init_beta0_TVS Initial value for beta0 (default: 1)
#' @param init_nu_TVS Initial value for nu (default: 1)
#' @param init_gamma_TVS Initial value for gamma (default: 1)
#' @param init_theta_TVS Initial value for theta (default: 0.5)
#' @param SS_t0_TVS Initial temperature for simulated annealing (default: 10.0)
#' @param SS_t1_TVS Final temperature for simulated annealing (default: 1.0)
#' @param hyper_mu_beta0_TVS Hyperparameter for beta0 mean (default: 0)
#' @param hyper_sigma_beta0_TVS Hyperparameter for beta0 variance (default: 1e6)
#' @param hyper_mu_nu_TVS Hyperparameter for nu mean (default: 1)
#' @param hyper_sigma_nu_TVS Hyperparameter for nu variance (default: 1)
#' @param hyper_c_gamma_TVS Hyperparameter for gamma shape (default: 0.0001)
#' @param hyper_d_gamma_TVS Hyperparameter for gamma rate (default: 0.0001)
#' @param hyper_a_theta_TVS Hyperparameter for theta shape (default: 1)
#' @param hyper_b_theta_TVS Hyperparameter for theta shape (default: -1.0)
#' @param max_iter_TVS Maximum number of iterations (default: 100)
#' @param tol_TVS Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 0.001)
#' @return List of multi-stage TVS results
#' @export
TVS_multi_stage <- function(wrapper_beta,
                            wrapper_beta0,
                            wrapper_nu,
                            wrapper_gamma,
                            dataXY,
                            init_beta_TVS = rep(1,8),
                            group_B = 20,
                            indiv_B = 20,
                            B_final = 300,
                            group_cutoff = 3.0/20.0,
                            indiv_cutoff = 3.0/20.0,
                            sig_cutoff = 0.05,
                            group_size = 4,
                            init_beta0_TVS = 1,
                            init_nu_TVS = 1,
                            init_gamma_TVS = 1,
                            init_theta_TVS = 0.5,
                            SS_t0_TVS = 10.0,
                            SS_t1_TVS = 1.0,
                            hyper_mu_beta0_TVS = 0,
                            hyper_sigma_beta0_TVS = 1e6,
                            hyper_mu_nu_TVS = 1,
                            hyper_sigma_nu_TVS = 1,
                            hyper_c_gamma_TVS = 0.0001,
                            hyper_d_gamma_TVS = 0.0001,
                            hyper_a_theta_TVS = 1,
                            hyper_b_theta_TVS = -1.0,
                            max_iter_TVS = 100,
                            tol_TVS = 1e-6,
                            add_correc_CiS = 0.001) {

  TVS_multi_stage_cpp(
    wrapper_beta = wrapper_beta,
    wrapper_beta0 = wrapper_beta0,
    wrapper_nu = wrapper_nu,
    wrapper_gamma = wrapper_gamma,
    dataXY = dataXY,
    init_beta_TVS = init_beta_TVS,
    group_B = group_B,
    indiv_B = indiv_B,
    B_final = B_final,
    group_cutoff = group_cutoff,
    indiv_cutoff = indiv_cutoff,
    sig_cutoff = sig_cutoff,
    group_size = group_size,
    init_beta0_TVS = init_beta0_TVS,
    init_nu_TVS = init_nu_TVS,
    init_gamma_TVS = init_gamma_TVS,
    init_theta_TVS = init_theta_TVS,
    SS_t0_TVS = SS_t0_TVS,
    SS_t1_TVS = SS_t1_TVS,
    hyper_mu_beta0_TVS = hyper_mu_beta0_TVS,
    hyper_sigma_beta0_TVS = hyper_sigma_beta0_TVS,
    hyper_mu_nu_TVS = hyper_mu_nu_TVS,
    hyper_sigma_nu_TVS = hyper_sigma_nu_TVS,
    hyper_c_gamma_TVS = hyper_c_gamma_TVS,
    hyper_d_gamma_TVS = hyper_d_gamma_TVS,
    hyper_a_theta_TVS = hyper_a_theta_TVS,
    hyper_b_theta_TVS = hyper_b_theta_TVS,
    max_iter_TVS = max_iter_TVS,
    tol_TVS = tol_TVS,
    add_correc_CiS = add_correc_CiS
  )
}
