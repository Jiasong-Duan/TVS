#' Dispatcher for beta optimization
#'
#' @param method Estimation function for beta (default: Quasi-Newton method implemented by R optim() function)
#' @return Optimized beta parameter
#' @export
wrapper_beta <- function(beta, beta0, nu, gamma, betaPRE, t0, t1, Y_lk, X_lk, theta_lk, method = getOption("TVS_beta_method", "maxLik")) {
  method <- match.arg(method, c("optim", "maxLik", "cd_maxLik", "cd_maxLik_R"))
  if (method == "optim") {
    wrapper_beta_optim(beta, beta0, nu, gamma, betaPRE, t0, t1, Y_lk, X_lk, theta_lk)
  } else if (method == "maxLik") {
    wrapper_beta_maxLik(beta, beta0, nu, gamma, betaPRE, t0, t1, Y_lk, X_lk, theta_lk)
  } else if (method == "cd_maxLik"){
    beta_coordinate_descent_cpp_maxLik(beta_cd=beta, beta0_cd=beta0, nu_cd=nu, ga_cd=gamma,
                                       betaPRE, t0, t1, Y_cd=Y_lk, X_cd=X_lk, theta_cd=theta_lk,
                                       maX_cd_iter=2, tol=1e-6)
  } else if (method == "cd_maxLik_R"){
    beta_coordinate_descent_R_maxLik(beta_cd=beta, beta0_cd=beta0, nu_cd=nu, ga_cd=gamma,
                                     betaPRE, t0, t1, Y_cd=Y_lk, X_cd=X_lk, theta_cd=theta_lk,
                                     maX_cd_iter=2, tol=1e-6)
  }
}


#' Wrapper function for beta optimization using quasi-Newton method
#'
#' @param beta Initial beta parameter
#' @param beta0 Beta0 parameter
#' @param nu Nu parameter
#' @param gamma Gamma parameter
#' @param betaPRE Previous beta value
#' @param t0 the "spike" hyperparameter in beta's spike-and-slab prior
#' @param t1 the "slab" hyperparameter in beta's spike-and-slab prior
#' @param Y_lk Response vector
#' @param X_lk Predictor matrix
#' @param theta_lk Theta parameter (default: number of predictors)
#' @return Optimized beta parameter
#' @export
wrapper_beta_optim <- function(beta, beta0, nu, gamma, betaPRE, t0, t1, Y_lk, X_lk, theta_lk) {
  func_lk <- function(x) beta_neg_lk_cpp(beta_lk=x, beta0_lk=beta0, nu_lk=nu, ga_lk=gamma, betaPRE=betaPRE,
                                         t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk)
  func_gradient <- function(x) beta_neg_gradient_cpp(beta_lk=x, beta0_lk=beta0, nu_lk=nu, ga_lk=gamma, betaPRE=betaPRE,
                                                     t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk)
  result <- tryCatch(
    optim(par = beta, fn = func_lk, gr= func_gradient, method = "BFGS"),
    error = function(e) {
      warning("beta optimization failed: ", conditionMessage(e))
      list(par = NA)
    }
  )
  return(result$par)
}

#' Wrapper function for beta optimization using Newton–Raphson method
#'
#' @param beta Initial beta parameter
#' @param beta0 Beta0 parameter
#' @param nu Nu parameter
#' @param gamma Gamma parameter
#' @param betaPRE Previous beta value
#' @param t0 the "spike" hyperparameter in beta's spike-and-slab prior
#' @param t1 the "slab" hyperparameter in beta's spike-and-slab prior
#' @param Y_lk Response vector
#' @param X_lk Predictor matrix
#' @param theta_lk Theta parameter (default: number of predictors)
#' @return Optimized beta parameter
#' @export
wrapper_beta_maxLik <- function(beta, beta0, nu, gamma, betaPRE, t0, t1, Y_lk, X_lk, theta_lk) {
  func_lk <- function(x) -beta_neg_lk_cpp(beta_lk=x, beta0_lk=beta0, nu_lk=nu, ga_lk=gamma, betaPRE=betaPRE,
                                          t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk)
  func_gradient <- function(x) -t(beta_neg_gradient_cpp(beta_lk=x, beta0_lk=beta0, nu_lk=nu, ga_lk=gamma, betaPRE=betaPRE,
                                                        t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk))
  func_Hessian <- function(x) -beta_neg_hessian_cpp(beta_lk=x, beta0_lk=beta0, nu_lk=nu, ga_lk=gamma, betaPRE=betaPRE,
                                                        t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk)
  result <- tryCatch(
    maxLik::maxLik(logLik = func_lk, grad = func_gradient, hess = func_Hessian, start = beta),
    error = function(e) {
      warning("maxLik() failed: ", conditionMessage(e))
      list(estimate = NA)
    }
  )
  return(result$estimate)
}

#' Wrapper function for beta optimization using coordinate ascend based on Newton–Raphson update
#'
#' @param beta_cd Initial beta parameter
#' @param beta0_cd Beta0 parameter
#' @param nu_cd Nu parameter
#' @param ga_cd Gamma parameter
#' @param betaPRE Previous beta value
#' @param t0 the "spike" hyperparameter in beta's spike-and-slab prior
#' @param t1 the "slab" hyperparameter in beta's spike-and-slab prior
#' @param Y_cd Response vector
#' @param X_cd Predictor matrix
#' @param theta_cd Theta parameter (default: number of predictors)
#' @return Optimized beta parameter
#' @export
beta_coordinate_descent_R_maxLik <- function(beta_cd, beta0_cd, nu_cd, ga_cd, betaPRE, t0, t1,
                                             Y_cd, X_cd, theta_cd, maX_cd_iter, tol) {
  # Get p
  p <- ncol(X_cd)
  # Update beta
  for (iter in 1: maX_cd_iter) {
    beta_old <- beta_cd
    # Update beta sequentially
    for (j in 1:p) {
      # Newton update using maxLik
      beta_cd[j] <- wrapper_beta_cd(j, beta_cd, beta0_cd, nu_cd, ga_cd, betaPRE, t0, t1, Y_cd, X_cd, theta_cd)
    }
    # Check convergence
    if (sqrt(sum((beta_cd - beta_old)^2)) < tol)
      break
  }
  # Return updated beta
  return(beta_cd);
}

#' Wrapper function for beta optimization using coordinate ascend Newton–Raphson method
#'
#' @param j the index of beta to update
#' @param beta Initial beta parameter
#' @param beta0 Beta0 parameter
#' @param nu Nu parameter
#' @param gamma Gamma parameter
#' @param betaPRE Previous beta value
#' @param t0 the "spike" hyperparameter in beta's spike-and-slab prior
#' @param t1 the "slab" hyperparameter in beta's spike-and-slab prior
#' @param Y_lk Response vector
#' @param X_lk Predictor matrix
#' @param theta_lk Theta parameter (default: number of predictors)
#' @return Optimized beta parameter
#' @export
wrapper_beta_cd <- function(j, beta, beta0, nu, gamma, betaPRE, t0, t1, Y_lk, X_lk, theta_lk){
  func_lk <- function(x) -jbeta_neg_lk_cpp_maxLik(beta_j=x, j_index = j, beta_noj = beta[-j],
                                                  beta0_lk=beta0, ga_lk=gamma, nu_lk=nu, betaPRE=betaPRE,
                                                  t0=t0, t1=t1, Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk)
  func_gradient <- function(x) -jbeta_neg_gradient_cpp_maxLik(beta_j=x, j_index = j, beta_noj = beta[-j],
                                                              beta0_lk = beta0, ga_lk = gamma, nu_lk = nu,
                                                              betaPRE = betaPRE, t0 = t0, t1 = t1,
                                                              Y_lk=Y_lk, X_lk=X_lk, theta_lk=theta_lk)
  func_Hessian <- function(x) -jbeta_neg_hessian_cpp_maxLik(beta_j=x, j_index = j, beta_noj = beta[-j],
                                                            beta0_lk = beta0, ga_lk = gamma, nu_lk = nu,
                                                            Y_lk=Y_lk, X_lk=X_lk)
  result <- tryCatch(
    maxLik::maxLik(logLik = func_lk, grad = func_gradient, hess = func_Hessian, start = beta[j]),
    error = function(e) {
      warning("maxLik() failed: ", conditionMessage(e))
      list(estimate = NA)
    }
  )
  return(result$estimate)
}

#' Wrapper function for beta0 optimization
#'
#' @param beta0.lk Initial beta0 parameter
#' @param beta.lk Beta parameter
#' @param nu.lk Nu parameter
#' @param gamma.lk Gamma parameter
#' @param hyper.mu.beta0 the mean hyperparameter in beta0's Normal prior
#' @param hyper.sigma.beta0 the variance hyperparameter in beta0's Normal prior
#' @param Y.lk Response vector
#' @param X.lk Predictor matrix
#' @return Optimized beta0 parameter
#' @export
wrapper_beta0 <- function(beta0.lk, beta.lk, nu.lk, gamma.lk, hyper.mu.beta0, hyper.sigma.beta0, Y.lk, X.lk){
  func <- function(x) beta0_neg_lk_cpp(beta0_lk=x, beta_lk=beta.lk, nu_lk=nu.lk, gamma_lk=gamma.lk,
                                       hyper_mu_beta0= hyper.mu.beta0, hyper_sigma_beta0= hyper.sigma.beta0, Y_lk=Y.lk, X_lk=X.lk)
  result <- tryCatch(
    optim(par = beta0.lk, fn = func, method = "BFGS"),
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
#' @param hyper_mu the location hyperparameter in nu's log-normal prior
#' @param hyper_sigma the scale hyperparameter in nu's log-normal prior
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
#' @param hyper_c the shape hyperparameter in gamma's gamma prior
#' @param hyper_d the rate hyperparameter in gamma's gamma prior
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
#' @param dataXY List containing data X and Y
#' @param init_beta Initial vector for beta (default: rep(1, number of predictors))
#' @param init_beta0 Initial value for beta0 (default: 1)
#' @param init_nu Initial value for nu (default: 1)
#' @param init_gamma Initial value for gamma (default: 1)
#' @param init_theta Initial value for theta (default: 0.5)
#' @param SS_t0 the "spike" hyperparameter in beta's spike-and-slab prior (default: 10)
#' @param SS_t1 the "slab" hyperparameter in beta's spike-and-slab prior (default: 1)
#' @param hyper_mu_beta0 the mean hyperparameter in beta0's Normal prior (default: 0)
#' @param hyper_sigma_beta0 the variance hyperparameter in beta0's Normal prior (default: 1e6)
#' @param hyper_mu_nu the location hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_sigma_nu the scale hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_c_gamma the shape hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_d_gamma the rate hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_a_theta the shape hyperparameter a in theta's beta prior (default: 1)
#' @param hyper_b_theta the shape hyperparameter b in theta's beta prior (default: the number of predictors)
#' @param max_iter Maximum number of iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param conv_type Convergence type (default: "param")
#' @param beta_method Optimization method for beta updates ("maxLik", "optim", or "cd_maxLik"). Defaults to "maxLik".
#' @return List of results from the EM algorithm
#' @export
TVS_EM <- function(
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
    conv_type = "param",
    beta_method = c("maxLik", "optim", "cd_maxLik", "cd_maxLik_R")) {

  #Pick the estimation method for beta estimation
  beta.method <- match.arg(beta_method)
  # set temporarily for dispatcher
  old <- options(TVS_beta_method = beta.method)
  on.exit(options(old), add = TRUE)

  TVS_EM_cpp(
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
#' @return Vector of slope of errors
#' @export
slope_error <- function(err, nu_slo, ga_slo) {
  slope_error_cpp(err = err, nu_slo = nu_slo, ga_slo = ga_slo)
}

#' Calculate curvature error
#'
#' @param err Error vector
#' @param nu_cur Nu parameter for curvature
#' @param ga_cur Gamma parameter for curvature
#' @return Vector of curvature of errors
#' @export
curvature_error <- function(err, nu_cur, ga_cur) {
  curvature_error_cpp(err = err, nu_cur = nu_cur, ga_cur = ga_cur)
}

#' Calculate CiS for a single index
#'
#' @param test_index Index of predictor to test
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

#' Calculate permutation value of CiS for a single index
#'
#' @param j_index Index of predictor for permutation
#' @param dataXY List containing data X and Y
#' @param init_beta_per Initial vector for beta (default: rep(1, number of predictors))
#' @param init_beta0_per Initial value for beta0 (default: 1)
#' @param init_nu_per Initial value for nu (default: 1)
#' @param init_gamma_per Initial value for gamma (default: 1)
#' @param init_theta_per Initial value for theta (default: 0.5)
#' @param SS_t0_per the "spike" hyperparameter in beta's spike-and-slab prior (default: 10)
#' @param SS_t1_per the "slab" hyperparameter in beta's spike-and-slab prior (default: 1)
#' @param hyper_mu_beta0_per the mean hyperparameter in beta0's Normal prior (default: 0)
#' @param hyper_sigma_beta0_per the variance hyperparameter in beta0's Normal prior (default: 1e6)
#' @param hyper_mu_nu_per the location hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_sigma_nu_per the scale hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_c_gamma_per the shape hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_d_gamma_per the rate hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_a_theta_per the shape hyperparameter a in theta's beta prior (default: 1)
#' @param hyper_b_theta_per the shape hyperparameter b in theta's beta prior (default: the number of predictors)
#' @param max_iter_per Maximum number of iterations (default: 100)
#' @param tol_per Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 1e-3)
#' @param beta_method Optimization method for beta updates ("maxLik", "optim", or "cd_maxLik"). Defaults to "maxLik".
#' @return Permutation function value
#' @export
per_fun <- function(j_index,
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
                    add_correc_CiS = 0.001,
                    beta_method = c("maxLik", "optim", "cd_maxLik", "cd_maxLik_R")) {

  #Pick the estimation method for beta estimation
  beta.method <- match.arg(beta_method)
  # set temporarily for dispatcher
  old <- options(TVS_beta_method = beta.method)
  on.exit(options(old), add = TRUE)

  per_fun_cpp(
    j_index = j_index,
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
#' @param test_indices Vector of indices of predictors to test
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

#' Calculate permutation value of CiS for a group of indices
#'
#' @param j_indices Vector of indices of predictors for permutation
#' @param dataXY List containing data X and Y
#' @param init_beta_per Initial vector for beta (default: rep(1, number of predictors))
#' @param init_beta0_per Initial value for beta0 (default: 1)
#' @param init_nu_per Initial value for nu (default: 1)
#' @param init_gamma_per Initial value for gamma (default: 1)
#' @param init_theta_per Initial value for theta (default: 0.5)
#' @param SS_t0_per the "spike" hyperparameter in beta's spike-and-slab prior (default: 10)
#' @param SS_t1_per the "slab" hyperparameter in beta's spike-and-slab prior (default: 1)
#' @param hyper_mu_beta0_per the mean hyperparameter in beta0's Normal prior (default: 0)
#' @param hyper_sigma_beta0_per the variance hyperparameter in beta0's Normal prior (default: 1e6)
#' @param hyper_mu_nu_per the location hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_sigma_nu_per the scale hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_c_gamma_per the shape hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_d_gamma_per the rate hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_a_theta_per the shape hyperparameter a in theta's beta prior (default: 1)
#' @param hyper_b_theta_per the shape hyperparameter b in theta's beta prior (default: the number of predictors)
#' @param max_iter_per Maximum number of iterations (default: 100)
#' @param tol_per Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 1e-3)
#' @param beta_method Optimization method for beta updates ("maxLik", "optim", or "cd_maxLik"). Defaults to "maxLik".
#' @return Group permutation function value
#' @export
per_group_fun <- function(j_indices,
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
                          add_correc_CiS = 0.001,
                          beta_method = c("maxLik", "optim", "cd_maxLik", "cd_maxLik_R")) {

  #Pick the estimation method for beta estimation
  beta.method <- match.arg(beta_method)
  # set temporarily for dispatcher
  old <- options(TVS_beta_method = beta.method)
  on.exit(options(old), add = TRUE)

  per_group_fun_cpp(
    j_indices = j_indices,
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

#' TVS variable selection without Multi-stage pre-screening process
#'
#' @param dataXY List containing data X and Y
#' @param init_beta_TVS Initial vector for beta (default: rep(1, number of predictors))
#' @param B Number of iterations (default: 300)
#' @param sig_cutoff Threshold used for significance cutoff in variable selection (default: 0.05)
#' @param init_beta0_TVS Initial value for beta0 (default: 1)
#' @param init_nu_TVS Initial value for nu (default: 1)
#' @param init_gamma_TVS Initial value for gamma (default: 1)
#' @param init_theta_TVS Initial value for theta (default: 0.5)
#' @param SS_t0_TVS the "spike" hyperparameter in beta's spike-and-slab prior (default: 10)
#' @param SS_t1_TVS the "slab" hyperparameter in beta's spike-and-slab prior (default: 1)
#' @param hyper_mu_beta0_TVS the mean hyperparameter in beta0's Normal prior (default: 0)
#' @param hyper_sigma_beta0_TVS the variance hyperparameter in beta0's Normal prior (default: 1e6)
#' @param hyper_mu_nu_TVS the location hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_sigma_nu_TVS the scale hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_c_gamma_TVS the shape hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_d_gamma_TVS the rate hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_a_theta_TVS the shape hyperparameter a in theta's beta prior (default: 1)
#' @param hyper_b_theta_TVS the shape hyperparameter b in theta's beta prior (default: the number of predictors)
#' @param max_iter_TVS Maximum number of iterations (default: 100)
#' @param tol_TVS Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 1e-3)
#' @param beta_method Optimization method for beta updates ("maxLik", "optim", or "cd_maxLik"). Defaults to "maxLik".
#' @return List containing variable selection results
#' @export
TVS <- function(
    dataXY,
    init_beta_TVS = rep(1,8),
    B = 300,
    sig_cutoff = 0.05,
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
    add_correc_CiS = 0.001,
    beta_method = c("maxLik", "optim", "cd_maxLik", "cd_maxLik_R")) {

  #Pick the estimation method for beta estimation
  beta.method <- match.arg(beta_method)
  # set temporarily for dispatcher
  old <- options(TVS_beta_method = beta.method)
  on.exit(options(old), add = TRUE)

  TVS_cpp(
    dataXY = dataXY,
    init_beta_TVS = init_beta_TVS,
    B = B,
    sig_cutoff = sig_cutoff,
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

#' TVS variable selection for a single predictor j
#'
#' @param test_index Index of predictor to test
#' @param dataXY List containing data X and Y
#' @param init_beta_TVS Initial vector for beta (default: rep(1, number of predictors))
#' @param B Number of iterations (default: 300)
#' @param sig_cutoff Threshold used for significance cutoff in variable selection (default: 0.05)
#' @param init_beta0_TVS Initial value for beta0 (default: 1)
#' @param init_nu_TVS Initial value for nu (default: 1)
#' @param init_gamma_TVS Initial value for gamma (default: 1)
#' @param init_theta_TVS Initial value for theta (default: 0.5)
#' @param SS_t0_TVS the "spike" hyperparameter in beta's spike-and-slab prior (default: 10)
#' @param SS_t1_TVS the "slab" hyperparameter in beta's spike-and-slab prior (default: 1)
#' @param hyper_mu_beta0_TVS the mean hyperparameter in beta0's Normal prior (default: 0)
#' @param hyper_sigma_beta0_TVS the variance hyperparameter in beta0's Normal prior (default: 1e6)
#' @param hyper_mu_nu_TVS the location hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_sigma_nu_TVS the scale hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_c_gamma_TVS the shape hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_d_gamma_TVS the rate hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_a_theta_TVS the shape hyperparameter a in theta's beta prior (default: 1)
#' @param hyper_b_theta_TVS the shape hyperparameter b in theta's beta prior (default: the number of predictors)
#' @param max_iter_TVS Maximum number of iterations (default: 100)
#' @param tol_TVS Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 1e-3)
#' @param beta_method Optimization method for beta updates ("maxLik", "optim", or "cd_maxLik"). Defaults to "maxLik".
#' @return List containing variable selection results
#' @export
TVS_j <- function(
    test_index,
    dataXY,
    init_beta_TVS = rep(1,8),
    B = 300,
    sig_cutoff = 0.05,
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
    add_correc_CiS = 0.001,
    beta_method = c("maxLik", "optim", "cd_maxLik", "cd_maxLik_R")) {

  #Pick the estimation method for beta estimation
  beta.method <- match.arg(beta_method)
  # set temporarily for dispatcher
  old <- options(TVS_beta_method = beta.method)
  on.exit(options(old), add = TRUE)

  TVS_j_cpp(
    test_index = test_index,
    dataXY = dataXY,
    init_beta_TVS = init_beta_TVS,
    B = B,
    sig_cutoff = sig_cutoff,
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

#' TVS variable selection for a group of predictors
#'
#' @param test_indices Indices of predictors to test
#' @param dataXY List containing data X and Y
#' @param init_beta_TVS Initial vector for beta (default: rep(1, number of predictors))
#' @param B Number of iterations (default: 300)
#' @param sig_cutoff Threshold used for significance cutoff in variable selection (default: 0.05)
#' @param init_beta0_TVS Initial value for beta0 (default: 1)
#' @param init_nu_TVS Initial value for nu (default: 1)
#' @param init_gamma_TVS Initial value for gamma (default: 1)
#' @param init_theta_TVS Initial value for theta (default: 0.5)
#' @param SS_t0_TVS the "spike" hyperparameter in beta's spike-and-slab prior (default: 10)
#' @param SS_t1_TVS the "slab" hyperparameter in beta's spike-and-slab prior (default: 1)
#' @param hyper_mu_beta0_TVS the mean hyperparameter in beta0's Normal prior (default: 0)
#' @param hyper_sigma_beta0_TVS the variance hyperparameter in beta0's Normal prior (default: 1e6)
#' @param hyper_mu_nu_TVS the location hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_sigma_nu_TVS the scale hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_c_gamma_TVS the shape hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_d_gamma_TVS the rate hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_a_theta_TVS the shape hyperparameter a in theta's beta prior (default: 1)
#' @param hyper_b_theta_TVS the shape hyperparameter b in theta's beta prior (default: the number of predictors)
#' @param max_iter_TVS Maximum number of iterations (default: 100)
#' @param tol_TVS Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 1e-3)
#' @param beta_method Optimization method for beta updates ("maxLik", "optim", or "cd_maxLik"). Defaults to "maxLik".
#' @return List containing variable selection results
#' @export
TVS_group <- function(
    test_indices,
    dataXY,
    init_beta_TVS = rep(1,8),
    B = 300,
    sig_cutoff = 0.05,
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
    add_correc_CiS = 0.001,
    beta_method = c("maxLik", "optim", "cd_maxLik", "cd_maxLik_R")) {

  #Pick the estimation method for beta estimation
  beta.method <- match.arg(beta_method)
  # set temporarily for dispatcher
  old <- options(TVS_beta_method = beta.method)
  on.exit(options(old), add = TRUE)

  TVS_group_cpp(
    test_indices = test_indices,
    dataXY = dataXY,
    init_beta_TVS = init_beta_TVS,
    B = B,
    sig_cutoff = sig_cutoff,
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

#' TVS variable selection with Multi-stage pre-screening process
#'
#' @param dataXY List containing data X and Y
#' @param init_beta_TVS Initial vector for beta (default: rep(1, number of predictors))
#' @param group_B Number of iterations for group stage (default: 20)
#' @param indiv_B Number of iterations for individual stage (default: 20)
#' @param B_final Number of iterations for final stage (default: 300)
#' @param group_cutoff Cutoff value for group stage (default: 3.0/20.0)
#' @param indiv_cutoff Cutoff value for individual stage (default: 3.0/20.0)
#' @param sig_cutoff Threshold used for significance cutoff in variable selection (default: 0.05)
#' @param group_size Size of groups (default: 4)
#' @param init_beta0_TVS Initial value for beta0 (default: 1)
#' @param init_nu_TVS Initial value for nu (default: 1)
#' @param init_gamma_TVS Initial value for gamma (default: 1)
#' @param init_theta_TVS Initial value for theta (default: 0.5)
#' @param SS_t0_TVS the "spike" hyperparameter in beta's spike-and-slab prior (default: 10)
#' @param SS_t1_TVS the "slab" hyperparameter in beta's spike-and-slab prior (default: 1)
#' @param hyper_mu_beta0_TVS the mean hyperparameter in beta0's Normal prior (default: 0)
#' @param hyper_sigma_beta0_TVS the variance hyperparameter in beta0's Normal prior (default: 1e6)
#' @param hyper_mu_nu_TVS the location hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_sigma_nu_TVS the scale hyperparameter in nu's log-normal prior (default: 1)
#' @param hyper_c_gamma_TVS the shape hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_d_gamma_TVS the rate hyperparameter in gamma's gamma prior (default: 1e-4)
#' @param hyper_a_theta_TVS the shape hyperparameter a in theta's beta prior (default: 1)
#' @param hyper_b_theta_TVS the shape hyperparameter b in theta's beta prior (default: the number of predictors)
#' @param max_iter_TVS Maximum number of iterations (default: 100)
#' @param tol_TVS Convergence tolerance (default: 1e-6)
#' @param add_correc_CiS Correction factor for CiS (default: 1e-3)
#' @param beta_method Optimization method for beta updates ("maxLik", "optim", or "cd_maxLik"). Defaults to "maxLik".
#' @return List of multi-stage TVS results
#' @export
TVS_multi_stage <- function(
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
    add_correc_CiS = 0.001,
    beta_method = c("maxLik", "optim", "cd_maxLik", "cd_maxLik_R")) {

  #Pick the estimation method for beta estimation
  beta.method <- match.arg(beta_method)
  # set temporarily for dispatcher
  old <- options(TVS_beta_method = beta.method)
  on.exit(options(old), add = TRUE)

  TVS_multi_stage_cpp(
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
