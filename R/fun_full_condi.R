###########################
# def of sampling functions
###########################

MH_accept=function(Z_prop, Z, log_a){
  if (is.na(log_a)) 
    Z
  else
    if (log(runif(1))<log_a)
      Z_prop
  else
    Z
}

### 1. Z

sample_Z=function(Z, j, sigma_Z, M, K, data_freq){
  # proposal
  Z_prop = mvtnorm::rmvnorm(1, Z, .02*sigma_Z^2*K, method = "svd")
  V_prop = Z_to_V(Z_prop, sigma_Z, M)
  V      = Z_to_V(Z, sigma_Z, M)
  # acceptance proba
  log_a = log_like_V(V_prop, j, data_freq) - log_like_V(V, j, data_freq) +
    log_prior_Z(Z_prop, sigma_Z, K) - log_prior_Z(Z, sigma_Z, K)
  # choice
  MH = MH_accept(Z_prop, Z, log_a)
  list(Z = MH, log_a = log_a)
}

### 2. sigma_Z

log_prior_sigma_Z = function(sigma_Z, a_Z = 1, b_Z = 1, IG = FALSE){
  # Inverse gamma prior on sigma_Z^2
  if (IG)
    dgamma(sigma_Z^(-2), a_Z, b_Z, log = TRUE) - 3*sigma_Z + log(2)
  else
    dgamma(sigma_Z, a_Z, b_Z, log = TRUE)
}

sample_sigma_Z = function(Z_matrix, sigma_Z, M, K, data_freq, sigma_Z_max,
                          a_Z, b_Z){
  # symmetric proposal
  sigma_Z_prop = runif(1,0.9*sigma_Z, 1.1*sigma_Z)
  # uniform proposal
  #   sigma_Z_prop = runif(1,0,sigma_Z_max)  
  llZ = log_like_Z_total(Z_matrix, sigma_Z, M, data_freq)
  llZ_prop = log_like_Z_total(Z_matrix, sigma_Z_prop, M, data_freq)
  lpZ = log_prior_Z_total(Z_matrix, sigma_Z, K)
  lpZ_prop = log_prior_Z_total(Z_matrix, sigma_Z_prop, K)
  log_a = llZ_prop + lpZ_prop - llZ - lpZ +
    log_prior_sigma_Z(sigma_Z_prop, a_Z, b_Z) - 
    log_prior_sigma_Z(sigma_Z, a_Z, b_Z)
  MH_accept(sigma_Z_prop, sigma_Z, log_a)  
}

### 3.1 lambda

log_prior_lambda = function(lambda, a_lambda = 10, b_lambda = 1/2, IG = TRUE){
  # Inverse gamma prior on lambda
  if (IG)
    dgamma(lambda^(-1), a_lambda, b_lambda, log = TRUE) - 2*lambda
  # gamma prior on lambda
  else
    dgamma(lambda, a_lambda, b_lambda, log = TRUE)  
}


sample_lambda = function(Z_matrix, sigma_Z, M, lambda, K, data_freq, lambda_min, lambda_max,
                         dist_X, dist2_X, GP, a_lambda, b_lambda, sigma_n){
  # symmetric proposal
  lambda_prop = runif(1, .9*lambda, 1.1*lambda)
  # uniform proposal
  #   lambda_prop = runif(1, lambda_min, lambda_max)
  K_prop = K_fun(lambda_prop, dist_X, dist2_X, GP, sigma_n)
  # accept rate
  lpZ = log_prior_Z_total(Z_matrix, sigma_Z, K)
  lpZ_prop = log_prior_Z_total(Z_matrix, sigma_Z, K_prop)
  log_a = 
    lpZ_prop - lpZ + 
    log_prior_lambda(lambda_prop, a_lambda, b_lambda) - 
    log_prior_lambda(lambda, a_lambda, b_lambda)
  MH_accept(lambda_prop, lambda, log_a)  
}

### 3.2 d

log_prior_d = function(d, a_d = 1, b_d = 1){
  # gamma prior on lambda
  dgamma(d, a_d, b_d, log = TRUE)  
}


sample_d = function(Z_matrix, sigma_Z, M, lambda, K, data_freq, lambda_min, lambda_max,
                    dist_X, dist2_X, GP, d){
  
  # symmetric proposal
  d_prop = runif(1, .9*d, 1.1*d)
  # uniform proposal
  #   d_prop = runif(1, d_min, d_max)
  K_prop = K_fun(lambda, dist_X, dist2_X, GP, d)
  # accept rate
  lpZ = log_prior_Z_total(Z_matrix, sigma_Z, K)
  lpZ_prop = log_prior_Z_total(Z_matrix, sigma_Z, K_prop)
  log_a = 
    lpZ_prop - lpZ + 
    log_prior_d(d_prop) - 
    log_prior_d(d)
  MH_accept(d_prop, d, log_a)  
}

### 4. M

log_prior_M = function(M, a_M, b_M){
  # gamma prior on M
  dgamma(M, a_M, b_M, log = TRUE)
}

sample_M = function(Z_matrix, sigma_Z, M, data_freq, M_min, M_max, a_M, b_M){
  # symmetric prop
  M_prop = runif(1, .9*M, 1.1*M)
  # unif prop
  #   M_prop = runif(1, M_min, M_max)
  # accept rate
  llZ = log_like_Z_total(Z_matrix, sigma_Z, M, data_freq)
  llZ_prop = log_like_Z_total(Z_matrix, sigma_Z, M_prop, data_freq)
  log_a = llZ_prop - llZ + log_prior_M(M_prop, a_M, b_M) - log_prior_M(M, a_M, b_M)
  MH_accept(M_prop, M, log_a)
}




