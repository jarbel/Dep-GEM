
dissim_matrix = function(X){
  I=length(X)
  abs(X%o%rep(1,I)-rep(1,I)%o%X)
}

#   data frequency table
compute_data_freq=function(Y){
  I=dim(Y)[1]
  J=dim(Y)[2]
  data_freq=array(0,c(I,J,2))
  data_freq[,,1]=Y
  for (i in 1:I){
    data_freq[i,,2]=c(rev(cumsum(rev(Y[i,2:J]))),0)
  }
  data_freq
}


# ##################################################################
# GIBBS Sampler and sampling functions
# ##################################################################

#' @export
gibbs=function(n.iter=10, Y, X, cat="", GP="SE", show_pb=TRUE, burnin_coef = 4/5,
               sigma_Z_max = 5, sigma_Z_0 = 1, a_Z = 1, b_Z = 1,
               lambda_min = .01, lambda_max = .3, 
               lambda_0 = .05, a_lambda = 4, b_lambda = 20,
               M_min = 0, M_max = 50, a_M = 2, b_M = .2,
               d_0 = 1){
  Y <- as.matrix(Y)
  I=dim(Y)[1]
  J=dim(Y)[2]
  
  dist_X=dissim_matrix(X)
  dist2_X=dist_X^2
  
  data_freq=compute_data_freq(Y)
  
  Z_store=array(0,c(I,J,n.iter))
  V_store=array(0,c(I,J,n.iter))
  sigma_Z_store=rep(sigma_Z_0,n.iter)
  lambda_store=rep(lambda_0,n.iter)
  d_store=rep(d_0,n.iter)
  M_store=rep(15,n.iter)
  p_store=array(0,c(I,J,n.iter))
  D_store=array(0,c(I,n.iter))
  log_a_store=array(0,c(J,n.iter))
  
  sigma_Z = sigma_Z_store[1]
  lambda = lambda_store[1]
  M = M_store[1]
  d = d_store[1]
  K = K_fun(lambda, dist_X, dist2_X, GP, d)
  
  if (show_pb)  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  
  for(iter in 2:n.iter){
    if (show_pb)  setTxtProgressBar(pb, iter)
    
    for (j in 1:J){
      # update
      Z_list = sample_Z(Z_store[,j,iter-1], j, sigma_Z, M, K, data_freq)
      Z = Z_list$Z
      log_a_store[j,iter] = Z_list$log_a
      # storage
      Z_store[,j,iter]=Z
      V=Z_to_V(Z, sigma_Z=sigma_Z, M)
      V_store[,j,iter]=V
    }
    # update
    Z_matrix = Z_store[,,iter]
    sigma_Z = sample_sigma_Z(Z_matrix, sigma_Z, M, K, data_freq, 
                             sigma_Z_max, a_Z, b_Z)
    lambda = sample_lambda(Z_matrix, sigma_Z, M, lambda, K, data_freq, lambda_min, lambda_max,
                           dist_X, dist2_X, GP, a_lambda, b_lambda, d)
    K = K_fun(lambda, dist_X, dist2_X, GP, diag_off = FALSE, d)
    d = sample_d(Z_matrix, sigma_Z, M, lambda, K, data_freq, lambda_min, lambda_max,
                 dist_X, dist2_X, GP, d = d)
    K = K_fun(lambda, dist_X, dist2_X, GP, diag_off = FALSE, d)
    M = sample_M(Z_matrix, sigma_Z, M, data_freq, M_min, M_max, a_M, b_M)
    # storage
    sigma_Z_store[iter]=sigma_Z
    lambda_store[iter]=lambda
    M_store[iter]=M
    d_store[iter]=d
    p_store[,,iter]=t(apply(V_store[,,iter],1,V_to_p)) 
    D_store[,iter]=apply(p_store[,,iter],1,D_fun)
  }
  burnin = (round(n.iter*burnin_coef,0)+1):n.iter
  Z_store = Z_store[,,burnin]
  V_store = V_store[,,burnin]
  sigma_Z_store = sigma_Z_store[burnin]
  lambda_store = lambda_store[burnin]
  M_store = M_store[burnin]
  d_store = d_store[burnin]
  p_store = p_store[,,burnin]
  D_store = D_store[,burnin]
  log_a_store = log_a_store[,burnin]
  list("Z_store" = Z_store,
       "V_store" = V_store,
       "sigma_Z_store" = sigma_Z_store,
       "lambda_store" = lambda_store,
       "M_store" = M_store,
       "d_store" = d_store,
       "p_store" = p_store,
       "D_store" = D_store,
       "log_a_store" = log_a_store)
}


