###########################################################
# Predictive functions
###########################################################

# 1. useful functions
symmetrize=function(M){(M+t(M))/2}
test_sym=function(M){prod(M==symmetrize(M))}
dissim_matrix2 = function(X, Xs){
  I=length(X)
  Is=length(Xs)
  abs(Xs%o%rep(1,I)-rep(1,Is)%o%X)
}

q_025 = function(D){quantile(D,probs = .025)}
q_975 = function(D){quantile(D,probs = .975)}
q_1 = function(D){quantile(D,probs = .1)}
q_9 = function(D){quantile(D,probs = .9)}

# 2. main function
#' @export
predictive=function(Xs, X, cat="", GP = "SE",
                    Z_store, lambda_store, sigma_Z_store, M_store, d_store,
                    D_store, 
                    by_burn = 4, d_bar = .0001, d_bar_s = 1){
  subset = seq(1,length(M_store), by = by_burn)
  Z_store = Z_store[,,subset]
  sigma_Z_store = sigma_Z_store[subset]
  lambda_store = lambda_store[subset]
  M_store = M_store[subset]
  D_store = D_store[,subset]
  
  n.iter = dim(Z_store)[3]
  I=dim(Z_store)[1]
  Is = length(Xs)
  J=dim(Z_store)[2]
  
  dist_X=dissim_matrix(X)
  dist2_X=dist_X^2
  
  dist_Xs=dissim_matrix(Xs)
  dist2_Xs=dist_Xs^2
  
  dist_Xs_X=dissim_matrix2(X,Xs)
  dist2_Xs_X=dist_Xs_X^2
  
  Z_star=array(0,dim=c(Is,J,n.iter))
  V_star=array(0,dim=c(Is,J,n.iter))
  p_star=array(0,dim=c(Is,J,n.iter))
  D_star=array(50,dim=c(Is,n.iter))
  
  pb <- txtProgressBar(min = 0, max = n.iter, style = 3)
  for (iter in 1:n.iter){
    # update progress bar
    setTxtProgressBar(pb, iter)
    # compute parameters of predictive distribution, mean and covariance matrix
    lambda = lambda_store[iter]
    sigma_Z = sigma_Z_store[iter]
    M = M_store[iter]
    d = d_store[iter]
    K_X = K_fun(lambda = lambda, R = dist_X, R2 = dist2_X, GP=GP, d=d*d_bar, diag_off = FALSE)
    K_Xs = K_fun(lambda = lambda, R = dist_Xs, R2 = dist2_Xs, GP=GP, d=d*d_bar_s, diag_off = FALSE) 
    K_Xs_X = K_fun(lambda = lambda, R = dist_Xs_X, R2 = dist2_Xs_X, GP=GP, d=d, diag_off = TRUE)
    K_inv_mat = solve(K_X)
    K_inv_mat = symmetrize(K_inv_mat)
    K_Xs_K_inv = K_Xs_X%*%K_inv_mat
    K_star = K_Xs-K_Xs_K_inv%*%t(K_Xs_X)
    K_star=symmetrize(K_star)
    m_star_fun=function(Z){
      K_Xs_K_inv%*%Z
    }
    for (j in 1:J){
      Zs = mvtnorm::rmvnorm(1, mean = c(m_star_fun(Z_store[,j,iter])), 
                   sigma = sigma_Z^2*K_star, method="svd")
      Z_star[,j,iter] = Zs
      Vs = Z_to_V(Zs, sigma_Z = sigma_Z, M)
      V_star[,j,iter] = Vs
    }
    p_star[,,iter]=t(apply(V_star[,,iter],1,V_to_p))
    D_star[,iter] = apply(p_star[,,iter],1,D_fun)
  }
  list("Xs" = Xs, 
       "Z_star" = Z_star, 
       "V_star" = V_star, 
       "p_star" = p_star, 
       "D_star" = D_star)
}



