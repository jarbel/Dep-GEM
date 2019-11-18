


##################
# log-density a priori 
###################

# pour un vecteur Z_j
log_prior_Z = function(Z, sigma_Z, K){
  mvtnorm::dmvnorm(Z,rep(0, length(Z)), sigma_Z^2*K, log=TRUE)
}

# pour la matrix de tous les Z_j
log_prior_Z_total = function(Z_matrix, sigma_Z, K){
  I = dim(Z_matrix)[1]
  J = dim(Z_matrix)[2]
  res = 0
  for (j in 1:J){
    res = res + mvtnorm::dmvnorm(Z_matrix[,j], rep(0, I), sigma_Z^2*K, log=TRUE)
  }
  res
}

##################
# likelihood in terms of V
# log likelihood of each V_j of V: it's a vector of size I_X
log_like_V_total=function(V, data_freq){
  matrice=data_freq[,,1]*log(V)+data_freq[,,2]*log(1-V)
  apply(X=matrice,MARGIN=2,FUN=sum)
}

log_like_Z_total = function(Z_matrix, sigma_Z, M, data_freq){
  I = dim(Z_matrix)[1]
  J = dim(Z_matrix)[2]
  res = 0
  for (j in 1:J){
    V = Z_to_V(Z_matrix[,j], sigma_Z, M)
    res = res + log_like_V(V, j, data_freq)
  }
  res
}

# compute for V_j only, res is a real nb
log_like_V=function(V, j, data_freq){
  V[V==1]=.999
  V[V==0]=.001
  sum(data_freq[,j,1]*log(V)+data_freq[,j,2]*log(1-V))
}

F_inv_V = function(U,M){
  1-(1-U)^(1/M)
}

Z_to_V = function(Z, sigma_Z, M){
  F_inv_V(pnorm(Z,0,sd=sigma_Z), M)
}