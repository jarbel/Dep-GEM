# GP covariance functions



K_fun = function(lambda, R, R2, GP="SE", diag_off = FALSE, d = 1){
  # squared-exponential
  if (GP=="SE" && !diag_off)
    exp(-R2/(2*lambda^2)) + d*diag(dim(R2)[1])
  else if (GP=="SE" && diag_off)
    exp(-R2/(2*lambda^2))
  # Ornstein-Ulenbeck
  else if (GP=="OU")
    exp(-R/lambda) + diag(dim(R)[1])
  # rational quadratic, alpha=1
  else if (GP=="RQ"){
    alpha = 1
    (1+R2/(2*alpha*lambda^2))^(-alpha) + diag(dim(R2)[1])}
}