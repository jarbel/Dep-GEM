# fun on the species weights

p_log_p = function(p){if (p==0) 0
                      else abs(p)*log(abs(p))
                      # in case some p is < 0
}
p_log_p = Vectorize(p_log_p)

D_fun = function(p){
  -sum(p_log_p(p))
}

V_to_p = function(V){
  J=length(V)
  p=rep(0,J)
  product=1
  for (j in 1:(J-1)){
    p[j]=product*V[j]
    product=product*(1-V[j])
  }
  p[J]=1-sum(p)
  p
}

p_to_V = function(p){
  J=length(p)
  sum_p = c(0,cumsum(p[1:(J-1)]))
  V = p/(1-sum_p)
  V[is.nan(V)] = 0
  V
}



