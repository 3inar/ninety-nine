# Simulate non-identical
indep_nonid_pmf <- function(n, theta_vec, m, rep){   
  
  min_fail = numeric(rep) # min number of failures 

  for (ell in 1:rep){
    x_fail = rbinom(m,n,1-theta_vec) # number of failures for classifier j=1, .., m
    min_fail[ell] = min(x_fail)
  }
  X = list(min_fail = min_fail, x_fail = x_fail)
  
  return(X)
}

# Analytical result
nonid_cdf <- function(n, theta_vec, m){   
  Fz = numeric(n+1)
  for (z in 0:n){
    P = numeric(m)
    term = numeric(m)
    
    for (j in 1:m){
      P[j] = pbinom(z,n,(1-theta_vec[j])) 
      term[j] = 1-P[j]
    }
    i = z+1
    Fz[i] = 1 - prod(term)
  }
  return(Fz)
}

nonid_pmf <- function(n, theta_vec, m){   
  fz = numeric(n+1)
  Fz = nonid_cdf(n, theta_vec, m)
  fz[1] = Fz[1]
  fz[2:(n+1)] = Fz[2:(n+1)]-Fz[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)
  
  return(fz)
}