# Simulate non-identical
indep_nonid_pmf <- function(n, theta_vec, m, rep){   
  
  min_nonid = numeric(rep) # min number of failures 

  for (ell in 1:rep){
    x_nonid = rbinom(m,n,1-theta_vec) # number of failures for classifier j=1, .., m
    min_nonid[ell] = min(x_nonid)
  }
  X = list(min_nonid = min_nonid, x_nonid = x_nonid)
  
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