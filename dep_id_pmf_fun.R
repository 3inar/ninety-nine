# Simulate dependent, identical

# n - number of trials/size of test set
# m - number of experiments/classifiers
# rho - correlation coefficient
# rep - number of repetitions
# theta - theta (probability of correct prediction)
# theta_0 = theta - probability of correct prediction for leading classifier

# The simulated pmf of dependent, identical classifiers
# parameters: n, theta, m, rho, rep, fixed = F is a variant
dep_id_pmf <- function(n, theta, m, rho, rep, fixed = F){   
  
  # Set-up from Boland et al (1989) 'Modelling dependence in simple and indirect majority systems',
  # where we have a leading classifier with classifications Y_0, and then the m classifiers with 
  # correlation rho = corr(Y_0, Y_j). The m classifiers are independent of each other given Y_0.
  
  # the probabilities of Y_j being the opposite of Y_0
  p_flip1 = 1-theta - rho*(1-theta) # P(Y_j = 0|Y_0 = 1), same as 1-p_dep
  p_flip0 = theta - rho*theta# P(Y_j = 1|Y_0 = 0), same as (mu/(n-mu))*(1-p_dep)
  
  # Simulations is the only way
  min_dep = numeric(rep) # min number of failures with dependency
  min_indep = numeric(rep) # for independent, as a check
  
  theta_y0 = numeric(rep)
  x_dep_hist = numeric(rep)
  
  if (fixed){
    y0 = numeric(n) # vector of zeros of length n
    y0[1:mu] = 1 # exactly \mu of them are correct classifications
  } 
  
  for (ell in 1:rep){
    
    x_dep = numeric(m)  # number of failures for m experiments
    
    if (!fixed){
      y0 = rbinom(n,1,theta)
    }
    theta_y0[ell] = sum(y0) # keeping this mu #s
    
    flip1 = rbinom(m,theta_y0[ell],p_flip1) # flipping correct predictions
    flip0 = rbinom(m,n-theta_y0[ell],p_flip0) # flipping incorrect predictions
    
    x_dep = n-(theta_y0[ell]-flip1+flip0) # number of wrong predictions for each classifier
    x_dep_hist[ell] = x_dep[1] # keeping this as an example
    
    min_dep[ell] = min(x_dep) # minimum number of wrong predictions for each rep
    
    x_indep = rbinom(m,n,1-theta) # independent classifiers for reference
    min_indep[ell] = min(x_indep)
  }
  
  X = list(theta_y0 = theta_y0, x_dep_hist = x_dep_hist, min_dep = min_dep, min_indep = min_indep, x_dep = x_dep, x_indep = x_indep)
  
  return(X)
}

# parallel version

dep_id_pmf_parallel <- function(n, theta, m, rho, rep, fixed = F){   
  
  # Set-up from Boland et al (1989) 'Modelling dependence in simple and indirect majority systems',
  # where we have a leading classifier with classifications Y_0, and then the m classifiers with 
  # correlation rho = corr(Y_0, Y_j). The m classifiers are independent of each other given Y_0.
  
  # the probabilities of Y_j being the opposite of Y_0
  p_flip1 = 1-theta - rho*(1-theta) # P(Y_j = 0|Y_0 = 1), same as 1-p_dep
  p_flip0 = theta - rho*theta# P(Y_j = 1|Y_0 = 0), same as (mu/(n-mu))*(1-p_dep)
  
  # Simulations is the only way
  min_dep = numeric(rep) # min number of failures with dependency
  min_indep = numeric(rep) # for independent, as a check
  
  theta_y0 = numeric(rep)
  x_dep_hist = numeric(rep)
  
  Xb <- furrr::future_map(1:rep, function (x) { 
    
    x_dep = numeric(m)  # number of failures for m experiments
    
    y0 = rbinom(n,1,theta)
    
    theta_y0 = sum(y0) # keeping this mu #s
    
    flip1 = rbinom(m,theta_y0,p_flip1) # flipping correct predictions
    flip0 = rbinom(m,n-theta_y0,p_flip0) # flipping incorrect predictions
    
    x_fail = n-(theta_y0-flip1+flip0) # number of wrong predictions for each classifier
    
    min_fail = min(x_fail) # minimum number of wrong predictions for each rep
    
    return(min_fail)
    
  }, .progress=T, .options= furrr::furrr_options(seed=T))
  
  # Not very elegant with the for-loop, but that's ok for now
  Xmin_fail = numeric(rep)
  for (r in 1:rep){ 
    Xmin_fail[r] = Xb[[r]]
  }
  
  return(Xmin_fail)
}
