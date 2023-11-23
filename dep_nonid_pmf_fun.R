# Simulate dependent, nonidentical
dep_nonid_pmf <- function(n, m, rho, rep, theta_vec, theta_0){   
  
  # Set-up from Boland et al (1989) 'Modelling dependence in simple and indirect majority systems',
  # where we have a leading classifier with classifications Y_0, and then the m classifiers with 
  # correlation rho = corr(Y_0, Y_j). The m classifiers are independent of each other given Y_0.
  
  sigma_0 = theta_0*(1-theta_0)
  sigma_vec = theta_vec*(1-theta_vec)
  
  p_flip1 = (theta_0 - rho*sqrt(sigma_0*sigma_vec) - theta_0*theta_vec)/theta_0 # P(Y_j = 0|Y_0 = 1)
  p_flip0 = (-rho*sqrt(sigma_0*sigma_vec)+(1-theta_0)*theta_vec)/(1-theta_0) # P(Y_j = 1|Y_0 = 0)
  
  min_dep = numeric(rep) # min number of failures with dependency
  min_indep = numeric(rep) # for independent, as a check
  teamsSOTA = numeric(rep) # for independent, as a check
  
  for (ell in 1:rep){
    
    x_dep = numeric(m)  # number of failures for m experiments
    
    y0 = rbinom(n,1,theta_0)
    theta_y0 = sum(y0) 
    
    flip1 = rbinom(m,theta_y0,p_flip1) # flipping correct predictions
    flip0 = rbinom(m,n-theta_y0,p_flip0) # flipping incorrect predictions
    
    x_dep = n-(theta_y0-flip1+flip0) # number of wrong predictions for each classifier
    hat_theta = (n-x_dep)/n
    teamsSOTA[ell] = length(hat_theta[hat_theta > theta_0])
    
    min_dep[ell] = min(x_dep)
  }
  
  X = list(min_dep = min_dep, x_dep = x_dep, teamsSOTA = teamsSOTA)
  
  return(X)
}
