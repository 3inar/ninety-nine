# Simulate dependent, nonidentical

# n - number of trials/size of test set
# m - number of experiments/classifiers
# rho - correlation coefficient
# rep - number of repetitions
# theta_vec - vector of thetas (probabilities of correct prediction)
# theta_0 - probability of correct prediction for leading classifier

# output: list X
# min_fail - vector of length rep with min number of failures for each rep
# teamsSOTA - vector of length rep with number of teams performing above theta_0 for each rep
# x_fail - vector of length m with number of failures for each team, as example


dep_nonid_pmf <- function(n, m, rho, rep, theta_vec, theta_0){   
  
  # Set-up from Boland et al (1989) 'Modelling dependence in simple and indirect majority systems',
  # where we have a leading classifier with classifications Y_0, and then the m classifiers with 
  # correlation rho = corr(Y_0, Y_j). The m classifiers are independent of each other given Y_0.
  
  sigma_0 = theta_0*(1-theta_0)
  sigma_vec = theta_vec*(1-theta_vec)
  
  p_flip1 = (theta_0 - rho*sqrt(sigma_0*sigma_vec) - theta_0*theta_vec)/theta_0 # P(Y_j = 0|Y_0 = 1)
  p_flip0 = (-rho*sqrt(sigma_0*sigma_vec)+(1-theta_0)*theta_vec)/(1-theta_0) # P(Y_j = 1|Y_0 = 0)
  
  #if (any(p_flip1<0)){
  #  print(p_flip1)
  #  print(rho)
  #  readline(prompt="Press [enter] to continue")
  #}
  
  
  min_fail = numeric(rep) # min number of failures with dependency
  teamsSOTA = numeric(rep) # number of teams above theta_0
  
  for (ell in 1:rep){
    
    x_fail = numeric(m)  # number of failures for m experiments
    
    y0 = rbinom(n,1,theta_0) # leading classifier outcome
    theta_y0 = sum(y0) # observed theta for leading classifier
    
    flip1 = rbinom(m,theta_y0,p_flip1) # flipping correct predictions
    flip0 = rbinom(m,n-theta_y0,p_flip0) # flipping incorrect predictions
    
    x_fail = n-(theta_y0-flip1+flip0) # number of wrong predictions for each classifier
    hat_theta = (n-x_fail)/n
    teamsSOTA[ell] = length(hat_theta[hat_theta > theta_0])
    
    min_fail[ell] = min(x_fail)
  }
  
  X = list(min_fail = min_fail, x_fail = x_fail, teamsSOTA = teamsSOTA)
  
  return(X)
}
