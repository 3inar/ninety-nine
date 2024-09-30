# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) February 2023 updated March 2023

###############################################################################
############################ 3.1 Nomenclature #############################
###############################################################################

# n     - number of trials/size of test set
# m     - number of experiments/number of classifiers
# theta     - prob of success/prob of correct prediction

# Y     - r.v., indicates if prediction is correct

# k     - number of successes
# X     - r.v., number of failures; n-k

# P_    - probability, specified when used

# px    - probability of x failures in one experiment
# Px    - probability of at most x failures in one experiment

# Cx    - r.v., number of experiments with at most x failures

# theta_hat - accuracy of a classifier: p_hat = (n-x)/n

# theta_hat_SOTA - max_j p_hat_j

# Z     - r.v., number of failures on at least one classifier

# Fz    - cdf of Z: P(C_z > 0|m,n,theta) prob of at least one classifier having at most 
#         z failures (identical to at least n-z successes)
# fz    - pmf of Z

# alpha - significance level

# SOTA - state-of-the-art


library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing
library(latex2exp)  # mathematical notation
library(e1071)      # skewness

# Contains the function 

# simulated pmf - dep_id_pmf

###############################################################################
############################ 4.1 Dependent, identical classifiers #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is denoted by $theta_hat$.

source("Parameters_PublicCompetition.R") # n, theta, m, alpha

mu = n*theta # the expected number of correct predictions

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

tic()
X = dep_id_pmf(n, theta, m, rho, rep, fixed = T)
toc()

# Example histogram of the number of failures for m classifiers.
hist(X$x_dep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = 10, ylim = c(0,m/3))

# This is the main outcome.
# Histograms of the minimum number of failures for m dependent classifiers, in rep repetitions. 
# histbreaks = seq(min(c(X$min_dep,X$min_indep)), max(c(X$min_dep,X$min_indep))+8,3)
hist(X$min_dep, xlab = 'minimum number of failures', ylab = 'number of classifiers', 
     breaks = 20, ylim = c(0,rep/3))

source("ProbDistr_thetaSOTA.R")
# source("Parameters_PublicCompetition.R") # n, theta, m, alpha

# The confidence interval
min_dep_alpha2 = sim_ci(alpha, X$min_dep)
sprintf("The simulated dependent upper bound of the %s confidence interval is %.7f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)

# The expected value
Esota = sim_mean(X$min_dep)
sprintf("The simulated expected value is %.7f, with %s repetitions.",  
        (n-Esota)/n, rep)

# The variance
Vsota = sim_var(X$min_dep)
sprintf("The simulated standard deviation is %.7f, with %s repetitions.",  
        sqrt(Vsota)/n, rep)

##################### Check-ups  ################

# Example histogram of the number of failures for m classifiers.
histbreaks = seq(min(c(X$x_dep,X$x_indep)), max(c(X$x_dep,X$x_indep))+8,10)
hist(X$x_indep, xlab = 'number of failures', ylab = 'number of classifiers', 
     #breaks = histbreaks, 
     ylim = c(0,m/3))

# Displays the distribution of the observed theta_y0's (we know it's binomial, so it's more of a check-up)
hist(X$theta_y0/n, xlab = mean(round(X$theta_y0)/n)) # only make sense if fixed = F

# Displays the distribution of one observed x_dep per rep (we know it's binomial, same as 1-theta_y0)
hist(X$x_dep_hist,xlab = mean((n-X$x_dep_hist)/n))

# Displays the independent counterpart
histbreaks = seq(min(c(X$min_dep,X$min_indep)), max(c(X$min_dep,X$min_indep))+8,3)
hist(X$min_indep, xlab = 'minimum number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,rep/3))

# The upper bound for the independent counterpart
min_indep_alpha2 = sim_ci(alpha, X$min_indep)
sprintf("The simulated independent upper bound of the %s confidence interval is %.7f, with %s repetitions.",  
        1-alpha, (n-min_indep_alpha2)/n, rep)

###########################################################################
### Expected value/bias/variance as a function of rho ################
###########################################################################

rho_vec = seq(0.0, 1.0, by=0.01)

ylm = c(0.0,0.035)

Esota_theta_vec = numeric(length(rho_vec))

for (j in 1:length(rho_vec)){
  rho = rho_vec[j] # for the correlated
  
  X = dep_id_pmf(n, theta, m, rho, rep, fixed = T)
  
  # Expected value
  Esota = sim_mean(X$min_dep)
  Esota_theta_vec[j] = (1-Esota/n)-theta
  print(j)
}
plot(rho_vec, Esota_theta_vec,"l", lty = "solid", col = "darkgreen", ylim = ylm,
     main = TeX(r'(Bias as a function of correlation)'), xlab = TeX(r'(${rho}_0$)'), ylab = TeX(r'($E \hat{theta}_{SOTA} - {theta}_{SOTA}$)'))
abline(v=0.6, col="gray")

