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


###############################################################################
############################ 4.3 Dependent non-identical #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is denoted by $theta_hat$.

alpha = 0.05

m = 1000
n = 3000

theta_min = 0.875
theta_max = 0.9
step = (theta_max-theta_min)/m
theta_vec = seq(theta_min, theta_max, step)
mu_theta = mean(theta_vec)

mu = floor(n*mu_theta)

rho = 0.6 # correlation coefficient

# Simulate dependency

# Set-up from Boland et al (1989) 'Modelling dependence in simple and indirect majority systems',
# where we have a leading classifier with classifications Y_0, and then the m classifiers with 
# correlation rho = corr(Y_0, Y_j). The m classifiers are independent of each other given Y_0.

# For simplicity, let hat{\theta}_0 = mu_theta
y0 = numeric(n) # vector of zeros of length n
y0[1:mu] = 1 # exactly \mu of them are correct classifications




rep = 1000         # 60 sec for a thousand, 9,000 sec for 100,000

min_dep = numeric(rep) # min number of failures with dependency
min_indep = numeric(rep) # for independent, as a check

tic()
for (ell in 1:rep){
  
  x_dep = numeric(m)  # number of failures for m experiments

  for (j in 1:m){ 
    
    # the probabilities of Y_j being the opposite of Y_0
    p_flip1 = 1-theta_vec[j] - rho*(1-theta_vec[j]) # P(Y_j = 0|Y_0 = 1)
    p_flip0 = theta_vec[j] - rho*theta_vec[j] # P(Y_j = 1|Y_0 = 0)
    
    
    # vectors of 0s and 1s indicating a flip relative to y0
    flip1 = rbinom(mu,1,p_flip1) # flipping correct predictions
    flip0 = rbinom(n-mu,1,p_flip0) # flipping incorrect predictions
    flip = c(flip1,flip0)
    
    y = abs(y0-flip) # correct predictions
    
    x_dep[j] = n-sum(y)
    
    
  }
  
  min_dep[ell] = min(x_dep)
  x_indep = rbinom(m,n,1-theta_vec)
  
  min_indep[ell] = min(x_indep)
  
}
toc()

# Histograms of the minimum number of failures for m classifiers, in rep repetitions.

histbreaks = seq(min(c(x_dep,x_indep)), max(c(x_dep,x_indep))+8,10)

hist(x_dep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,300))
hist(x_indep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,300))

# The upper bound of the 95% confidence interval
sort_min_dep = sort(min_dep) # sort the minimum number of failures
min_dep_alpha2 = sort_min_dep[(alpha/2)*rep] # find the alpha/2 bound

sprintf("The simulated dependent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)

sort_min_indep = sort(min_indep)
min_indep_alpha2 = sort_min_indep[(alpha/2)*rep]

sprintf("The simulated independent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_indep_alpha2)/n, rep)

sprintf("The true SOTA is %s.",  
        theta_max)


