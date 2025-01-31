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
library(future)     # parallelisation
library(furrr)

# Contains the function 

# simulated pmf - dep_id_pmf

###############################################################################
############################ Dependent, identical classifiers #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is denoted by $theta_hat$.

future::plan("multisession", workers=6) # let's engage 6 out of 8 cores

source("Parameters_PublicCompetition.R") # n, theta, m, alpha

source("ProbDistr_thetaSOTA.R") 
# for functions 
# sim_ci()
# sim_mean()
# sim_var()

source("dep_id_pmf_fun.R") 
# for function dep_id_pmf
# returns X = list(theta_y0, x_dep_hist, min_dep, min_indep, x_dep, x_indep)
# parallel version dep_id_pmf_parallel
# returns Xmin_fail

# mu = n*theta_SOTA # the expected number of correct predictions for fixed = T


tic()
X = dep_id_pmf(n, theta_SOTA, m, rho, rep, fixed = F) # 20 sec, rep = 100,000
toc()

# Example histogram of the number of failures for m classifiers.
hist(X$x_dep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = 10, ylim = c(0,m/3))

# This is the main outcome.
# Histograms of the minimum number of failures for m dependent classifiers, in rep repetitions. 
# histbreaks = seq(min(c(X$min_dep,X$min_indep)), max(c(X$min_dep,X$min_indep))+8,3)
hist(X$min_dep, xlab = 'minimum number of failures', ylab = 'number of classifiers', 
     breaks = 20, ylim = c(0,rep/3))

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

# parallel version

tic()
Xmin_fail = dep_id_pmf_parallel(n, theta_SOTA, m, rho, rep, fixed = F) # 20 sec, rep = 100,000
toc()

# The confidence interval
min_dep_alpha2 = sim_ci(alpha, Xmin_fail)
sprintf("The simulated dependent upper bound of the %s confidence interval is %.7f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)

# The expected value
Esota = sim_mean(Xmin_fail)
sprintf("The simulated expected value is %.7f, with %s repetitions.",  
        (n-Esota)/n, rep)

# The variance
Vsota = sim_var(Xmin_fail)
sprintf("The simulated standard deviation is %.7f, with %s repetitions.",  
        sqrt(Vsota)/n, rep)





##################### Check-ups  ################
check = 0
if (check){
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
}


### Bias  a function of rho ################

fig = 0
if (fig){
# This figure corresponds to the red line in `bias_thetamin_rho`, found in 
# `dependent_nonidentical.R`

rho_vec = seq(0.0, 1.0, by=0.1) # adjust until smooth

Esota_theta_vec = numeric(length(rho_vec))

for (j in 1:length(rho_vec)){
  rho = rho_vec[j] # for the correlated
  
  X = dep_id_pmf(n, theta, m, rho, rep, fixed = T)
  
  # Expected value
  Esota = sim_mean(X$min_dep)
  Esota_theta_vec[j] = (1-Esota/n)-theta
  print(c(j,length(rho_vec)-j))
}
plot(rho_vec, Esota_theta_vec,"l", lty = "solid", col = "red", ylim = ylm_bias,
     main = '', xlab = TeX(r'(${rho}_0$)'), ylab = ylab_bias)
abline(v=0.6, col="gray")
}
