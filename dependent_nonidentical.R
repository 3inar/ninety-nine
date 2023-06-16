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

# Functions

# dep_nonid_pmf - simulated pmf

source("Parameters_PublicCompetition.R") # n, theta, m, alpha, rho, theta_min, theta_max, theta_vec, 

source("dep_nonid_pmf_fun.R")



tic()
X = dep_nonid_pmf(n, theta, m, rho, rep, theta_vec, theta_0 = theta)
toc()

# Histograms of the minimum number of failures for m classifiers, in rep repetitions.

hist(X$x_dep, xlab = 'number of failures', ylab = 'number of classifiers', 
     ylim = c(0,250))


source("ProbDistr_thetaSOTA.R")
# The upper bound of the 95% confidence interval
min_dep_alpha2 = sim_ci(alpha, X$min_dep)
sprintf("The simulated dependent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)

# The expected value
Esota = mean(X$min_dep)

# The standard deviation
Vsota = mean(X$min_dep*X$min_dep) - Esota*Esota

sprintf("The simulated expected value is %.7f and a standard deviation is %.7f, with %s repetitions.",  
        (n-Esota)/n, sqrt(Vsota)/n, rep)
