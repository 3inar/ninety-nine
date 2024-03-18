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

# functions

# simulated pmf - sim_nonid_pmf
# cdf - nonid_cdf
# cdf based pmf - nonid_pmf


###############################################################################
############################ 4.2 Non-identical aka Poisson Binomial #############################
###############################################################################

# Consider a classification problem with a test set of size $n$, and 
# $m$ classifiers with varying probability of correct classification 

source("Parameters_PublicCompetition.R") # n, theta, m, alpha, theta_min, theta_max, theta_vec

source("indep_nonid_pmf_fun.R") # for the function 'indep_nonid_pmf' - simulated pmf
# nonid_cdf, nonid_pmf - analytical pmf

tic()
X = indep_nonid_pmf(n, theta_vec, m, rep)
toc()

# Histograms of the minimum number of failures for m classifiers, in rep repetitions.

histbreaks = seq(200,300,1)
hist(X$min_fail, xlab = 'number of failures', ylab = 'number of classifiers', breaks = 100, xlim = c(200,300)) 
#  ylim = c(0,300), 

# The upper bound of the 95% confidence interval
sort_min_fail = sort(X$min_fail) # sort the minimum number of failures
min_fail_alpha2 = sort_min_fail[(alpha/2)*rep] # find the alpha/2 bound

sprintf("The simulated non-identical upper bound of the %s confidence interval is %.5f, with %s repetitions. Distance to SOTA: %s.",  
        1-alpha, (n-min_fail_alpha2)/n, rep, (n-min_fail_alpha2)/n-theta_max)

# If there are $m$ classifiers, what must be the prob, $Palpha2$, of each 
# classifier having at most $x$ failures, for the probability of at least Cx = 1
# classifier having at most $x$ failures to be P(Cx > 0) = alpha/2?
# We can calculate $x$ from $Palpha2$.

Cx = 1

# This is easily calculated:
# P(Cx > 0) = 1 - P(Cx = 0) = alpha/2
# binom coeff is 1 since cx = 0, Palpha2^0 = 1
# and we are left with (1-Palpha2)^m = 1 - alpha\2

Palpha2 = 1 - (1-alpha/2)^(1/m)

sprintf("For the probability of at least %s classifier to achieve at most $x$ failures to be alpha/2=%s, the prob for each classifier achieving at most $x$ failures must be %.8f.",  
        Cx, alpha/2, Palpha2)

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
P = pbinom(Cx-1,m,Palpha2, lower.tail = F) # should be 0.025
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

# If the probability of wrongly predicting at most x out of n data points 
# is Palpha2. What is x when theta = 0.9?

theta_mean = (theta_max+theta_min)/2

k_alpha2 = qbinom(Palpha2, n, theta_mean, lower.tail = F)
x_alpha2 = n-k_alpha2
theta_alpha2 = (n-x_alpha2)/n
sprintf("With a probablitiy of alpha/2 = %s, at least %s team will achieve an accuracy of at least %.4f.",
        alpha/2, Cx, theta_alpha2)

# Analytical results:

# # # # # # # # # # # # # # # # # Figure 6 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Fz = nonid_cdf(n, theta_vec, m)
# the whole range, not very much information
plot(0:n,Fz, type = 'l', xlab = 'number of failures', 
     ylab = 'probability of at least one team')

# zooming in, and it gets more interesting, discreet curve 
z = 200:300
plot(z,Fz[z+1], type = 's', xlab = 'number of failures/accuracy', 
     ylab = 'F(z)', axes = F)

xax = seq(z[1],tail(z,1), 10)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
axis(2, las = 2)
axis(3, las = 2, at=xax, labels = as.character(klab))

# Does this align with the simulations?

print(c(Fz[min_fail_alpha2], Fz[min_fail_alpha2+1])) # ok 
print(c(which(Fz>alpha/2)[1]-1,min_fail_alpha2)) # ok


################# probability mass function ##########################

# # # # # # # # # # # # # # # # # Figure 7 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fz = nonid_pmf(n, theta_vec, m)

plot(0:n,fz, type = 's')

# need to zoom in
z = 200:300
plot(z,fz[z], type = 'h', xlab = 'number of failures/accuracy', 
     ylab = 'f(z)', axes = F)
xax = seq(z[1],tail(z,1), 10)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
yax = seq(0,max(fz)+0.01,0.02)
axis(2, las = 1, at=yax, labels = as.character(yax))
axis(3, las = 2, at=xax, labels = as.character(klab))
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

####################### Expected value and variance ###################################

# Expected value
Eterm = numeric(n+1)
for (z in 0:n){
  i = z+1
  Eterm[i] = z*fz[i]
}

Esota = sum(Eterm)
Esota_theta = 1-Esota/n

# Variance
vterm = numeric(n+1)
for (z in 0:n){
  i = z+1
  vterm[i] = z^2*fz[i]
}
esquare = sum(vterm)

Vsota = esquare - Esota^2

sprintf("The expected number of failures is %.4f, with a variance of %.4f.",
        Esota, Vsota)
sprintf("The expected theta_hat_SOTA is %.6f, with standard deviation of %.6f.",
        (n-Esota)/n, sqrt(Vsota)/n)






