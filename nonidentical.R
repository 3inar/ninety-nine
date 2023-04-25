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
############################ 4.2 Non-identical aka Poisson Binomial #############################
###############################################################################

# Consider a classification problem with a test set of size $n =3,000$, and 
# $m = 1,000$ classifiers with varying probability of correct classification 

alpha = 0.05

m = 1000
n = 3000
theta_min = 0.875
theta_max = 0.9
theta_mean = (theta_max+theta_min)/2
step = (theta_max-theta_min)/m
theta_vec = seq(theta_min, theta_max, step)

# Simulate non-identical

rep = 100000          # 6 sec for 100,000, 3 min for 1,000,000

min_nonid = numeric(rep) # min number of failures 
min_id = numeric(rep)


tic()
for (ell in 1:rep){
  
  x_nonid = rbinom(m,n,1-theta_vec) # number of failures for classifier j=1, .., m
  x_id = rbinom(m,n,1-theta_mean)
    
  min_nonid[ell] = min(x_nonid)
  min_id[ell] = min(x_id)
    
}
toc()

# Histograms of the minimum number of failures for m classifiers, in rep repetitions.

histbreaks = seq(200,300,1)


hist(min_nonid, xlab = 'number of failures', ylab = 'number of classifiers', breaks = 100, xlim = c(200,300)) 
#  ylim = c(0,300), 

# The upper bound of the 95% confidence interval
sort_min_nonid = sort(min_nonid) # sort the minimum number of failures
min_nonid_alpha2 = sort_min_nonid[(alpha/2)*rep] # find the alpha/2 bound

sprintf("The simulated non-identical upper bound of the %s confidence interval is %.5f, with %s repetitions. Bias: %s.",  
        1-alpha, (n-min_nonid_alpha2)/n, rep, (n-min_nonid_alpha2)/n-theta_max)

# The upper bound of the 95% confidence interval
sort_min_id = sort(min_id) # sort the minimum number of failures
min_id_alpha2 = sort_min_id[(alpha/2)*rep] # find the alpha/2 bound

sprintf("The simulated iid upper bound of the %s confidence interval is %.5f, with %s repetitions.True SOTA: %s",  
        1-alpha, (n-min_id_alpha2)/n, rep, (n-min_id_alpha2)/n-theta_mean)

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

k_alpha2 = qbinom(Palpha2, n, theta_mean, lower.tail = F)
x_alpha2 = n-k_alpha2
theta_alpha2 = (n-x_alpha2)/n
sprintf("With a probablitiy of alpha/2 = %s, at least %s team will achieve an accuracy of at least %.4f.",
        alpha/2, Cx, theta_alpha2)



# Analytical result

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

# # # # # # # # # # # # # # # # # Figure 2 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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


################# probability mass function ##########################

fz = numeric(n)
fz[1:n] = Fz[2:(n+1)]-Fz[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)

# # # # # # # # # # # # # # # # # Figure 3 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

plot(1:n,fz, type = 's')

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
# We do not have access to f(0), because f(0) = F(0)-F(-1), which does not exist. 
# This will not influence the expectation, since 0*f(0) = 0
Eterm = numeric(n+1)
for (z in 1:n){
  i = z
  Eterm[z] = z*fz[i]
}

Esota = sum(Eterm)
Esota_p = 1-Esota/n

# Variance
vterm = numeric(n+1)
for (z in 1:n){
  i = z
  vterm[i] = z^2*fz[i]
}
esquare = sum(vterm)

Vsota = esquare - Esota^2

sprintf("The expected number of failures is %.4f, with a variance of %.4f.",
        Esota, Vsota)
sprintf("The expected theta_hat_SOTA is %.4f, with standard deviation of %.4f.",
        (n-Esota)/n, sqrt(Vsota)/n)






