# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) 

library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing
library(latex2exp)  # mathematical notation

###############################################################################
############################ 3.1 Nomenclature #############################
###############################################################################

# n         - number of trials/size of test set
# theta     - prob of success/prob of correct prediction
# k         - number of successes
# X         - r.v., number of failures; n-k
# theta_hat = (n-X)/n - estimator of theta
# px    - probability of x failures in one experiment
# Px    - probability of at most x failures in one experiment

# m     - number of experiments/number of classifiers
# Cx    - r.v., number of experiments with at most x failures
# Z     - r.v., number of failures on at least one classifier

# theta_SOTA_hat - estimator of theta_SOTA 

# P_    - probability, specified when used

# Fz    - cdf of Z: P(C_z > 0|m,n,theta) prob of at least one classifier having at most 
#         z failures (identical to at least n-z successes)
# fz    - pmf of Z

# alpha - significance level

# SOTA - state-of-the-art

###############################################################################
############################ 3.2 Coin flipping example #############################
###############################################################################

# see CoinFlip.R

###############################################################################
########### 3.3 The probability distribution of p_SOTA ########################
###############################################################################

n = 3000
theta = 0.9
m = 1000
mu = n*theta # the expected number of correct predictions
alpha = 0.05

################################# Cumulative distribution function SOTA ############################################

# Z is the number of failures in at least one classifier.
# Let's have a look at probabilities for at least one classifier having at most 
# z failures. This gives the cumulative distribution function

Fz = numeric(n+1)

# use i for counting, since z = 0, ..., n

for (z in 0:n){
  Pz = pbinom(z,n,(1-theta)) 
  i = z+1
  Fz[i] = pbinom(0,m,Pz,lower.tail = F) #  P[C>0]
}

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
x_sota_hat = which(Fz>0.025)[1]-1 # should be same as x_alpha2 = 236
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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

# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# sum(fz)
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Fz = numeric(n+1)

# for (z in 0:n){
#   pz = pbinom(z,n,1-p)
#   Fz[z+1] = 1 - (1 - pz)^m
# }
# fz = Fz[2:(n+1)]-Fz[1:n]

# z = 200:300
# plot(z,Fz[z+1], type = 's')
# plot(z,fz[z], type = 'h', xlab = '', ylab = 'probability mass', axes = F)
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


###############################################################################
################# 3.4 A simulated public competition example ##################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $p$.
# The classifier's performance on the test set, denoted by $p_hat$ and referred 
# to as the accuracy is $(n-x_j)/n$, where $x_j$ is the observed number of failures
# for classifier $j$.

n = 3000
theta = 0.9
mu = n*theta # the expected number of correct predictions

# 95\% confidence interval
alpha = 0.05
ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %s is (%.4f,%.4f).",  
        (1-alpha)*100, theta, ci_binom["lower"], ci_binom["upper"])

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# The probability of exceeding CI_{up} should be close to alpha/2 = 0.025 # # # 

k_up = floor(ci_binom[["upper"]]*n) # number of successes exceeding the CI
# flooring the CI bound, so P_up > alpha/2
P_up = pbinom(k_up,n,theta, lower.tail = F) # P[X>x]
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

# If there are $m$ teams, each with a classifier with $\theta$, what is
# the probability of at least one team achieving at least $ci_binom["upper"]$ 
# correct predictions?
# Following the logic of the coin-flip multiple experiments, we have a binomial 
# distribution with $m$ trials and probability of success (exceeding 95% CI) 
# is by definition $alpha/2$. 

Cx = 1 # number of teams
Palpha2 = pbinom(Cx-1,m,alpha/2, lower.tail = F) # P[X>x].
sprintf("The probability of at least %s out of %s teams exceeding the upper limit of the CI is %s.",  
        Cx, m, Palpha2)
sprintf("Notice that the accuracy, theta, is not part of the equation.")

# # # # # # # # # # # # # # # # # Simulations # # # # # # # # # # # # # # # # # # # # # # # 

rep = 1000000 # one million repetitions takes about 80 seconds
# I still don't expect any failures

k_up = floor(ci_binom[["upper"]]*n) # number of successes exceeding the CI

n_success = numeric(rep) # pre-allocate
tic()
for (i in 0: rep){
  sim_success = rbinom(m, n, theta) # the number of successes for each team
  n_success[i] = length(which(sim_success>k_up)) # how many are above CI_{up}?
}

which(n_success == 0) # any repetitions with all below CI_{up}?
toc()
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

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

k_alpha2 = qbinom(Palpha2, n, theta, lower.tail = F)
x_alpha2 = n-k_alpha2
theta_alpha2 = (n-x_alpha2)/n
sprintf("With a probablitiy of alpha/2 = %s, at least %s team will achieve an accuracy of at least %.4f.",
        alpha/2, Cx, theta_alpha2)

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
Palpha2_discr = pbinom(k_alpha2,n,theta,lower.tail = F) # < Palpha2
P_discr = pbinom(Cx-1,m,Palpha2_discr, lower.tail = F) # = 0.02475 < P
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

sprintf("In summary, if there are %s teams, each with %s probability of correct prediction, and a test set of size %s, there is a probability of %s that at least %s team will have at most %s incorrect predictions, corresponding to an estimated accuracy of %.4f.",  
        m, theta, n, alpha/2, Cx, x_alpha2, theta_alpha2)

# # # # # # # # # # # # # # # # # Simulations # # # # # # # # # # # # # # # # # # # # # # # 
# I hope to recreate alpha/2
rep = 1000000 # one million repetitions takes about 80 seconds

n_success = numeric(rep)
tic()
for (i in 1: rep){
  sim_success = rbinom(m, n, theta) # the number of successes for each team
  n_success[i] = length(which(sim_success>k_alpha2)) # how many are above k_alpha2?
}

sum(n_success)/rep # fraction of successes, should around alpha/2 = 0.025
# there is no continuity correction, so I thought I would get < 0.025
toc()
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

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

# When the probability of correct prediction is $0.9$, there is a non-negligible 
# (at least in the common statistical significance sense) probability that the 
# top-ranked team will have an estimated accuracy of at least 0.9213, which is 
# often referred to as state-of-the-art (SOTA) performance. 

# What is the probability of beating the SOTA for a method with significantly better
# accuracy?

P_beat_sota = pbinom(k_alpha2, n, ci_binom[["upper"]], lower.tail = F) # P[X>x]

sprintf("The probability of achieving an estimated accuracy better than the upper bound of the CI for theta_hat_SOTA, %.4f, for a classifier with significantly better probability of correct prediction, p=%.4f, is just %.4f.",
        theta_alpha2, ci_binom[["upper"]], P_beat_sota)

# # # # # # # # # # # # # # # # # Figure 1 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

x = (mu-60):(mu+90) # this is the plot range, adjust to your liking

y = dbinom(x, n, theta) # probability of x successes in n trials
plot(x, y, type='h', xlab = '', ylab = '', axes=F) # histogram

par(new=TRUE) # new plot in same window

y = dbinom(x, n, ci_binom[["upper"]]) # significantly better classifier
plot(x, y, type = 'l', col = "red", xlab = '', ylab = '', axes=F)

# area under curve for expected SOTA performance

polygon(c(n-Esota, x[x>=n-Esota], max(x)), c(0,y[x>=n-Esota], 0), col="blue")

# area under curve for SOTA performance
polygon(c(k_alpha2, x[x>=k_alpha2], max(x)), c(0,y[x>=k_alpha2], 0), col="red")


axis(1 , las = 2, at=c(ci_binom[["lower"]]*n, n*theta, 
                       ci_binom[["upper"]]*n, n-Esota, k_alpha2), 
     labels=c(TeX('$\\theta_{alpha/2}$'), TeX('$\\theta$'), TeX('$\\theta_{1-alpha/2}$'), 
              TeX(r'($E \hat{theta}_{SOTA}$)'), TeX('$\\theta^m_{1-alpha/2}$')))

# # # # # # # # # # # # # # # # # # end figure # # # # # # # # # 

# beat esota?

beat_esota = pbinom(n-Esota, n, ci_binom[["upper"]], lower.tail = F) 

sprintf("The probability of achieving an accuracy better than E(SOTA)=%.4f, for a classifier with significantly better proability of correct prediction, %.4f, is %.4f.",
        Esota_p, ci_binom[["upper"]], beat_esota)


######################### varying m, n, p ####################################

# Expectation and standard deviation for p_SOTA

ESp_SOTA <- function(m,n,p){
  
  # Cumulative density function
  Fz = numeric(n+1)
  # z is the number of failures
  # use i for counting, since z = 0, ..., n
  for (z in 0:n){
    Pz = pbinom(z,n,(1-p)) 
    i = z+1
    Fz[i] = pbinom(0,m,Pz,lower.tail = F) #  P[C>0]
  }
  
  # Probability mass function
  fz = numeric(n)
  fz[1:n] = Fz[2:(n+1)]-Fz[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)
  
  # Expectation
  
  # Calculate each term z*f(z)
  Eterm = numeric(n)
  for (z in 1:n){
    Eterm[z] = z*fz[z]
  }
  
  Esota = sum(Eterm) # the sum, makes the expected number of failures
  Ep_SOTA = (n-Esota)/n # translated to accuracy p_SOTA = (n-z)/n
  
  # Variance
  vterm = numeric(n+1)
  for (z in 1:n){
    vterm[z] = z^2*fz[z]
  }
  esquare = sum(vterm)
  
  Vp_SOTA = esquare - Esota^2
  Sp_SOTA = sqrt(Vp_SOTA)/n
  
  return(c(Ep_SOTA, Sp_SOTA))
   
}

m = 1000
n = 3000
theta = 0.9
estd =  ESp_SOTA(m, n, theta)
sprintf("The expected theta_sota is %.4f, with a standard deviation of %.6f, for m=%s, n=%s, theta=%s.",
        estd[1], estd[2], m, n, theta)









