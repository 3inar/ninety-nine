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


# comes with the functions 'cdf', 'pmf', 'expect' and 'variance'
source("ProbDistr_thetaSOTA.R") 

# figures at the end

###############################################################################
################# 3.4 A simulated public competition example ##################
###############################################################################

# Consider a classification problem with a test set of size $n=3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The classifier's performance on the test set, denoted by $theta_hat$ and referred 
# to as the accuracy is $(n-x_j)/n$, where $x_j$ is the observed number of failures
# for classifier $j$.

n = 3000
theta = 0.9
mu = n*theta # the expected number of correct predictions
m = 1000

# 95% confidence interval for hat{theta} = mu, single classifier
alpha = 0.05
ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %s is (%.4f,%.4f).",  
        (1-alpha)*100, theta, ci_binom["lower"], ci_binom["upper"])

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# The probability of exceeding the upper limit of the CI should be close to alpha/2 = 0.025 
k_up = floor(ci_binom[["upper"]]*n) # number of successes exceeding the CI
# flooring the CI bound, so P_up > alpha/2
P_up = pbinom(k_up,n,theta, lower.tail = F) # P[X>x]
# # # # # # # # # # # # # # # # # P_up = 0.02619 OK # # # # # # # # # # # # # # # # # # # # # # # 

# If there are $m=1,000$ teams, each with a classifier with $\theta$, what is
# the probability of at least one team achieving at least $ci_binom["upper"]$ 
# correct predictions?
# Following the logic of the coin-flip multiple experiments, we have a binomial 
# distribution with $m$ trials and probability of success (exceeding 95% CI) 
# is by definition $alpha/2$. 

# The probability of at least one team exceeding the upper CI bound:
Cx = 1 # number of teams
Px = alpha/2 # by definition
Palpha2 = 1 - dbinom(Cx-1,m,Px)

sprintf("The probability of at least %s out of %s teams exceeding the upper limit of the CI is %s.",  
        Cx, m, Palpha2)

# If there are $m$ classifiers, what must be the prob, $Px$, of each 
# classifier having at most $x$ failures, for the probability of at least Cx = 1
# classifier having at most $x$ failures to be P(Cx > 0) = alpha/2?
# We can calculate $x$ from $Px$.

Cx = 1 # number of teams

# This is easily calculated:
# P(Cx > 0) = 1 - P(Cx = 0) = alpha/2 
# binom coeff is 1 since cx = 0, Palpha2^0 = 1
# and we are left with (1-Palpha2)^m = 1 - alpha\2

Px = 1 - (1-alpha/2)^(1/m)

sprintf("For the probability of at least %s out of %s classifiers to achieve at most $x$ failures to be alpha/2=%s, the prob for each classifier achieving at most $x$ failures must be %.8f.",  
        Cx, m, alpha/2, Px)

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
P = 1-pbinom(Cx-1,m,Px) # should be 0.025
# # # # # # # # # # # # # # # # # P = 0.025000 OK # # # # # # # # # # # # # # # # # # # # # # # 

# If the probability of wrongly predicting at most x out of n data points is Px, 
# what is x when theta = 0.9?

#"The quantile is defined as the smallest value x such that F(x)≥Px, where F is the distribution function."

x_alpha2 = qbinom(Px, n, 1-theta) 
theta_alpha2 = (n-x_alpha2)/n
sprintf("With a probablitiy of alpha/2 = %s, at least %s team will achieve an accuracy of at least %.4f, corresponding to x = %s failures.",
        alpha/2, Cx, theta_alpha2, x_alpha2)

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Fz = cdf(n,theta,m) # the cdf
x = which(Fz>alpha/2)[1]-1 # should be same as x_alpha2 = 236
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# Since x_alpha2 is discreet, but Px is not, we'll calculate for x_alpha2 and x_alpha2-1, 
# and then the two results should be on each side of alpha/2

Palpha2_discr = pbinom(x_alpha2-1,n,1-theta) # 
P_low = 1-pbinom(Cx-1,m,Palpha2_discr) #

Palpha2_discr = pbinom(x_alpha2,n,1-theta) # 
P_up = 1-pbinom(Cx-1,m,Palpha2_discr) #
# # # # # # # # # # # # # # # # # ok 0.02475 0.03243# # # # # # # # # # # # # # # # # # # # # # # 

sprintf("In summary, if there are %s teams, each with %s probability of correct prediction, and a test set of size %s, there is a probability of %s that at least %s team will have at most %s incorrect predictions, corresponding to an estimated accuracy of %.4f.",  
        m, theta, n, alpha/2, Cx, x_alpha2, theta_alpha2)

# Expected value
Esota = expect(n,theta,m)
Esota_theta = 1-Esota/n
sprintf("The expected value is %.4f",
        Esota_theta)

# When the probability of correct prediction is $0.9$, there is a non-negligible probability that the 
# top-ranked team will have an estimated accuracy of at least theta^m_alpha2 = 0.9213. 

# What is the probability of beating the SOTA for a method with significantly better
# accuracy?

P_beat_sota[1] = pbinom(x_alpha2, n, 1-ci_binom[["upper"]]) # 
P_beat_sota[2] = pbinom(x_alpha2-1, n, 1-ci_binom[["upper"]]) # 
P_beat_sota

sprintf("The probability of achieving an estimated accuracy better than the upper bound of the CI for theta_hat_SOTA, %.4f, for a classifier with significantly better probability of correct prediction, p=%.4f, is just %.4f.",
        theta_alpha2, ci_binom[["upper"]], P_beat_sota[1])

# beat esota?

beat_esota = pbinom(Esota, n, 1-ci_binom[["upper"]]) 

sprintf("The probability of achieving an accuracy better than E(SOTA)=%.4f, for a classifier with significantly better proability of correct prediction, %.4f, is %.4f.",
        Esota_theta, ci_binom[["upper"]], beat_esota)

######################### varying m, n, p ####################################

n = c(3000,1000,10000)
m = c(1000,100,500)
theta = c(0.85,0.9,0.95)

for (i in 1:3){
  for (j in 1:3){
    for (k in 1:3){
      Esota =  expect(n[i], theta[k], m[j])
      Esota_theta = 1-Esota/n[i]
      Vsota = variance(n[i], theta[k], m[j])
      Std_sota_theta = sqrt(Vsota)/n[i]
      sprintf("The expected theta_sota is %.4f, with a standard deviation of %.6f, for m=%s, n=%s, theta=%s.",
              Esota_theta, Std_sota_theta, m[j], n[i], theta[k])
    }
  }
}






# # # # # # # # # # # # # # # # # Figure 1 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# let k be the number of successes, = n-x

k = (mu-60):(mu+90) # this is the plot range, adjust to your liking

y = dbinom(k, n, theta) # probability of x successes in n trials
plot(k, y, type='h', xlab = '', ylab = '', axes=F) # histogram

par(new=TRUE) # new plot in same window

y = dbinom(k, n, ci_binom[["upper"]]) # significantly better classifier
plot(k, y, type = 'l', col = "red", xlab = '', ylab = '', axes=F)

# area under curve for expected SOTA performance
polygon(c(n-Esota, k[k>=n-Esota], max(k)), c(0,y[k>=n-Esota], 0), col="blue")

# area under curve for SOTA performance
k_alpha2 = n-x_alpha2
polygon(c(k_alpha2, k[k>=k_alpha2], max(k)), c(0,y[k>=k_alpha2], 0), col="red")

axis(1 , las = 2, at=c(ci_binom[["lower"]]*n, n*theta, 
                       ci_binom[["upper"]]*n, n-Esota, k_alpha2), 
     labels=c(TeX('$\\theta_{alpha/2}$'), TeX('$\\theta$'), TeX('$\\theta_{1-alpha/2}$'), 
              TeX(r'($E \hat{theta}_{SOTA}$)'), TeX('$\\theta^m_{1-alpha/2}$')))

# # # # # # # # # # # # # # # # # # end figure 1 # # # # # # # # # 


# # # # # # # # # # # # # # # # # Figure 2 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Fz = cdf(n,theta,m)

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

# # # # # # # # # # # # # # # # # end figure 2 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # Figure 3 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

fz = pmf(n,theta,m,f0 = T)

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
# # # # # # # # # # # # # # # # # end figure 3 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


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














