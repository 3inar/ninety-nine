# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) 

library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing
library(latex2exp)  # mathematical notation

###############################################################################
############################ 3.1 Notation #############################
###############################################################################

# n         - number of trials/size of test set
# theta     - prob of success/prob of correct prediction
# k         - number of successes
# X         - r.v., number of failures; n-k
# theta_hat = (n-X)/n - estimator of theta
# px        - probability of x failures in one experiment
# Px        - probability of at most x failures in one experiment

# m         - number of experiments/number of classifiers
# Cx        - r.v., number of experiments with at most x failures
# Z         - r.v., the most number of failures on at least one classifier

# theta_SOTA_hat - estimator of theta_SOTA 

# P_    - probability, specified when used

# Fz    - cdf of Z: P(C_z > 0|m,n,theta) prob of at least one classifier having at most z failures 
# the probability that the most wrong predictions in at least one classifier is less than or equal to a number $z$

# fz    - pmf of Z

# alpha - significance level

# SOTA - state-of-the-art

###############################################################################
############################ 3.2 Two coin-flip examples #############################
###############################################################################

# a) Consider $n$ flips of a fair coin, and the outcome is the number of heads, 
# referred to as the number of successes in a binomial distribution. 

n = 20        # number of trials
theta = 0.5   # probability of success
m = 100       # number of experiments

k = 15        # number of heads (successes)
x = n-k       # number of failures

# The probability of observing at least $n-x$ successes in one series, equivalent
# to observing at most $x$ failures P(X \leq x)

Px = pbinom(x,n,1-theta) # cdf of a binomial distribution
px = dbinom(x,n,1-theta) # pmf of a binomial distribution

sprintf("The probability of %s heads or more in %s trials is %.5f.",  
        k, n, Px)
sprintf("The probability of hat{\theta} = %s or above in %s trials is %.5f.",  
        (n-x)/n, n, Px)
sprintf("The probability of hat{\theta} = %s in %s trials is %.5f.",  
        (n-x)/n, n, px)

# The probability of observing at most $x$ failures at least once, in other words, more than 0
Cx = 1
P_Cx = pbinom(Cx-1,m,Px, lower.tail = F) # P[X>x] 
# equivalent to
P_Cx = 1 - dbinom(Cx-1,m,Px)

sprintf("The probability of at least %s out of %s experiments having at most %s out of %s failures is %.5f.",  
        Cx, m, x, n, P_Cx)
sprintf("The probability of at least %s out of %s experiments having hat{\theta} = %s or above in %s trials is %.5f.",  
        Cx, m, (n-x)/n, n, P_Cx)

# Equivalently, we can do the probability of observing $x$ failures at least once
Cx = 1
P_Cx = pbinom(Cx-1,m,px, lower.tail = F) # P[X>x]

sprintf("The probability of at least %s out of %s experiments having %s out of %s failures is %.5f.",  
        Cx, m, x, n, P_Cx)


# b) Consider $n$ flips of a fair coin, and the number of successes is the 
# number of times a classifier correctly predicts the outcome. 

n = 20      # number of trials
theta = 0.5 # probability of success
m = 1000    # number of experiments

k = 18      # number of correct predictions
x = n-k     # number of failures (wrong predictions)

# The probability of observing at most $x$ failures in one series 
Px = pbinom(x,n,1-theta) 

sprintf("The probability of at most %s failures in %s trials is %.5f.",  
        x, n, Px)

# If the experiment is performed $m$ times, we can calculate the 
# probability of at least one classifier achieving at most $x$ failures as 
# follows: The probability of at least one success in a binomial distribution
# with $m$ trials, and a probability of success, $Px$ is calculated 
# as a single experiment's probability of at most $x$ failures in $n$ trials.
Cx = 1
P_Cx = pbinom(Cx-1,m,Px, lower.tail = F) # P[X>x]


theta_hat = (n-x)/n 

sprintf("The probability of at least %s out of %s experiments having at most %s out of %s failures is %.4f.",  
        Cx, m, x, n, P_Cx)
# In other words
sprintf("The top-ranked accuracy is at least %.4f with a probability of %.4f. The true accuracy is %s.",  
        theta_hat, P_Cx, theta)

################################################### Simulations ###########################################

rep = 10000 # We recommend 10 million repetitions, because of the small Px = 0.00020 in b). It takes about 900 seconds

# a) Consider $n$ flips of a fair coin, and the outcome is the number of heads, 
# referred to as the number of successes in a binomial distribution. 

n = 20      # number of trials
theta = 0.5     # probability of success
m = 100    # number of experiments

k = 15      # number of heads (successes)
x = n-k     # number of failures


# Draw at random the number of failures in $m$ independent experiments, each having $n$ trials and probability $theta$
x_vec = rbinom(m, n, 1-theta)

# Have a quick look
hist(x_vec)


# 10 million repetitions takes 120 seconds 
n_success = numeric(rep) # pre-allocate empty vector


tic()

for (i in 0: rep){
  x_vec = rbinom(m, n, 1-theta) # the number of failures in each of the m experiment
  
  # The definition of success in the multiple experiment is for a single trial to have
  # X \leq x failures
  n_success[i] = length(which(x_vec<=x)) # how many experiments have at most x failures
}

sim_P_Cx = length(which(n_success > 0))/rep # the number of simulations with success outcome
toc()

# The probability of observing at most $x$ failures in one series 
Px = pbinom(x,n,1-theta) 
# The probability of at least one out of $m$ trials resulting in at most $5$ failures is 
Cx = 1
P_Cx = pbinom(Cx-1,m,Px, lower.tail = F) # P[X>x]

sprintf("The theoretical probability of at least %s out of %s trials resulting in at most %s failures is %.5f.",  
        Cx, m, x, P_Cx)

sprintf("The simulated result for %s simulations of at least %s out of %s trials resulting in at most %s failures is %.5f.",  
        rep, Cx, m, x, sim_P_Cx)


# b) Consider $n$ flips of a fair coin, and the number of successes is the 
# number of times a classifier correctly predicts the outcome. 

n = 20      # number of trials
theta = 0.5     # probability of success
m = 1000    # number of experiments

k = 18      # number of correct predictions
x = n-k     # number of failures

# Draw at random the number of failures in $m$ independent experiments, each having $n$ trials and probability $theta$
x_vec = rbinom(m, n, 1-theta)

# Have a quick look
hist(x_vec)

# 10 million repetitions takes 860 seconds 
n_success = numeric(rep) # pre-allocate empty vector


tic()

for (i in 0: rep){
  x_vec = rbinom(m, n, 1-theta) # the number of failures in each of the m experiment
  
  # The definition of success in the multiple experiment is for a single trial to have
  # X \leq x failures
  n_success[i] = length(which(x_vec<=x)) # how many experiments have at most x failures
}

sim_P_Cx = length(which(n_success > 0))/rep # the number of simulations with success outcome
toc()


# The probability of observing at most $x$ failures in one trial 
Px = pbinom(x,n,1-theta) 
# The probability of at least one out of $m$ trials resulting in at most $5$ failures is 
Cx = 1
P_Cx = pbinom(Cx-1,m,Px, lower.tail = F) # P[X>x]

sprintf("The theoretical probability of at least %s out of %s trials resulting in at most %s failures is %.5f.",  
        Cx, m, x, P_Cx)

sprintf("The simulated result for %s simulations of at least %s out of %s trials resulting in at most %s failures is %.5f.",  
        rep, Cx, m, x, sim_P_Cx)

############################# Probability functions #########################


# comes with the functions 'cdf', 'pmf', 'expect' and 'variance'
source("ProbDistr_thetaSOTA.R") 

# Parameters from coin-flip experiment a)
n = 20
theta = 0.5
m = 100

# In coin-flip experiment a), we investigated for z = 5, and got P = 0.87646
z = 5
Fz = cdf(n, theta, m) # the cdf of the coin-flip experiment a)
sprintf("The probability of at most %s out of %s failures in at least one out of %s classifiers is %.5f.",  
        z, n, m, Fz[z+1])


# Parameters from coin-flip experiment b)
n = 20
theta = 0.5
m = 1000

Fz = cdf(n, theta, m) # the cdf of the coin-flip experiment b)

# From coin-flip experiment b), we have that P(\hat{theta}_SOTA > 0.9) = 0.1823. 
z = n*(1-0.9) # r gets all confused, and as.integer(z) becomes 1. wtf?
z = 2
sprintf("The probability of at most %s out of %s failures in at least one out of %s classifiers is %.5f.",  
        z, n, m, Fz[z+1])



