# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) 

# Contains the following functions: 

# cdf       - cumulative distribution function
# pmf       - probability mass function
# expect    - expectation
# variance  - variance


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
########### 3.3 The probability distribution of p_SOTA ########################
###############################################################################

################################# Cumulative distribution function SOTA ############################################

# Z is the most failures in at least one classifier.
# If Z <= z in at least one classifier, then \hat{\theta}_j = (n-x_j)/n >= (n-z)/n,
# and we have that \hat{\theta}_SOTA >= (n-z)/n.
# Let's have a look at probabilities for at least one classifier having at most 
# z failures. This gives the cumulative distribution function


# The cumulative distribution function of Z
cdf <- function(n, theta, m){   
  Fz = numeric(n+1)             # pre-allocate

  for (z in 0:n){
    # The individual probabilities of each classifier P(X\leq z|n,1-theta) having at most z wrong prediction. 
    # This is the definition of 'success' for a Bernoulli r.v.
    Pz = pbinom(z,n,(1-theta)) 
    i = z+1                     # use i for counting, since z = 0, ..., n
    Fz[i] = 1 - dbinom(0,m,Pz)  #  P[C>0]
  }
  return(Fz)
}


################# Probability mass function ##########################

# The corresponding probability mass function of Z
pmf <- function(n, theta, m, f0 = F){   
  
  Fz = cdf(n, theta, m) # the cumulative distribution function
  
  fz = numeric(n) # pre-allocate
  fz[2:(n+1)] = Fz[2:(n+1)]-Fz[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)
  
  # We are missing out on f(z=0), and is therefore calculated "by hand", but only if m < 1000, because of
  # the insane binomial coefficient. It is possible to override this by f0=T.
  if (f0==T|| m<1000) {
  # special case of z = 0
  z=0
  P0 = dbinom(z,n,(1-theta)) 
  binom_term = numeric(m) # pre-allocate
  for (k in 1:m){
    binom_term[k] = (-1)^k*choose(m,k)*P0^k
  }
  
  fz[1] = -sum(binom_term)
}
  
  return(fz)
}


# Expectation of Z

expect <- function(m,n,theta){
  
  # Probability mass function
  fz = pmf(n, theta, m, f0 = F)
 
  # Expectation
  
  # Calculate each term z*f(z)
  Eterm = numeric(n)
  for (z in 0:n){
    Eterm[z] = z*fz[z]
  }
  
  EZ = sum(Eterm) # the sum
  
  return(EZ)
  
}

# Variance for Z

variance <- function(m,n,theta){
  
  # Probability mass function
  fz = pmf(n, theta, m, f0 = F)
  
  # Variance
  
  # Calculate the expectation of Z^2
  vterm = numeric(n+1)
  for (z in 1:n){
    vterm[z] = z^2*fz[z]
  }
  esquare = sum(vterm)
  
  # Get the expectation of Z
  EZ = expect(m,n,theta)
  
  VarZ = esquare - EZ^2
  
  return(VarZ)
  
}




###################################### Check-ups ############################

# Parameters from coin-flip experiment a)
n = 20
theta = 0.5
m = 100

Fz = numeric(n+1) # cumulative distribution
Fzz= numeric(n+1) # alternative equation for cumulative distribution, see manuscript

# use i for counting, since z = 0, ..., n

for (z in 0:n){
  # The individual probabilities of each classifier P(X\leq z|n,1-theta) having at most z wrong prediction. 
  # This is the definition of 'success' for a Bernoulli r.v.
  Pz = pbinom(z,n,(1-theta)) 
  i = z+1 
  Fz[i] = 1 - dbinom(0,m,Pz)  #  P[C>0]
  Fzz[i] = 1 - (1-Pz)^m       # they are identical
}

fz = pmf(n, theta, m)
sum(fz) # check that this is 1

# take a look
z = 0:n
plot(z,Fz[z+1], type = 's')
#par(new=TRUE) # new plot in same window
plot(z,fz[z+1], type = 'h', xlab = '', ylab = 'probability mass') #, axes = F

EZ = expect(m,n,theta)
VarZ = variance(m,n,theta)

sprintf("The expected value is %.4f, and the variance is %.4f.",
        EZ, VarZ)










n = 3000
theta = 0.9
m = 1000
mu = n*theta # the expected number of correct predictions
alpha = 0.05

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




