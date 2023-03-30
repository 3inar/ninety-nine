# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) February 2023 updated March 2023

# P     - probability, specified when used
# n     - number of trials, size of a test set
# p     - for prob of success in a binomial distr when this is a constant param
# k     - number of successes, n-k is the number of failures
# gk    - for P(Z=k|n,p) binomal distribution: probability of k successes
# gf    - for P(Z=f|n,1-p) binomal distribution: probability of f failures
# Gxk   - for P(X \leq x) binomial cum distr: prob of at most x successes
# Gxf   - for P(X \leq x) binomial cum distr: prob of at most x failures
# n_c   - number of experiments, number of classifiers
# p_hat - estimated accuracy of a classifier: p_hat = k/n
# Cx    - r.v., cx = 0,1,...,n_c. number of classifiers w p_hat>(n-x)/n
# Fx    - for P(C > 0|n_c,n,p) prob of at least one classifier having at most 
#         x failures (identical to at least n-x successes)
# sota  - state-of-the-art

library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing
library(latex2exp)  # mathematical notation

###############################################################################
############################ Coin flipping example #############################
###############################################################################

# a) Consider $n$ flips of a fair coin, and the outcome is the number of heads, 
# referred to as the number of successes in a binomial distribution. 
# b) Consider $n$ flips of a fair coin, and the number of successes is the 
# number of times a classifier correctly predicts the outcome. 

n = 20      # number of trials
p = 0.5     # probability of success

# k = 15      # a) number of heads
k = 18      # b) number of correct predictions

# The probability of observing at least $k$ successes in one series 
Pk = pbinom(k-1,n,p,lower.tail = F) # P[X>x].

sprintf("The probability of at least %s successes in %s trials is %.5f.",  
        k, n, Pk)

# If the experiment is performed $n_c$ times, we can calculate the 
# probability of at least one classifier achieving at least $k$ successes as 
# follows: The probability of at least one success in a binomial distribution
# with $n_c$ trials, and a probability of success, $gk$ is calculated 
# as a single experiment's probability of $k$ successes in $nl$ trials.

n_c = 1000
Fx = pbinom(0,n_c,Pk, lower.tail = F) # P[X>x]

p_hat = k/n 

sprintf("The probability of at least one out of %s experiments having at least %s out of %s successes is %.5f.",  
        n_c, k, n, Fx)
# In other words
sprintf("The top-ranked accuracy is at least %.4f with a probability of %.4f. The true accuracy is %s.",  
        p_hat, Fx, p)


###############################################################################
############################ Competition example #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $P(\hat{y}=y)$, where
# $y$ is the true label, and $\hat{y}$ is the label predicted by the classifier.
# We will refer to this as the classifier's accuracy and denote it by $p$. 
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is $\hat{P}(\hat{y}=y)$, denoted by $p_hat$.

n = 3000
p = 0.9
mu = n*p # the expected number of correct predictions

# 95\% confidence interval
alpha = 0.05
ci_binom = binom.confint(mu,n,conf.level=1-alpha) # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %s is (%.4f,%.4f).",  
        1-alpha, p, ci_binom[5,"lower"], ci_binom[5,"upper"])

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# The probability of exceeding CI_{up} should be close to alpha/2 = 0.025
k_up = floor(ci_binom[5,"upper"]*n) # number of successes exceeding the CI
P_up = pbinom(k_up,n,p, lower.tail = F) # P[X>x]
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

# If there are $n_c$ teams, each with a classifier with $P(\hat{y}=y)$, what is
# the probability of at least one team achieving at least $ci_binom[5,"upper"]$ 
# correct predictions?
# Following the logic of the coin-flip multiple experiments, we have a binomial 
# distribution with $n_c$ trials and probability of success (exceeding 95% CI) 
# is by definition $alpha/2$. 

Palpha2 = pbinom(0,n_c,alpha/2, lower.tail = F) # P[X>x].
sprintf("The probability of at least one out of %s teams exceeding the upper limit of the CI is 
        %s.",  
        n_c, Palpha2)
sprintf("Notice that the accuracy, p, is not part of the equation.")

# # # # # # # # # # # # # # # # # Simulations # # # # # # # # # # # # # # # # # # # # # # # 

rep = 1000 # one million repetitions takes about 80 seconds
# I still don't expect any failures

k_up = floor(ci_binom[5,"upper"]*n) # number of successes exceeding the CI

n_success = numeric(rep) # pre-allocate
tic()
for (i in 0: rep){
  sim_success = rbinom(n_c, n, p) # the number of successes for each team
  n_success[i] = length(which(sim_success>k)) # how many are above CI_{up}?
}

which(n_success == 0) # any repetitions with all below CI_{up}?
toc()
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 


# If there are $n_c$ classifiers, what must be the prob, $Palpha2$, of each 
# classifier exceeding $k-1$ successes, for the probability of at least one
# classifier exceeding $k-1$ successes to be P(C > 0) = alpha/2?
# We can calculate $k$ from $Palpha2$.

# This is easily calculated:
# P(C > 0) = 1 - P(C = 0) = alpha/2
# binom coeff is 1 since c = 0, Palpha2^0 = 1
# and we are left with (1-Palpha2)^(n_c) = 1 - alpha\2

Palpha2 = 1 - (1-alpha/2)^(1/n_c)

sprintf("For the probability of at least one classifier to achieve at least $k$ correct predictions to be 
        alpha/2=%s, the prob for each classifier achieving at least $k$ correct predictions must be %.8f.",  
        alpha/2, Palpha2)

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
P = pbinom(0,n_c,Palpha2, lower.tail = F) # P[C>0] = 0.025
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

# If the probability of correctly predicting at least k out of n data points 
# is Palpha2. What is k when p = 0.9?

k_alpha2 = qbinom(Palpha2, n, p, lower.tail = F)

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
Palpha2_discr = pbinom(k_alpha2,n,p,lower.tail = F) # < Palpha2
P_discr = pbinom(0,n_c,Palpha2_discr, lower.tail = F) # = 0.02475 < P

# k_alpha2_cont < k_alpha2
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 

p_hat = k_alpha2/n # this is \hat{p}

sprintf("In summary, if there are %s teams, each with %s accuracy, and a test set of size %s, there is a      
        probability of %s that at least one team will achieve %s correct predictions, corresponding to an        
        estimated accuracy of %.4f.",  
        n_c, p, n, alpha/2, k_alpha2, p_hat)

# # # # # # # # # # # # # # # # # Simulations # # # # # # # # # # # # # # # # # # # # # # # 
# I hope to recreate alpha/2
rep = 10000 # one million repetitions takes about 80 seconds

n_success = numeric(rep)
tic()
for (i in 1: rep){
  sim_success = rbinom(n_c, n, p) # the number of successes for each team
  n_success[i] = length(which(sim_success>k_alpha2)) # how many are above k_alpha2?
}

sum(n_success)/rep # fraction of successes, should around alpha/2 = 0.025
# there is no continuity correction, so I thought I would get < 0.025
toc()
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 


# When the true accuracy is $0.9$, there is a non-negligible (at least in the 
# common statistical significance sense) probability that the top-ranked team
# will have an estimated accuracy of at least 0.9213, which is often referred to 
# as state-of-the-art (SOTA) performance. 

# What is the probability of beating the SOTA for a method with significantly better
# accuracy?

P_beat_sota = pbinom(k_alpha2, n, ci_binom[5,"upper"], lower.tail = F) # P[X>x]

sprintf("The probability of achieving an estimated accuracy better than SOTA, %.4f, for a classifier with         
        significantly better accuracy, p=%.4f, is just %.4f.",
        p_hat, ci_binom[5,"upper"], P_beat_sota)

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

x = (mu-60):(mu+90) # this is the plot range, adjust to your liking

y = dbinom(x, n, p) # probability of x successes in n trials
plot(x, y, type='h', xlab = '', ylab = '', axes=F) # histogram

par(new=TRUE) # new plot in same window

y = dbinom(x, n, ci_binom[5,"upper"]) # significantly better classifier
plot(x, y, type = 'l', col = "red", xlab = 'estimated probility of success', ylab = '', axes=F)

# area under curve for SOTA performance
polygon(c(k_alpha2, x[x>=k_alpha2], max(x)), c(0,y[x>=k_alpha2], 0), col="red")

axis(1 , las = 2, at=c(x[1], ci_binom[5,"lower"]*n, n*p, 
                       ci_binom[5,"upper"]*n, k_alpha2, tail(x,1)), 
     labels=c(as.character(x[1]/n), TeX('$\\p_{alpha/2}$'), as.character(p), TeX('$\\p_{1-alpha/2}$'), 
              TeX('$\\p^m_{1-alpha/2}$'), as.character(tail(x,1)/n)))

# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


###################################### Multiple classifiers ####################################################

# First, let's have a look at probabilities for at least one team having exactly 
# k successes out of n = 3000.

Psota = numeric(n+1)

for (k in 0: n){
  Pk = dbinom(k,n,p)
  i = k+1
  Psota[i] = pbinom(0,n_c,Pk,lower.tail = F) # P[C>0]
}

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

x = (mu-70):(mu+70) # this is the plot range, adjust to your liking

plot(x,Psota[x+1], type = 'h', axes=F, 
     xlab = 'estimated accuracy/number of successes', 
     ylab = 'prob of at least one classifier > hat{p}')
xax = seq(2650,2750, 10)
klab = xax 
plab = round(1000*xax/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
axis(2, las = 2)
axis(3, las = 2, at=xax, labels = as.character(klab))

# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# What does this figure illustrate, and why isn't it a pmf?
# Looking at a single classifier, the events x = k and x = k', k \neq k' are mutually exclusive.
# Looking at all classifiers, the events P(C>0|k) and P(C>0|k') are not, because we can simultaneously have at
# least one classifier with k successes and another one with k' successes

# What about exactly one classifier with exactly k successes?

Psota = numeric(n+1)

for (k in 0: n){
  Pk = dbinom(k,n,p)
  i = k+1
  Psota[i] = dbinom(1,n_c,Pk) 
}
# This is also not a pmf, because we can simultaneously have exactly one classifier with k successes and
# exactly one classifier with k' \neq k successes, out of n_c classifiers. 

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

x = (mu-80):(mu+80) # this is the plot range, adjust to your liking

plot(x,Psota[x+1], type = 'h', axes=F, 
     xlab = 'estimated accuracy/number of successes', 
     ylab = 'prob of at least one classifier > hat{p}')
xax = seq(x[1],tail(x,1), 10)
klab = xax 
plab = round(1000*xax/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
axis(2, las = 2)
axis(3, las = 2, at=xax, labels = as.character(klab))
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Fun figure. Interestingly, the left part has the same shape as the pmf. I have to think about that a bit. 


################################# Cumulative distribution function SOTA ############################################

# Now, let's have a look at probabilities for at least one classifier having at least 
# k successes out of n = 3000. This gives the lower limit of Psota with a
# certain probability

Fx = numeric(n+1)

# use i for counting, since x = 0, ..., n

# k is the number of successes, then x = n-k is the number of failures
# probability of failure = 1-p

for (x in 0:n){
  Gxf = pbinom(x,n,(1-p)) #  P[X \leq x]
  i = x+1
  Fx[i] = pbinom(0,n_c,Gxf,lower.tail = F) #  P[C>0]
}

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# should be same as n-k_alpha2 = 236
k_sota_est = which(Fx>0.025)[1]-1
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# the whole range, not very much information
plot(0:n,Fx, type = 'l', xlab = 'number of failures', 
     ylab = 'probability of at least one team')

# zooming in, and it gets more interesting, discreet curve 
x = 200:300
plot(x,Fx[x+1], type = 's', xlab = 'number of failures/estimated accuracy', 
     ylab = 'F(z)', axes = F)

xax = seq(x[1],tail(x,1), 10)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
axis(2, las = 2)
axis(3, las = 2, at=xax, labels = as.character(klab))

# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# The window where Fx goes from 0 to 1 is very narrow. I don't know what it means

k_50 = x[1]+which(Fx<0.5)[1]
CI_50 = binom.confint(n-k_50,n,conf.level=1-alpha)[5,] # CI for binomial
CI_width = CI_50[1,"upper"]-CI_50[1,"lower"]
k_0 = which(Fx>0.00001)[1]
k_1 = which(Fx>0.99999)[1]
range01 = (k_1-k_0)/n

sprintf("The width of a 95 CI is %.2f, which is approx same range as Fx=1 - Fx=0 = %.2f.",  
        CI_width, range01)

################# probability mass function ##########################

fx = numeric(n)
fx[1:n] = Fx[2:(n+1)]-Fx[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

plot(1:n,fx, type = 's')

# need to zoom in
x = 200:300
plot(x,fx[x], type = 'h', xlab = 'number of failures/estimated accuracy', 
     ylab = 'f(z)', axes = F)
xax = seq(x[1],tail(x,1), 10)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
yax = seq(0,max(fx)+0.01,0.02)
axis(2, las = 1, at=yax, labels = as.character(yax))
axis(3, las = 2, at=xax, labels = as.character(klab))
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
sum(fx)
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Fx = numeric(n+1)

for (x in 0:n){
  px = pbinom(x,n,1-p)
  Fx[x+1] = 1 - (1 - px)^n_c
}
fx = Fx[2:(n+1)]-Fx[1:n]

x = 200:300
plot(x,Fx[x+1], type = 's')
plot(x,fx[x], type = 'h', xlab = '', ylab = 'probability mass', axes = F)
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


####################### Expected value and variance ###################################

# expected value
Eterm = numeric(n+1)
for (x in 1:n){
  i = x
  Eterm[i] = x*fx[i]
}

Esota = sum(Eterm)
Esota_p = 1-Esota/n

# variance
vterm = numeric(n+1)
for (x in 1:n){
  i = x
  vterm[i] = x^2*fx[i]
}
esquare = sum(vterm)

Vsota = esquare - Esota^2

sprintf("The expected number of failures is %.4f, with a variance of %.4f.",
        Esota, Vsota)

# beat esota?

beat_esota = pbinom(n-Esota, n, ci_binom[5,"upper"], lower.tail = F) # P[X>x]

sprintf("The probability of achieving an estimated accuracy better than E(SOTA)=%.4f, for a classifier with      
        significantly better accuracy, %.4f, is %.4f.",
        Esota_p, ci_binom[5,"upper"], beat_esota)

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

x = (mu-60):(mu+90) # this is the plot range, adjust to your liking

y = dbinom(x, n, p) # probability of x successes in n trials
plot(x, y, type='h', xlab = '', ylab = '', axes=F) # histogram

par(new=TRUE) # new plot in same window

y = dbinom(x, n, ci_binom[5,"upper"]) # significantly better classifier
plot(x, y, type = 'l', col = "red", xlab = '', ylab = '', axes=F)

# area under curve for expected SOTA performance

polygon(c(n-Esota, x[x>=n-Esota], max(x)), c(0,y[x>=n-Esota], 0), col="blue")

# area under curve for SOTA performance
polygon(c(k_alpha2, x[x>=k_alpha2], max(x)), c(0,y[x>=k_alpha2], 0), col="red")


axis(1 , las = 2, at=c(ci_binom[5,"lower"]*n, n*p, 
                       ci_binom[5,"upper"]*n, n-Esota, k_alpha2), 
     labels=c(TeX('$\\p_{alpha/2}$'), 'p', TeX('$\\p_{1-alpha/2}$'), 
              TeX('$\\Ep_{SOTA}$'), TeX('$\\p^m_{1-alpha/2}$')))

######################### varying m, n, p ####################################

m = 1000
n = 10000
p = 0.9


Fx = numeric(n+1)

# use i for counting, since x = 0, ..., n

# k is the number of successes, then x = n-k is the number of failures
# probability of failure = 1-p

for (x in 0:n){
  Gxf = pbinom(x,n,(1-p)) #  P[X \leq x]
  i = x+1
  Fx[i] = pbinom(0,m,Gxf,lower.tail = F) #  P[C>0]
}

fx = numeric(n)
fx[1:n] = Fx[2:(n+1)]-Fx[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)

####################### Expected value and variance ###################################

# expected value
Eterm = numeric(n+1)
for (x in 1:n){
  i = x
  Eterm[i] = x*fx[i]
}

Esota = sum(Eterm)
Esota_p = 1-Esota/n

# variance
vterm = numeric(n+1)
for (x in 1:n){
  i = x
  vterm[i] = x^2*fx[i]
}
esquare = sum(vterm)

Vsota = esquare - Esota^2

sprintf("The expected p_sota is %.4f, with a standard deviation of %.6f, for m=%s, n=%s, p=%s.",
        (n-Esota)/n, sqrt(Vsota)/n, m, n, p)




# need to zoom in
x = 50:700
plot(x,fx[x], type = 'l', col = 'blue', xlab = 'number of failures/estimated accuracy', 
     ylab = 'f(z)', axes = F)
xax = seq(x[1],tail(x,1), 10)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
yax = seq(0,max(fx)+0.01,0.02)
axis(2, las = 1, at=yax, labels = as.character(yax))
axis(3, las = 2, at=xax, labels = as.character(klab))
par(new=TRUE) # new plot in same window

