# Closed form 
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) March 2023

# Checking that our analytical results of the closed form are correct

###############################################################################
############################ 3.1 Nomenclature #############################
###############################################################################

# n     - number of trials/size of test set
# m     - number of experiments/number of classifiers
# p     - prob of success/prob of correct prediction
# k     - number of successes
# X     - r.v., number of failures; n-k

# P_    - probability, specified when used

# px    - probability of x failures in one experiment
# Px    - probability of at most x failures in one experiment

# Cx    - r.v., number of experiments with at most x failures

# p_hat - accuracy of a classifier: p_hat = (n-x)/n

# p_SOTA - max_j p_hat_j

# Z     - r.v., number of failures on at least one classifier

# Fz    - cdf of Z: P(C_z > 0|m,n,p) prob of at least one classifier having at most 
#         z failures (identical to at least n-z successes)
# fz    - pmf of Z

# alpha - significance level

# SOTA - state-of-the-art


library(tictoc)     # for timing

###############################################################################
############################ Coin flipping example #############################
###############################################################################

# Using the coin flipping example with low m, so that it is possible to calculate
# the binomial coefficient

# Consider $n$ flips of an unfair coin, and the number of successes is the 
# number of times a classifier correctly predicts the outcome. 

n = 20      # number of trials
p = 0.6     # probability of success
m = 20    # number of classifiers

# z is the number of failures, from 0 to n
# probability of failure = 1-p

Fz = numeric(n+1)

for (z in 0:n){
  Pz = pbinom(z,n,(1-p)) #  P[X \leq x]
  i = z+1 # cannot have 0th entry 
  Fz[i] = pbinom(0,m,Pz,lower.tail = F) #  P[C>0]
}

# figure
x = 0:n
plot(x,Fz, type = 's', xlab = 'number of failures/estimated accuracy', 
     ylab = 'probability of at least one classifier', axes = F)

xax = seq(x[1],tail(x,1), 2)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
axis(2, las = 2)
axis(3, las = 2, at=xax, labels = as.character(klab))

fz = numeric(n)
fz[1:n] = Fz[2:(n+1)]-Fz[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)


x = 1:n
plot(x,fz, type = 'h', col='red', xlab = '', 
     ylab = 'probability mass', axes = F)
xax = seq(x[1],tail(x,1), 2)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
yax = seq(0,max(fz)+0.01,0.02)
axis(2, las = 1, at=yax, labels = as.character(yax))
axis(3, las = 2, at=xax, labels = as.character(klab))

########################## closed form ##############################

f = numeric(n)

for (x in 1:n){
  term_c = numeric(m)
  for (c in 1:m){
    p_ell = numeric(x) 
    for (ell in 0:x-1){
      j = ell+1
      p_ell[j] = dbinom(ell, n, 1-p)
    }
    term_c[c] = (-1)^c*choose(m,c)*(1-sum(p_ell))^(m-c)*dbinom(x, n, 1-p)^c
  }
  i = x
  f[i] = -sum(term_c)
}

par(new=TRUE) # new plot in same window
x = 1:n
plot(x,f, type = 'h', xlab = '', 
     ylab = 'probability mass', axes = F)

sum(abs(fz-f)) # sum of absolute differences is ^-15. that's ok

# expected value (no closed'er form)

term_e = numeric(n)

for (x in 1:n){
  i = x
  term_e[i] = x*f[i] 
}

EX = sum(term_e)

# simulate expected value

n_sim = 1000000
x_sim = numeric(n_sim)

for (i in 1:n_sim){
  x_sim[i] = min(rbinom(m,n,1-p))
}
mean(x_sim)

###############################################################################
########### 3.3 The probability distribution of p_SOTA ########################
###############################################################################

# Here are some of the struggles I had before I was able to get it right

n = 3000
p = 0.9
m = 1000
mu = n*p # the expected number of correct predictions
alpha = 0.05

# First, let's have a look at probabilities for at least one team having exactly 
# x failures out of n.

Cx = 1 # one team

P = numeric(n+1)
for (x in 0: n){
  px = dbinom(x,n,1-p)
  i = x+1
  P[i] = pbinom(Cx-1,m,px,lower.tail = F) # P[C>0]
}

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

x = (n-mu-70):(n-mu+70) # this is the plot range, adjust to your liking

plot(x,P[x+1], type = 'h', axes=F, 
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

# What about exactly one classifier with exactly x failures?

P = numeric(n+1)

for (x in 0: n){
  px = dbinom(x,n,1-p)
  i = x+1
  P[i] = dbinom(1,m,px) 
}
# This is also not a pmf, because we can simultaneously have exactly one classifier with x failures and
# exactly one classifier with x' \neq x failures, out of m classifiers. 

# # # # # # # # # # # # # # # # # Figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

x = (n-mu-80):(n-mu+80) # this is the plot range, adjust to your liking

plot(x,P[x+1], type = 'h', axes=F, 
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



