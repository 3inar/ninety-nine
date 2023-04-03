# Closed form 
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) March 2023

# P     - probability, specified when used
# n     - number of trials, size of a test set
# p     - for prob of success in a binomial distr when this is a constant param
# n_c   - number of experiments, number of classifiers
# Fx    - for P(C > 0|n_c,n,p) prob of at least one classifier having at most 
#         x failures (identical to at least n-x successes)
# fx    - pmf of Fx

library(tictoc)     # for timing

###############################################################################
############################ Coin flipping example #############################
###############################################################################

# Consider $n$ flips of an unfair coin, and the number of successes is the 
# number of times a classifier correctly predicts the outcome. 

n = 20      # number of trials
p = 0.6     # probability of success
n_c = 20    # number of classifiers

# x is the number of failures, from 0 to n

Fx = numeric(n+1)

# if k is the number of successes, then ell = n-k is the number of failures
# probability of failure = 1-p

for (x in 0:n){
  Gxf = pbinom(x,n,(1-p)) #  P[X \leq x]
  i = x+1 # cannot have 0th entry 
  Fx[i] = pbinom(0,n_c,Gxf,lower.tail = F) #  P[C>0]
}

# figure
x = 0:n
plot(x,Fx, type = 's', xlab = 'number of failures/estimated accuracy', 
     ylab = 'probability of at least one classifier', axes = F)

xax = seq(x[1],tail(x,1), 2)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
axis(2, las = 2)
axis(3, las = 2, at=xax, labels = as.character(klab))

fx = numeric(n)
fx[1:n] = Fx[2:(n+1)]-Fx[1:n] # this is how pmf is defined: f(x) = F(x)-F(x-1)


x = 1:n
plot(x,fx, type = 'h', col='red', xlab = '', 
     ylab = 'probability mass', axes = F)
xax = seq(x[1],tail(x,1), 2)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
yax = seq(0,max(fx)+0.01,0.02)
axis(2, las = 1, at=yax, labels = as.character(yax))
axis(3, las = 2, at=xax, labels = as.character(klab))

########################## closed form ##############################

f = numeric(n)

for (x in 1:n){
  term_c = numeric(n_c)
  for (c in 1:n_c){
    p_ell = numeric(x) 
    for (ell in 0:x-1){
      j = ell+1
      p_ell[j] = dbinom(ell, n, 1-p)
    }
    term_c[c] = (-1)^c*choose(n_c,c)*(1-sum(p_ell))^(n_c-c)*dbinom(x, n, 1-p)^c
  }
  i = x
  f[i] = -sum(term_c)
}

par(new=TRUE) # new plot in same window
x = 1:n
plot(x,f, type = 'h', xlab = '', 
     ylab = 'probability mass', axes = F)

sum(abs(fx-f)) # sum of absolute differences is ^-15. that's ok

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
  x_sim[i] = min(rbinom(n_c,n,1-p))
}
mean(x_sim)

