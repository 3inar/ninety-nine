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
############################ Competition example - dependency #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $P(\hat{y}=y)$, where
# $y$ is the true label, and $\hat{y}$ is the label predicted by the classifier.
# We will refer to this as the classifier's accuracy and denote it by $p$. 
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is $\hat{P}(\hat{y}=y)$, denoted by $p_hat$.

m = 1000
n = 3000
p = 0.9
mu = n*p # the expected number of correct predictions

# Simulate dependency

q0 = numeric(n)
q0[1:mu] = 1

rho = 0.6

p_dep = p + rho*(1-p)

p1 = 1-p_dep
p0 = (mu/(n-mu))*(1-p_dep)



rep = 100000          # 60 sec for a thousand, 9,000 sec for 100,000

min_dep = numeric(rep)
min_indep = numeric(rep)
tic()
for (ell in 1:rep){
  diff = numeric(m)
  x_dep = numeric(m)  # number of failures

  for (k in 1:m){
    flip1 = rbinom(mu,1,p1)
    flip0 = rbinom(n-mu,1,p0)
  
    x_dep[k] = n-sum(abs(q0-c(flip1,flip0)))
  
  #  diff[k] = sum(flip1)-sum(flip0)
  }

  min_dep[ell] = min(x_dep)
  

  #plot(1:m,sort(x_dep), col = 'red', type = 's', axes = F, xlab = '', ylab ='')
  

  x_indep = rbinom(m,n,1-p)
  min_indep[ell] = min(x_indep)
  
  
  #par(new=TRUE) # new plot in same window

  #plot(1:m,sort(x_indep), type = 's', xlab = '', ylab ='number of failures')
  
}

toc()

# print(min(x_dep))
# print(min(x_indep))

histbreaks = seq(min(c(x_dep,x_indep)), max(c(x_dep,x_indep))+8,10)

hist(x_dep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,300))
hist(x_indep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,300))

sort_min_dep = sort(min_dep)
x_dep_fails = sort_min_dep[0.025*rep]

(n-x_dep_fails)/n


sort_min_indep = sort(min_indep)
x_indep_fails = sort_min_indep[0.025*rep]
(n-x_indep_fails)/n

histbreaks = seq(210,273,3)

hist(min_dep, xlab = 'minimum number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,30000))

hist(min_indep, xlab = 'minimum number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,30000))

