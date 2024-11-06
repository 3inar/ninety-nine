# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) May 2024 - mainly copied from SOTA_bootstrap.R

# thetas - probabilities of correct prediction

library(dplyr)
library(magrittr)
library(rvest)
library(purrr)
library(readr)
library(ggplot2)

library(plotly)
library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing

###############################################################################
############################ \section{A kaggle challenge example} #############################
###############################################################################

# We have two examples from Kaggle; Cassava Leaf Disease Classification and Multi-Class Prediction of Obesity Risk. 

casava = T # to distinguish the two data sets

# Data comes from this website:
# https://www.kaggle.com/competitions/cassava-leaf-disease-classification/leaderboard
# https://www.kaggle.com/competitions/playground-series-s4e2/leaderboard
# It has been scraped and saved locally 

if (casava) {
  comb_data <- readRDS('/Users/kajsam/Documents/Casava Webscrape/Casava_kaggle_leadboard_data.RDS')
} else {
  comb_data <- readRDS('/Users/kajsam/Documents/Obesity Webscrape/Obesity_kaggle_leadboard_data.RDS')
}

head(comb_data) # have a look

# Only interested in the private scores, i.e., the independent test set results
dat <- data.frame(theta = comb_data$prv_score, dataset = "test")

# Both are multiclass problems. I will exclude teams with \hat{\theta} <= 1/#classes.
# Information about the number of classes is found on the websites.
c = 5
exl_thresh = 1/c

theta_obs = dat$theta 
m_tot = length(theta_obs)
theta_obs = theta_obs[theta_obs>1/c]
m = length(theta_obs)

# The size of the test set. This information comes from the website, or the test data set. 
if (casava){
  n = 15000 # approximate number, test data set not available
} else {
  n = 13840 # size of test data set
}


# Have a quick look at the histogram
hist(theta_obs, breaks=200, main = paste("Histogram of thetas"), xlab = "probablility of correct class", ylab = "number of teams")

# Figure settings - these might need to be adjusted manually
kslim = c(0.88, 0.916)
n_breks = 175         # adjust for pleasant graphics
whylim = c(0,120)


#################################################################################
####################### \subsection{} #################
###############################################################################

# Calculating the CI for hat{\theta}^max without multiplicity correction
alpha = 0.05
mu = floor(n*max(theta_obs))
ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial

print(sprintf("The %s confidence interval for an estimated accuracy of %.4f with n = %s is (%.4f,%.4f), width = %.4f.",  
              (1-alpha)*100, max(theta_obs), n, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] ))

teams_single95 = length(theta_obs[theta_obs > ci_binom["lower"][1,1]])
print(sprintf("%s teams have accuracies above the lower 95 CI.", 
              teams_single95))


# The SOTA estimation is done by cropping
# the kagge observations, and then a simulation is performed with correlation. 
# The parameters for cropping are adjusted according to the wanted 
# outcome: either that the expected value of the maximum simulated theta is equal 
# to the maximum kaggle theta, or that the upper limit of the 95% CI of the 
# maximum simulated theta is equal to the maximum kaggle theta. CIs are estimated by bootstrapping

rho = 0.6 # correlation coefficient, the number is calculated from Mania (2019)

option = 2

step = (kslim[2]-kslim[1])/n_breks
breks = c(kslim[1], seq(kslim[1]+step, kslim[2], step))

if (option == 1){
  theta_SOTA = max(theta_obs) # Demonstrating that max(theta_obs) is a biased estimate for theta_SOTA 0.91157
  main_title = 'uncropped'
} else if (option == 2){ # crop
  theta_SOTA = breks[which(breks>0.90615)[1]] # parameter adjusted until E(theta_sota) = max(theta_obs). 0.9067->0.9122, 0.906 -> 0.91148, 0.9063->0.91173, 0.90615-> 0.91156
  main_title = paste('Cropped at ', round(theta_SOTA*10000)/10000)
} else if (option == 3){ # crop
  if (casava){
    theta_SOTA = 0.9131
  } else {
    theta_SOTA = 0.90275 # parameter until upper limit of 95% CI = max(theta_obs) 0.90 -> 0.90870, 0.905 -> 0.91384, 0.9025 -> 0.91132, 0.903->0.91176, 0.90275 -> 0.91158
  }
   main_title = 'crop upper CI'
}



hist(theta_obs[theta_obs>kslim[1]], breaks=breks, freq = T,
     main = paste("Test set performance"), xlab = "probablility of correct class", 
     ylab = "number of teams", xlim = kslim, ylim = whylim, col = "gray20", border = "gray20")
par(new = T)
plot(c(ci_binom["lower"][1,1], ci_binom["upper"][1,1]), c( -1,-1), "l", col = "red", xlim = kslim, ylim = whylim,
     ylab = "number of teams", xlab = "probablility of correct class", )

# Creating the ghost effect for the cropped out observations
# Because of the deceiving nature of the histogram ("poor man's density function"), the bins need to be exactly the same


hist(theta_obs[(theta_obs>kslim[1])&(theta_obs<=theta_SOTA)], breaks=breks[breks<=theta_SOTA], freq = T,
     main = "", xlab = "probablility of correct class", 
     ylab = "number of teams", xlim = kslim, ylim = whylim, col = "gray20", border = "gray20")
par(new = T)
hist(theta_obs[(theta_obs>theta_SOTA)], breaks=breks[breks>=theta_SOTA], freq = T,
     main = main_title, xlab = "probablility of correct class", 
     ylab = "number of teams", xlim = kslim, ylim = whylim, col = "gray50", border = "gray50")


# Simulate dependency
theta_0 = theta_SOTA # theta_0 is the probability of correct prediction for the leading classifier. adjust to E(theta_SOTA)? 

# The theta_obs must be truncated
# with a dependency of rho, the minimum theta_j is 
trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))

# the two solutions for the quadratic equation (does not influence the lower cut-off)
# a = -rho^2*theta_0*(1-theta_0)-theta_0^2
# b = rho^2*theta_0*(1-theta_0)+2*theta_0
# c = -theta_0^2
# x1 = (-b+sqrt(b^2-4*a*c))/(2*a)
# x2 = (-b-sqrt(b^2-4*a*c))/(2*a)

print(sprintf("With an estimated SOTA of %s and a correlation coefficient of %s, the minimum value for theta_j is %.4f.",
              theta_0, rho, trunc_min))

m = length(theta_obs[(theta_obs > trunc_min)]) # number of teams left
print(sprintf("The number of teams above the lower threshold is %s.",
              m))

# Have a quick look at the truncated histogram
trunc_dat = theta_obs[theta_obs>trunc_min]
# hist(trunc_dat, breaks=n_breks, main = paste(m, "theta's truncated from below"), xlab = "theta's", ylab = "number of teams")


# bootstrap sampling from kaggle data lower than theta_SOTA, and higher that the lower cut-off
trunc_dat = theta_obs[(theta_obs > trunc_min)&(theta_obs <= theta_SOTA)]

source("dep_nonid_pmf_fun.R") # for the function 'dep_nonid_pmf' - simulated pmf

rep = 10000 # 10 000 number of repetitions. high number gives low variation. 
if (option == 1){
  B = 1
} else {
  B = 3 #1000 # number of bootstraps. high number gives stable error estimation
}

lowerCI = numeric(B) # lower limit of confidence interval
upperCI = numeric(B) # upper limit of confidence interval

E_SOTA = numeric(B) # the mean (over rep) minimum number of failures among all classifiers
V_SOTA = numeric(B) # the variance of the expected value of minimum number of failures among all classifiers
teamsSOTA = numeric(B) # mean number of teams above SOTA
qq_x = matrix(nrow = B, ncol = m) # keep some for the qq-plot
boot_x = matrix(nrow = B, ncol = rep) # keep some for bootstrap illustration

tic()
for (b in 1:B){
  if (B == 1){
    theta_vec = trunc_dat
  } else {
    theta_vec = sample(x=trunc_dat, size=m, replace=TRUE) # bootstrap from truncated kaggle observations
  }

  # minimum number of wrong classifications among all classifiers, vectors of length rep
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
  
  # The bounds of the 95% confidence interval
  sort_min_dep = sort(X$min_fail) # sort the minimum number of failures
  lowerCI[b] = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound
  upperCI[b] = sort_min_dep[(alpha/2)*rep] # find the alpha/2 upper bound
  
  E_SOTA[b] = (n-mean(X$min_fail))/n # The expected value
  V_SOTA[b] = mean(X$min_fail*X$min_fail) - mean(X$min_fail)*mean(X$min_fail) # The variance
  
  teamsSOTA[b] = mean(X$teamsSOTA)
  
  qq_x[b,] = X$x_fail
  boot_x[b,] = X$min_fail  
  
}

if (B>1){
  hist(V_SOTA)
}


qq_x = as.vector(qq_x)
boot_x = as.vector(boot_x)

tid = toc()
print(sprintf("%s repetitions, %s bootstraps took %s",
              rep, B, tid))

# have a look at the histograms
if (B>10){
  par(mfrow = c(2,2))
  # a bootstrap
  hist(theta_vec[theta_vec>kslim[1]], breaks=n_breks, xlim = kslim, ylim = whylim, 
       main = 'one bootstrap', xlab = 'bootstrapped accuracies', ylab = 'number of teams')
  # the kaggle 
  hist(theta_obs[theta_obs>kslim[1]], breaks=n_breks, xlim = kslim, ylim = whylim, 
       main = 'kaggle observations', xlab = 'accuracies', ylab = 'number of teams')
  # one realisation 
  theta_real = (n-X$x_fail)/n
  hist(theta_real[theta_real>kslim[1]], breaks=n_breks, xlim = kslim, ylim = whylim, 
       main = paste('one realisation'), xlab = 'simulated accuracies', ylab = 'number of teams')
  # many realisations
  hat_theta = (n-qq_x)/n
  hist(hat_theta[hat_theta>kslim[1]], breaks=breks, xlim = kslim, ylim = whylim, 
       main = paste(B,'realisations'), xlab = 'simulated accuracies', ylab = 'number of teams')
} else {
  # one realisation 
  theta_real = (n-X$x_fail)/n
  hist(theta_real[theta_real>kslim[1]], breaks=breks, xlim = kslim, ylim = whylim, freq = T,
     main = paste('one realisation'), xlab = 'simulated accuracies', ylab = 'number of teams', col="gray20", border = "gray20")
  par(new = T)
  plot(c(ci_binom["lower"][1,1], ci_binom["upper"][1,1]), c( -1,-1), "l", col = "red", xlim = kslim, ylim = whylim,
       ylab = "number of teams", xlab = 'simulated accuracies')
}

# Expected value and its standard deviation
print(sprintf("The expected value of max(theta_obs) is %.5f. The standard deviations is %.6f.", 
              mean(E_SOTA), sqrt(var(E_SOTA))))

# print(sprintf("The mean number of simulated teams have accuracies above the true sota, %s, is %.1f.", theta_SOTA, mean(teamsSOTA)))
teamSOTA = length(theta_obs[theta_obs > theta_SOTA])
print(sprintf("%s teams have accuracies above the true sota, %s.", 
              teamSOTA, theta_SOTA))

# Mean and standard deviation of the 95 CI upper and lower bound
mean_lowerCI = mean(lowerCI)
mean_upperCI = mean(upperCI)
V_CI = mean(upperCI*upperCI) - mean_upperCI*mean_upperCI

print(sprintf("The mean bootstrapped simulated %s confidence interval is (%.5f,%.5f) with %s repetitions and %s bootstraps. The standard deviation of the upper CI is %.7f  %.7f",  
              1-alpha, (n-mean_lowerCI)/n, (n-mean_upperCI)/n, rep, B,  sqrt(V_CI)/n, sqrt(var(upperCI))/n))
teams95CI = length(theta_obs[theta_obs > ((n-mean_lowerCI)/n)])
print(sprintf("%s teams have accuracies above the lower 95 CI.", 
              teams95CI))


# the distribution of SOTA
hist((n-X$min_fail)/n, breaks = 30, xlim=c((n-mean_lowerCI)/n-0.01, (n-mean_upperCI)/n+0.01), freq = F, main = 'max accuracies one bootstrap')
if (B>1){
  hist((n-boot_x)/n, breaks = 50, xlim=c((n-mean_lowerCI)/n-0.01, (n-mean_upperCI)/n+0.01), freq = F, main = paste('max accuracies', B,'bootstraps'))
}
# checking if this corresponds to mean CI
sort_min_dep = sort((n-boot_x[-1])/n)
lowerCIboot = sort_min_dep[(1-alpha/2)*rep*B] # find the alpha/2 lower bound
upperCIboot = sort_min_dep[(alpha/2)*rep*B] # find the alpha/2 upper bound
# it does




