# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) September 2023 

# pretending that the AUCs from kaggle are thetas - probabilities of correct prediction

library(dplyr)
library(magrittr)
library(rvest)
library(purrr)
library(readr)
library(ggplot2)

library(plotly)

library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing


# https://www.kaggle.com/c/siim-isic-melanoma-classification/leaderboard

comb_data <- readRDS('/Users/kajsam/Documents/kaggle-leaderboard-scrape/SIIM-ISIC_Melanoma_kaggle_leadboard_data.RDS')
head(comb_data) # have a look

dat <- data.frame(AUC = comb_data$prv_score, dataset = "test")

# have a look at the histogram
hist(dat$AUC, breaks=200)

###############################################################################
############################ 4.4 A kaggle challenge example #############################
###############################################################################

# With balanced accuracy and binary prediction, we have that accuracy = AUC
theta_obs = dat$AUC 

# We don't have access to the test set labels, so we'll have to estimate based on assumptions. 
# Here is what we do know:

malignant_rate = 584/33126 # in training set ( from https://arxiv.org/ftp/arxiv/papers/2008/2008.07360.pdf)
n_val_test = 10982 # size of test and validation set combined
test_prop = 0.7 # 30/70 split

# Approximations and estimates
n_test = round(test_prop*n_val_test) # approximated test set size - this should be fairly accurate
n_mal = round(test_prop*n_val_test*malignant_rate) # estimated number of malignant cases
# based on equal proportions in validation and test set. this is a strong assumption, but the most reasonable one
n_ben = n_test-n_mal # estimated number of benign cases

sprintf("Estimated number of malignant cases: %s. Estimated test set size: %s. Number of teams: %s.",
        n_mal, n_test, length(theta_obs))

########################### Dependent, non-identical #########################

rho = 0.6 # correlation coefficient, the number is calculated from Mania (2019)
lambda = 4 # shrinking parameter, I adjusted it so that number of teams above SOTA is same for sim and kaggle

option = 4

if (option == 1){
  theta_SOTA = max(theta_obs) # Demonstrating that max(theta_obs) is a biased estimate for theta_SOTA
} else if (option == 2){
  theta_SOTA = 0.9378 # This is the parameter that needs to be adjusted until E(theta_sota) = max(theta_obs)
} else if (option == 3){
  maud = 0.9335 # Still E(theta_sota), shrink instead of crop, 0.923 for lambda = 1.5, 0.93 for lambda = 2, 0.9325 for lambda = 3
  theta_shrink = (theta_obs-maud)/lambda + maud 
  theta_SOTA = max(theta_shrink)
} else if (option == 4){
  theta_SOTA = 0.9295 # This is the parameter that needs to be adjusted until upper limit of 95% CI = max(theta_obs)
} else if (option == 5){ #Still 95% CI = max(theta_obs), shrink 
  maud = 0.92415
  theta_shrink = (theta_obs-maud)/lambda + maud 
  theta_SOTA = max(theta_shrink)
}
theta_SOTA
teamSOTA = length(theta_obs[theta_obs > theta_SOTA])
sprintf("%s teams have accuracies above the true sota, %s.", 
        teamSOTA, theta_SOTA)
# Simulate dependency
theta_0 = theta_SOTA # theta_0 is the probability of correct prediction for the leading classifier. adjust to E(theta_SOTA)? 

# The theta_obs must be truncated
# with a dependency of rho, the minimum theta_j is 
trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))

# the two solutions for the quadratic equation (does not influence the lower cut-off)
a = -rho^2*theta_0*(1-theta_0)-theta_0^2
b = rho^2*theta_0*(1-theta_0)+2*theta_0
c = -theta_0^2

trunc_min1 = (-b+sqrt(b^2-4*a*c))/(2*a)
trunc_min2 = (-b-sqrt(b^2-4*a*c))/(2*a)

# bootstrap sampling from kaggle data lower than theta_SOTA, and higher that the lower cut-off
if ((option == 1)|(option == 2)|(option == 4)){
  trunc_dat = theta_obs[(theta_obs > trunc_min)&(theta_obs < theta_SOTA)]
} else {
  # shrink instead of crop
  trunc_dat = theta_shrink[(theta_shrink > trunc_min)]
}

m = length(theta_obs[(theta_obs > trunc_min)])

# The class imbalance gives a false sense of stability. I'm adjusting n so that the width of the 95% CI 
# corresponds to the AUC CI = 0.0230 (from fig 8, middle panel). Only if we pretend AUC to be theta. 
if (option == 1){
  n_adj = 1481 # for theta_SOTA = max(theta_obs)
} else if ((option == 2)|(option == 3)){
  n_adj = 1750 # for E(theta_SOTA) = max(theta_obs)
} else if ((option == 4)|(option == 5)){
  n_adj= 1972 # for 95 CI = max(theta_obs)
}

mu = floor(n_adj*theta_SOTA)
alpha = 0.05 # 95\% confidence interval
ci_binom = binom.confint(mu,n_adj,conf.level=1-alpha, methods = "exact") # CI for binomial
sprintf("The %s confidence interval for an estimated accuracy of %.4f with n = %s is (%.4f,%.4f), width = %.4f.",  
        (1-alpha)*100, theta_SOTA, n_adj, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] )
n = n_adj # Number of images, n

sprintf("Adjusted test set size: %s. Number of teams with accuracy above %s is %s.",
        n, trunc_min, m)




source("dep_nonid_pmf_fun.R") # for the function 'dep_nonid_pmf' - simulated pmf
source("indep_nonid_pmf_fun.R") # for the function 'indep_nonid_pmf' - simulated pmf
# nonid_cdf, nonid_pmf - analytical pmf


rep = 100000

# Bootstrap a theta-vector of length m from the kaggle observations truncated at theta_trunc
B = 1000

lowerCI = numeric(B) 
upperCI = numeric(B) 

lowerCIindep = numeric(B) 
upperCIindep = numeric(B) 

E_SOTA = numeric(B) # the mean (over the bootstraps) minimum number of failures among all classifiers
max_SOTA = numeric(B) # the max for each bootstraps
teamsSOTA = numeric(B) #number of teams above SOTA

for (b in 1:B){
  theta_vec = sample(x=trunc_dat, size=m, replace=TRUE) # bootstrap from truncated kaggle observations
  
  # minimum number of wrong classifications among all classifiers, vector of length rep
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
  
  # Xindep = indep_nonid_pmf(n, theta_vec, m, rep)
  
  # The bounds of the 95% confidence interval
  sort_min_dep = sort(X$min_dep) # sort the minimum number of failures
  min_dep_alpha2 = sort_min_dep[(alpha/2)*rep] # find the alpha/2 upper bound
  min_dep_alpha2_low = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound
  
  lowerCI[b] = min_dep_alpha2_low
  upperCI[b] = min_dep_alpha2
  
  # The bounds of the 95% confidence interval, indep
  # sort_min_indep = sort(Xindep$min_nonid) # sort the minimum number of failures
  # min_indep_alpha2 = sort_min_indep[(alpha/2)*rep] # find the alpha/2 upper bound
  # min_indep_alpha2_low = sort_min_indep[(1-alpha/2)*rep] # find the alpha/2 lower bound
  # lowerCIindep[b] = min_indep_alpha2_low
  # upperCIindep[b] = min_indep_alpha2
  
  E_SOTA[b] =  (n-mean(X$min_dep))/n # The expected value
  max_SOTA[b] = sort_min_dep[1]
  
  print(c(b, sort_min_dep[(1/2)*rep], max_SOTA[b]))
  print((n-upperCI[b])/n)
  
  teamsSOTA[b] = mean(X$teamsSOTA)
  
  print(teamsSOTA[b])
}

# have a look at the histograms
hist(theta_vec, breaks=250, xlim = c(0.84, 0.96),
     main = 'one bootstrap', xlab = 'bootstrapped accuracies', ylab = 'number of teams')

hist(theta_obs[theta_obs > trunc_min], breaks=200, xlim = c(0.84, 0.96),
     main = 'kaggle realisations', xlab = 'accuracies', ylab = 'number of teams')

hat_theta = (n-X$x_dep)/n
hist(hat_theta, breaks=200, ylim = c(0,100), xlim = c(0.84, 0.96),
     main = 'one realisation', xlab = 'simulated accuracies', ylab = 'number of teams')



sprintf("%s teams have accuracies above the true sota, %s.", 
        mean(teamsSOTA), theta_SOTA)
teamSOTA = length(theta_obs[theta_obs > theta_SOTA])
sprintf("%s teams have accuracies above the true sota, %s.", 
        teamSOTA, theta_SOTA)




sprintf("The expected value of max(theta_obs) is %.4f. The standard deviations is %.4f.", 
        mean(E_SOTA), sqrt(var(E_SOTA)))


# Mean and standard deviation of the 95 CI upper and lower bound
mean_lowerCI = mean(lowerCI)
mean_upperCI = mean(upperCI)
V_CI = mean(upperCI*upperCI) - mean_upperCI*mean_upperCI

sprintf("The mean bootstrapped simulated %s confidence interval is (%.5f,%.5f) with %s repetitions and %s bootstraps. The standard deviation of the upper CI is %.7f  %.7f",  
        1-alpha, (n-mean_lowerCI)/n, (n-mean_upperCI)/n, rep, B,  sqrt(V_CI)/n, sqrt(var(upperCI))/n)

hist((n-X$min_dep)/n, breaks = 30, xlim=c(trunc_min,0.99), freq = F, main = 'max accuracies one bootstrap')
hist((n-max_SOTA)/n, xlim=c(trunc_min,0.99), freq = F, main = 'max accuracies all bootstraps')



teams95CI = length(theta_obs[theta_obs > ((n-mean_lowerCI)/n)])
sprintf("%s teams have accuracies above the lower 95 CI.", 
        teams95CI)








