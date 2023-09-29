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

sprintf("Estimated test set size: %s.",
        n_test)

n_mal = round(test_prop*n_val_test*malignant_rate) # estimated number of malignant cases
# based on equal proportions in validation and test set. this is a strong assumption, but the most reasonable one

sprintf("Estimated number of malignant cases: %s.",
        n_mal)

n_ben = n_test-n_mal # estimated number of benign cases

m = length(theta_obs)

########################### Dependent, non-identical #########################

rho = 0.6

theta_SOTA = 0.9255 # This is the parameter that needs to be adjusted until, say upper limit of 99% CI = max(theta_obs)

# Simulate dependency

theta_0 = theta_SOTA # theta_0 is the probability of correct prediction for the leading classifier. 

# with a dependency of rho, the minimum theta_j is 
trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))

# the two solutions for the quadratic equation (does not influence the lower cut-off)
a = -rho^2*theta_0*(1-theta_0)-theta_0^2
b = rho^2*theta_0*(1-theta_0)+2*theta_0
c = -theta_0^2

trunc_min1 = (-b+sqrt(b^2-4*a*c))/(2*a)
trunc_min2 = (-b-sqrt(b^2-4*a*c))/(2*a)

# bootstrap sampling from kaggle data lower than theta_SOTA, and higher that the lower cut-off
trunc_dat = theta_obs[(theta_obs < theta_SOTA)]
trunc_dat = trunc_dat[trunc_dat > trunc_min]

m = length(theta_obs[(theta_obs > trunc_min)])

# The class imbalance gives a false sense of stability. I'm adjusting n so that the width of the 95% CI 
# corresponds to the AUC CI = 0.0240. Only if we pretend AUC to be theta. 
n_adj= round(n_test/4+5)
mu = floor(n_adj*theta_trunc)
# 95\% confidence interval
alpha = 0.05
ci_binom = binom.confint(mu,n_adj,conf.level=1-alpha, methods = "exact") # CI for binomial
sprintf("The %s confidence interval for an estimated accuracy of %.4f with n = %s is (%.4f,%.4f), width = %.4f.",  
        (1-alpha)*100, mean(trunc_dat), n_adj, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] )
n = n_adj # Number of images, n




source("dep_nonid_pmf_fun.R") # for the function 'dep_nonid_pmf' - simulated pmf

rep = 10000

# Bootstrap a theta-vector of length m from the kaggle observations truncated at theta_trunc
B = 10

theta_vec = sample(x=trunc_dat, size=m, replace=TRUE)

# have a look at the histogram
hist(theta_vec, breaks=200)
hist(theta_obs[theta_obs > trunc_min], breaks=200)

lowerCI = numeric(B) 
upperCI = numeric(B) 

lowerCI99 = numeric(B) 
upperCI99 = numeric(B) 
mdn = numeric(B)

tic
for (b in 1:B){
  theta_vec = sample(x=trunc_dat, size=m, replace=TRUE)
  
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0)
  
  # The bounds of the 95% confidence interval
  sort_min_dep = sort(X$min_dep) # sort the minimum number of failures
  min_dep_alpha2 = sort_min_dep[(alpha/2)*rep] # find the alpha/2 lower bound
  min_dep_alpha2_low = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound
  
  lowerCI[b] = min_dep_alpha2_low
  upperCI[b] = min_dep_alpha2
  
  # The bounds of the 99% confidence interval
  sort_min_dep = sort(X$min_dep) # sort the minimum number of failures
  min_dep_alpha2 = sort_min_dep[(0.01/2)*rep] # find the alpha/2 lower bound
  min_dep_alpha2_low = sort_min_dep[(1-0.01/2)*rep] # find the alpha/2 lower bound
  
  lowerCI99[b] = min_dep_alpha2_low
  upperCI99[b] = min_dep_alpha2
  mdn[b] =  sort_min_dep[(1/2)*rep]
  
  print(c(b,mdn[b]))
  print((n-upperCI99[b])/n)
  
}
toc

mdn

hat_theta = (n-X$x_dep)/n
hist(hat_theta, breaks=200, ylim = c(0,100), 
     main = 'one bootstrap repetition', xlab = 'simulated accuracies', ylab = 'number of teams')

#hist(comb_data$prv_score[(comb_data$prv_score > trunc_min)], breaks = 200, ylim = c(0,100), 
#     main = 'kaggle data', xlab = 'observed accuracies', ylab = 'number of teams')

hist(theta_obs[(theta_obs > trunc_min)], breaks = 200, ylim = c(0,100), 
     main = 'kaggle data', xlab = 'observed accuracies', ylab = 'number of teams')


# Mean and standard deviation of the 95 CI upper and lower bound
mean_lowerCI = mean(lowerCI)
mean_upperCI = mean(upperCI)
Vsota = mean(upperCI*upperCI) - mean_upperCI*mean_upperCI

sprintf("The mean bootstrapped simulated %s confidence interval is (%.5f,%.5f) with %s repetitions and %s bootstraps. The standard deviation of the upper CI is %.7f",  
        1-alpha, (n-mean_lowerCI)/n, (n-mean_upperCI)/n, rep, B,  sqrt(Vsota)/n)

# Mean and standard deviation of the 99 CI upper and lower bound
mean_lowerCI99 = mean(lowerCI99)
mean_upperCI99 = mean(upperCI99)
Vsota99 = mean(upperCI99*upperCI99) - mean_upperCI99*mean_upperCI99

sprintf("The mean bootstrapped simulated 99 confidence interval is (%.5f,%.5f) with %.0f repetitions and %s bootstraps. The standard deviation of the upper CI is %.7f",   
        (n-mean_lowerCI99)/n, (n-mean_upperCI99)/n, rep, B,  sqrt(Vsota99)/n)


hat_theta = (n-X$x_dep)/n
hist(hat_theta[hat_theta > trunc_min], xlim=c(trunc_min,0.98), breaks = 30, ylim = c(0,200), main = 'example from one repetition')

hist(theta_obs[theta_obs > trunc_min], breaks = 50, 
     xlim=c(trunc_min,0.98), ylim = c(0,200), main = 'kaggle data', xlab = 'pretend it is accurcies')

hist((n-X$min_dep)/n, breaks = 30, xlim=c(trunc_min,0.98), freq = F, main = 'distribution of max accuracies')

sprintf("There are %s teams with accuracies above %.4f. With true SOTA of %.4f, the upper 99 CI is %.4f, whereas max kaggle accuracy is %.4f.",  
        m, trunc_min, theta_trunc, (n-mean_upperCI99)/n, max(theta_obs))


teams99 = length(theta_obs[theta_obs > (n-mean_upperCI99)/n])
sprintf("%s out of %s teams have accuracies above 99 CI, %.1f percent.",  
        teams99, m, 100*teams99/m)
teams95 = length(theta_obs[theta_obs > (n-mean_upperCI)/n])
sprintf("%s out of %s teams have accuracies above 95 CI, %.1f percent.", 
        teams95, m, 100*teams95/m)

teamsSOTA = length(theta_obs[theta_obs > theta_trunc])
sprintf("%s out of %s teams have accuracies above the true sota, %.1f percent.", 
        teamsSOTA, m, 100*teamsSOTA/m)









########################### Independent, non-identical #########################

# A quick estimate for AUC_SOTA, see AUROC.R -> 0.9216

theta_trunc = 0.9405 # This is the parameter that needs to be adjusted until, say upper limit of 99% CI = max(theta_obs)

teamsSOTA = length(theta_obs[theta_obs > theta_trunc])
sprintf("%s out of %s teams have accuracies above the true sota, %.1f percent.", 
        teamsSOTA, m, 100*teamsSOTA/m)

source("indep_nonid_pmf_fun.R") # for the function 'indep_nonid_pmf' - simulated pmf
                                # nonid_cdf, nonid_pmf - analytical pmf

# Bootstrap a theta-vector of length m from the kaggle observations truncated at theta_trunc
B = 100000

trunc_dat = theta_obs[(theta_obs < theta_trunc)]

theta_vec = sample(x=trunc_dat, size=m, replace=TRUE)

# have a look at the histogram
hist(theta_vec, breaks=200)
hist(theta_obs, breaks=200)

upperCI95 = numeric(B)
upperCI99 = numeric(B) 

rep = 1
for (b in 1:B){
  theta_vec = sample(x=trunc_dat, size=m, replace=TRUE)
  
  X = indep_nonid_pmf(n, theta_vec, m, rep)
  
  
  upperCI99[b] = X$min_nonid
  
  n99 = length(theta_obs[theta_obs > (n-upperCI99[b])/n])
  
  #print(c(b, (n-upperCI99[b])/n, n99))
  
}

length(upperCI99[upperCI99>(n-max(theta_obs)*n)])







for (b in 1:B){
  theta_vec = sample(x=trunc_dat, size=m, replace=TRUE)
  theta_vec = sort(theta_vec)
  
  # Fz = nonid_cdf(n, theta_vec, m)
  fz = nonid_pmf(n, theta_vec, m)
  
  # Expected value
  Eterm = numeric(n+1)
  for (z in 0:n){
    i = z+1
    Eterm[i] = z*fz[i]
  }
  
  Esota = sum(Eterm)
  Esota_theta = 1-Esota/n
  
  # Variance
  vterm = numeric(n+1)
  for (z in 0:n){
    i = z+1
    vterm[i] = z^2*fz[i]
  }
  esquare = sum(vterm)
  
  Vsota = esquare - Esota^2
  
  sprintf("The expected number of failures is %.4f, with a variance of %.4f.",
          Esota, Vsota)
  sprintf("The expected theta_hat_SOTA is %.6f, with standard deviation of %.6f.",
          (n-Esota)/n, sqrt(Vsota)/n)
  
  print(c(b, (n-Esota)/n, sqrt(Vsota)/n, (n-which(Fz>0.01/2)[1]-1)/n))

  
}







hist(theta_obs, breaks = 200, ylim = c(0,300), 
     main = 'kaggle data', xlab = 'observed accuracies', ylab = 'number of teams')


# Mean and standard deviation of the 95 CI upper bound
mean_upperCI95 = mean(upperCI95)
Vsota = mean(upperCI95*upperCI95) - mean_upperCI95*mean_upperCI95

sprintf("The mean bootstrapped upper bound of the 95 confidence interval is %.5f with %s bootstraps. The standard deviation of the upper CI is %.7f",  
        (n-mean_upperCI95)/n, B,  sqrt(Vsota)/n)

# Mean and standard deviation of the 99 CI upper bound
mean_upperCI99 = mean(upperCI99)
Vsota99 = mean(upperCI99*upperCI99) - mean_upperCI99*mean_upperCI99

sprintf("The mean bootstrapped upper bound of the 99 confidence interval is %.5f with %s bootstraps. The standard deviation of the upper CI is %.7f",   
        (n-mean_upperCI99)/n, B,  sqrt(Vsota99)/n)

hist(theta_obs, breaks = 50, 
     xlim=c(0.3,0.98), ylim = c(0,200), main = 'kaggle data', xlab = 'pretend it is accurcies')

sprintf("There are %s teams. With true SOTA of %.4f, the upper 99 CI is %.4f, whereas max kaggle accuracy is %.4f.",  
        m, theta_trunc, (n-mean_upperCI99)/n, max(theta_obs))


teams99 = length(theta_obs[theta_obs > (n-mean_upperCI99)/n])
sprintf("%s out of %s teams have accuracies above 99 CI, %.1f percent.",  
        teams99, m, 100*teams99/m)
teams95 = length(theta_obs[theta_obs > (n-mean_upperCI95)/n])
sprintf("%s out of %s teams have accuracies above 95 CI, %.1f percent.", 
        teams95, m, 100*teams95/m)

teamsSOTA = length(theta_obs[theta_obs > theta_trunc])
sprintf("%s out of %s teams have accuracies above the true sota, %.1f percent.", 
        teamsSOTA, m, 100*teamsSOTA/m)











