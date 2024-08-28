# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) September 2023 

# approximating the AUCs from kaggle with thetas - probabilities of correct prediction

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

# Data comes from this website:
# https://www.kaggle.com/c/siim-isic-melanoma-classification/leaderboard
# It has been scraped and saved locally 

comb_data <- readRDS('/Users/kajsam/Documents/kaggle-leaderboard-scrape/SIIM-ISIC_Melanoma_kaggle_leadboard_data.RDS')
head(comb_data) # have a look

# Only interested in the private scores, i.e., the independent test set results
dat <- data.frame(AUC = comb_data$prv_score, dataset = "test")

# Have a quick look at the histogram
hist(dat$AUC, breaks=200, main = paste("Histogram of AUCs"), xlab = "AUCs", ylab = "number of teams")

# We don't have access to the test set labels, so we'll have to estimate based on assumptions. 
# Numbers from https://arxiv.org/ftp/arxiv/papers/2008/2008.07360.pdf
malignant_rate = 584/33126 # in training set 
n_val_test = 10982 # size of test and validation set combined
test_prop = 0.7 # 30/70 split

# Approximations and estimates
n_test = round(test_prop*n_val_test) # approximated test set size - this should be fairly accurate
n_mal = round(test_prop*n_val_test*malignant_rate) # estimated number of malignant cases

sprintf("Estimated number of malignant cases: %s. Estimated test set size: %s. Number of teams: %s. Maximum kaggle AUC: %s",
        n_mal, n_test, length(dat$AUC), max(dat$AUC))

#################################################################################
####################### \subsection{Strategies for estimating SOTA} #################
###############################################################################

# With balanced accuracy and binary prediction, we have that accuracy = AUC
theta_obs = dat$AUC 

# The SOTA estimation is done by bootstrapping from a cropped/shrinked version 
# of the kagge observations, and then a simulation is performed with correlation. 
# The parameters for cropping/shrinking are adjusted according to the wanted 
# outcome: either that the expected value of the maximum simulated theta is equal 
# to the maximum kaggle AUC, or that the upper limit of the 95% CI of the 
# maximum simulated theta is equal to the maximum kaggle AUC

rho = 0.6 # correlation coefficient, the number is calculated from Mania (2019)
lambda = 4 # shrinking parameter, I adjusted it so that number of teams above SOTA is same for sim and kaggle

option = 2

if (option == 1){
  theta_SOTA = max(theta_obs) # Demonstrating that max(theta_obs) is a biased estimate for theta_SOTA
  main_title = 'direct bootstrap'
} else if (option == 2){ # crop
  theta_SOTA = 0.9378 # parameter adjusted until E(theta_sota) = max(theta_obs)
  main_title = 'crop expected value'
} else if (option == 3){ # shrink
  maud = 0.9335 # parameter adjusted until E(theta_sota) = max(theta_obs)
  theta_shrink = (theta_obs-maud)/lambda + maud 
  theta_SOTA = max(theta_shrink)
  main_title = 'shrink expected value'
} else if (option == 4){ # crop
  theta_SOTA = 0.9295 # parameter until upper limit of 95% CI = max(theta_obs)
  main_title = 'crop upper CI'
} else if (option == 5){ # shrink
  maud = 0.92415
  theta_shrink = (theta_obs-maud)/lambda + maud 
  theta_SOTA = max(theta_shrink)
  main_title = 'shrink upper CI'
}

# The class imbalance gives a false sense of stability. I'm adjusting n so that the width of the 95% CI 
# corresponds to the AUC CI = 0.0230 (from fig 8, middle panel). Only if we approximate AUC by theta. 
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
print(sprintf("The %s confidence interval for an estimated accuracy of %.4f with n = %s is (%.4f,%.4f), width = %.4f.",  
        (1-alpha)*100, theta_SOTA, n_adj, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] ))
n = n_adj # Number of images, n

print(sprintf("Adjusted test set size: %s.",
        n))


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
hist(trunc_dat, breaks=200, main = paste("Test set performance"), xlab = "AUCs", ylab = "number of teams")



# bootstrap sampling from kaggle data lower than theta_SOTA, and higher that the lower cut-off
if ((option == 1)|(option == 2)|(option == 4)){
  trunc_dat = theta_obs[(theta_obs > trunc_min)&(theta_obs < theta_SOTA)]
} else {
  # shrink instead of crop
  trunc_dat = theta_shrink[(theta_shrink > trunc_min)]
}



source("dep_nonid_pmf_fun.R") # for the function 'dep_nonid_pmf' - simulated pmf

rep = 1000 # 100 000 number of repetitions. high number gives low variation. 
B = 100 #1000 number of bootstraps. high number gives stable error estimation


lowerCI = numeric(B) # lower limit of confidence interval
upperCI = numeric(B) # upper limit of confidence interval

E_SOTA = numeric(B) # the mean (over rep) minimum number of failures among all classifiers
V_SOTA = numeric(B) # the variance of the expected value of minimum number of failures among all classifiers
teamsSOTA = numeric(B) # mean number of teams above SOTA
qq_x = matrix(nrow = B, ncol = m) # keep some for the qq-plot
boot_x = matrix(nrow = B, ncol = rep) # keep some for bootstrap illustration

tic()
for (b in 1:B){
  theta_vec = sample(x=trunc_dat, size=m, replace=TRUE) # bootstrap from truncated kaggle observations
  
  # minimum number of wrong classifications among all classifiers, vectors of length rep
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
  
  # The bounds of the 95% confidence interval
  sort_min_dep = sort(X$min_fail) # sort the minimum number of failures
  lowerCI[b] = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound
  upperCI[b] = sort_min_dep[(alpha/2)*rep] # find the alpha/2 upper bound
  
  E_SOTA[b] =  (n-mean(X$min_fail))/n # The expected value
  V_SOTA[b] = mean(X$min_fail*X$min_fail) - mean(X$min_fail)*mean(X$min_fail) # The variance
  
  teamsSOTA[b] = mean(X$teamsSOTA)
  
  qq_x[b,] = X$x_fail
  boot_x[b,] = X$min_fail  
  
}

hist(V_SOTA)

qq_x = as.vector(qq_x)
boot_x = as.vector(boot_x)
tid = toc()
print(sprintf("%s repetitions, %s bootstraps took %s",
              rep, B, tid))

# have a look at the histograms

par(mfrow = c(2,2))
# a bootstrap
hist(theta_vec, breaks=250, xlim = c(0.84, 0.96), ylim = c(0,100), 
     main = 'one bootstrap', xlab = 'bootstrapped accuracies', ylab = 'number of teams')
# the kaggle 
hist(theta_obs[theta_obs > trunc_min], breaks=250, xlim = c(0.84, 0.96), ylim = c(0,100), 
     main = 'kaggle observations', xlab = 'accuracies', ylab = 'number of teams')
# one realisation 
hist((n-X$x_fail)/n, breaks=250, xlim = c(0.84, 0.96), ylim = c(0,100), 
     main = paste('one realisation'), xlab = 'simulated accuracies', ylab = 'number of teams')
# many realisations
hat_theta = (n-qq_x)/n
hist(hat_theta, breaks=250, xlim = c(0.84, 0.96), #ylim = c(0,100), 
     main = paste(B,'realisations'), xlab = 'simulated accuracies', ylab = 'number of teams')

par(mfrow = c(1,1))
# one realisation 
hist((n-X$x_fail)/n, breaks=250, xlim = c(0.84, 0.96), ylim = c(0,100), col = 'red',
     main = main_title, xlab = 'simulated accuracies', ylab = 'number of teams')
par(new=TRUE)
hist(theta_obs[theta_obs > trunc_min], breaks=250, xlim = c(0.84, 0.96), ylim = c(0,100), 
     main = main_title, xlab = 'accuracies', ylab = 'number of teams')
par(new=FALSE)
hist(theta_obs[theta_obs > trunc_min], breaks=250, xlim = c(0.84, 0.96), ylim = c(0,100), 
     main = main_title, xlab = 'accuracies', ylab = 'number of teams')
par(new=TRUE)
# one realisation 
hist((n-X$x_fail)/n, breaks=250, xlim = c(0.84, 0.96), ylim = c(0,100), col = 'red',
     main = main_title, xlab = 'simulated accuracies', ylab = 'number of teams')

# have a look at the qq-plot

qqplot(theta_obs[theta_obs > 0.88], hat_theta[hat_theta>0.88], cex = 0.5, xlab = 'kaggle observations',
       ylab = paste('option=',option,' size=',length(qq_x)), conf.level = 0.95,
       xaxt = "n", yaxt = "n", main = main_title,
       conf.args = list(exact = NULL, simulate.p.value = TRUE, B = 2000, col = NA, border = NULL))
axis(1,at=c(0.88, 0.90, 0.92, 0.93, 0.94,0.945, 0.95), tck = 1, lty = 2, col = "gray")
axis(2,at=c(0.88, 0.90, 0.92, 0.93, 0.94,0.945, 0.95), tck = 1, lty = 2, col = "gray")
abline(0,1, col = "red", lwd = 2, lty = 2) #lm(sort(hat_theta) ~ sort(kaggle_line))



print(sprintf("The mean number of simulated teams have accuracies above the true sota, %s, is %.1f.", 
        theta_SOTA, mean(teamsSOTA)))
teamSOTA = length(theta_obs[theta_obs > theta_SOTA])
print(sprintf("%s teams have accuracies above the true sota, %s.", 
        teamSOTA, theta_SOTA))

print(sprintf("The expected value of max(theta_obs) is %.4f. The standard deviations is %.6f.", 
        mean(E_SOTA), sqrt(var(E_SOTA))))


# Mean and standard deviation of the 95 CI upper and lower bound
mean_lowerCI = mean(lowerCI)
mean_upperCI = mean(upperCI)
V_CI = mean(upperCI*upperCI) - mean_upperCI*mean_upperCI

print(sprintf("The mean bootstrapped simulated %s confidence interval is (%.5f,%.5f) with %s repetitions and %s bootstraps. The standard deviation of the upper CI is %.7f  %.7f",  
        1-alpha, (n-mean_lowerCI)/n, (n-mean_upperCI)/n, rep, B,  sqrt(V_CI)/n, sqrt(var(upperCI))/n))

# the distribution of SOTA
hist((n-X$min_fail)/n, breaks = 30, xlim=c(trunc_min,0.99), freq = F, main = 'max accuracies one bootstrap')
hist((n-boot_x)/n, breaks = 50, xlim=c(trunc_min,0.99), freq = F, main = paste('max accuracies', B,'bootstraps'))

# checking if this corresponds to mean CI
sort_min_dep = sort((n-boot_x[-1])/n)
lowerCIboot = sort_min_dep[(1-alpha/2)*rep*B] # find the alpha/2 lower bound
upperCIboot = sort_min_dep[(alpha/2)*rep*B] # find the alpha/2 upper bound
# it does


teams95CI = length(theta_obs[theta_obs > ((n-mean_lowerCI)/n)])
print(sprintf("%s teams have accuracies above the lower 95 CI.", 
        teams95CI))

teamSOTA = length(theta_obs[theta_obs > theta_SOTA])
print(sprintf("%s teams have accuracies above the estimated sota, %s.", 
        teamSOTA, theta_SOTA))

