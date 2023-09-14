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

# all teams
dat <- data.frame(AUC = comb_data$prv_score, dataset = "test")

# violin plot of validation and test

p = ggplot(dat, aes(x = dataset, y = AUC)) +
  geom_violin(trim = T)

p + scale_x_discrete(limits=c("test"))

# best performing teams
# C2 <- data.frame(AUC = comb_data$prv_score[comb_data$prv_score > 0.8], dataset = "test")

# have a look at the histogram
hist(dat$AUC, breaks=200)

hist(comb_data$prv_score[comb_data$prv_score > 0.9], breaks = 30, xlim=c(0.9,0.95), ylim = c(0,200))

###############################################################################
############################ 4.4 A kaggle challenge example #############################
###############################################################################

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

# Parameters
source("Parameters_PublicCompetition.R") # n, theta, m, alpha, rho, theta_min, theta_max, theta_vec, theta_trunc

# The class imbalance gives a false sense of stability. I'm adjusting n so that the width of the 95% CI 
# corresponds to the AUC CI = 0.0240

n_adj= round(n_test/4+300)
mu = floor(n_adj*theta_trunc)

# 95\% confidence interval
ci_binom = binom.confint(mu,n_adj,conf.level=1-alpha, methods = "exact") # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %.4f with n = %s is (%.4f,%.4f), width = %.4f.",  
        (1-alpha)*100, mean(theta_vec), n_adj, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] )

# Parameters
# update the overwritten (sure there is a better way)
n = n_adj # Number of images, n

# Simulate dependency
theta_0 = theta_trunc # theta_0 is the probability of correct prediction for the leading classifier. 

# with a dependency of rho, the minimum theta_j is 
trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))

a = -rho^2*theta_0*(1-theta_0)-theta_0^2
b = rho^2*theta_0*(1-theta_0)+2*theta_0
c = -theta_0^2

trunc_min1 = (-b+sqrt(b^2-4*a*c))/(2*a)
trunc_min2 = (-b-sqrt(b^2-4*a*c))/(2*a)

trunc_dat = comb_data$prv_score[(comb_data$prv_score < theta_trunc)]
trunc_dat = trunc_dat[trunc_dat > trunc_min]

m = length(comb_data$prv_score[(comb_data$prv_score > trunc_min)])




source("dep_nonid_pmf_fun.R") # for the function 'dep_nonid_pmf' - simulated pmf




# rep = 1000

# Bootstrap a theta-vector of length m from the kaggle observations truncated at theta_trunc
B = 100

theta_vec = sample(x=trunc_dat, size=m, replace=TRUE)

# have a look at the histogram
hist(theta_vec, breaks=200)

lowerCI = numeric(B) 
upperCI = numeric(B) 

lowerCI99 = numeric(B) 
upperCI99 = numeric(B) 

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
  
  print(b)

}
toc

hat_theta = (n-X$x_dep)/n
hist(hat_theta, breaks=200)

hist(comb_data$prv_score[(comb_data$prv_score > trunc_min)], breaks = 200)



sprintf("The mean bootstrapped simulated  %s confidence interval is (%.5f,%.5f) with %s repetitions and %s bootstraps.",  
        1-alpha, (n-mean(lowerCI))/n, (n-mean(upperCI))/n, rep, B)

# The standard deviation
Vsota = mean(upperCI*upperCI) - mean(upperCI)*mean(upperCI)
sprintf("The bootstrapped simulated standard deviation of the upper CI is %.7f, with %s repetitions and %s bootstraps.",  
        sqrt(Vsota)/n, rep, B)

sprintf("The mean bootstrapped simulated 99 confidence interval is (%.5f,%.5f) with %s repetitions and %s bootstraps.",  
        (n-mean(lowerCI99))/n, (n-mean(upperCI99))/n, rep, B)

# The standard deviation
Vsota = mean(upperCI99*upperCI99) - mean(upperCI99)*mean(upperCI99)
sprintf("The bootstrapped simulated standard deviation of the upper 99 CI is %.7f, with %s repetitions and %s bootstraps.",  
        sqrt(Vsota)/n, rep, B)





# The expected value
Esota = mean(X$min_dep)

# The standard deviation
Vsota = mean(X$min_dep*X$min_dep) - Esota*Esota

sprintf("The simulated expected value is %.7f and standard deviation is %.7f, with %s repetitions.",  
        (n-Esota)/n, sqrt(Vsota)/n, rep)

hat_theta = (n-X$x_dep)/n
hist(hat_theta, breaks = 30, xlim=c(0.89,0.95), ylim = c(0,500), main = 'example from one repetition')

hist(comb_data$prv_score[comb_data$prv_score > 0.89], breaks = 30, 
     xlim=c(0.89,0.95), ylim = c(0,500), main = 'kaggle data', xlab = 'pretend it is accurcies')

hist((n-X$min_dep)/n, breaks = 30, xlim=c(0.89,0.95), freq = F, main = 'simulated distribution')




hat_theta_up = (n-min_dep_alpha2)/n
hat_theta_low = (n-min_dep_alpha2_low)/n


sprintf("%s/2 percent of the %s teams: %s",  
        alpha, m, m*alpha/2)

sprintf("Number of teams above the upper bound: %s",  
        length(C2$AUC[C2$AUC>hat_theta_up]))

sprintf("The simulated dependent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-0.01, (n-min_dep_one)/n, rep)
sprintf("%s/2 percent of the %s teams: %s",  
        0.01, m, m*0.01/2)

sprintf("Number of teams above the upper bound: %s",  
        length(C2$AUC[C2$AUC>(n-min_dep_one)/n]))


sprintf("The true SOTA is %s.",  
        theta_max)







