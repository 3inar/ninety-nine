# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) April 2023 


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
C1 <- data.frame(AUC = comb_data$pub_score, dataset = "validation")
C2 <- data.frame(AUC = comb_data$prv_score, dataset = "test")

# violin plot of validation and test
dat <- rbind(C1, C2)

p = ggplot(dat, aes(x = dataset, y = AUC)) +
  geom_violin(trim = T)

p + scale_x_discrete(limits=c("validation", "test"))

# best performing teams
C1 <- data.frame(AUC = comb_data$pub_score[comb_data$prv_score > 0.8], dataset = "validation")
C2 <- data.frame(AUC = comb_data$prv_score[comb_data$prv_score > 0.8], dataset = "test")

# have a look
hist(C2$AUC, breaks=200)

hist(comb_data$prv_score[comb_data$prv_score > 0.9], breaks = 30, xlim=c(0.9,0.95), ylim = c(0,200))

###############################################################################
############################ 4.4 Estimating $\theta_{SOTA}$ from accuracies #############################
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
alpha = 0.05 # significance level
rho = 0.6 # correlation coefficient

# These parameters are used to simulate a nonidentical and dependent distribution
# that is comparable to the observed distribution. We do not want it to be as similar 
# possible, just comparable

theta_min = 0.9 
theta_max = 0.9223 # this is the true theta_SOTA

# Number of classifiers, m

m = length(comb_data$prv_score > theta_min)
sprintf("Number of classifiers above %s, m = %s.",
        theta_min, m)


# So, from the observed accuracies, it seems like we have a mixture of a Gaussian 
# and a Uniform, with a couple of spikes. Obviously too complicated to serve any 
# illustrative purpose. Instead, I'll create a theta vector with equal step 
# size

# theta_mid = theta_min + 2*(theta_max-theta_min)/3
# step0 = (theta_max-theta_min)/m
# step1 = (theta_max-theta_min)/m
# theta_vec0 = seq(theta_min, theta_mid-step0, step0)
# theta_vec1 = seq(theta_mid, theta_max, step1)
# theta_vec = c(theta_vec0, theta_vec1)

step = (theta_max-theta_min)/m
theta_vec = seq(theta_min, theta_max, step)


# The class imbalance gives a false sense of stability. I'm adjusting n so that the width of the 95% CI 
# corresponds to the AUC CI = 0.0240

n_adj= round(n_test/4+300)
mu = floor(n_adj*mean(theta_vec))

# 95\% confidence interval
ci_binom = binom.confint(mu,n_adj,conf.level=1-alpha, methods = "exact") # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %.4f with n = %s is (%.4f,%.4f), width = %.4f.",  
        (1-alpha)*100, mean(theta_vec), n_adj, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] )

n = n_adj

source("dep_nonid_pmf_fun.R") # for the function 'dep_nonid_pmf' - simulated pmf

# Simulate dependency
theta_0 = theta_max # theta_0 is the probability of correct prediction for the leading classifier. 
# I think it should be theta_SOTA, but I'm not entirely sure
rep = 10
tic()
X = dep_nonid_pmf(n, theta, m, rho, rep, theta_vec, theta_0)
toc()

# The bounds of the 95% confidence interval
sort_min_dep = sort(X$min_dep) # sort the minimum number of failures
min_dep_alpha2 = sort_min_dep[(alpha/2)*rep] # find the alpha/2 lower bound
min_dep_alpha2_low = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound
min_dep_one = sort_min_dep[(0.01/2)*rep] # find the alpha/2 upper bound

sprintf("The simulated  %s confidence interval is (%.5f,%.5f) with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2_low)/n, (n-min_dep_alpha2)/n, rep)

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







