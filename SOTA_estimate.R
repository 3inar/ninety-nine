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

# https://www.kaggle.com/c/siim-isic-melanoma-classification/leaderboard

comb_data <- readRDS('/Users/kajsam/Documents/kaggle-leaderboard-scrape/SIIM-ISIC_Melanoma_kaggle_leadboard_data.RDS')

head(comb_data) # have a look

# all teams
C1 <- data.frame(AUC = comb_data$pub_score, dataset = "validation")
C2 <- data.frame(AUC = comb_data$prv_score, dataset = "test")

dat <- rbind(C1, C2)

p = ggplot(dat, aes(x = dataset, y = AUC)) +
  geom_violin(trim = T)

p + scale_x_discrete(limits=c("validation", "test"))

# best performing teams
C1 <- data.frame(AUC = comb_data$pub_score[comb_data$prv_score > 0.8], dataset = "validation")
C2 <- data.frame(AUC = comb_data$prv_score[comb_data$prv_score > 0.8], dataset = "test")

# have a look
hist(C2$AUC, breaks=200)

hist(comb_data$prv_score[comb_data$prv_score > 0.9], breaks = 30, xlim=c(0.9,0.95), ylim = c(0,50))

###############################################################################
############################ 4.3 Dependent non-identical #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is denoted by $theta_hat$.

alpha = 0.05



n_val_test = 10982 # size of test and validation set
test_prop = 0.7 # 30/70 split
n_test = round(test_prop*n_val_test) # test set size

malignant_rate = 584/33126 # in training set ( from https://arxiv.org/ftp/arxiv/papers/2008/2008.07360.pdf)

n_mal = round(test_prop*n_val_test*malignant_rate) # estimated number of malignant cases
sprintf("Estimated number of malignant cases: %s.",
  n_mal)

n_ben = n_test-n_mal # estimated number of benign cases

sprintf("Estimated test set size: %s.",
        n_test)

theta_ref = 0.905


rate = 2
shape = 2
xgamma = rgamma(m/2, shape, scale = 1/rate)
xgam = max(xgamma)-xgamma
#hist(xgam, breaks = 200)
width = 0.025
xgam = width*xgam/max(xgam)
#hist(xgam, breaks = 200)

xunif = runif(m/2, 0.9, 0.93)

theta_vec = c(theta_ref+xgam,xunif)
hist(theta_vec,breaks = 30, xlim=c(0.9,0.95), ylim = c(0,500))

mu_theta = mean(theta_vec)
theta_max = max(theta_vec)
theta_min = min(theta_vec)

m = length(comb_data$prv_score > theta_min)
sprintf("Number of classifiers above %s, m = %s.",
        theta_min, m)



# The class imbalance gives a false sense of stability. I'm adjusting n so that the width of the 95% CI 
# corresponds to the AUC CI = 0.0240

n_adj= round(n_test/4)
mu = floor(n_adj*mu_theta)

# 95\% confidence interval
alpha = 0.05
ci_binom = binom.confint(mu,n_adj,conf.level=1-alpha, methods = "exact") # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %s with n = %s is (%.4f,%.4f), width = %.4f.",  
        (1-alpha)*100, mu_theta, n_adj, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] )

n = n_adj



rho = 0.6 # correlation coefficient

# Simulate dependency

# Set-up from Boland et al (1989) 'Modelling dependence in simple and indirect majority systems',
# where we have a leading classifier with classifications Y_0, and then the m classifiers with 
# correlation rho = corr(Y_0, Y_j). The m classifiers are independent of each other given Y_0.

# For simplicity, let hat{\theta}_0 = mu_theta
y0 = numeric(n) # vector of zeros of length n
mu0 = round(n*0.93)
y0[1:mu0] = 1 # exactly \mu of them are correct classifications




rep = 1000         # 60 sec for a thousand, 10,000 sec for 100,000

min_dep = numeric(rep) # min number of failures with dependency
min_indep = numeric(rep) # for independent, as a check

tic()
for (ell in 1:rep){
  
  x_dep = numeric(m)  # number of failures for m experiments
  
  for (j in 1:m){ 
    
    # the probabilities of Y_j being the opposite of Y_0
    p_flip1 = 1-theta_vec[j] - rho*(1-theta_vec[j]) # P(Y_j = 0|Y_0 = 1)
    p_flip0 = theta_vec[j] - rho*theta_vec[j] # P(Y_j = 1|Y_0 = 0)
    
    # vectors of 0s and 1s indicating a flip relative to y0
    flip1 = rbinom(mu0,1,p_flip1) # flipping correct predictions
    flip0 = rbinom(n-mu0,1,p_flip0) # flipping incorrect predictions
    flip = c(flip1,flip0)
    
    y = abs(y0-flip) # correct predictions
    
    x_dep[j] = n-sum(y)
    
    
  }
  
  
  min_dep[ell] = min(x_dep)
  x_indep = rbinom(m,n,1-theta_vec)
  
  min_indep[ell] = min(x_indep)
  
}
toc()

hat_theta = (n-x_dep)/n

hist(hat_theta, breaks = 30, xlim=c(0.9,0.95), ylim = c(0,500))


# Histograms of the minimum number of failures for m classifiers, in rep repetitions.

histbreaks = seq(min(c(x_dep,x_indep)), max(c(x_dep,x_indep))+9,10)

hist(x_dep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,250))
hist(x_indep, xlab = 'number of failures', ylab = 'number of classifiers', 
     breaks = histbreaks, ylim = c(0,250))

# The upper bound of the 95% confidence interval
sort_min_dep = sort(min_dep) # sort the minimum number of failures
min_dep_alpha2 = sort_min_dep[(alpha/2)*rep] # find the alpha/2 upper bound
min_dep_one = sort_min_dep[(0.01/2)*rep] # find the alpha/2 upper bound

min_dep_alpha2_low = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound

sprintf("The simulated dependent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)
sprintf("The simulated dependent lower bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2_low)/n, rep)

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






sort_min_indep = sort(min_indep)
min_indep_alpha2 = sort_min_indep[(alpha/2)*rep]

#sprintf("The simulated independent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
#        1-alpha, (n-min_indep_alpha2)/n, rep)



sprintf("The true SOTA is %s.",  
        theta_max)







