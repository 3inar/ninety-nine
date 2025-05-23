# All parameters used in the manuscript, except for the coin-flip example

n = 3000 # size of test set
theta_SOTA = 0.9 # probability of correct prediction
m = 1000 # number of classifiers
alpha = 0.05 # significance level

mu = n*theta_SOTA # the expected number of correct predictions for a single classifier

# the non-identical

d = 0.025

# theta_min = 0.875 # for the non-identical
# theta_max = theta 
#step = (theta_max-theta_min)/(m-1)
# theta_vec = seq(theta_min, theta_max, step)

malignant_rate = 584/33126 # rate of positives in melanoma data training set 
# ( from https://arxiv.org/ftp/arxiv/papers/2008/2008.07360.pdf)

rho = 0.6 # correlation coefficient, the number is calculated from Mania (2019)

rep = 100000 #  number of repetitions, 100,000 gives nice and smooth figures, 1 million is doable

B = 1000 # the outer loop - parameter sampling

theta_0 =  theta_SOTA

# Figure parameters

ylm_bias = c(0.0,0.035)
ylm_sd = c(0.0,0.005)
z_range = 200:300

ylab_bias = TeX(r'($E \hat{theta}_{\max} - {theta}_{SOTA}$)')
ylab_sd = TeX(r'($\sigma_{\hat{\theta}_{\max}}$)')