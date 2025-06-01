# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) May 2025 - copied from SOTA_bootstrap_accuracy.R

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
library(latex2exp)  # mathematical notation
library(future)     # parallelisation
library(furrr)

future::plan("multisession", workers=6) # let's engage 6 out of 8 cores

###############################################################################
############################ Estimating $\theta_{SOTA}$ #############################
###############################################################################

source('Parameters_PublicCompetition.R')
# correlation coefficient, the number is calculated from Mania (2019)
# some figure parameters

source("dep_nonid_pmf_fun.R") # for the function 'dep_nonid_pmf' - simulated pmf


# We have two examples from Kaggle; Multi-Class Prediction of Obesity Risk and Cassava Leaf Disease Classification
# The file SOTA_bootstrap_accuracy.R gives estimates for theta_SOTA by cropping. 
# It turned out quite ugly, so here we do shrinking instead, after we decided on 
# an 'intuitive' shrinking point: 1/c, where c is the number of classes. The 
# results are, as expected, more or less the same; for cropping we had 
# theta_SOTA=0.9063, for shrinking we have 0.9070 for Obesity


# Figures for the manuscipt at the end

# Figure settings - these might need to be adjusted manually
kslim = c(0.88, 0.92)
whylim = c(0,120)
kslab = TeX(r'($\theta$)')
whylab = '          m'

n_breks = 175         # adjust for pleasant graphics
step = (kslim[2]-kslim[1])/n_breks
breks = c(kslim[1], seq(kslim[1]+step, kslim[2], step))

casava = F # to distinguish the two data sets

# Data comes from this website:
# https://www.kaggle.com/competitions/cassava-leaf-disease-classification/leaderboard
# https://www.kaggle.com/competitions/playground-series-s4e2/leaderboard
# It has been scraped and saved locally 

# I will exclude teams with \hat{\theta} <= 1/#classes.
# Information about the number of classes is found on the websites/in train data.
if (casava) {
  comb_data <- readRDS('/Users/kajsam/Documents/Casava Webscrape/Casava_kaggle_leadboard_data.RDS')
  c = 5
  maintitle = 'cassava'
} else {
  comb_data <- readRDS('/Users/kajsam/Documents/Obesity Webscrape/Obesity_kaggle_leadboard_data.RDS')
  c = 7
  maintitle = 'obesity'
}

head(comb_data) # have a look

# Only interested in the private scores, i.e., the independent test set results
dat <- data.frame(theta = comb_data$prv_score, dataset = "test")

theta_obs = dat$theta # the sample estimates (observed) probability of success
m_tot = length(theta_obs) # total number of participating teams
theta_obs = theta_obs[theta_obs>1/c] # exclude performance below chance
m = length(theta_obs)
max_theta_hat = max(theta_obs)

# The size of the test set. This information comes from the website, or the test data set. 
if (casava){
  n = 15000 # approximate number, test data set not available
} else {
  n = 13840 # size of test data set
}

# Have a quick look at the histogram
hist(theta_obs, breaks=200, main = maintitle, xlab = kslab, ylab = whylab)

# Calculating the CI for hat{\theta}^max without multiplicity correction
alpha = 0.05
mu = floor(n*max_theta_hat)
ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial

print(sprintf("The %s confidence interval for an estimated accuracy of %.4f with n = %s is (%.4f,%.4f), width = %.4f.",  
              (1-alpha)*100, max_theta_hat, n, ci_binom["lower"], ci_binom["upper"],ci_binom["upper"]-ci_binom["lower"] ))

teams_single95 = length(theta_obs[theta_obs > ci_binom["lower"][1,1]])
print(sprintf("%s teams have accuracies above the lower 95 CI.", 
              teams_single95))


shrink = F # shrink or not. not shrinking is to produce figures

if (shrink){
  # The SOTA estimation is done by shrinking the kaggle observations, and then a 
  # simulation is performed with correlation. The parameters for shrinking are 
  # adjusted according to the wanted outcome: either that the expected value of 
  # the maximum simulated theta is equal to the maximum kaggle theta, or that the 
  # upper limit of the 95% CI of the maximum simulated theta is equal to the 
  # maximum kaggle theta. 
  # search is done with lower B and rep, and intermediate numbers are not accurate

  shrink_point = 1/c
  # for expected value:
  weight <- 0.9941 # parameter adjusted until E(theta_sota) = max_theta_hat. 
  # 0.99  -> 0.90847
  # 0.9925-> 0.91034
  # 0.9935-> 0.91109
  # 0.994 -> 0.91146
  # 0.9942-> 0.91161
  # 0.9945-> 0.91183
  # 0.995 -> 0.91220 

  # for upper CI:
  # weight <- 0.9896 # parameter adjusted until upper CI = max_theta_hat. 
  # 0.99  -> 0.91179
  # 0.9896-> 0.91150 
  # 0.9895-> 0.91142
  # 0.9875-> 0.90996
  # 0.9885-> 0.91069
  # 0.985 -> 0.90811

  # max_theta_hat = 0.91157 # the aim

  theta_shrunk <- weight*theta_obs + (1-weight)*shrink_point

  # Have a quick look at the shrunk data
  hist(theta_shrunk, breaks=200, main = maintitle, xlab = kslab, ylab = whylab)

  theta_SOTA = max(theta_shrunk)
}

# Simulate dependency
theta_0 = theta_SOTA # theta_0 is the probability of correct prediction for the leading classifier. adjust to E(theta_SOTA)? 

# The theta_obs must be truncated from below because of the a dependency 
# the minimum theta_j is 
trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))

# the two solutions for the quadratic equation (does not influence the lower cut-off)
# a = -rho^2*theta_0*(1-theta_0)-theta_0^2
# b = rho^2*theta_0*(1-theta_0)+2*theta_0
# c = -theta_0^2
# x1 = (-b+sqrt(b^2-4*a*c))/(2*a)
# x2 = (-b-sqrt(b^2-4*a*c))/(2*a)

print(sprintf("With an estimated SOTA of %.4f and a correlation coefficient of %s, the minimum value for theta_j is %.4f.",
              theta_0, rho, trunc_min))

if (shrink){
  shrunk_dat = theta_shrunk[theta_shrunk>trunc_min]
  m = length(shrunk_dat) # number of teams left
  print(sprintf("The number of teams above the lower threshold is %s.",
              m)) 

  # Have a quick look at the shrunk histograms
  hist(shrunk_dat, breaks=n_breks, main = c(maintitle, 'shrunk > min_trunc'), xlab = kslab, ylab = whylab)
}

rep = 1000# 0 # 10 000 number of repetitions. high number gives low variation. 
B = 50 #1000 # number of parameter samples. high number gives stable error estimation

if (B < 50){
  samp_x = matrix(nrow = B, ncol = rep) # keep for illustration
}

tic()
if (shrink){
  EVsotab <- furrr::future_map(1:B, function (x) { 
    theta_vec = sample(x=shrunk_dat, size=m, replace=TRUE) # sample from shrunk kaggle observations
  
    # minimum number of wrong classifications among all classifiers, vectors of length rep
    X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
  
    # The bounds of the 95% confidence interval
    sort_min_dep = sort(X$min_fail) # sort the minimum number of failures
    lowerCI = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound
    upperCI = sort_min_dep[(alpha/2)*rep] # find the alpha/2 upper bound
  
    E_SOTA = mean(X$min_fail) # The expected value
    V_SOTA = mean(X$min_fail*X$min_fail) - mean(X$min_fail)*mean(X$min_fail) # The variance
  
    teamsSOTA = mean(X$teamsSOTA)
  
    return(cbind(lowerCI, upperCI, E_SOTA, V_SOTA, teamsSOTA))
  
  }, .progress=T, .options= furrr::furrr_options(seed=T))
}else{ # not shrinking
  theta_SOTA = max_theta_hat
  theta_0 = theta_SOTA
  trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))
  obs_dat = theta_obs[theta_obs>trunc_min]
  m = length(obs_dat) # number of teams left
  
  EVsotab <- furrr::future_map(1:B, function (x) { 
    theta_vec = sample(x=obs_dat, size=m, replace=TRUE) # sample from shrunk kaggle observations
    
    # minimum number of wrong classifications among all classifiers, vectors of length rep
    X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
    
    # The bounds of the 95% confidence interval
    sort_min_dep = sort(X$min_fail) # sort the minimum number of failures
    lowerCI = sort_min_dep[(1-alpha/2)*rep] # find the alpha/2 lower bound
    upperCI = sort_min_dep[(alpha/2)*rep] # find the alpha/2 upper bound
    
    E_SOTA = mean(X$min_fail) # The expected value
    V_SOTA = mean(X$min_fail*X$min_fail) - mean(X$min_fail)*mean(X$min_fail) # The variance
    
    teamsSOTA = mean(X$teamsSOTA)
    
    return(cbind(lowerCI, upperCI, E_SOTA, V_SOTA, teamsSOTA))
    
  }, .progress=T, .options= furrr::furrr_options(seed=T))
}

tid = toc()
print(sprintf("%s repetitions, %s param samples took %s",
              rep, B, tid))

# Not very elegant with the for-loop, but that's ok for now
lowerCI = numeric(B) # lower limit of confidence interval
upperCI = numeric(B) # upper limit of confidence interval
E_SOTA = numeric(B) # the mean (over rep) minimum number of failures among all classifiers
V_SOTA = numeric(B) # the variance of the expected value of minimum number of failures among all classifiers
teamsSOTA = numeric(B) # mean number of teams above SOTA

for (b in 1:B){ 
  lowerCI[b] = EVsotab[[b]][1]
  upperCI[b] = EVsotab[[b]][2]
  E_SOTA[b] = EVsotab[[b]][3]
  V_SOTA[b] = EVsotab[[b]][4]
  teamsSOTA[b] = EVsotab[[b]][5]
}

if (B>1){
  hist(sqrt(V_SOTA)/n)
}

# Expected value and the standard deviation
print(sprintf("The expected value of max_theta_hat is %.5f. The standard deviations is %.5f. From the expected variance: %.5f", 
              (n-mean(E_SOTA))/n, sqrt(mean(V_SOTA)+var(E_SOTA))/n, sqrt(mean(V_SOTA))/n))

# print(sprintf("The mean number of simulated teams have accuracies above the true sota, %s, is %.1f.", theta_SOTA, mean(teamsSOTA)))
teamSOTA = length(theta_obs[theta_obs > theta_SOTA])
print(sprintf("%s teams have accuracies above the true sota, %s.", 
              teamSOTA, theta_SOTA))

# Mean and standard deviation of the 95 CI upper and lower bound
mean_lowerCI = mean(lowerCI)
mean_upperCI = mean(upperCI)
V_CI = mean(upperCI*upperCI) - mean_upperCI*mean_upperCI

print(sprintf("The mean simulated %s confidence interval is (%.5f,%.5f) with %s repetitions and %s param samples. The standard deviation of the upper CI is %.7f  %.7f",  
              1-alpha, (n-mean_lowerCI)/n, (n-mean_upperCI)/n, rep, B,  sqrt(V_CI)/n, sqrt(var(upperCI))/n))
teams95CI = length(theta_obs[theta_obs > ((n-mean_lowerCI)/n)])
print(sprintf("%s teams have accuracies above the lower 95 CI.", 
              teams95CI))

# the distribution of SOTA, just to have a look
if (shrink){
  theta_vec = sample(x=shrunk_dat, size=m, replace=TRUE) # sample from shrunk kaggle observations
  # minimum number of wrong classifications among all classifiers, vectors of length rep
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
  hist((n-X$min_fail)/n, breaks = 30, xlim=c((n-mean_lowerCI)/n-0.01, (n-mean_upperCI)/n+0.01), 
       freq = F, main = 'max accuracies one param sample')
}else{
  theta_vec = sample(x=obs_dat, size=m, replace=TRUE) # sample from the kaggle observations
  # minimum number of wrong classifications among all classifiers, vectors of length rep
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
  hist((n-X$min_fail)/n, breaks = 30, xlim=c((n-mean_lowerCI)/n-0.01, (n-mean_upperCI)/n+0.01), 
       freq = F, main = 'max accuracies one param sample')
}



################################################################################
####################### Figures ################################################
################################################################################



if (shrink == F){
  
  ##################### obesity_kaggle / casava_kaggle ##################################
  
  hist(theta_obs[theta_obs>kslim[1]], breaks=breks, freq = T,
       main = '', xlab = '',  ylab = '', xlim = kslim, ylim = whylim, 
       col = "lightgray", border = "lightgray",axes=F)
  
  clr = 'magenta'
  par(new = T) # plot the confidence interval
  plot(c(ci_binom["lower"][1,1], ci_binom["upper"][1,1]), c( -1,-1), "l", lwd = 2, 
       col = clr, xlim = kslim, ylim = whylim, ylab = '', xlab = '', axes=F)
  par(new=TRUE) 
  plot(c(ci_binom[["lower"]], ci_binom[["lower"]]), c(-4,2),"l", lwd = 2, 
       col=clr, xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)
  par(new=TRUE) 
  plot(c(ci_binom[["upper"]], ci_binom[["upper"]]), c(-4,2),"l", lwd = 2,
       col=clr, xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)
  
  par(new=TRUE) # dottet vertical line for max_theta_hat
  plot(c(max_theta_hat, max_theta_hat), c(0,max(whylim)),"l", lty = 5, col='magenta', xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  
  # axis, ticks and labels
  axis(1 , cex.axis=1.2, las = 2, at=c(kslim[1],max_theta_hat,kslim[2]), 
       labels=c(kslim[1], TeX(r'($\hat{\theta}_{max}$)'), kslim[2]))
  axis(2 , cex.axis=1.2, las = 2)
  
  title(main = "", xlab = kslab, ylab = whylab, line = 2, cex.lab=1.2)
  legend(0.88, 100, legend=c(TeX(r'($\hat{\Theta}$)')), col=c("lightgray"), pch=c(15), cex=0.8, bty="n")
  
  ####################### end figure ##########################################
  
  ############ obesity_direct_sampling / casava_direct_sampling ###############
  
  # Sampling from the rho-truncated data, once to show a histogram
  trunc_dat = theta_obs[theta_obs>trunc_min]
  theta_vec = sample(x=trunc_dat, size=m, replace=TRUE) # sample from truncated kaggle observations
  
  # minimum number of wrong classifications among all classifiers, vectors of length rep
  tic()
  rep = 10000 
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_SOTA) # one realisation
  toc() # 110 sec for rep = 100,000
  
  theta_real = (n-X$x_fail)/n
  hist(theta_real[theta_real>kslim[1]], breaks=breks, xlim = kslim, ylim = whylim, freq = T,
       main = '', xlab = '', ylab = '', col="gray20", border = "gray20",axes=F)
  
  par(new = T) #
  hist(theta_obs[theta_obs>kslim[1]], breaks=breks, freq = T,
       main = '', xlab = '',  ylab = '', xlim = kslim, ylim = whylim, 
       col = NULL, border = "lightgray",axes=F)
  
  par(new = T) # plot the confidence interval
  plot(c((n-mean_lowerCI)/n, (n-mean_upperCI)/n), c( -1,-1), "l", lwd = 2, 
       col = "blue", xlim = kslim, ylim = whylim, ylab = '', xlab = '', axes=F)
  par(new=TRUE) 
  plot(c((n-mean_lowerCI)/n, (n-mean_lowerCI)/n), c(-4,2),"l", lwd = 2, 
       col="blue", xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)
  par(new=TRUE) 
  plot(c((n-mean_upperCI)/n, (n-mean_upperCI)/n), c(-4,2),"l", lwd = 2,
       col="blue", xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)
  
  par(new=TRUE) # dottet vertical line for expected value
  plot(c((n-mean(E_SOTA))/n, (n-mean(E_SOTA))/n), c(0,max(whylim)),"l", lty = 5, col="blue", xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  
  par(new=TRUE) # dottet vertical line for max_theta_hat
  plot(c(max_theta_hat, max_theta_hat), c(0,max(whylim)),"l", lty = 5, col='magenta', xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  
  # axis, ticks and labels
  axis(1 , cex.axis=1.2, las = 2, at=c(kslim[1],(n-mean(E_SOTA))/n,kslim[2]), 
       labels=c(kslim[1], TeX(r'($E \max \hat{\Theta}'$)'), kslim[2]))
  axis(2 , cex.axis=1.2, las = 2)
  
  title(main = "", xlab = kslab, ylab = whylab, line = 2, cex.lab=1.2)
  legend(0.88, 100, legend=c(TeX(r'($\hat{\Theta}$)'), TeX(r'($\hat{\theta}|{\Theta}' = \hat{\Theta}$)')), col=c("lightgray", "gray20"), pch=c(0,15), cex=0.8, bty="n")
  # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # #
}

if (shrink==T){
  
  ##################### obesity_shrunk_for_expect ##################################
  hist(theta_obs[theta_obs>kslim[1]], breaks=breks, freq = T,
       main = '', xlab = '',  ylab = '', xlim = kslim, ylim = whylim, 
       col = NULL, border = "lightgray", axes=F)
  
  
  #hist(theta_obs[(theta_obs>kslim[1])&(theta_obs<=theta_SOTA)], breaks=breks[breks<=theta_SOTA], freq = T,
  #     main = "", xlab = "", ylab = "", xlim = kslim, ylim = whylim, col = "gray20", border = "gray20")
  par(new = T)
  
  hist(theta_shrunk[theta_shrunk>kslim[1]], breaks=breks, freq = T,
       main = '', xlab = "", ylab = "", xlim = kslim, ylim = whylim, 
       col = "gray20", border = "gray20", axes=F)
  
  par(new=TRUE) # dottet vertical line for theta_SOTA
  plot(c(theta_SOTA, theta_SOTA), c(0,max(whylim)),"l", lty = 5, col="green", xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  par(new=TRUE) # dottet vertical line for max_theta_hat
  plot(c(max_theta_hat, max_theta_hat), c(0,max(whylim)),"l", lty = 5, col='magenta', xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  
  # axis, ticks and labels
  axis(1 , cex.axis=1.2, las = 2, at=c(kslim[1],theta_SOTA, max_theta_hat,kslim[2]), 
       labels=c(kslim[1], TeX(r'(${\theta_{SOTA}}$)'), TeX(r'($\hat{\theta}_{max}$)'), kslim[2]))
  axis(2 , cex.axis=1.2, las = 2)
  
  title(main = "", xlab = kslab, ylab = whylab, line = 2, cex.lab=1.2)
  legend(0.88, 100, legend=c(TeX(r'($\hat{\Theta}$)'), TeX(r'(${\Theta}'$)')), col=c("lightgray", "gray20"), pch=c(0,15), cex=0.8, bty="n")
  
  ##################### end figure ########################################
  
  
  ######################## obesity_shrunk_for_expect_realisation #######################
  
  hist(theta_obs[theta_obs>kslim[1]], breaks=breks, freq = T,
       main = '', xlab = '',  ylab = '', xlim = kslim, ylim = whylim, 
       col = NULL, border = "lightgray", axes=F)
  
  # Sampling from the shrunk data, once to show a histogram
  
  theta_vec = sample(x=shrunk_dat, size=m, replace=TRUE) # sample from shrunk kaggle observations
  
  # minimum number of wrong classifications among all classifiers, vectors of length rep
  tic()
  rep = 10000
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_0) 
  toc()
  
  # one realisation 
  par(new = T)
  theta_real = (n-X$x_fail)/n
  hist(theta_real[theta_real>kslim[1]], breaks=breks, xlim = kslim, ylim = whylim, freq = T,
       main = '', xlab = '', ylab = '', col="gray20", border = "gray20", axes=F)
  
  par(new = T) # plot the confidence interval
  plot(c((n-mean_lowerCI)/n, (n-mean_upperCI)/n), c( -1,-1), "l", lwd = 2, 
       col = "red", xlim = kslim, ylim = whylim, ylab = '', xlab = '', axes=F)
  par(new=TRUE) 
  plot(c((n-mean_lowerCI)/n, (n-mean_lowerCI)/n), c(-4,2),"l", lwd = 2, 
       col="red", xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)
  par(new=TRUE) 
  plot(c((n-mean_upperCI)/n, (n-mean_upperCI)/n), c(-4,2),"l", lwd = 2,
       col="red", xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)
  
  par(new=TRUE) # dottet vertical line for expected value
  plot(c((n-mean(E_SOTA))/n,(n-mean(E_SOTA))/n), c(0,max(whylim)),"l", lty = 5, col="red", xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  
  par(new=TRUE) # dottet vertical line for theta_SOTA
  plot(c(theta_SOTA, theta_SOTA), c(0,max(whylim)),"l", lty = 5, col="green", xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  
  par(new=TRUE) # dottet vertical line for max_theta_hat
  plot(c(max_theta_hat, max_theta_hat), c(0,max(whylim)),"l", lty = 5, col='magenta', xlab = '', ylab = '', 
       xlim = kslim, ylim = whylim, axes=F)
  
  # axis, ticks and labels
  axis(1 , cex.axis=1.2, las = 2, at=c(kslim[1],theta_SOTA, max_theta_hat,kslim[2]), 
       labels=c(kslim[1], TeX(r'(${\theta_{SOTA}}$)'), TeX(r'($E \max \hat{\Theta}'$)'), kslim[2]))
  axis(2 , cex.axis=1.2, las = 2)
  
  title(main = "", xlab = kslab, ylab = whylab, line = 2, cex.lab=1.2)
  legend(0.88, 100, legend=c(TeX(r'($\hat{\Theta}$)'), TeX(r'($\hat{\theta}|\Theta'$)')), col=c("lightgray", "gray20"), pch=c(0,15), cex=0.8, bty="n")
  
}


