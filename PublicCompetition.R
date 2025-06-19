# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) 

library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing
library(latex2exp)  # mathematical notation

###############################################################################
############################ 3.1 Nomenclature #############################
###############################################################################

# n         - number of trials/size of test set
# theta     - prob of success/prob of correct prediction
# k         - number of successes
# X         - r.v., number of failures; n-k
# theta_hat = (n-X)/n - estimator of theta
# px    - probability of x failures in one experiment
# Px    - probability of at most x failures in one experiment

# m     - number of experiments/number of classifiers
# Cx    - r.v., number of experiments with at most x failures
# Z     - r.v., number of failures on at least one classifier

# theta_SOTA_hat - estimator of theta_SOTA 

# P_    - probability, specified when used

# Fz    - cdf of Z: P(C_z > 0|m,n,theta) prob of at least one classifier having at most 
#         z failures (identical to at least n-z successes)
# fz    - pmf of Z

# alpha - significance level

# SOTA - state-of-the-art



source("ProbDistr_thetaSOTA.R") 
# for functions 'cdf', 'pmf', 'expect' and 'variance'

# figures at the end

###############################################################################
################# 3.4 A simulated public competition example ##################
###############################################################################

# Consider a classification problem with a test set of size $n=3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The classifier's performance on the test set, denoted by $theta_hat$ and referred 
# to as the accuracy is $(n-x_j)/n$, where $x_j$ is the observed number of failures
# for classifier $j$.

source("Parameters_PublicCompetition.R") # n, theta_SOTA, m, alpha, mu = n*theta, rep

# 95% confidence interval for hat{theta} = mu, single classifier

ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %s is (%.4f,%.4f).",  
        (1-alpha)*100, theta_SOTA, ci_binom["lower"], ci_binom["upper"])

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# The probability of exceeding the upper limit of the CI should be close to alpha/2 = 0.025 
k_up = floor(ci_binom[["upper"]]*n) # number of successes exceeding the CI
# flooring the CI bound, so P_up > alpha/2
P_up = pbinom(k_up,n,theta_SOTA, lower.tail = F) # P[X>x]
# # # # # # # # # # # # # # # # # P_up = 0.02619 OK # # # # # # # # # # # # # # # # # # # # # # # 

# If there are $m=1,000$ teams, each with a classifier with $\theta$, what is
# the probability of at least one team achieving at least $ci_binom["upper"]$ 
# correct predictions?
# Following the logic of the coin-flip multiple experiments, we have a binomial 
# distribution with $m$ trials and probability of success (exceeding 95% CI) 
# is by definition $alpha/2$. 

# The probability of at least one team exceeding the upper CI bound:
Cx = 1 # number of teams
Px = alpha/2 # by definition
Palpha2 = 1 - dbinom(Cx-1,m,Px)

sprintf("The probability of at least %s out of %s teams exceeding the upper limit of the CI is %s.",  
        Cx, m, Palpha2)

# If there are $m$ classifiers, what must be the prob, $Px$, of each 
# classifier having at most $x$ failures, for the probability of at least Cx = 1
# classifier having at most $x$ failures to be P(Cx > 0) = alpha/2?
# We can calculate $x$ from $Px$.

Cx = 1 # number of teams

# This is easily calculated:
# P(Cx > 0) = 1 - P(Cx = 0) = alpha/2 
# binom coeff is 1 since cx = 0, Palpha2^0 = 1
# and we are left with (1-Palpha2)^m = 1 - alpha\2

Px = 1 - (1-alpha/2)^(1/m)

sprintf("For the probability of at least %s out of %s classifiers to achieve at most $x$ failures to be alpha/2=%s, the prob for each classifier achieving at most $x$ failures must be %.8f.",  
        Cx, m, alpha/2, Px)

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
P = 1-pbinom(Cx-1,m,Px) # should be 0.025
# # # # # # # # # # # # # # # # # P = 0.025000 OK # # # # # # # # # # # # # # # # # # # # # # # 

# If the probability of wrongly predicting at most x out of n data points is Px, 
# what is x when theta = 0.9?

#"The quantile is defined as the smallest value x such that F(x)≥Px, where F is the distribution function."

x_alpha2 = qbinom(Px, n, 1-theta_SOTA) 
theta_alpha2 = (n-x_alpha2)/n
sprintf("With a probablitiy of alpha/2 = %s, at least %s team will achieve an accuracy of at least %.4f, corresponding to x = %s failures.",
        alpha/2, Cx, theta_alpha2, x_alpha2)

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Fz = cdf(n,theta_SOTA,m) # the cdf
x = which(Fz>alpha/2)[1]-1 # should be same as x_alpha2 = 236
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# Since x_alpha2 is discreet, but Px is not, we'll calculate for x_alpha2 and x_alpha2-1, 
# and then the two results should be on each side of alpha/2

Palpha2_discr = pbinom(x_alpha2-1,n,1-theta_SOTA) # 
P_low = 1-pbinom(Cx-1,m,Palpha2_discr) #

Palpha2_discr = pbinom(x_alpha2,n,1-theta_SOTA) # 
P_up = 1-pbinom(Cx-1,m,Palpha2_discr) #
# # # # # # # # # # # # # # # # # ok 0.02475 0.03243# # # # # # # # # # # # # # # # # # # # # # # 

sprintf("In summary, if there are %s teams, each with %s probability of correct prediction, and a test set of size %s, there is a probability of %s that at least %s team will have at most %s incorrect predictions, corresponding to an estimated accuracy of %.4f.",  
        m, theta_SOTA, n, alpha/2, Cx, x_alpha2, theta_alpha2)

# Expected value
Esota = expect(n,theta_SOTA,m)
Esota_theta = 1-Esota/n
sprintf("The expected value is %.4f",
        Esota_theta)



# When the probability of correct prediction is $0.9$, there is a non-negligible probability that the 
# top-ranked team will have an estimated accuracy of at least theta^m_alpha2 = 0.9213. 

# What is the probability of beating the SOTA for a method with significantly better
# accuracy?
P_beat_sota = numeric(2)
P_beat_sota[1] = pbinom(x_alpha2, n, 1-ci_binom[["upper"]]) # 
P_beat_sota[2] = pbinom(x_alpha2-1, n, 1-ci_binom[["upper"]]) # 
P_beat_sota

sprintf("The probability of achieving an estimated accuracy better than the upper bound of the CI for theta_hat_SOTA, %.4f, for a classifier with significantly better probability of correct prediction, p=%.4f, is just %.4f.",
        theta_alpha2, ci_binom[["upper"]], P_beat_sota[1])

# Beat esota?

beat_esota = pbinom(Esota, n, 1-ci_binom[["upper"]]) 

sprintf("The probability of achieving an accuracy better than E(SOTA)=%.4f, for a classifier with significantly better proability of correct prediction, %.4f, is %.4f.",
        Esota_theta, ci_binom[["upper"]], beat_esota)

######################### Expected values and standard deviations for hat theta_max (X) ####################################

# varying m, n, p
# these params also used in plots below
n_vec = c(3000,1000,10000)
m_vec = c(1000,100,5000)
theta_vec = c(0.9,0.85,0.95)

table = matrix(0,27,5)
t = 1 # counting table rows

for (j in 1:length(m_vec)){
  for (i in 1:length(n_vec)){
    for (k in 1:length(theta_vec)){
      
      # parameter values
      table[t,1:3] = c( m_vec[j], n_vec[i], theta_vec[k])
      
      Esota =  expect(n_vec[i], theta_vec[k], m_vec[j])
      Esota_theta = 1-Esota/n_vec[i]
      table[t,4] = Esota_theta
      Vsota = variance(n_vec[i], theta_vec[k], m_vec[j])
      Std_sota_theta = sqrt(Vsota)/n_vec[i]
      table[t,5] = Std_sota_theta
      # sprintf("The expected theta_sota is %.4f, with a standard deviation of %.6f, for m=%s, n=%s, theta=%s.", Esota_theta, Std_sota_theta, m_vec[j], n_vec[i], theta_vec[k])
      
      t = t+1
    }
  }
}

table

##############################################################################################
################################### Figures ##################################################
##############################################################################################

source("plotting_params.R")

################################### Figure multi_ci ##########################################
# The pmfs of two single $\hat{\theta}(X)$ compared to $\hat{\theta}_{\max}$ statistics#######

Esota = expect(n, theta_SOTA, m) # updating the value

# let k be the number of successes = n-x
k = (mu-60):(mu+90) # this is the plot range, adjust to your liking
kslim = c(k[1],k[length(k)])
whylim = c(0,0.03)

# the pmf of $hat theta (X)$ with \theta = \theta

{
  new_png("multi_ci.png", n_figures=1)

  mars <- par()$mar
  mars[2] <- mars[2] - 1 # shrinks margin a litte, nothing there anyway
  mars[1] <- mars[1] + .25  # need some room for the long labels
  par(mar=mars)

  y = dbinom(k, n, theta_SOTA) # probability of x successes in n trials
  plot(k, y, type='l', col = "blue", xlab = '', ylab = '', xlim = kslim, ylim =
       whylim, axes=F)

  # the confidence interval of $hat theta (x) = \theta$ with bars
  par(new=TRUE) 
  plot(c(ci_binom[["lower"]]*n, ci_binom[["upper"]]*n), c( -0.0005,-0.0005),
       "l", col = "blue",  xlab = '', ylab = '', xlim = kslim, ylim = whylim,
       axes=F)

  par(new=TRUE) 
  plot(c(ci_binom[["lower"]]*n, ci_binom[["lower"]]*n), c(-0.002,0.001),"l",
       col="blue", xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)

  par(new=TRUE) 
  plot(c(ci_binom[["upper"]]*n, ci_binom[["upper"]]*n), c(-0.002,0.001),"l",
       col="blue", xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)

  # dottet vertical line for \theta
  par(new=TRUE)
  plot(c(theta_SOTA*n, theta_SOTA*n), c(0,max(y)),"l", lty = 5,
       col=rgb(0,0,1,0.25), xlab = '', ylab = '', xlim = kslim, ylim = whylim,
       axes=F)

  # the pmf of $hat theta (X)$ with \theta = upper CI
  par(new=TRUE) 
  y = dbinom(k, n, ci_binom[["upper"]]) # significantly better classifier
  plot(k, y, type = 'l', col = "red", xlab = '', ylab = '', axes=F, xlim =
       kslim, ylim = whylim)

  # dottet vertical line for \theta
  par(new=TRUE) # 
  plot(c(ci_binom[["upper"]]*n, ci_binom[["upper"]]*n), c(0,max(y)),"l", lty =
       5, col=rgb(1, 0, 0,0.25), xlab = '', ylab = '', xlim = kslim, ylim =
       whylim, axes=F)

  # area under curve for expected SOTA performance
  polygon(c(n-Esota, k[k>=n-Esota], max(k)), c(0,y[k>=n-Esota], 0),
          col=rgb(0,1,0,0.25)) 

  # area under curve for SOTA performance
  k_alpha2 = n-x_alpha2 
  polygon(c(k_alpha2, k[k>=k_alpha2], max(k)),
                                c(0,y[k>=k_alpha2], 0), col=rgb(0,1,0,0.25)) 

  # axis, ticks and labels
  axis(1, las = 3, at = c(kslim[1], ci_binom[["lower"]] * n,
                        n * theta_SOTA, ci_binom[["upper"]] * n,
                        n - Esota, k_alpha2, kslim[2]),
     labels = FALSE) 

  mtext(c('', TeX('$\\theta_{alpha/2}$'), TeX('$\\theta$'),
        TeX('$\\theta_{1-alpha/2}$'), TeX(r'($E \hat{\theta}_{max}$)'),
        TeX('$\\theta^m_{1-alpha/2}$'), ''),
      side = 1, at = c(kslim[1], ci_binom[["lower"]] * n,
                       n * theta_SOTA, ci_binom[["upper"]] * n,
                       n - Esota, k_alpha2, kslim[2]),
      line = 1, las = 2) # Adjust 'line' to move labels up or down

  dev.off()
}
# # # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################### Figure cumul_fail ##########################################################

{
  new_png("cumul_fail.png", n_figures=2)

  mars <- par()$mar
  mars[3] = mars[3] + 1.75
  par(mar=mars)

  Fz = cdf(n,theta_SOTA,m) # updating the value

  # Plotting the cfd
  #  plot(0:n,Fz, type = 'l', xlab = 'number of failures', 
  #       ylab = 'probability of at least one team')
  # the whole range, not very much information

  # zooming in
  z = z_range # parameter
  plot(z,Fz[z+1], type = 's', xlab ='', 
       ylab = '', axes = F)

  # axis, ticks and labels
  xax = seq(z[1],tail(z,1), 20)
  klab = xax 
  plab = round(1000*(n-xax)/n)/1000
  axis(1, cex.axis=1, las = 1, at=xax, labels = as.character(plab))
  axis(2, cex.axis=1, las = 2)
  axis(3, cex.axis=1, las = 1, at=xax, labels = as.character(klab))

  mtext(TeX(r'($z$)'), side=3, line=1.5) 
  mtext(TeX(r'($\hat{\theta}_{max}$)'), side=1, line=2)

  dev.off()
}
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################### Figure pmf_fail ###################################################

{
  new_png("pmf_fail.png", n_figures=2)

  mars <- par()$mar
  mars[3] = mars[3] + 1.75
  par(mar=mars)

  fz = pmf(n,theta_SOTA,m,f0 = T) # updating the value

  # Plotting the pmf
  plot(0:n,fz, type = 's')

  # zooming in, same values z as for the cdf
  z = z_range # parameter
  plot(z,fz[z], type = 'h', xlab ='', ylab = '', axes = F)

  xax = seq(z[1],tail(z,1), 20)
  klab = xax 
  plab = round(1000*(n-xax)/n)/1000

  # axis, ticks and labels, same values as for the cdf
  axis(1, cex.axis=1, las = 1, at=xax, labels = as.character(plab))
  axis(2, cex.axis=1, las = 2)
  axis(3, cex.axis=1, las = 1, at=xax, labels = as.character(klab))

  mtext(TeX(r'($z$)'), side=3, line=1.5) 
  mtext(TeX(r'($\hat{\theta}_{max}$)'), side=1, line=2)
  dev.off()
}
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #! open 


################################### Figure bias_sd_m_n_theta ###################################################

# There are six subfigures, each with 5 curves and a vertical line

# Three bias subfigures m,n,theta, and corresponding sd figures

col_vec = c("green","red","blue") # param = m,n,theta 
line_vec = c("solid","dotted","longdash") #subparam

# these computations done v. often
Esth <- function(n,theta,m) {
  Esota =  expect(n, theta, m)
  (1-Esota/n) - theta
}

SDsth <- function(n,theta,m) {
  Vsota = variance(n, theta, m)
  sqrt(Vsota)/n
}

# I coudn't figure out how to switch between graphics and save them at the end,
# so I'll just manually do bias and sd for now

########################### Figure bias_m and sd_m ##########################################

# plots the "identical" curve parts of the bias and the SD plot
bsm_curve <- function(bias) {
  if (bias) {
    stat_fn = Esth
    ylm = ylm_bias
  } else {
    stat_fn = SDsth
    ylm = ylm_sd
  }

  # middle black line
  stat_curve = sapply(m_x, \(xx) stat_fn(n, theta_SOTA, xx))

  plot(m_x, stat_curve, "l", lty = "solid", col = "black", ylim=ylm,
       xlab = "", ylab ="", axes=F)
  abline(v=m, col="gray") # intersection corresponding to upper row in table

  # two red lines
  #################### param = n -- means that n is changed from the black line
  param = 2
  for (k in 2:3){
    stat_curve = sapply(m_x, \(xx) stat_fn(n_vec[k], theta_SOTA, xx))
    lines(m_x, stat_curve, lty = line_vec[k], col = col_vec[param])
  }

  # two blue lines
  #################### param = theta
  param = 3
  for (k in 2:3){
    stat_curve = sapply(m_x, \(xx) stat_fn(n, theta_vec[k], xx))
    lines(m_x, stat_curve, lty = line_vec[k], col = col_vec[param])
  }
}

# common for both figures
lgnds = c(TeX(r'($n=1000$)'),TeX(r'(${theta}=0.85$)'), TeX(r'($n = 3000, {theta}=0.90$)'), 
          TeX(r'(${theta}=0.95$)'), TeX(r'($n=10000$)'))
cls = c("red","blue","black","blue", "red")
m_x = seq(1, 5000, by=10) # m is on the x-axis

# bias subfigure
bsm_curve(bias=T)

axis(1, cex.axis=1.2, las = 1, at=c(1, 1000, 2000, 3000, 4000, 5000), # ticks
     labels=c('1','1000','2000', '3000', '4000', '5000'))
axis(2, cex.axis=1.2, las = 1)

title(main = "", xlab = "m", ylab = ylab_bias, line = 2, cex.lab=1.2)
legend(2900, 0.035, legend=lgnds, col=cls, lty=c(3,3,1,5,5), cex=0.8)

# sd subfigure
bsm_curve(bias=F)

axis(1, cex.axis=1.2, las = 1, at=c(1, 1000, 2000, 3000, 4000, 5000), # ticks
     labels=c('1','1000','2000', '3000', '4000', '5000'))
axis(2, cex.axis=1.2, las = 1)

title(main = "", xlab = "m", ylab = ylab_sd, line = 2, cex.lab=1.2)
legend(2900, 0.005, legend=lgnds_sd, col=cls_sd, lty=c(3,3,1,5,5), cex=0.8)


##################### Figure bias_n and sd_n ###################################################

n_x = seq(1000, 10000, by=10) # n is on the x-axis

Esota_theta_vec = numeric(length(n_x)) # pre-allocate for bias
SDsota_theta_vec = numeric(length(n_x)) # for standard deviation
for (i in 1:length(n_x)){
  Esota =  expect(n_x[i], theta_SOTA, m)
  Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta_SOTA
  
  Vsota = variance(n_x[i], theta_SOTA, m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
}

bias = 1

# black line
if (bias){
  plot(n_x, Esota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm_bias,
       xlab = "", ylab ="")
  abline(v=n, col="gray") # intersection corresponding to upper row in table
} else {
# x11()
plot(n_x, SDsota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm_sd,
     xlab = "", ylab ="")
abline(v=n, col="gray") # intersection corresponding to upper row in table
}

# two green lines
#################### param = m
param = 1
for (k in 2:3){
  for (i in 1:length(n_x)){
    Esota =  expect(n_x[i], theta_SOTA, m_vec[k])
    Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta_SOTA
    
    Vsota =  variance(n_x[i], theta_SOTA, m_vec[k])
    SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
  }
  if (bias){
  # dev.set(dev.prev())
  par(new=TRUE)   
  plot(n_x, Esota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], ylim=ylm_bias,
       xlab = "", ylab ="")
  } else {
  #dev.set(dev.next())
  par(new=TRUE)   
  plot(n_x, SDsota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], ylim=ylm_sd,
       xlab = "", ylab ="")
  }
}

#################### param = theta
param = 3

for (k in 2:3){
  for (i in 1:length(n_x)){
    Esota =  expect(n_x[i], theta_vec[k], m)
    Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta_vec[k]
    
    Vsota =  variance(n_x[i], theta_vec[k], m)
    SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
  }
  if (bias){
  # dev.set(dev.prev())
  par(new=TRUE)
  plot(n_x, Esota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], 
       ylim=ylm_bias, xlab = "", ylab ="")
  } else {
  # dev.set(dev.next())
  par(new=TRUE)
  plot(n_x, SDsota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], ylim=ylm_sd,
       xlab = "", ylab ="")
  }
}

# title and legends
lgnds_bias = c(TeX(r'(${theta}=0.85$)'),TeX(r'($m = 5000$)'), TeX(r'($m = 1000, {theta}=0.90$)'), 
               TeX(r'($m = 100$)'), TeX(r'(${theta}=0.95$)'))
cls_bias = c("blue","green","black","green","blue")

if (bias){
#dev.set(dev.prev())
title(main = "", xlab = "n", ylab = ylab_bias, line = 2, cex.lab=1.2)
legend(6000, 0.035, legend=lgnds_bias, col=cls_bias, lty=c(3,5,1,3,5), cex=0.8)

lgnds_sd = c(TeX(r'($m = 100$)'),TeX(r'(${theta}=0.85$)'),TeX(r'($n = 1000, {theta}=0.90$)'), 
             TeX(r'($m = 5000$)'), TeX(r'(${theta}=0.95$)'))
cls_sd = c("green","blue","black","green","blue")
} else {
# dev.set(dev.next())
title(main = "", xlab = "n", ylab = ylab_sd, line = 2, cex.lab=1.2)
legend(6000, 0.005, legend=lgnds_sd, col=cls_sd, lty=c(3,3,1,5,5), cex=0.8)
}

################################ Figure bias_theta and sd_theta ###################################################

theta_x = seq(0.85, 0.95, by=0.0001) # theta is on the x-axis

Esota_theta_vec = numeric(length(theta_x)) # pre-allocate
SDsota_theta_vec = numeric(length(theta_x)) # for standard deviation
for (i in 1:length(theta_x)){
  Esota =  expect(n, theta_x[i], m)
  Esota_theta_vec[i] = (1-Esota/n) - theta_x[i]
  
  Vsota = variance(n, theta_x[i], m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
bias = 0
if (bias){
#x11()
plot(theta_x, Esota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm_bias,
     xlab = "", ylab ="")
abline(v=theta_SOTA, col="gray") #intersection corresponding to upper row in table
} else {
#x11()
plot(theta_x, SDsota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm_sd,
     xlab = "", ylab ="")
abline(v=theta_SOTA, col="gray") # intersection corresponding to upper row in table
}
########### param = n
param = 2

for (k in 2:3){
  for (i in 1:length(theta_x)){
    Esota =  expect(n_vec[k], theta_x[i], m)
    Esota_theta_vec[i] = (1-Esota/n_vec[k]) - theta_x[i]
    
    Vsota = variance(n_vec[k], theta_x[i], m)
    SDsota_theta_vec[i] = sqrt(Vsota)/n_vec[k]
  }
  if (bias){
  #dev.set(dev.prev())
  par(new=TRUE)   
  plot(theta_x, Esota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], 
       ylim=ylm_bias,xlab = "", ylab ="")
  } else{
  #dev.set(dev.next())
  par(new=TRUE)   
  plot(theta_x, SDsota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], 
       ylim=ylm_sd, xlab = "", ylab ="")
  }
}

############## param = m
param = 1

for (k in 2:3){
  for (i in 1:length(theta_x)){
    Esota =  expect(n, theta_x[i], m_vec[k])
    Esota_theta_vec[i] = (1-Esota/n) - theta_x[i]
    
    Vsota =  variance(n, theta_x[i], m_vec[k])
    SDsota_theta_vec[i] = sqrt(Vsota)/n
  }
  if (bias){
  # dev.set(dev.prev())
  par(new=TRUE) 
  plot(theta_x, Esota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], 
       ylim=ylm_bias,xlab = "", ylab ="")
  } else {
  # dev.set(dev.next())
  par(new=TRUE)
  plot(theta_x, SDsota_theta_vec,"l", lty = line_vec[k], col = col_vec[param], 
       ylim=ylm_sd, xlab = "", ylab ="")
  }
}

# title and legends
lgnds_bias = c(TeX(r'($n=1000$)'),TeX(r'($m=5000$)'), TeX(r'($n = 3000, m=1000$)'), 
               TeX(r'($m=100$)'), TeX(r'($n=10000$)'))
cls_bias = c("red","green","black","green", "red")

if (bias){
#dev.set(dev.prev())
title(main = "", xlab =  TeX(r'(${theta}$)'), ylab = ylab_bias, line = 2, cex.lab=1.2)
legend(0.905, 0.035, legend=lgnds_bias, col=cls_bias, lty=c(3,5,1,3,5), cex=0.8)
} else {

lgnds_sd = c(TeX(r'($n=1000$)'),TeX(r'($m=100$)'), TeX(r'($n = 3000, m=1000$)'), 
             TeX(r'($m=5000$)'), TeX(r'($n=10000$)'))
cls_sd = c("red","green","black","green", "red")

# dev.set(dev.next())
title(main = "", xlab =  TeX(r'(${theta}$)'), ylab = ylab_sd, line = 2, cex.lab=1.2)
legend(0.905, 0.005, legend=lgnds_sd, col=cls_sd, lty=c(3,3,1,5,5), cex=0.8)
}



##############################################################################################
################################### Simulations ##############################################
##############################################################################################

simul  = 0
if (simul){

# I will here recreate the numbers from above. 

# I will not simulate the binomial distribution and its conf.int. 
ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial
ci_up = ci_binom["upper"][[1]] 

sprintf('Simulate the multiplicity adjusted upper limit of the (1-alpha) confidence interval, %s',
        x_alpha2)

# Draw at random the number of failures in $m$ independent experiments, each having $n$ trials and probability $theta$
x_vec = rbinom(m, n, 1-theta_SOTA)

# Have a quick look
hist(x_vec, xlab = mean(x_vec))

# fz_ma = matrix(0,n+1,rep) # requires too much vector memory rep = 1,000,000
fz_sim = numeric(n+1)

tic()
for (ell in 1: rep){
  fz_rep = numeric(n+1) # clear this for each repetition
  x_vec = rbinom(m, n, 1-theta_SOTA) # the number of failures in each of the m classifiers
  
  # Here, we want to find out for how many repetitions did at least one classifier have i failures, 
  # but none had fewer than i failures
  
  for (i in 0: n){
    j = i+1 # for counting
    if(any(x_vec==i)){ # at least one with i failures
      fz_rep[j] = 1 # then the rest is not 'none had fewer'
      break
    } 
  }
  fz_sim = fz_sim+fz_rep
}

toc() # rowSums(fz_ma)
fz_sim = fz_sim/rep

plot(0:n,fz_sim, type = 's')

# need to zoom in
z = 200:300
plot(z,fz_sim[z], type = 'h', xlab = 'number of failures/accuracy', 
     ylab = 'f(z)')

# Let's see if these add up
Palpha2_discr = pbinom(x_alpha2-1,n,1-theta_SOTA) # 
P_low = 1-pbinom(Cx-1,m,Palpha2_discr) #

Palpha2_discr = pbinom(x_alpha2,n,1-theta_SOTA) # 
P_up = 1-pbinom(Cx-1,m,Palpha2_discr) #
print(c(P_low,P_up))

# Note that x_alpha2 is shifted by 1 because P_up is 1 - P or something: need to verify that
P_simlow = sum(fz_sim[1:x_alpha2])
P_simup = sum(fz_sim[1:x_alpha2+1])
print(c(P_simlow,P_simup))

# Expectation

Eterm = numeric(n+1)
for (z in 0:n){
  i = z+1
  Eterm[z] = z*fz_sim[i]
}

Esota_s = sum(Eterm)
Esota_sim = 1-Esota_s/n

Esota = expect(n, theta_SOTA, m)
Esota_theta = 1-Esota/n

print(c(Esota_sim,Esota_theta))

# Variance

vterm = numeric(n+1)
for (z in 0:n){
  i = z+1
  vterm[i] = z^2*fz_sim[i]
}
esquare = sum(vterm)

Vsota_sim = esquare - Esota_s^2

Vsota = variance(n, theta_SOTA, m)

print(c(Vsota_sim,Vsota))
print(c(sqrt(Vsota_sim)/n, sqrt(Vsota)/n))


# Cumulative distribution function:
fz_cumsum = cumsum(fz_sim)

plot(fz_cumsum, type = 's')

# need to zoom in
z = 200:300
plot(z,fz_cumsum[z], type = 's', xlab = 'number of failures/accuracy', 
     ylab = 'F(z)')

Fz_sim = numeric(n+1)

tic()
for (ell in 1: rep){
  Fz_rep = numeric(n+1)
  x_vec = rbinom(m, n, 1-theta_SOTA) # the number of failures in each of the m classifiers
  
  # Here, we want to find out for how many repetitions did at least one classifier have i or fewer failures, 
  
  for (i in 0: n){
    j = i+1 # for counting
    if(any(x_vec==i)){ # at least one with i failures
      Fz_rep[j:n+1] = 1 # then the rest has 'or fewer'
      break
    } 
  }
  Fz_sim = Fz_sim+Fz_rep
}
toc()
Fz_sim = Fz_sim/rep

plot(Fz_sim, type = 's')

# need to zoom in
z = 200:300
plot(z,Fz_sim[z], type = 's', xlab = 'number of failures/accuracy', 
     ylab = 'F(z)')

# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # 
}

