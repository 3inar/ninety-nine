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


# comes with the functions 'cdf', 'pmf', 'expect' and 'variance'
source("ProbDistr_thetaSOTA.R") 

# figures at the end

###############################################################################
################# 3.4 A simulated public competition example ##################
###############################################################################

# Consider a classification problem with a test set of size $n=3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The classifier's performance on the test set, denoted by $theta_hat$ and referred 
# to as the accuracy is $(n-x_j)/n$, where $x_j$ is the observed number of failures
# for classifier $j$.

source("Parameters_PublicCompetition.R") # n, theta, m, alpha, mu = n*theta, rep

# 95% confidence interval for hat{theta} = mu, single classifier

ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial

sprintf("The %s confidence interval for an estimated accuracy of %s is (%.4f,%.4f).",  
        (1-alpha)*100, theta, ci_binom["lower"], ci_binom["upper"])

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# The probability of exceeding the upper limit of the CI should be close to alpha/2 = 0.025 
k_up = floor(ci_binom[["upper"]]*n) # number of successes exceeding the CI
# flooring the CI bound, so P_up > alpha/2
P_up = pbinom(k_up,n,theta, lower.tail = F) # P[X>x]
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

x_alpha2 = qbinom(Px, n, 1-theta) 
theta_alpha2 = (n-x_alpha2)/n
sprintf("With a probablitiy of alpha/2 = %s, at least %s team will achieve an accuracy of at least %.4f, corresponding to x = %s failures.",
        alpha/2, Cx, theta_alpha2, x_alpha2)

# # # # # # # # # # # # # # # # # check-up # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Fz = cdf(n,theta,m) # the cdf
x = which(Fz>alpha/2)[1]-1 # should be same as x_alpha2 = 236
# # # # # # # # # # # # # # # # # ok # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # Check-up # # # # # # # # # # # # # # # # # # # # # # # 
# Since x_alpha2 is discreet, but Px is not, we'll calculate for x_alpha2 and x_alpha2-1, 
# and then the two results should be on each side of alpha/2

Palpha2_discr = pbinom(x_alpha2-1,n,1-theta) # 
P_low = 1-pbinom(Cx-1,m,Palpha2_discr) #

Palpha2_discr = pbinom(x_alpha2,n,1-theta) # 
P_up = 1-pbinom(Cx-1,m,Palpha2_discr) #
# # # # # # # # # # # # # # # # # ok 0.02475 0.03243# # # # # # # # # # # # # # # # # # # # # # # 

sprintf("In summary, if there are %s teams, each with %s probability of correct prediction, and a test set of size %s, there is a probability of %s that at least %s team will have at most %s incorrect predictions, corresponding to an estimated accuracy of %.4f.",  
        m, theta, n, alpha/2, Cx, x_alpha2, theta_alpha2)

# Expected value
Esota = expect(n,theta,m)
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

######################### varying m, n, p ####################################

n_vec = c(3000,1000,10000)
m_vec = c(1000,100,5000)
theta_vec = c(0.9,0.85,0.95)

for (i in 1:3){
  for (j in 1:3){
    for (k in 1:3){
      Esota =  expect(n_vec[i], theta_vec[k], m_vec[j])
      Esota_theta = 1-Esota/n_vec[i]
      Vsota = variance(n_vec[i], theta_vec[k], m_vec[j])
      Std_sota_theta = sqrt(Vsota)/n_vec[i]
      sprintf("The expected theta_sota is %.4f, with a standard deviation of %.6f, for m=%s, n=%s, theta=%s.",
              Esota_theta, Std_sota_theta, m_vec[j], n_vec[i], theta_vec[k])
    }
  }
}

j = 1; i = 1; k = 3
Esota = expect(n_vec[i], theta_vec[k], m_vec[j])
Esota_theta = 1-Esota/n_vec[i]
Vsota = variance(n_vec[i], theta_vec[k], m_vec[j])
Std_sota_theta = sqrt(Vsota)/n_vec[i]
sprintf("The expected theta_sota is %.4f, with a standard deviation of %.6f, for m=%s, n=%s, theta=%s.",
        Esota_theta, Std_sota_theta, m_vec[j], n_vec[i], theta_vec[k])



##############################################################################################
################################### Figures ##################################################
##############################################################################################

################################### Figure multi_ci #################################################
Esota = expect(n, theta, m)

# let k be the number of successes, = n-x

k = (mu-60):(mu+90) # this is the plot range, adjust to your liking
kslim = c(k[1],k[length(k)])
whylim = c(0,0.03)

y = dbinom(k, n, theta) # probability of x successes in n trials
plot(k, y, type='l', col = "blue", xlab = '', ylab = '', xlim = kslim, ylim = whylim, axes=F)
par(new=TRUE) # new plot in same window
plot(c(ci_binom[["lower"]]*n, ci_binom[["upper"]]*n), c( -0.0005,-0.0005), "l", col = "blue",  xlab = '', ylab = '', 
     xlim = kslim, ylim = whylim, axes=F)
par(new=TRUE) # new plot in same window
plot(c(ci_binom[["lower"]]*n, ci_binom[["lower"]]*n), c(-0.002,0.001),"l", col="blue", xlab = '', ylab = '', 
     xlim = kslim, ylim = whylim, axes=F)
par(new=TRUE) # new plot in same window
plot(c(ci_binom[["upper"]]*n, ci_binom[["upper"]]*n), c(-0.002,0.001),"l", col="blue", xlab = '', ylab = '', 
     xlim = kslim, ylim = whylim, axes=F)
par(new=TRUE) # new plot in same window
plot(c(theta*n, theta*n), c(0,max(y)),"l", lty = 5, col=rgb(0, 0, 1,0.25), xlab = '', ylab = '', 
     xlim = kslim, ylim = whylim, axes=F)

par(new=TRUE) # new plot in same window

y = dbinom(k, n, ci_binom[["upper"]]) # significantly better classifier
plot(k, y, type = 'l', col = "red", xlab = '', ylab = '', axes=F, xlim = kslim, ylim = whylim)
par(new=TRUE) # new plot in same window
plot(c(ci_binom[["upper"]]*n, ci_binom[["upper"]]*n), c(0,max(y)),"l", lty = 5, col=rgb(1, 0, 0,0.25), xlab = '', ylab = '', 
     xlim = kslim, ylim = whylim, axes=F)

# area under curve for expected SOTA performance
polygon(c(n-Esota, k[k>=n-Esota], max(k)), c(0,y[k>=n-Esota], 0), col=rgb(0, 1, 0,0.25)) 

# area under curve for SOTA performance
k_alpha2 = n-x_alpha2
polygon(c(k_alpha2, k[k>=k_alpha2], max(k)), c(0,y[k>=k_alpha2], 0), col=rgb(0, 1, 0,0.25)) 

axis(1 , cex.axis=1.2, las = 2, at=c(kslim[1],ci_binom[["lower"]]*n, n*theta, 
                       ci_binom[["upper"]]*n, n-Esota, k_alpha2,kslim[2]), 
     labels=c('',TeX('$\\theta_{alpha/2}$'), TeX('$\\theta$'), TeX('$\\theta_{1-alpha/2}$'), 
              TeX(r'($E \hat{\theta}_{max}$)'), TeX('$\\theta^m_{1-alpha/2}$'),''))

title(ylab = '', line=2, cex.lab=1.2, xlab = '')
# # # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################### Figure cumul_fail ##########################################################
Fz = cdf(n,theta,m)

# the whole range, not very much information
plot(0:n,Fz, type = 'l', xlab = 'number of failures', 
     ylab = 'probability of at least one team')

# zooming in, and it gets more interesting, discreet curve 
z = 200:300
plot(z,Fz[z+1], type = 's', xlab ='', 
     ylab = '', axes = F)

xax = seq(z[1],tail(z,1), 20)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, cex.axis=1, las = 2, at=xax, labels = as.character(plab))
axis(2, cex.axis=1, las = 2)
axis(3, cex.axis=1, las = 2, at=xax, labels = as.character(klab))
title(main = list(TeX(r'($z$)'), cex = 1.2,
                  col = "black"), sub = list(TeX(r'($\hat{\theta}_{max}$)'),cex = 1.2))
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################### Figure pmf_sota ###################################################

fz = pmf(n,theta,m,f0 = T)

plot(0:n,fz, type = 's')

# need to zoom in
z = 200:300
plot(z,fz[z], type = 'h', xlab ='', 
     ylab = '', axes = F)
# xax = seq(z[1],tail(z,1), 10)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
yax = seq(0,max(fz)+0.01,0.02)
axis(2, las = 1, at=yax, labels = as.character(yax))
axis(3, las = 2, at=xax, labels = as.character(klab))
title(main = list(TeX(r'($z$)'), cex = 1.2,
                  col = "black"), sub = list(TeX(r'($\hat{\theta}_{max}$)'),cex = 1.2))
# # # # # # # # # # # # # # # # # end figure 3 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


################################### Figure bias_m ###################################################

ylm = c(0.0,0.035)

m_x = seq(1, 5000, by=10)
Esota_theta_vec = numeric(length(m_x))

for (i in 1:length(m_x)){
  Esota =  expect(n, theta, m_x[i])
  Esota_theta_vec[i] = (1-Esota/n) - theta
}
plot(m_x, Esota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ns

for (i in 1:length(m_x)){
  Esota =  expect(n_vec[2], theta, m_x[i])
  Esota_theta_vec[i] = (1-Esota/n_vec[2]) - theta
}
plot(m_x, Esota_theta_vec,"l", lty = "dotted", col = "red", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(m_x)){
  Esota =  expect(n_vec[3], theta, m_x[i])
  Esota_theta_vec[i] = (1-Esota/n_vec[3]) - theta
}
plot(m_x, Esota_theta_vec,"l", lty = "longdash", col = "red", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three thetas

for (i in 1:length(m_x)){
  Esota =  expect(n, theta_vec[2], m_x[i])
  Esota_theta_vec[i] = (1-Esota/n) - theta_vec[2]
}
plot(m_x, Esota_theta_vec,"l", lty = "dotted", col = "blue", ylim=ylm, 
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(m_x)){
  Esota =  expect(n, theta_vec[3], m_x[i])
  Esota_theta_vec[i] = (1-Esota/n) - theta_vec[3]
}
plot(m_x, Esota_theta_vec,"l", lty = "longdash", col = "blue", ylim=ylm, ylab = '', xlab = '')
title(main = "", xlab = "m", ylab = TeX(r'($E \hat{theta}_{\max} - {theta}_{SOTA}$)'), line = 2, cex.lab=1.2)

legend(2900, 0.035, legend=c(TeX(r'($n=1000$)'),TeX(r'(${theta}=0.85$)'), TeX(r'($n = 3000, {theta}=0.90$)'), 
                             TeX(r'(${theta}=0.95$)'), TeX(r'($n=10000$)')),
       col=c("red","blue","black","blue", "red"), lty=c(3,3,1,5,5), cex=0.8)


abline(v=1000, col="gray")

################################### Figure bias_n ###################################################

ylm = c(0.0,0.035)

n_x = seq(1000, 10000, by=10)
Esota_theta_vec = numeric(length(n_x))

for (i in 1:length(n_x)){
  Esota =  expect(n_x[i], theta, m)
  Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta
}
plot(n_x, Esota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ms

for (i in 1:length(n_x)){
  Esota =  expect(n_x[i], theta, m_vec[2])
  Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta
}
plot(n_x, Esota_theta_vec,"l", lty = "dotted", col = "green", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(n_x)){
  Esota =  expect(n_x[i], theta, m_vec[3])
  Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta
}
plot(n_x, Esota_theta_vec,"l", lty = "longdash", col = "green", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three thetas

for (i in 1:length(n_x)){
  Esota =  expect(n_x[i], theta_vec[2], m)
  Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta_vec[2]
}
plot(n_x, Esota_theta_vec,"l", lty = "dotted", col = "blue", ylim=ylm, 
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(n_x)){
  Esota =  expect(n_x[i], theta_vec[3], m)
  Esota_theta_vec[i] = (1-Esota/n_x[i]) - theta_vec[3]
}
plot(n_x, Esota_theta_vec,"l", lty = "longdash", col = "blue", ylim=ylm,ylab = '', xlab = '')
title(main = "", xlab = "n", ylab = TeX(r'($E \hat{theta}_{\max} - {theta}_{SOTA}$)'), line = 2, cex.lab=1.2)
     

abline(v=3000, col="gray")

legend(6000, 0.035, legend=c(TeX(r'(${theta}=0.85$)'),TeX(r'($m = 5000$)'), TeX(r'($n = 1000, {theta}=0.90$)'), 
                             TeX(r'($m = 100$)'), TeX(r'(${theta}=0.95$)')),
       col=c("blue","green","black","green","blue"), lty=c(3,5,1,3,5), cex=0.8)




################################### Figure bias_theta ###################################################

ylm = c(0.0,0.035)

theta_x = seq(0.85, 0.95, by=0.0001)
Esota_theta_vec = numeric(length(theta_x))

for (i in 1:length(theta_x)){
  Esota =  expect(n, theta_x[i], m)
  Esota_theta_vec[i] = (1-Esota/n) - theta_x[i]
}
plot(theta_x, Esota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ns

for (i in 1:length(theta_x)){
  Esota =  expect(n_vec[2], theta_x[i], m)
  Esota_theta_vec[i] = (1-Esota/n_vec[2]) - theta_x[i]
}
plot(theta_x, Esota_theta_vec,"l", lty = "dotted", col = "red", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(theta_x)){
  Esota =  expect(n_vec[3], theta_x[i], m)
  Esota_theta_vec[i] = (1-Esota/n_vec[3]) - theta_x[i]
}
plot(theta_x, Esota_theta_vec,"l", lty = "longdash", col = "red", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ms

for (i in 1:length(theta_x)){
  Esota =  expect(n, theta_x[i], m_vec[2])
  Esota_theta_vec[i] = (1-Esota/n) - theta_x[i]
}
plot(theta_x, Esota_theta_vec,"l", lty = "dotted", col = "green", ylim=ylm, 
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(theta_x)){
  Esota =  expect(n, theta_x[i], m_vec[3])
  Esota_theta_vec[i] = (1-Esota/n) - theta_x[i]
}
plot(theta_x, Esota_theta_vec,"l", lty = "longdash", col = "green", ylim=ylm,ylab = '', xlab = '')
title(main = "", xlab = TeX(r'(${theta}$)'), ylab = TeX(r'($E \hat{theta}_{\max} - {theta}_{SOTA}$)'), line = 2, cex.lab=1.2)

legend(0.905, 0.035, legend=c(TeX(r'($n=1000$)'),TeX(r'($m=5000$)'), TeX(r'($n = 3000, m=1000$)'), 
                              TeX(r'($m=100$)'), TeX(r'($n=10000$)')),
       col=c("red","green","black","green", "red"), lty=c(3,5,1,3,5), cex=0.8)

abline(v=0.9, col="gray")

################################### Figure sd_m ###################################################

ylm = c(0.0,0.005)

m_x = seq(1, 5000, by=10)
SDsota_theta_vec = numeric(length(m_x))

for (i in 1:length(m_x)){
  Vsota = variance(n, theta, m_x[i])
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(m_x, SDsota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ns

for (i in 1:length(m_x)){
  Vsota =  variance(n_vec[2], theta, m_x[i])
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(m_x, SDsota_theta_vec,"l", lty = "dotted", col = "red", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(m_x)){
  Vsota =  variance(n_vec[3], theta, m_x[i])
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(m_x, SDsota_theta_vec,"l", lty = "longdash", col = "red", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three thetas

for (i in 1:length(m_x)){
  Vsota =  variance(n, theta_vec[2], m_x[i])
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(m_x, SDsota_theta_vec,"l", lty = "dotted", col = "blue", ylim=ylm, 
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(m_x)){
  Vsota =  variance(n, theta_vec[3], m_x[i])
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(m_x, SDsota_theta_vec,"l", lty = "longdash", col = "blue", ylim=ylm, ylab = '', xlab = '')
title(main = "", xlab = "m", ylab = TeX(r'($\sigma_{\hat{\theta}_{\max}}$)'), line = 2, cex.lab=1.2)

legend(2900, 0.005, legend=c(TeX(r'(${theta}=0.85$)'),TeX(r'($m = 5000$)'), TeX(r'($n = 1000, {theta}=0.90$)'), 
                             TeX(r'($m = 100$)'), TeX(r'(${theta}=0.95$)')),
       col=c("blue","green","black","green","blue"), lty=c(3,5,1,3,5), cex=0.8)


abline(v=1000, col="gray")

################################### Figure sd_n ###################################################

ylm = c(0.0,0.005)

n_x = seq(1000, 10000, by=10)
SDsota_theta_vec = numeric(length(n_x))

for (i in 1:length(n_x)){
  Vsota = variance(n_x[i], theta, m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
}
plot(n_x, SDsota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ms

for (i in 1:length(n_x)){
  Vsota =  variance(n_x[i], theta, m_vec[2])
  SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
}
plot(n_x, SDsota_theta_vec,"l", lty = "dotted", col = "green", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(n_x)){
  Vsota =  variance(n_x[i], theta, m_vec[3])
  SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
}
plot(n_x, SDsota_theta_vec,"l", lty = "longdash", col = "green", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three thetas

for (i in 1:length(n_x)){
  Vsota =  variance(n_x[i], theta_vec[2], m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
}
plot(n_x, SDsota_theta_vec,"l", lty = "dotted", col = "blue", ylim=ylm, 
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(n_x)){
  Vsota =  variance(n_x[i], theta_vec[3], m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n_x[i]
}
plot(n_x, SDsota_theta_vec,"l", lty = "longdash", col = "blue", ylim=ylm, ylab = '', xlab = '')
title(main = "", xlab = "n", ylab = TeX(r'($\sigma_{\hat{\theta}_{\max}}$)'), line = 2, cex.lab=1.2)

abline(v=3000, col="gray")

legend(6000, 0.005, legend=c(TeX(r'(${theta}=0.85$)'),TeX(r'($m = 5000$)'), TeX(r'($n = 1000, {theta}=0.90$)'), 
                             TeX(r'($m = 100$)'), TeX(r'(${theta}=0.95$)')),
       col=c("blue","green","black","green","blue"), lty=c(3,5,1,3,5), cex=0.8)


################################### Figure sd_theta ###################################################

ylm = c(0.0,0.005)

theta_x = seq(0.85, 0.95, by=0.0001)
SDsota_theta_vec = numeric(length(n_x))

for (i in 1:length(theta_x)){
  Vsota = variance(n, theta_x[i], m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(theta_x, SDsota_theta_vec,"l", lty = "solid", col = "black", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ms

for (i in 1:length(theta_x)){
  Vsota =  variance(n, theta_x[i], m_vec[2])
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(theta_x, SDsota_theta_vec,"l", lty = "dotted", col = "green", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(theta_x)){
  Vsota =  variance(n, theta_x[i], m_vec[3])
  SDsota_theta_vec[i] = sqrt(Vsota)/n
}
plot(theta_x, SDsota_theta_vec,"l", lty = "longdash", col = "green", ylim=ylm,
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

# with the three ns

for (i in 1:length(theta_x)){
  Vsota =  variance(n_vec[2], theta_x[i], m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n_vec[2]
}
plot(theta_x, SDsota_theta_vec,"l", lty = "dotted", col = "red", ylim=ylm, 
     xlab = "", ylab ="")
par(new=TRUE) # new plot in same window

for (i in 1:length(theta_x)){
  Vsota =  variance(n_vec[3], theta_x[i], m)
  SDsota_theta_vec[i] = sqrt(Vsota)/n_vec[3]
}
plot(theta_x, SDsota_theta_vec,"l", lty = "longdash", col = "red", ylim=ylm, ylab = '', xlab = '')
title(main = "", xlab = TeX(r'(${theta}$)'), ylab = TeX(r'($\sigma_{\hat{\theta}_{\max}}$)'), line = 2, cex.lab=1.2)


legend(0.905, 0.005, legend=c(TeX(r'($n=1000$)'),TeX(r'($m=100$)'), TeX(r'($n = 3000, m=1000$)'), 
                              TeX(r'($m=5000$)'), TeX(r'($n=10000$)')),
       col=c("red","green","black","green", "red"), lty=c(3,3,1,5,5), cex=0.8)

abline(v=0.9, col="gray")






##############################################################################################
################################### Simulations ##############################################
##############################################################################################

# I will here recreate the numbers from above. 

# I will not simulate the binomial distribution and its conf.int. 
ci_binom = binom.confint(mu,n,conf.level=1-alpha, methods = "exact") # CI for binomial
ci_up = ci_binom["upper"][[1]] 

sprintf('Simulate the multiplicity adjusted upper limit of the (1-alpha) confidence interval, %s',
        x_alpha2)

# Draw at random the number of failures in $m$ independent experiments, each having $n$ trials and probability $theta$
x_vec = rbinom(m, n, 1-theta)

# Have a quick look
hist(x_vec, xlab = mean(x_vec))

# fz_ma = matrix(0,n+1,rep) # requires too much vector memory rep = 1,000,000
fz_sim = numeric(n+1)

tic()
for (ell in 1: rep){
  fz_rep = numeric(n+1) # clear this for each repetition
  x_vec = rbinom(m, n, 1-theta) # the number of failures in each of the m classifiers
  
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
Palpha2_discr = pbinom(x_alpha2-1,n,1-theta) # 
P_low = 1-pbinom(Cx-1,m,Palpha2_discr) #

Palpha2_discr = pbinom(x_alpha2,n,1-theta) # 
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

Esota = expect(n, theta, m)
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

Vsota = variance(n, theta, m)

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
  x_vec = rbinom(m, n, 1-theta) # the number of failures in each of the m classifiers
  
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

