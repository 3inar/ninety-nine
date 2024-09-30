# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) February 2023 updated March 2023

###############################################################################
############################ Nomenclature #############################
###############################################################################

# n     - number of trials/size of test set
# m     - number of experiments/number of classifiers
# theta     - prob of success/prob of correct prediction

# Y     - r.v., indicates if prediction is correct

# k     - number of successes
# X     - r.v., number of failures; n-k

# P_    - probability, specified when used

# px    - probability of x failures in one experiment
# Px    - probability of at most x failures in one experiment

# Cx    - r.v., number of experiments with at most x failures

# theta_hat - accuracy of a classifier: p_hat = (n-x)/n

# theta_hat_SOTA - max_j p_hat_j

# Z     - r.v., number of failures on at least one classifier

# Fz    - cdf of Z: P(C_z > 0|m,n,theta) prob of at least one classifier having at most 
#         z failures (identical to at least n-z successes)
# fz    - pmf of Z

# alpha - significance level

# SOTA - state-of-the-art


library(binom)      # confidence interval for binomial distribution
library(tictoc)     # for timing
library(latex2exp)  # mathematical notation


###############################################################################
############################ Dependent non-identical #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is denoted by $theta_hat$.

# Functions

# dep_nonid_pmf - simulated pmf

source("Parameters_PublicCompetition.R") # n, theta, m, alpha, rho, theta_min, theta_max, theta_vec

length(theta_vec)

source("dep_nonid_pmf_fun.R") # for function dep_nonid_pmf
# returns X = list(min_fail, x_fail, teamsSOTA)


tic()
X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta)
toc()

# Histograms of the minimum number of failures for m classifiers, in rep repetitions.
source("plotting_params.R")

new_png("test.png", n_figures=1)
hist(X$x_fail, xlab = 'number of failures', ylab = 'number of classifiers', 
     ylim = c(0,250), main=NULL)
dev.off()


hist(X$x_fail, xlab = 'number of failures', ylab = 'number of classifiers', 
     ylim = c(0,250))


source("ProbDistr_thetaSOTA.R")
# The upper bound of the 95% confidence interval
min_dep_alpha2 = sim_ci(alpha, X$min_fail)
sprintf("The simulated dependent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)

# The expected value
Esota = mean(X$min_fail)

# The standard deviation
Vsota = mean(X$min_fail*X$min_fail) - Esota*Esota

sprintf("The simulated expected value is %.7f and a standard deviation is %.7f, with %s repetitions.",  
        (n-Esota)/n, sqrt(Vsota)/n, rep)

###########################################################################
### Expected value/bias/variance as a function of rho ################
###########################################################################

# theta_min can be only so small
trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))

rho_vec = seq(0.0, 1.0, by=0.01)

ylm = c(0.0,0.035)

theta_min_vec = seq(0.8,theta, by=0.025) 


Esota_theta_vec = matrix(0,length(theta_min_vec),length(rho_vec))
SDsota_theta_vec = matrix(0,length(theta_min_vec),length(rho_vec))

rep = 1000
theta_max = theta 
for (i in 1:length(theta_min_vec)){
  theta_min = theta_min_vec[i] # for the non-identical
  
  # theta_min = max(trunc_min,theta_min_vec[i])
  step = (theta_max-theta_min)/(m-1)
  theta_vec = seq(theta_min, theta_max, step)
  
  for (j in 1:length(rho_vec)){ #
    rho = rho_vec[j] # for the correlated
    # trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))
    
    
    
 #    print(trunc_min)
#     print(theta_min_vec[i])
  
    X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta)
    #X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta)
  
    # Expected value
    Esota = mean(X$min_fail)
    Esota_theta_vec[i,j] = (1-Esota/n)-theta
    
    # Variance
    Vsota = mean(X$min_fail*X$min_fail) - Esota*Esota
    SDsota_theta_vec[i,j] = sqrt(Vsota)/n
    
    
    print(c(i,j))
  }
}

new_png("test2.png", n_figures=2)

# plotting for theta_min = 0.875
plot(rho_vec, Esota_theta_vec[4,],"l", lty = 1, col = "black", ylim = ylm,
     xlab = "", ylab = "")
par(new=TRUE) # new plot in same window

for (i in 1:3){ 
  plot(rho_vec, Esota_theta_vec[i,],"l", lty = i+2, col = "black", ylim = ylm,
       xlab = "", ylab = "")
  par(new=TRUE) # new plot in same window
}

# plotting for theta_min = theta
plot(rho_vec, Esota_theta_vec[5,],"l", lty = 1, col = "red", ylim = ylm,
     main = '', xlab = '', ylab = '')
par(new=TRUE) # new plot in same window
title(ylab = TeX(r'($E \hat{\theta}_{\max} - {theta}_{SOTA}$)'), line=2, cex.lab=1.2, xlab = TeX(r'(${rho}_0$)'))

abline(v=0.6, col="gray")

legend(0.55, 0.035, legend=c(TeX(r'(${min}{(Theta)}=0.9$)'),TeX(r'(${min}{(Theta)}=0.875$)'), TeX(r'(${min}{(Theta)}=0.850$)'), 
                           TeX(r'(${min}{(Theta)}=0.825$)'), TeX(r'(${min}{(Theta)}=0.800$)')),
       col=c("red","black","black","black","black"), lty=c(1,1,5,4,3), cex=0.8)

dev.off()

######################################### standard deviation
# plotting for theta_min = 0.875
plot(rho_vec, SDsota_theta_vec[4,],"l", lty = 1, col = "black", ylim = c(0,0.0055),
     xlab = "", ylab = "")
par(new=TRUE) # new plot in same window

for (i in 1:3){ 
  plot(rho_vec, SDsota_theta_vec[i,],"l", lty = i+2, col = "black", ylim = c(0,0.0055),
       xlab = "", ylab = "")
  par(new=TRUE) # new plot in same window
}

# plotting for theta_min = theta
plot(rho_vec, SDsota_theta_vec[5,],"l", lty = 1, col = "red", ylim = c(0,0.0055),
     main = '', xlab = '', ylab = '')
par(new=TRUE) # new plot in same window
title(ylab = TeX(r'($\sigma_{\hat{\theta}_{\max}}$)'), line=2, cex.lab=1.2, xlab = TeX(r'(${rho}_0$)'))

abline(v=0.6, col="gray")

legend(0.55, 0.0026, legend=c(TeX(r'(${min}{(Theta)}=0.9$)'),TeX(r'(${min}{(Theta)}=0.875$)'), TeX(r'(${min}{(Theta)}=0.850$)'), 
                             TeX(r'(${min}{(Theta)}=0.825$)'), TeX(r'(${min}{(Theta)}=0.800$)')),
       col=c("red","black","black","black","black"), lty=c(1,1,5,4,3), cex=0.8)



