# Multiple classifiers and biased state-of-the-art estimation
# https://www.overleaf.com/project/63c8012bf045548a94e2d140
# by Kajsa Møllersen (kajsa.mollersen@uit.no) February 2023 updated March 2023

###############################################################################
############################ 3.1 Nomenclature #############################
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

# functions

# simulated pmf - sim_nonid_pmf
# cdf - nonid_cdf
# cdf based pmf - nonid_pmf


###############################################################################
################### Non-identical, independent classifiers ####################
###############################################################################

# Consider a classification problem with a test set of size $n$, and 
# $m$ classifiers with varying probability of correct classification 

source("Parameters_PublicCompetition.R") # n, theta, m, alpha, theta_min, theta_max, theta_vec

source("indep_nonid_pmf_fun.R") # for the functions 
# indep_nonid_pmf' - simulated x
# nonid_pmf - analytical pmf
# nonid_cdf - analytical cdf

# theta_vec is sampled from uni(theta_min,theta_max)
theta_vec = runif(m, min=theta_min, max=theta_max)

tic()
X = indep_nonid_pmf(n, theta_vec, m, rep) # 10 sec for rep = 100,000
toc()

# Histograms of the minimum number of failures for m classifiers, in rep repetitions.
x11()
histbreaks = z_range
hist(X$min_nonid, xlab = 'number of failures', ylab = 'm', breaks = 100, 
     xlim = c(min(z_range),max(z_range))) 

# The upper bound of the 95% confidence interval
sort_min_nonid = sort(X$min_nonid) # sort the minimum number of failures
min_nonid_alpha2 = sort_min_nonid[(alpha/2)*rep] # find the alpha/2 bound

sprintf("The simulated non-identical upper bound of the %s confidence interval is %.5f, with %s repetitions. Distance to SOTA: %s.",  
        1-alpha, (n-min_nonid_alpha2)/n, rep, (n-min_nonid_alpha2)/n-theta_max)


# Analytical results:

####################### Expected value and variance ###################################

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

################################################################################
################################### Figures ####################################
################################################################################

# # # # # # # # # # # # # # # # # Figure noniid_cdf.png # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
x11()
Fz = nonid_cdf(n, theta_vec, m)
# the whole range, not very much information
plot(0:n,Fz, type = 'l', xlab = 'number of failures', 
     ylab = 'probability of at least one team')

# zooming in, and it gets more interesting, discreet curve 
z = z_range
plot(z,Fz[z+1], type = 's', xlab = '', ylab = '', axes = F)

# axis, ticks and labels
xax = seq(z[1],tail(z,1), 20)
klab = xax 
plab = round(1000*(n-xax)/n)/1000
axis(1, cex.axis=1, las = 2, at=xax, labels = as.character(plab))
axis(2, cex.axis=1, las = 2)
axis(3, cex.axis=1, las = 2, at=xax, labels = as.character(klab))
title(main = list(TeX(r'($z$)'), cex = 1.2,
                  col = "black"), sub = list(TeX(r'($\hat{\theta}_{max}$)'),cex = 1.2))
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Does this align with the simulations?

print(c(Fz[min_nonid_alpha2], Fz[min_nonid_alpha2+1])) # ok 
print(c(which(Fz>alpha/2)[1]-1,min_nonid_alpha2)) # ok


################# probability mass function ##########################

# # # # # # # # # # # # # # # # # Figure noniid_pmf # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
fz = nonid_pmf(n, theta_vec, m)

plot(0:n,fz, type = 's')

# need to zoom in
z = z_range
plot(z,fz[z], type = 'h', xlab = '', ylab = '', axes = F)

# axis, ticks and labels, same values as for the cdf
# plab = round(1000*(n-xax)/n)/1000
axis(1, las = 2, at=xax, labels = as.character(plab))
yax = seq(0,max(fz)+0.01,0.02)
axis(2, las = 1, at=yax, labels = as.character(yax))
axis(3, las = 2, at=xax, labels = as.character(klab))
title(main = list(TeX(r'($z$)'), cex = 1.2,
                  col = "black"), sub = list(TeX(r'($\hat{\theta}_{max}$)'),cex = 1.2))
# # # # # # # # # # # # # # # # # end figure # # # # # # # # # # # # # # # # #


############################ bias_sd_thetamin #################################


theta_min_vec = seq(0.5, 0.9, by=0.005) # adjust until smooth
print(length(theta_min_vec))

Bias_theta_vec = numeric(length(theta_min_vec))
SD_theta_vec = numeric(length(theta_min_vec))

B = 500

fzb = matrix(B,n+1) #fz returns vector of length n+1


for (j in 1:length(theta_min_vec)){
  theta_min = theta_min_vec[j] # for the non-identical
  theta_max = theta 
  step = (theta_max-theta_min)/(m-1)
  theta_vec = seq(theta_min, theta_max, step)
  
  Esotab = numeric(B)
  Vsotab = numeric(B)
  for (b in 1:B){
    theta_vec = runif(m, min=theta_min, max=theta_max)
    fz = nonid_pmf(n, theta_vec, m)
  
    # Expected value and variance
    Eterm = numeric(n+1)
    Vterm = numeric(n+1)
    for (z in 0:n){
      i = z+1
      Eterm[i] = z*fz[i]
      Vterm[i] = z^2*fz[i]
    }
  
    Esotab[b] = sum(Eterm)
    Vsotab[b] = sum(Vterm)-Esotab[b]^2
  }
  Esota = mean(Esotab)
  Vsota = mean(Vsotab) + var(Esotab)
  
  Bias_theta_vec[j] = (1-Esota/n)-theta
  SD_theta_vec[j] = sqrt(Vsota)/n
  
  print(c(j,length(theta_min_vec)-j))
}
cols = c("black","darkgreen")
plot(theta_min_vec, Bias_theta_vec,"l", lty = "solid", col = cols[1], 
     ylim = ylm_bias, xlab = "", ylab = "")
par(new=TRUE)
plot(theta_min_vec, SD_theta_vec,"l", lty = "solid", col = cols[2], 
     ylim = ylm_sd,  xlab = "", ylab = "", axes = FALSE)
axis(4)

abline(v=0.875, col="gray")
abline(v=0.85, col="gray",lty = 5)
abline(v=0.825, col="gray",lty = 4)
abline(v=0.8, col="gray",lty = 3)

title(main = "", xlab = TeX(r'(${min}{(Theta)}$)'), ylab = "", line = 2, cex.lab=1.2)
legend(0.6, 0.005, legend=c(ylab_bias,ylab_sd), col=cols, lty=c(1,1), cex=0.8)



