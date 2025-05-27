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

library(future)     # paralellisation
library(furrr)


# functions

# simulated pmf - sim_nonid_pmf
# cdf - nonid_cdf
# cdf based pmf - nonid_pmf


###############################################################################
################### Non-identical, independent classifiers ####################
###############################################################################


future::plan("multisession", workers=6) # let's engage 6 out of 8 cores


# Consider a classification problem with a test set of size $n$, and 
# $m$ classifiers with varying probability of correct classification 

source("Parameters_PublicCompetition.R") # n, theta, m, alpha, theta_SOTA, d = b-a

source("ProbDistr_thetaSOTA.R") # for function sim_ci()

source("indep_nonid_pmf_fun.R") # for the functions 
# indep_nonid_pmf' - simulated x
# nonid_pmf - analytical pmf
# nonid_cdf - analytical cdf, Fz = numeric(n+1)

# theta_vec is sampled from uni(theta_min,theta_max)
theta_max = theta_SOTA + d/(m+1)
theta_min = theta_max - d

B = 200

tic()
fz = furrr::future_map(1:B, function (x) { 

  theta_vec = runif(m, min=theta_min, max=theta_max)
  X = indep_nonid_pmf(n, theta_vec, m, rep) # 10 sec for rep = 100,000
  
  return(X$min_nonid)
}, .progress=T, .options= furrr::furrr_options(seed=T))

# Not very elegant with the for-loop, but that's ok for now
fzb = matrix(data = NA, nrow = B, ncol = rep) 
Esotab = numeric(B)
Vsotab = numeric(B)
alpha2b = numeric(B)
for (b in 1:B){ 
  fzb[b,] = fz[[b]]
  Esotab[b] = mean(fz[[b]])
  Vsotab[b] = mean(fz[[b]]*fz[[b]]) - Esotab[b]*Esotab[b]
  alpha2b[b] = sim_ci(alpha, fz[[b]])
}
toc() # B = 25, loop 256 sec, parallel 70 sec. B = 200, parallel 501 sec

# An example histograms of the minimum number of failures for m classifiers, in rep repetitions.
# x11()
histbreaks = z_range
hist(fz[[1]], xlab = 'number of failures', ylab = 'm', breaks = 100, 
     xlim = c(min(z_range),max(z_range))) 

Esota = mean(Esotab)
Vsota = mean(Vsotab) + var(Esotab)
min_dep_alpha2 = mean(alpha2b) # not sure about this one

sprintf("The simulated nonidentical upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)
# B = 100: 0.91781 
# B = 200: 0.91782 
# B = 200: 0.91780

sprintf("The simulated expected value is %.7f and a standard deviation is %.7f, with %s repetitions.",  
        (n-Esota)/n, sqrt(Vsota)/n, rep)
# B = 100: 0.9129823 and 0.0021333
# B = 200: 0.9129668 and 0.0021338
# B = 200: 0.9129637 and 0.0021353



# "Analytical" results:

####################### Expected value and variance ###################################

Esotab = numeric(B)
Vsotab = numeric(B)
tic()
EVsotab <- furrr::future_map(1:B, function (x) { 
  
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

  Esota = sum(Eterm)
  Vsota = sum(Vterm)-Esota^2
  
  return(cbind(Esota, Vsota))
  
}, .progress=T, .options= furrr::furrr_options(seed=T))
toc()

# Not very elegant with the for-loop, but that's ok for now
Esotab = numeric(B)
Vsotab = numeric(B)
for (b in 1:B){ 
  Esotab[b] = EVsotab[[b]][1]
  Vsotab[b] = EVsotab[[b]][2]
}

Esota = mean(Esotab)
Vsota = mean(Vsotab) + var(Esotab)

sprintf("The expected number of failures is %.4f, with a variance of %.4f.",
        Esota, Vsota)
sprintf("The expected theta_hat_SOTA is %.6f, with standard deviation of %.6f.",
        (n-Esota)/n, sqrt(Vsota)/n)

# theta-vector: 0.9130   0.0021
# hierarchical: 
# B = 200:      0.912963 0.002135

################################################################################
################################### Figures ####################################
################################################################################

# # # # # # # # # # # # # # # # # Figure noniid_cdf.png # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

B = 1000

# paralellised version: takes 5-10 minutes

Fb <- furrr::future_map(1:B, function (x) { 
  theta_vec = runif(m, min=theta_min, max=theta_max)
  Fb = nonid_cdf(n, theta_vec, m)
  
  return(Fb)
    
}, .progress=T, .options= furrr::furrr_options(seed=T))

# Not very elegant with the for-loop, but that's ok for now
Fzb = matrix(data = NA, nrow = B, ncol = n+1) 
for (b in 1:B){ 
  Fzb[b,] = Fb[[b]]
}

Fz = colMeans(Fzb)
  
# saving
# naming convention: fig name_vector_name
saveRDS(Fz, file = "noniid_cdf_Fz.rds")

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

fb <- furrr::future_map(1:B, function (x) { 
  theta_vec = runif(m, min=theta_min, max=theta_max)
  fb = nonid_pmf(n, theta_vec, m)
  
  return(fb)
  
}, .progress=T, .options= furrr::furrr_options(seed=T))

# Not very elegant with the for-loop, but that's ok for now
fzb = matrix(data = NA, nrow = B, ncol = n+1) 
for (b in 1:B){ 
  fzb[b,] = fb[[b]]
}

fz = colMeans(fzb)

# saving
# naming convention: fig name_vector_name
saveRDS(fz, file = "noniid_pmf_fz.rds")

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


############################ bias_sd_d #################################

d_vec = seq(0, 0.4, by=0.005) # adjust until smooth, d = 0.4 gives \theta = 0.5
print(length(d_vec))

Bias_theta_vec = numeric(length(d_vec))
SD_theta_vec = numeric(length(d_vec))


# paralellised version: this takes several hours, like 6 maybe or more

for (j in 1:length(d_vec)){ 
  
  theta_max = theta_SOTA + d_vec[j]/(m+1)
  theta_min = theta_max - d_vec[j]

  EVsotab <- furrr::future_map(1:B, function (x) { 
    theta_vec = runif(m, min=theta_min, max=theta_max)
    #source("indep_nonid_pmf_fun.R") 
    fz = nonid_pmf(n, theta_vec, m)
    
    # Expected value and variance
    Eterm = numeric(n+1)
    Vterm = numeric(n+1)
    for (z in 0:n){
      i = z+1
      Eterm[i] = z*fz[i]
      Vterm[i] = z^2*fz[i]
    }
    
    Esotab = sum(Eterm)
    Vsotab = sum(Vterm)-Esotab^2
    
    return(cbind(Esotab, Vsotab))
  
  }, .progress=T, .options= furrr::furrr_options(seed=T))
  
  # Not very elegant with the for-loop, but that's ok for now
  Esotab = numeric(B)
  Vsotab = numeric(B)
  for (b in 1:B){ 
    Esotab[b] = EVsotab[[b]][1]
    Vsotab[b] = EVsotab[[b]][2]
  }
  
  
  Esota = mean(Esotab)
  Vsota = mean(Vsotab) + var(Esotab)
  
  Bias_theta_vec[j] = (1-Esota/n)-theta_SOTA
  SD_theta_vec[j] = sqrt(Vsota)/n
  
  print(c(j,length(d_vec)-j))

}

# saving
# naming convention: fig_name_vector_name
saveRDS(d_vec, file = "bias_sd_d_d_vec.rds")
saveRDS(Bias_theta_vec, file = "bias_sd_d_Bias_theta_vec.rds")
saveRDS(SD_theta_vec, file = "bias_sd_d_SD_theta_vec.rds")

cols = c("black","darkgreen")
plot(d_vec, Bias_theta_vec,"l", lty = "solid", col = cols[1], 
     ylim = ylm_bias, xlab = "", ylab = "")
par(new=TRUE)
plot(d_vec, SD_theta_vec,"l", lty = "solid", col = cols[2], 
     ylim = ylm_sd,  xlab = "", ylab = "", axes = FALSE)
axis(4)

abline(v=0.025, col="gray")
abline(v=0.05, col="gray",lty = 5)
abline(v=0.075, col="gray",lty = 4)
abline(v=0.1, col="gray",lty = 3)

title(main = "", xlab = TeX(r'($d = b-a$)'), ylab = "", line = 2, cex.lab=1.2)
legend(0.125, 0.005, legend=c(ylab_bias,ylab_sd), col=cols, lty=c(1,1), cex=0.8)
