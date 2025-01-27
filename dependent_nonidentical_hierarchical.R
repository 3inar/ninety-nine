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
library(future)     # parallelisation
library(furrr)


###############################################################################
############################ Dependent non-identical #############################
###############################################################################

# Consider a classification problem with a test set of size $3,000$, and a 
# classifier with probability of correct classification is $theta$.
# The estimated accuracy, that is, the classifier's performance on the test set, 
# is denoted by $theta_hat$.

future::plan("multisession", workers=8) # let's engage 6 out of 8 cores

# Functions

# dep_nonid_pmf - simulated pmf

source("Parameters_PublicCompetition.R") # n, theta, m, alpha, rho, theta_min, theta_max, theta_vec

source("dep_nonid_pmf_fun.R") # for function dep_nonid_pmf
# returns X = list(min_fail, x_fail, teamsSOTA)

source("ProbDistr_thetaSOTA.R") # for function sim_ci()

# theta_vec is sampled from uni(theta_min,theta_max)
theta_max = theta_SOTA + d/(m+1)
theta_min = theta_max - d

tic()
EVsotab <- furrr::future_map(1:B, function (x) { 
  
  theta_vec = runif(m, min=theta_min, max=theta_max)
  X = dep_nonid_pmf(n, m, rho, rep, theta_vec, theta_0 = theta_SOTA) 
  
  Esota = mean(X$min_fail)
  Vsota = mean(X$min_fail*X$min_fail) - Esota*Esota
  min_dep_alpha2 = sim_ci(alpha, X$min_fail) # The upper bound of the 95% confidence interval
  
  return(cbind(Esota, Vsota, min_dep_alpha2))
  
}, .progress=T, .options= furrr::furrr_options(seed=T))
toc() # 4228 sec

# Not very elegant with the for-loop, but that's ok for now
Esotab = numeric(B)
Vsotab = numeric(B)
alpha2b = numeric(B)
for (b in 1:B){ 
  Esotab[b] = EVsotab[[b]][1]
  Vsotab[b] = EVsotab[[b]][2]
  alpha2b[b] = EVsotab[[b]][3]
}

Esota = mean(Esotab)
Vsota = mean(Vsotab) + var(Esotab)
min_dep_alpha2 = mean(alpha2b) # not sure about this one

sprintf("The simulated dependent upper bound of the %s confidence interval is %.5f, with %s repetitions.",  
        1-alpha, (n-min_dep_alpha2)/n, rep)
#"The simulated dependent upper bound of the 0.95 confidence interval is 0.91733, with 1e+05 repetitions."

sprintf("The simulated expected value is %.7f and a standard deviation is %.7f, with %s repetitions.",  
        (n-Esota)/n, sqrt(Vsota)/n, rep)
# "The simulated expected value is 0.9101173 and a standard deviation is 0.0036547, with 1e+05 repetitions."

###########################################################################
########################### Figures ##################### ################
###########################################################################

########################## bias_thetamin_rho and sd_thetamin_rho ############

# theta_min can be only so small, see Eq.`rho_min`
trunc_min = ((rho*rho)*theta_0/(1-theta_0))/(1+(rho*rho)*theta_0/(1-theta_0))

rho_step = 0.01 # adjust for smoothness, x-axis
rho_vec = seq(0.0, 1.0, by=rho_step) 

d_vec = seq(0,0.1,by=0.025) # five curves

#pre-allocate
Esota_theta_vec = matrix(NA,length(d_vec),length(rho_vec))
SDsota_theta_vec = matrix(NA,length(d_vec),length(rho_vec))

# rep = 200000 # increase for smoothness for the standard deviation

tic()
for (i in 1:length(d_vec)){ 
  # create theta-vector
  theta_max = theta_SOTA + d_vec[i]/(m+1)
  theta_min = theta_max - d_vec[i]
  
  # truncate rho
  rho_trunc = sqrt((theta_min*(1-theta_0))/(theta_0*(1-theta_min)))
  rho_vec = seq(0.0, rho_trunc, by=rho_step) 
  
  for (j in 1:length(rho_vec)){ 
    
    EVsotab <- furrr::future_map(1:B, function (x) { 
      theta_vec = runif(m, min=theta_min, max=theta_max)
    
      X = dep_nonid_pmf(n, m, rho_vec[j], rep, theta_vec, theta_0 = theta_SOTA)
    
      # Expected value and variance
      Esota = mean(X$min_fail)
      Vsota = mean(X$min_fail*X$min_fail) - Esota*Esota
      
      return(cbind(Esota, Vsota))
      
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
    
    Esota_theta_vec[i,j] = (1-Esota/n)-theta_SOTA
    SDsota_theta_vec[i,j] = sqrt(Vsota)/n
      
    print(c(i,j,length(d_vec),length(rho_vec)))
  }
}
toc()

lgnds = c(TeX(r'($d=0$)'),TeX(r'($d=0.025$)'), TeX(r'($d=0.050$)'), 
          TeX(r'($d=0.075$)'), TeX(r'($d=0.100$)'))

############# bias_thetamin_rho ##############
rho_vec = seq(0.0, 1.0, by=rho_step) 

# plotting for the d's
plot(rho_vec, Esota_theta_vec[1,],"l", lty = 1, col = "red", ylim = ylm_bias,
     main = '', xlab = '', ylab = '')
par(new=TRUE) 
plot(rho_vec, Esota_theta_vec[2,],"l", lty = 1, col = "black", ylim = ylm_bias,
     xlab = "", ylab = "")
par(new=TRUE) 
plot(rho_vec, Esota_theta_vec[3,],"l", lty = 5, col = "black", ylim = ylm_bias,
     xlab = "", ylab = "")
par(new=TRUE) 
plot(rho_vec, Esota_theta_vec[4,],"l", lty = 4, col = "black", ylim = ylm_bias,
     xlab = "", ylab = "")
par(new=TRUE) 
plot(rho_vec, Esota_theta_vec[5,],"l", lty = 3, col = "black", ylim = ylm_bias,
     xlab = "", ylab = "")

# intersections correspond to table `noniid` 
abline(v=rho, col="gray")

title(ylab = ylab_bias, line=2, cex.lab=1.2, xlab = TeX(r'(${rho}_0$)'))

legend(0.55, 0.035, legend=lgnds, col=c("red","black","black","black","black"),
       lty=c(1,1,5,4,3), cex=0.8)

# saving
# naming convention: fig name_vector_name
saveRDS(rho_vec, file = "bias_thetamin_rho_rho_vec.rds")
saveRDS(d_vec, file = "bias_thetamin_rho_d_vec.rds")
saveRDS(Esota_theta_vec, file = "bias_thetamin_rho_Esota_theta_vec.rds")
saveRDS(SDsota_theta_vec, file = "sd_thetamin_rho_SDsota_theta_vec.rds")

################################### sd_thetamin_rho ############################
# repeat all, only for standard deviation

plot(rho_vec, SDsota_theta_vec[1,],"l", lty = 1, col = "red", ylim = ylm_sd,
     main = '', xlab = '', ylab = '')
par(new=TRUE)
plot(rho_vec, SDsota_theta_vec[2,],"l", lty = 1, col = "black", ylim = ylm_sd,
     xlab = "", ylab = "")
par(new=TRUE) 
plot(rho_vec, SDsota_theta_vec[3,],"l", lty = 5, col = "black", ylim = ylm_sd,
     xlab = "", ylab = "")
par(new=TRUE) 
plot(rho_vec, SDsota_theta_vec[4,],"l", lty = 4, col = "black", ylim = ylm_sd,
     xlab = "", ylab = "")
par(new=TRUE) 
plot(rho_vec, SDsota_theta_vec[5,],"l", lty = 3, col = "black", ylim = ylm_sd,
     xlab = "", ylab = "")

abline(v=rho, col="gray")

title(ylab = ylab_sd, line=2, cex.lab=1.2, xlab = TeX(r'(${rho}_0$)'))

legend(0.55, 0.0026, legend=lgnds, col=c("red","black","black","black","black"),
       lty=c(1,1,5,4,3), cex=0.8)



