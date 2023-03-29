# Kommentar øverst i fila

library(pROC) # I use this below

# This function converts a mean and variance for a beta distribution to its 
# canonical parameters, alpha and beta. Found at:
# https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  
  return(list(alpha = alpha, beta = beta))
}

# If the score for a randomly chosen negative is represented by random variable X
# Y is that of a randomly chosen positive, the AUC is 
#         p(X < Y) = \int F_X(y)f_Y(y) dy.
# The below function calculates this integral numerically for two beta
# distributions
calculate_auc = function(pars, pars2) {
  integration = integrate(lower=0, upper=1,
    f = \(y) { pbeta(y, pars$alpha, pars$beta)*dbeta(y, pars2$alpha, pars2$beta) })
  integration$value
}

# A function that generates a pair of beta distributions with an AUC very close
# to some desired value
auc_generator = function(desired_auc, tolerance=0.0001, rands=5000) {
  
  # Generate many random betas 
  muu2_grid = runif(rands)
  random_betas = cbind(muu2_grid,  runif(rands, max=muu2_grid*(1-muu2_grid)))
  param_list = plyr::alply(random_betas, 1, 
                           function(prm) { beta_params(prm[1], prm[2]) } )
  
  # select some random pairs
  while (T) {
    pars = param_list[sample(rands, 2)]
    
    ret = tryCatch(calculate_auc(pars[[1]], pars[[2]]), error = function(e) -1)
    if (abs(ret -desired_auc) < tolerance) break 
  }
  
  # get AUCs comparing the first beta to big list of random ones
#   auc_list = plyr::laply(param_list, 
#     function(pars2) { 
#       ret = tryCatch(calculate_auc(pars, pars2), error = function(e) NA)
#       ret
#     })
    
  #pars2 = param_list[[which.min(abs(auc_list - .95))]]
  
  return(list(negatives=pars[[1]], positives=pars[[2]]))
}

plot_betas = function(betas) {
  curve(dbeta(x, betas$negatives$alpha, betas$negatives$beta), col="grey")
  curve(dbeta(x, betas$positives$alpha, betas$positives$beta), add=T)
}

means_bs = function(betas) {
  n_m = betas$negatives$alpa/(betas$negatives$alpa + betas$negatives$beta)
  p_m = betas$positives$alpa/(betas$positives$alpa + betas$positives$beta)
  c(negatives = n_m, positives=p_m)
}

vars_bs = function(betas) {
  b = betas$negatives$beta
  a = betas$negatives$alpha
  n_v = (b*a)/(((a + b)^2) + (a + b + 1))
  
  b = betas$positives$beta
  a = betas$positives$alpha
  p_v = (b*a)/(((a + b)^2) + (a + b + 1))
  c(negatives = n_v, positives=p_v)
}

roc_betas = function(betas) {
  zeros = rbeta(7350, betas$negatives$alpha, betas$negatives$beta)
  ones = rbeta(150, betas$positives$alpha, betas$positives$beta)
  
  true = c(rep(0, 7350), rep(1, 150))
  pred = c(zeros, ones)
  
  roc(true, pred, quiet=T)
}


# # find various examples of .95 AUC scorers and their variances
# curves = list(); N = 500
# for (i in 1:N) {
#   betas = auc_generator(.95)
#   
#   aucs = rep(-1, 250)
#   for (j in 1:250) {
#     # generate scores according to our chosen distributions
#     zeros = rbeta(7350, betas$negatives$alpha, betas$negatives$beta)
#     ones = rbeta(150, betas$positives$alpha, betas$positives$beta)
#     
#     # ground truth
#     true = c(rep(0, 7350), rep(1, 150))
#     pred = c(zeros, ones)
#     
#     
#     roc1 = roc(true, pred, quiet=T)
#     aucs[j] = auc(roc1)
#   }
#   
#   res = list(betas=betas, aucs=aucs)
#   curves[[i]] <- res
#   
#   cat(as.character(i))
#   cat("/")
#   cat(as.character(N))
#   cat("\r")
# }; cat("\n")
# 
# save(curves, file = "auc_gallery.rda")
load("auc_gallery.rda")

sdevs = plyr::laply(curves, function(x) {sd(x$aucs)} )
means = plyr::laply(curves, function(x) {mean(x$aucs)} )
curves = curves[order(sdevs)]

plot_roc_b = function(betas, N=100, mn="You forgot main") {
  roc_high = roc_betas(betas); 
  ci = ci.auc(roc_high, method="boot")
  mt = paste0(mn, " (95% CI: ", signif(ci[1], 3), "—", signif(ci[3], 3),")")
  plot(roc_high, lwd = 1, col="grey", main=mt)
  
  for (i in 1:N) {
    roc_high = roc_betas(betas); 
    plot(roc_high, lwd = 1, col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5), add=T)
  }
  
}

plot_samples = function(rocc, main="You forgot main") {
  plot(c(rocc$cases, rocc$controls), 
       c(rep(0, length(rocc$cases)), rep(.8, length(rocc$controls))) + rnorm(7500, sd=.1), 
       pch=20, cex=.5, col=c(rep("black", length(rocc$cases)), rep("grey", length(rocc$controls))), 
       yaxt = "n", ylab="", xlab= "Score", sub="cases in black", bty="n", main=main)
}

not_variable = curves[[6]]$betas
medium_variable = curves[[round(length(curves)/2) + 0]]$betas
highly_variable = curves[[length(curves) - 10]]$betas

oldp = par(mfrow=c(1,3))
plot_samples(roc_betas(not_variable), main="Low variance AUC")
plot_samples(roc_betas(medium_variable), main="Medium variance AUC")
plot_samples(roc_betas(highly_variable), main="High variance AUC")
par(oldp)

plot(sdevs, means)

# plot_betas(not_variable)

oldp = par(mfrow=c(1,3))
plot_roc_b(not_variable, 20, "AUC: .95")
plot_roc_b(medium_variable, 20, "AUC: .95")
plot_roc_b(highly_variable, 20, "AUC: .95")
par(oldp)



plot_betas(medium_variable)
plot_roc_b(medium_variable, 20)

plot_betas(highly_variable)
plot_roc_b(highly_variable, 20)


# negs - poss: we see that the variance for the negatives is much smaller
# compared to that of the positives for those score distributions that are less
# variable. It is the other way for the very variable ones: the positives have a
# tighter distribution than do the negatives. means what?
vars = plyr::laply(curves, function(b) {vars_bs(b$betas)})
plot(log(vars[,1]/vars[,2]))

# Kajsa's approach: flat distributions with a certain overlap. leads to a very
# distinct curve
zeros = runif(7000, 0, .6)
ones = runif(3000, .4, 1)

true = c(rep(0, 7000), rep(1, 3000))
pred = c(zeros, ones)

# can we make the "dichotome" curve?
# not like this:
zeros = rbinom(7000, 1, .1)
ones = rbinom(3000, 1, .7)

true = c(rep(0, 7000), rep(1, 3000))
pred = c(zeros, ones)

rocc = roc(true, pred)
plot(rocc)

# some random experiments suggest that the dichotome curve stems from a kind of
# strange situation where the scores are very confidently right and also very
# confidently wrong for the negative class specifically
