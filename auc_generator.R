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

betas = auc_generator(.9)
calculate_auc(betas$negatives, betas$positives)

curve(dbeta(x, betas$negatives$alpha, betas$negatives$beta), col="grey")
curve(dbeta(x, betas$positives$alpha, betas$positives$beta), add=T)

# generate scores according to our chosen distributions
zeros = rbeta(7000, betas$negatives$alpha, betas$negatives$beta)
ones = rbeta(3000, betas$positives$alpha, betas$positives$beta)

# ground truth
true = c(rep(0, 7000), rep(1, 3000))
pred = c(zeros, ones)


roc1 = roc(true, pred)
plot(roc1)
auc(roc1)

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

# some random experiments suggest that the dichotome curve stems from a kind of
# strange situation where the scores are very confidently right and also very
# confidently wrong for the negative class specifically