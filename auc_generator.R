# https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
beta_params <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  
  return(list(alpha = alpha, beta = beta))
}

muu = runif(1)
varr = runif(1, max=muu*(1-muu))

pars = beta_params(muu, varr)

curve(dbeta(x, pars$alpha, pars$beta))

muu2 = runif(1, min=muu)
varr2 = runif(1, max=muu2*(1-muu2))

pars2 = beta_params(muu2, varr2)

#px<y = int p(x<y|y)p(y) dy
p_x_lt_y = function(y, ax, bx, ay, by) {
  pbeta(y, ax, bx)*dbeta(y, ay, by)
}

integrate(p_x_lt_y, 0, 1, pars$alpha, pars$beta, 
                          pars2$alpha, pars2$beta)

curve(dbeta(x, pars$alpha, pars$beta), col="grey")
curve(dbeta(x, pars2$alpha, pars2$beta), add=T)

zeros = rbeta(7000, pars$alpha, pars$beta)
ones = rbeta(3000, pars2$alpha, pars2$beta)

true = c(rep(0, 7000), rep(1, 3000))
pred = c(zeros, ones)

library(pROC)

roc1 = roc(true, pred)
plot(roc1)
auc(roc1)

# Kajsa's approach:
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
