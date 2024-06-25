# Idea: AUC is the probabilty of ranking a random T higher than a random F
# (Mann-Whitney-Wilcoxon U statistic), hence if the scores for class x are
# distributed as X and the scores for class y distributed as Y we have AUC=0.95
# is equivalent to P(0 > Y - X) = 0.95

classifier <- setClass("classifier", slots = c(mu = "numeric", sd = "numeric"))

# Looking at Gaussians for simplicity, taking mu_x, sd_x, sd_y at random and
# solving for mu_y so that a given AUC obtains
make_classifier <- function(auc) {
  if (auc < 0 || auc > 1) stop("AUCs should be between 0 and 1")
  x <- 1  # counted as the Trues 
  y <- 2
  mu <- c(0,0)
  sd <- c(1, 1)
  sd_z <- sqrt(sd[x]^2 + sd[y]^2)
  
  # pnorm is the CDF of a Gaussian rv., we optimize wrt mu_y so that the above
  # probability is as close as possible to .95 or whatever
  f <- function(mu_y) pnorm(0, (mu_y - mu[x]), sd_z)
  
  soln = uniroot(\(mm) {f(mm) - auc}, interval=mu[x] + 5*c(-sd_z, sd_z))
  mu[y] = soln$root
  
  return(classifier(mu=mu, sd=sd))
}


# Generate predictions from a classifier
predict.classifier <- function(obj, true_size=50, false_size=100) {
  return(list(
    truth = c(rep(1, true_size), rep(0, false_size)),
    predicted = c(rnorm(true_size, obj@mu[1], obj@sd[1]),
                  rnorm(false_size, obj@mu[2], obj@sd[2])))
  )
}

# These functions extract the means and standard deviations for the two
# prediction distributions score|true, score|false
sds <- function(obj) {
  slot(obj, "sd")
}

mus <- function(obj) {
  slot(obj, "mu")
}

# Given reference values x1...xn, and parameters for a bivariate normal for 
# (x, y), generates predictions from y|xi
draw_correlated <- function(reference_values, mu1, mu2, sd1, sd2, corl) {
  # https://www2.stat.duke.edu/courses/Spring12/sta104.1/Lectures/Lec22.pdf
  Z1 <- reference_values
  Z1 <- (Z1 - mu1)/sd1
  Z2 <- rnorm(length(Z1))

  conditionals <- sd2*(corl*Z1 + sqrt(1 - corl^2)*Z2) + mu2

  return(conditionals)
}

# Given some prediction from a "leader" classifier, generate orrelated
# predictions for a "follower" classifier. The requested correlation will be
# fulfilled for the conditional distributions score|true, score|false, but for
# the overall distrbution of scores it will be somewhat off.
correlated_predict <- function(follower, leader, leading_predictions, corl) {
  predictions <- leading_predictions$predicted
  truth <- leading_predictions$truth

  mu_l <- mus(leader)
  mu_f <- mus(follower)

  sd_l <- sds(leader)
  sd_f <- sds(follower)

  # trues
  pred_true <- draw_correlated(predictions[truth==1], 
                               mu_l[1], mu_f[1], 
                               sd_l[1], sd_f[1], corl)


  # falses
  pred_true <- draw_correlated(predictions[truth==0], 
                               mu_l[2], mu_f[2], 
                               sd_l[2], sd_f[2], corl)

  pred = numeric(length(predictions))
  pred[truth == 1] <- pred_true
  pred[truth == 0] <- pred_false

  return(list(truth = truth, predicted = pred))
}

# empirical_auc_slow <- function(predictions) {
#   wc = wilcox.test(predictions$predicted[predictions$truth==1], 
#                    predictions$predicted[predictions$truth==0])
#   wc$statistic/(sum(predictions$truth==1)*sum(predictions$truth == 0))
# }
# 
# # should be faster than wilcox.test or pROC https://blog.mbq.me/augh-roc/
# empirical_auc_newer <-function(predictions){
#   cls <- predictions$truth == 1
#   score <- predictions$predicted
# 
#   n1 <-sum(!cls); sum(cls)->n2;
#   U <-sum(rank(score)[!cls])-n1*(n1+1)/2;
# 
#   return(1-U/n1/n2);
# }
# 
# empirical_auc_newer_still <- function(predictions) {
#   bigstatsr::AUC(predictions$predicted, predictions$truth)
# }

empirical_auc<- function(predictions) {
  # need that third : because AUC2 is "not exported"
  bigstatsr:::AUC2(predictions$predicted, as.logical(predictions$truth))
}

# simulates a competition result based on "true"/expected aucs 
sim_competition <- function(aucs, n_true, n_false) {
  plyr::aaply(aucs, 1, \(x) { 
    classifier <- make_classifier(x)
    predictions <- predict(classifier, n_true, n_false)
    empirical_auc(predictions)
  })
}
