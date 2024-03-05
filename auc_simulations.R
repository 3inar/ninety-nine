library(pROC)

source("auc_functions.R")

aucs = seq(from=.89,to=.915, by=.0001)

# want CI to have upper limit .9490
experiment = plyr::aaply(aucs, 1, function(x) {
  xx = plyr::raply(10000, function() {
    cls_95 <- make_classifier(x)
    predicted <- predict(cls_95, true_size=50, false_size=2950)
    # I just use the default DeLong interval
    lim_95 = suppressMessages(ci.auc(predicted$truth, predicted$predicted,
                                 conf.level=.95)[3])
    lim_99 = suppressMessages(ci.auc(predicted$truth, predicted$predicted,
                                 conf.level=.99)[3])
    c(lim_95,lim_99)
  })
  
  colMeans(xx)
}, .progress = "text")

plot(aucs, experiment[,1], type="l")
lines(aucs, experiment[,2], type="l")
abline(h=.95, col="grey")

aucs[which.min(abs(experiment[,1] - 0.9490))]
aucs[which.min(abs(experiment[,2] - 0.9490))]

auc(predicted$truth, predicted$predicted)
plot(roc(predicted$truth, predicted$predicted))
boxplot(predicted$predicted, predicted$truth)

################################################################################
################################################################################
########### standard deviations for AUCs with values 0.9490, 0.9378 og 0.9295
################################################################################
################################################################################

# from kajsa.R file; visual intelligence seminar. Should perhaps just go with
# the ParametersPublicCompetition file

n_pos = round(n*malignant_rate)
n_neg = n - n_pos

sampling_auc <- function(true_auc) {
  cls <- make_classifier(true_auc)
  predicted <- predict(cls, true_size=n_pos, false_size=n_neg)
  
  empirical_auc(predicted)
}

sd(replicate(10000, sampling_auc(.9490)))
sd(replicate(10000, sampling_auc(.9378)))
sd(replicate(10000, sampling_auc(.9295)))


## Imports the following parameters (among others)
# n - size of test set
# m - number of classifiers
# alpha - significance level
# rho - correlation coefficient
# malignant_rate # proportion of positives in melanoma training set 
# rep - number of repetitions
# NB: TODO: the above overwrites the internal function rep and should be changed probably
source("Parameters_PublicCompetition.R")

# private leaderboard scores of the melanoma challenge
load("melanoma_private.rda")
melanoma_scores <- melanoma_private_scores[melanoma_private_scores > .5]

winner_auc = max(melanoma_scores)

m_vals <- c(1000, 1000, 1000, 100, 500, 1000, 1000)
n_vals <- c(3000, 3000, 3000, 3000, 3000, 1000, 10000)
t_vals <- c(.85, .90, .95, .90, .90, .90, .90)
params <- cbind(m_vals, n_vals, t_vals)


sim_competition <- function(aucs, n_true, n_false) {
  # this would look more beautiful as a piped function
  plyr::alply(aucs, 1, \(x) { 
    classifier <- make_classifier(x)
    predictions <- predict(classifier, n_true, n_false)
    oauc <- empirical_auc(predictions)
    list(classifier=classifier, predictions=predictions, 
         expected_auc=x, observed_auc=oauc)})
}

get_aucs <- function(obj) {
  plyr::laply(obj, \(x) { x$observed_auc })
}

get_predictions <- function(obj) {
  plyr::laply(obj, \(x) { x$predictions$predicted })
}

profvis({
  smallsim <- sim_competition(rep(.9, 10000), 900, 2000)
})


# row 2 of table 1 (used in fig3) recreated
ex_auc <- .9
num_classifiers <- 1000
test_size <- 3000
ntr <- floor(test_size*malignant_rate)
nfa <- test_size - ntr
system.time({
  auc_sims <- plyr::rlply(200, function () { 
                sim_competition(rep(ex_auc, num_classifiers), ntr, nfa)
              }, .progress="text")
})

mxx <- plyr::laply(auc_sims, \(x) max(get_aucs(x)))
# head(plyr::laply(auc_sims, \(x) { get_predictions(x) }))

aucc <- plyr::laply(auc_sims, \(x) get_aucs(x))

perm_predicts <- function(predicts) {
  predicts$truth <- sample(predicts$truth)
  predicts
}

sample_preds <- auc_sims[[1]][[1]]$predictions
permutation_aucs <- plyr::raply(100000, empirical_auc(perm_predicts(sample_preds)))
hist(permutation_aucs[permutation_aucs > 0.55], nclass=500)
abline(v=(1:52)/52)

hist(aucc[aucc > .90], nclass=500)


# The reference classifier is quite good but it shouldn't actually matter much
reference_classifier <- make_classifier(.99)
baby_classifier <- make_classifier(.95)

reference_predicts <- predict(reference_classifier)
corr_predictions <- correlated_predict(baby_classifier, reference_classifier,
                                       reference_predicts, 0.6)


corr_predictions

plot(reference_predicts$predicted, corr_predictions$predicted, pch=(1 + corr_predictions$truth))
empirical_auc(reference_predicts)

profvis({
  replicate(1000, empirical_auc(corr_predictions))
})



## this all seems to work
# TODO: 
# * Figure out what was Kajsas experiment
