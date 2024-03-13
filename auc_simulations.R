## Imports the following parameters (among others)
# n - size of test set
# m - number of classifiers
# alpha - significance level
# rho - correlation coefficient
# malignant_rate # proportion of positives in melanoma training set 
# rep - number of repetitions
# NB: TODO: the above has same name as the internal function rep and should be changed probably
source("Parameters_PublicCompetition.R")

## IMports functions to simulate classifiers w/ given AUCs
source("auc_functions.R")

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

    # used to return the original classifier and its predictions but it takes a
    # lot of memory space 
    list(classifier=NA, predictions=NA, 
         expected_auc=x, observed_auc=oauc)})
}

get_aucs <- function(obj) {
  plyr::laply(obj, \(x) { x$observed_auc })
}

get_predictions <- function(obj) {
  plyr::laply(obj, \(x) { x$predictions$predicted })
}

# row 2 of table 1 (used in fig3) recreated
ex_auc <- .9
num_classifiers <- 1000
test_size <- 3000
ntr <- floor(test_size*malignant_rate)
nfa <- test_size - ntr

## TODO:: future package + furrr to make it parallel
## I have 8 logical AND 8 physical cores
future::plan("multisession", workers=4)

system.time({
  auc_sims <- furrr::future_map(1:10000, function (x) { 
                source("auc_functions.R") # this might be slow but I don't think I care
                sim_competition(rep(ex_auc, num_classifiers), ntr, nfa)
              }, .progress=T, .options= furrr::furrr_options(seed=T))
})


mxx <- plyr::laply(auc_sims, \(x) max(get_aucs(x)))
aucc <- plyr::laply(auc_sims, \(x) get_aucs(x))

hist(aucc, nclass=100); abline(v=mean(aucc), lwd=3)
hist(mxx, nclass=100)

perm_predicts <- function(predicts) {
  predicts$truth <- sample(predicts$truth)
  predicts
}


# The reference classifier is quite good but it shouldn't actually matter much
reference_classifier <- make_classifier(.99)
baby_classifier <- make_classifier(.95)

reference_predicts <- predict(reference_classifier)
corr_predictions <- correlated_predict(baby_classifier, reference_classifier,
                                       reference_predicts, 0.6)

plot(reference_predicts$predicted, corr_predictions$predicted, pch=(1 + corr_predictions$truth))
empirical_auc(reference_predicts)


