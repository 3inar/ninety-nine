source("plotting_params.R")
load("auc_exp_1.rda")

get_aucs <- function(obj) {
  plyr::laply(obj, \(x) { x$observed_auc })
}

get_predictions <- function(obj) {
  plyr::laply(obj, \(x) { x$predictions$predicted })
}

mxx <- plyr::laply(first_experiment, \(x) max(get_aucs(x)))
aucc <- plyr::laply(first_experiment, \(x) get_aucs(x))

mean(mxx); sd(mxx)
quantile(mxx, prob=c(.025, .975))

mean(aucc); sd(aucc)
quantile(aucc, prob=c(.025, .975))

new_png("auc_figure_1.png", n_figures=1)
hist(mxx, nclass=200, prob=T,
     main=NULL, col="black",
     xlab="Sample maximum AUC among 1000 models.")
dev.off()


