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

png("auc_figure_1.png", width=6, height=4, units="in", res=300)
hist(mxx, nclass=200, prob=T,
     main="Simulation distribution of apparent SOTA AUC,\nall models have true AUC of .9",
     xlab="Observed largest AUC among 1000 models.")
dev.off()


