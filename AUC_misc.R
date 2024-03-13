# just need to put this somewhere it won't clutter up main simulations
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
# from params file
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
