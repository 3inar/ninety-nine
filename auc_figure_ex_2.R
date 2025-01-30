source("plotting_params.R")
source("auc_functions.R")
source("Parameters_PublicCompetition.R")

load("exp_auc_09491.rda")
load("melanoma_private.rda")

num_classifiers <- 1000
test_size <- 3000
ntr <- floor(test_size*malignant_rate)
nfa <- test_size - ntr

# file mxx loaded with the file exp_auc... contains sample maxima from a
# situation engineerd to have mean .949
max_mean <- mean(mxx)
# [1] 0.9490963
ci_max <- quantile(mxx, prob=c(.05/2, 1- .05/2))
#      2.5%     97.5% 
# 0.9409052 0.9599991

ci_single <- auc_ci(max(melanoma_scores), conf_level=.95, true_size=ntr,
                    false_size=nfa, n_sim=10000)

weight <- 0.8890625
truncated <- weight*melanoma_scores + (1-weight)*0.5

set.seed(2024-05-15)
boots <- sample(truncated, length(truncated), replace=T)
realization <- sim_competition(boots, ntr, nfa)


breaks <- hist(melanoma_scores, nclass=200)$breaks
xl=c(.6, 1)
yl=c(-1,31)

# new_png("auc_truncated_boot.png", n_figures=2)
# hist(melanoma_scores, nclass=200, prob=T, xlab="AUC",
#      xlim=xl, ylim=yl, col=colors, border=colors, main=NULL)
# #abline(v=threshold, col="green", lty="dashed")
# lines(ci_single, rep(yl[1], 2), col="red", lwd=2)
# abline(v=max(melanoma_scores), col="red", lty="dashed")
# dev.off()
# 
#yl=c(-0.25, 8.75)

new_png("auc_shrink_boot.png", n_figures=2)
hist(melanoma_scores, nclass=200, prob=T, xlab="AUC",
     xlim=xl, ylim=yl, col="lightgrey", border="lightgrey", main=NULL)
hist(truncated, nclass=200, prob=T, xlab="AUC",
     xlim=xl, ylim=yl, col="black", border="black", main=NULL, add=T)
#abline(v=threshold, col="green", lty="dashed")
lines(ci_single, rep(yl[1], 2), col="red", lwd=2)
abline(v=max(melanoma_scores), col="red", lty="dashed")
dev.off()

max(truncated)
# [1] 0.8992448
sum(melanoma_scores > max(truncated))
# [1] 2013

new_png("auc_truncated_simulation.png", n_figures=2)
hist(realization, nclass=200, prob=T, 
     xlab=latex2exp::TeX(r'(\widehat{AUC})'), xlim=xl, ylim=yl,
     col="black", border="black", main=NULL)
lines(ci_max, rep(yl[1], 2), col="red", lwd=2)
abline(v=mean(mxx), col="red", lty="dashed")
dev.off()

# kajsa figures go to .96
