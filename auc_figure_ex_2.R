source("plotting_params.R")
source("auc_functions.R")
source("Parameters_PublicCompetition.R")

load("exp_auc_09491.rda")
load("melanoma_private.rda")

num_classifiers <- 1000
test_size <- 3000
ntr <- floor(test_size*malignant_rate)
nfa <- test_size - ntr

threshold <- .8975
truncated <- melanoma_scores[melanoma_scores < threshold]

set.seed(2024-05-15)
boots <- sample(truncated, length(melanoma_scores), replace=T)
realization <- sim_competition(boots, ntr, nfa)


breaks <- hist(melanoma_scores, nclass=200)$breaks
colors <- rep("black", length(breaks))
colors[breaks > threshold] = "darkgrey"
xl=c(.4, 1)

new_png("auc_truncated_boot.png", n_figures=2)
hist(melanoma_scores, nclass=200, prob=T, xlab="AUC",
     xlim=xl, col=colors, border=colors, main=NULL)
dev.off()

new_png("auc_truncated_simulation.png", n_figures=2)
hist(realization, nclass=200, prob=T, xlab="AUC",
     xlim=xl, col="black", border="black", main=NULL)
dev.off()

# kajsa figures go to .96
