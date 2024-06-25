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

xl=c(.4, 1)
png("auc_truncated_boot.png", width=6, height=4, units="in", res=300)
hist(boots, nclass=200, prob=T, xlab="bootstrapped expected AUC",
     xlim=xl,
     main="AUCs bootstrap-sampled from truncated ( < .8975)\nkaggle results")
dev.off()

png("auc_truncated_simulation.png", width=6, height=4, units="in", res=300)
hist(realization, nclass=200, prob=T, xlab="simulated observed AUC",
     xlim=xl,
     main="AUCs simulated from bootstrapped distribution")
dev.off()

# kajsa figures go to .96
