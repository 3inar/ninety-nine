library(pROC)
source("Parameters_PublicCompetition.R")
source("auc_functions.R")

# row 2 of table 1 (used in fig3) recreated
ex_auc <- .9
num_classifiers <- 1000
test_size <- 3000
ntr <- floor(test_size*malignant_rate)
nfa <- test_size - ntr

clasr <- make_classifier(ex_auc)

ms <- mus(clasr)
ss <- sds(clasr)

png("auc_figure_rocB.png", width=6, height=4, units="in", res=300)
curve(dnorm(x, ms[1], ss[1]), xlim=c(-5+ms[2], 5), col="black", lwd=2,
        main="Density of S- and S+ for AUC = .9",
        xlab="Score",
        ylab="Density")
curve(dnorm(x, ms[2], ss[2]), add=T, col="grey", lwd=2)

set.seed(2024-05-24)
rg_f <- rnorm(nfa, ms[2], ss[2])
rg_t <- rnorm(ntr, ms[1], ss[2])
rug(rg_f, col="grey", lwd=1.5)
rug(rg_t, col="black", lwd=1.5)
dev.off()

png("auc_figure_rocA.png", width=6, height=4, units="in", res=300)
preds <- predict(clasr, ntr, nfa)
roccc <- roc(response=preds$truth, predictor=preds$predicted)
plot(roccc, lwd=1, main="ROC curves for 100 simulated competitions, AUC = .9")

for (i in 1:99) {
        preds <- predict(clasr, ntr, nfa)
        roccc <- roc(response=preds$truth, predictor=preds$predicted)
        plot(roccc, lwd=1, add=T)
}
dev.off()
