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

set.seed(2024-05-24)
png("auc_figure_rocB.png", width=6, height=4, units="in", res=300)
curve(dnorm(x, ms[1], ss[1]), xlim=c(-5+ms[2], 5), col="black", lwd=2,
        main="Density of S- and S+ for AUC = .9",
        xlab="Score",
        ylab="Density")
curve(dnorm(x, ms[2], ss[2]), add=T, col="grey", lwd=2)

rg_f <- rnorm(nfa, ms[2], ss[2])
rg_t <- rnorm(ntr, ms[1], ss[2])
rug(rg_f, col="grey", lwd=1.5)
rug(rg_t, col="black", lwd=1.5)
dev.off()

set.seed(2024-10-24)
png("auc_figure_rocA.png", width=6, height=4, units="in", res=300)
  # simulate 1000 competitors
  lst <- list()
  for (i in 1:1000) {
          preds <- predict(clasr, ntr, nfa)
          lst[[length(lst)+1]] <- roc(response=preds$truth, predictor=preds$predicted)
  }

  # what are their observed aucs? what are the inner quartiles and max of these>
  acs <- sapply(lst, auc)
  inner_quart <- quantile(acs, probs=c(1/4, 3/4))
  mxx <- which.max(acs)

  # pull out the maximum curve
  maxline <- lst[[mxx]]
  lst[[mxx]] <- NULL
  acs <- acs[-mxx]

  # which curves belong to the inner quartile region?
  inner <- acs > inner_quart[1] & acs < inner_quart[2]

  # plot the middle curves in dark grey, the outer ones in light grey, the max auc one in black
  plot(lst[[1]], type="n")
  for (roccc in lst[!inner]) {
    plot(roccc, lwd=1, add=T, col="#dddddd")
  }
  for (roccc in lst[inner]) {
    plot(roccc, lwd=1, add=T, col="#999999")
  }

  plot(maxline, lwd=2, add=T, col="black")
dev.off()
