# ROC curves from AUCs
# by Kajsa Møllersen (kajsa.mollersen@uit.no) February 2023

library(dplyr)
library(magrittr)
library(rvest)
library(purrr)
library(readr)
library(ggplot2)

library(plotly)


# All the confidence intervals ############################################################

library(pROC) # to find the conf int of auc

# These numbers are from https://www.kaggle.com/c/siim-isic-melanoma-classification/data?select=train.csv 
n_val_test = 10982 # size of test and validation set
test_prop = 0.7 # 30/70 split
n = round(test_prop*n_val_test)

malignant_rate = 584/33126 # in training set ( from https://arxiv.org/ftp/arxiv/papers/2008/2008.07360.pdf)
teams = 3308
AUC_top = 0.9490 # best performing team

n_mal = round(n*malignant_rate) # estimated number of malignant cases
n_ben = n-n_mal # estimated number of benign cases

sprintf("The test set size is %s, with %s malignant cases",
        n, n_mal)

class = rep(0, n)
class[1:n_mal] = 1 # true class label

# Create toy example with upper bound CI AUC = 0.9490

false1 = 0.376 # this parameter is adjusted until required AUC is produced
# 0.376 gives AUC=0.9308, and produces 95% CI upper bound = AUC_top=0.9490
predict1 = class

# malignant prediction
n_mal_05 = round(false1*n_mal)
n_mal_pred = seq(from = (1/n_mal_05), to = 1, by = (1/n_mal_05))
predict1[1:n_mal_05] = n_mal_pred

# benign prediction
n_ben_05 = round(false1*n_ben)
n_ben_pred = seq(from = (1/n_ben_05), to = 1, by = (1/n_ben_05))
predict1[(n_mal+1):(n_mal+n_ben_05)] = n_ben_pred

roc1 = roc(class,predict1)

auc_sota = auc(class,predict1)
sprintf("The estimated SOTA is %.4f",
        auc_sota)
auc_round = round(auc_sota*10000)/10000
plot(roc1, main = paste('Area under curve', auc_round))
# ci.auc(class, predict1) # deLong doesn't work that well

CI = ci.auc(class, predict1, method = "bootstrap", boot.n=50000, conf.level=0.95)

print(CI)

# kaggle data

comb_data <- readRDS('/Users/kajsam/Documents/kaggle-leaderboard-scrape/SIIM-ISIC_Melanoma_kaggle_leadboard_data.RDS')
head(comb_data) # have a look

dat <- data.frame(AUC = comb_data$prv_score, dataset = "test")

# number of teams performing above AUC_SOTA (and CI)

AUC = sort(dat$AUC) # all observed performances, sorted (why not)
m = length(AUC)
maxAUC = max(AUC) # best observed perfomance

m_above= length(AUC[AUC>CI[[1]]])
m_above = length(AUC[AUC>CI[[3]]])

m_above = length(AUC[AUC>auc_sota])
sprintf("Number of teams performing above the estimated SOTA is %s",
        m_above)
