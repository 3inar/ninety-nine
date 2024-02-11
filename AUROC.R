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
print(n)

malignant_rate = 584/33126 # in training set ( from https://arxiv.org/ftp/arxiv/papers/2008/2008.07360.pdf)
teams = 3308
AUC_1 = 0.9490 # best performing team

n_mal = round(n*malignant_rate) # estimated number of malignant cases
print(n_mal)
n_ben = n-n_mal # estimated number of benign cases

class = rep(0, n)
class[1:n_mal] = 1 # true class label

# Create toy example with AUC = 0.9490

false1 = 0.376 # this parameter is adjusted until required AUC is produced
predict1 = class


# param.  AUC.  upper 95% CI
# 0.2 0.9808
# 0.3 0.9558
# 0.35 0.9395
# 0.375 0.931
# 0.376 0.9308          0.949
# 0.378 0.9305
# 0.3785 0.9304
# 0.3786 0.9303
# 0.379 0.9289
# 0.38 0.9287
# 0.4 0.922






# malignant prediction
n_mal_05 = round(false1*n_mal)
n_mal_pred = seq(from = (1/n_mal_05), to = 1, by = (1/n_mal_05))
predict1[1:n_mal_05] = n_mal_pred

# benign prediction
n_ben_05 = round(false1*n_ben)
n_ben_pred = seq(from = (1/n_ben_05), to = 1, by = (1/n_ben_05))
predict1[(n_mal+1):(n_mal+n_ben_05)] = n_ben_pred

roc1 = roc(class,predict1)
plot(roc1)
auc_sota = auc(class,predict1)
print(auc_sota)
# ci.auc(class, predict1) # deLong doesn't work that well
CI = ci.auc(class, predict1, method = "bootstrap", boot.n=50000, conf.level=0.95)

print(CI)

print(CI[3]-CI[1])

# kaggle data

comb_data <- readRDS('/Users/kajsam/Documents/kaggle-leaderboard-scrape/SIIM-ISIC_Melanoma_kaggle_leadboard_data.RDS')
head(comb_data) # have a look

dat <- data.frame(AUC = comb_data$prv_score, dataset = "test")

# Just trying different AUC_SOTA until I find the one that gives the correct CI upper bound of max(dat$AUC)

AUC = sort(dat$AUC) # all observed performances, sorted (why not)
m = length(AUC)
maxAUC = max(AUC) # best observed perfomance

m_top = length(AUC[AUC>CI[[1]]])
print(m_top)
m_sota = length(AUC[AUC>CI[[2]]])
print(m_sota)
m_above = length(AUC[AUC>CI[[3]]])
print(m_above)
m_above = length(AUC[AUC>auc_sota])
print(m_above)
print(100*m_above/m_top)
