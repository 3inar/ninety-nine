# The variance of empirical cdf aka auc
# Kajsa Møllersen (kajsa.mollersen@uit.no) 16th of April 2024

# The AUC is equivalent to the cdf of S_{neg}-S_{pos} at value 0.
# The variance of an empirical cdf is 1/n {F(x) - F(x)^2}, where F(x) is the cdf.

# In our case, F(x) = 0.9, since we decided for it to be that way.
# The negative and the positive class are not balanced, but I would assume that 
# n = min(n_{pos},n_{neg}) for the purpose of variance calculations. 

n = 3000
# n_pos = round((584/33126)*n) 

n_pos = n/2 
n_neg = n-n_pos

# n_neg = 53
# n_pos = n-n_neg

class = rep(0, n)
class[1:n_pos] = 1 # true class label

mu_pos = 0
mu_neg = -1.812413
sigma_pos = 1
sigma_neg = 1

mu_neg = -3.624775
sigma_pos = 2
sigma_neg = 2

mu_neg = -12.87944
sigma_pos = 1
sigma_neg = 10

t = -(mu_neg - mu_pos)/(sqrt(sigma_pos^2 + sigma_neg^2))



th_auc = pnorm(t) # should be 0.9

th_varauc = (th_auc - th_auc^2)/n


m = 10000#0
auc = numeric(m)
oldw <- getOption("warn")
options(warn = -1)

for (j in 1:m) {
  predict_pos = rnorm(n_pos, mu_pos, sigma_pos)
  predict_neg = rnorm(n_neg, mu_neg, sigma_neg)
  auc[j] = auc(class,c(predict_pos, predict_neg))
}

options(warn = oldw)

print(n)
var(auc)
th_varauc

