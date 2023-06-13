# parameters

n = 3000 # size of test set
theta = 0.9 # probability of correct prediction
m = 1000 # number of classifiers
alpha = 0.05 # significance level

mu = n*theta # the expected number of correct predictions for a single classifier

theta_min = 0.875 # for the non-identical
theta_max = theta 
step = (theta_max-theta_min)/(m-1)
theta_vec = seq(theta_min, theta_max, step)

rho = 0.6 # correlation coefficient, the number is calculated from Mania (2019)

rep = 10000000 #  number of repetitions, 100,000 gives nice and smooth histograms, 1 million is doable

theta_0 =  theta