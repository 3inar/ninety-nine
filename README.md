# ninety-nine

Code for the manuscript *Accounting for multiplicity in machine learning benchmark performance* 	
https://doi.org/10.48550/arXiv.2303.07272

Files needed to reproduce numbers/figures: \
`Parameters_PublicCompetition.R` - all parameters used in the manuscript

Figures are found locally.

## *Estimating state-of-the-art performance*

Bla-bla-bla

## *Related work*

Cited bla-bla-bla

## *Multiple classifiers and biased state-of-the-art estimation*

### *Notation*

Aligns with the (simplified) variable names found here. 

### *A coin-flip example*

`CoinFlip.R` 

Produces the numbers found in this subsection. There are two examples, and 
example b) is used in the manuscript. The output is re-produced by 
simulation, and then finally the output is re-produced using functions from

`ProbDistr_thetaSOTA.R`

No figures

### The probability distribution of $\hat{\theta}_{max}(X)$

`ProbDistr_thetaSOTA.R` with functions:

`cdf()`       - cumulative distribution function \
`pmf()`       - probability mass function \
`expect()`    - expectation \
`variance()`  - variance \
`sim_var()` - simulated variance \
`sim_mean()` - simulated mean \
`sim_ci()` - simulated confidence interval

Does not produce any output. As an option at the end, you can check that the functions 
produces the correct numbers from the coin-flip example.

### A simulated public competition example

`PublicCompetition.R` \
calls `ProbDistr_thetaSOTA.R` and `Parameters_PublicCompetition.R`

Output: \
The numbers in the text and some intermediate calculations that we did
not include in the manuscript. There are check-ups along the way. \

The entries for table *Expected values and standard deviations for $\hat\theta_{max}(X)$ *.
The output table displays all combinations of m, n, theta.

Figures:\

`multi_ci` *The pmfs of two single $\hat{\theta}(X)$ compared to $\hat{\theta}_{\max}$ statistics*\
Against a mass of identical classifiers, a single significantly better
classifier has little chance of beating the observed sample maximum 

`cumul_fail` *Cdf and pmf for at least one classifier having at most z failures* a)\
Displaying the function `cdf()` (ProbDistr_thetaSOTA.R) for the parameters in Parameters_PublicCompetition.R

`pdf_fail` *Cdf and pmf for at least one classifier having at most z failures* b)\
Displaying the function `pmf()` (ProbDistr_thetaSOTA.R) for the parameters in Parameters_PublicCompetition.R

`bias_sd_m_n_theta` *The bias and standard deviation of the sample maximum estimator*\
The bias and standard deviation of $\hat{\theta}_{\max}$ as a function of each of its three parameters; number of classifiers, test set size and probability of success.\
Consists of six subfigures\
`bias_m` and `sd_m`\
`bias_n` and `sd_n`\
`bias_theta` and `sd_theta`\

There are optional simulations at the end that confirms that the functions are correct. 


## Dependency and non-identical theta's

dependence.R
calls ProbDistr_thetaSOTA.R and Parameters_PublicCompetition.R
`dep_id_pmf` - simulated pmf for dependent, identical classifiers, fixed or random Y0

## Section 'Non-identical, independent classifiers':

nonidentical.R
calls Parameters_PublicCompetition.R

indep_nonid_pmf_fun.R with functions
`indep_nonid_pmf` - simulated pmf for non-identical, independent classifiers
`nonid_cdf` - cdf for non-identical, independent classifiers
`nonid_pmf` - pmf based on cdf

## Section 'Non-identical, dependent classifiers':

dependent_nonidentical.R
calls Parameters_PublicCompetition.R and ProbDistr_thetaSOTA.R

dep_nonid_pmf_fun.R with function
`dep_nonid_pmf` - simulated pmf for dependent, nonidentical classifiers


## Section 'Real world examples':

AUROC.R for quick estimate of AUC_SOTA

SOTA_bootstrap.R for better estimate of AUC_SOTA

auc_functions.R contains functions for working with and simulating AUCs:

- `make_classifier` makes an object of class `classifier` with true aus as specified
- `predict.classifier` makes predictions from a `classifier` object
- `sds` gets the standard deviations for undelying distributions of a `classifier` object
- `mus` gets the means for underlying distributions of a `classifier` object
- `draw_correlated` conditional draw from a bivariate normal distribution
- `correlated_predict` make correlated predictions where a `classifier` has a certain correlation to a "hidden" reference classifier
- `empirical_auc` fast calculation of AUC from a `list(predicted=..., truth=...)`
- `sim_competition` simulates an uncorrelated competition with given true AUCs

auc_simulations.R generates data for the figures and numbers in this section

auc_figure_ex_1.R makes the figure that shows the distribution of maximum AUC among 3000 classifiers all with true AUC 0.9

auc_figure_ex_2.R makes the figures that show a true AUC distribution based on
truncated kaggle scores and a simulation based on these

auc_figure_roc_curves.R makes the figure that shows example ROC curves for a classifier with true AUC of .9

## Section 'Discussion':

melanoma.R: investigates the Melanoma 2020 kaggle competition

## Additional files:

SOTA_estimate.R: Creates toy example from the Melanoma 2020 kaggle competition to estimate theta_SOTA from observed performance

closed_form.R: Extra simulations

## Outdated files:

99Club_multiplicity.R: Contains text only
auc_generator.R: creates AUC curves based on underlying beta distributions
