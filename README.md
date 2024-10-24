# ninety-nine

Code for the manuscript *What is the state of the art? Accounting for multiplicity in machine learning benchmark performance* 	
https://doi.org/10.48550/arXiv.2303.07272

Files needed to reproduce numbers/figures:

`Parameters_PublicCompetition.R` - all parameters used in the manuscript

## Section *Estimating state-of-the-art performance*

Bla-bla-bla

## Section *Related work*

Cited bla-bla-bla

## Section *Multiple classifiers and biased state-of-the-art estimation*

## Subsection *Notation*

Aligns with the (simplified) variable names found here. 

## Subsection *A coin-flip example*

> CoinFlip.R \
> optional ProbDistr_thetaSOTA.R

## Section 'The probability distribution of p_SOTA':

ProbDistr_thetaSOTA.R with functions:

`cdf`       - cumulative distribution function
`pmf`       - probability mass function
`expect`    - expectation
`variance`  - variance
`sim_var` - simulated variance
`sim_mean` - simulated mean
`sim_ci` - simulated confidence interval

## Section 'A simulated public competition example':

PublicCompetition.R
calls ProbDistr_thetaSOTA.R and Parameters_PublicCompetition.R


## Section 'Dependent, identical classifiers':

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
