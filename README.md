# ninety-nine

Code for the manuscript *Accounting for multiplicity in machine learning benchmark performance* 	
https://doi.org/10.48550/arXiv.2303.07272

Files needed to reproduce numbers/figures: \
`Parameters_PublicCompetition.R` - all parameters used in the manuscript

Figures are found locally.

## Estimating state-of-the-art performance

Bla-bla-bla

## *Related work*

Cited bla-bla-bla

## *Multiple classifiers and biased state-of-the-art estimation*

### *Notation*

Aligns with the (simplified) variable names found here. 

### A coin-flip example

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

The entries for table `mnp` *Expected values and standard deviations for $\hat\theta_{max}(X)$ *.
The output table displays all combinations of m, n, theta.

Figures:\

`multi_ci` *The pmfs of two single $\hat{\theta}(X)$ compared to $\hat{\theta}_{\max}$ statistics*\
Against a mass of identical classifiers, a single significantly better
classifier has little chance of beating the observed sample maximum 

`cumul_fail` *Cdf and pmf for at least one classifier having at most z failures* a)\
Displaying the function `cdf()` (ProbDistr_thetaSOTA.R) for the parameters in Parameters_PublicCompetition.R

`pmf_fail` *Cdf and pmf for at least one classifier having at most z failures* b)\
Displaying the function `pmf()` (ProbDistr_thetaSOTA.R) for the parameters in Parameters_PublicCompetition.R

`bias_sd_m_n_theta` *The bias and standard deviation of the sample maximum estimator*\
The bias and standard deviation of $\hat{\theta}_{\max}$ as a function of each of its three parameters; number of classifiers, test set size and probability of success.\
Consists of six subfigures\
`bias_m` and `sd_m`\
`bias_n` and `sd_n`\
`bias_theta` and `sd_theta`\

There are optional simulations at the end that confirms that the functions are correct. 


## Non-identical $\theta$s and classifier dependency

The entries for table `noniid` are found under their respective subsections.

Figure:\
`bias_sd_thetamin_rho` *Bias for non-identical and dependent classifiers*\
Bias of the sample maximum estimator as a function of $\min (\Theta)$ for non-identical classifiers and of $\rho_0$ for conditional independent classifiers.

Consists of three subfigures\
`bias_sd_thetamin`: generated in `nonidentical_hierarchical.R`\
`bias_thetamin_rho` and `sd_thetamin_rho`: generated in `dependent_nonidentical_hierarchical.R`,\
see the respective subsections below.

### Non-identical, independent classifiers

`nonidentical_hierarchical.R`\
calls `Parameters_PublicCompetition.R`\
calls indep_nonid_pmf_fun.R with functions \ 
`indep_nonid_pmf` - simulated x for non-identical, independent classifiers \
`nonid_cdf` - cdf for non-identical, independent classifiers, `eq:noniid_cdf`\
`nonid_pmf` - pmf based on cdf, `eq:noniid_pmf`

Output:

The entries for table `noniid`.

Figures:

`noniid_cdf`: *Cdf for non-identical classifiers*\
`noniid_pmf`: *Pmf for non-identical classifiers*\
Corresponds to figures `cumul_fail` and `pmf_fail` in `PublicCompetition.R`, 
section *A simulated public competition example*\
`bias_sd_d`: Bias and standard deviation as a function of $d = b-a$ in $\mathcal{U}(a,b)$,
subfigure for `bias_sd_thetamin_rho`

### Identical, dependent classifiers

`dependence.R`\
calls `Parameters_PublicCompetition.R`\
calls `ProbDistr_thetaSOTA.R`

Output:

The entries for table `noniid`.

### Non-identical, dependent classifiers

`dependent_nonidentical_hierarchical.R`\
calls `Parameters_PublicCompetition.R`\
calls `dep_nonid_pmf_fun.R`\
calls `ProbDistr_thetaSOTA.R`

Output:

The entries for table `noniid`.

Figures:

`bias_thetamin_rho`: Bias as a function of $\rho$ for various $d = b-a$ in $\mathcal{U}(a,b)$,
subfigure for `bias_sd_thetamin_rho`\
`sd_thetamin_rho`: Standard deviation as a function of $\rho$ for various $d = b-a$ in $\mathcal{U}(a,b)$,
subfigure for `bias_sd_thetamin_rho`\


## Real world examples

bla-bla-bla

### Estimating $\theta_{SOTA}$

blah

#### Multi-Class Prediction of Obesity Risk

`SOTA_accuracy_shrink.R` ** need some more work\
calls `Parameters_PublicCompetition.R`\
calls `dep_nonid_pmf_fun.R`

Output:

The entries in Table `obesity`

Figures:

`obesity`: Suggested $\theta'$s and the corresponding realisations\
with subfigures. Exact numbers in Table `obesity`

`obesity_kaggle` Sample estimates and single CI for $\hat{\theta}_{\max}$ \
`obesity_direct_bootstrap` A simulated realisation from `obesity_kaggle`, the expected sample maximum and its CI\
`obesity_cropped_for_expect` Sample estimates cropped at $0.9063$, the proposed $\theta_{SOTA}$\
`obesity_cropped_for_expect_realisation` A simulated realisation from c), the expected sample maximum and its CI

#### Cassava Leaf Disease

### Estimating $AUC_{SOTA}$

#### Melanoma Classification

#### Simulation of AUC

#### Simulation results for uncorrelated classifiers

#### Correlation between classifiers in the AUC simulations

## Discussion



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

`SOTA_estimate.R`: Creates toy example from the Melanoma 2020 kaggle competition to estimate theta_SOTA from observed performance

`closed_form.R`: Extra simulations

`var_auc.R`: Interesting only for balanced classes

## Outdated files:

`99Club_multiplicity.R`: Contains text only\
`auc_generator.R`: creates AUC curves based on underlying beta distributions

early exploration of different options, treating AUC as accuracy:\
`AUROC.R`: confidence interval \
`SOTA_bootstrap.R`: melanoma

## List of files and the section where they are described:

"99Club_multiplicity.R"  			    - Outdated files\
 "auc_figure_ex_1.R" 		    	        - Einar      \        
"auc_figure_ex_2.R" 			            - Einar   \
 "auc_figure_roc_curves.R"        		- Einar    \           
"auc_functions.R"                		- Einar   \
"auc_generator.R"               		- Einar   \
 "AUC_misc.R"                     		- Einar  \ 
"auc_simulations.R"             		- Einar   \
"AUROC.R" 					- Outdated files           \          
"closed_form.R" 				- Additional files      \          
"CoinFlip.R" 					- A coin-flip example      \              
"dep_nonid_pmf_fun.R" 		    	- Non-identical, dependent classifiers          \
"dependence.R" 				            - Identical, dependent classifiers             \     
"dependent_nonidentical_hierarchical.R" 		    - Non-identical, dependent classifiers       \
"indep_nonid_pmf_fun.R" 			    - Non-identical, independent classifiers      \  
"melanoma.R" 				              - Einar   \
"nonidentical_hierarchical.R" 				          - Non-identical, independent classifiers\
"Parameters_PublicCompetition.R" 	- ninety-nine\
"plotting_params.R" 			- Einar   \
"ProbDistr_thetaSOTA.R" 			- The probability distribution of $\hat{\theta}_{max}(X)$         \
"PublicCompetition.R" 			- A simulated public competition example           \
"SOTA_bootstrap_accuracy.R"      	- Kajsa\
"SOTA_bootstrap.R" 			- Outdated        \      
"SOTA_estimate.R"               		- Additional files\
"var_auc.R"                  			- Additional files\


