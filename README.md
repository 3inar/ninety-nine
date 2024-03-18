# ninety-nine

Code for the manuscript *What is the state of the art? Accounting for multiplicity in machine learning benchmark performance* 	
https://doi.org/10.48550/arXiv.2303.07272

Files needed to reproduce numbers/figures:

`Parameters_PublicCompetition.R` - all parameters used in the manuscript

## Section *3 Multiple classifiers and biased state-of-the-art estimation*

### Section *3.2 Two coin-flip examples* 

> CoinFlip.R  
> ProbDistr_thetaSOTA.R - optional 

### Section *3.3 The probability distribution of $\hat{\theta}^{\, max}_{SOTA}$*  

No output, only functions  

> ProbDistr_thetaSOTA.R 

with functions

`cdf`       - cumulative distribution function  
`pmf`       - probability mass function  
`expect`    - expectation  
`variance`  - variance  
`sim_var` - simulated variance  
`sim_mean` - simulated mean  
`sim_ci` - simulated confidence interval 
                                                  
### Section *3.4 A simulated public competition example* 

> PublicCompetition.R  

calls 

`Parameters_PublicCompetition.R` for parameters  
`ProbDistr_thetaSOTA.R` for functions 

## Section *4 Non-i.i.d.*
                                                  
                                                  
### Section *4.1 Dependent, identical classifiers* 

> dependence.R  

`dep_id_pmf` - simulated pmf for dependent, identical classifiers, fixed or random Y0

calls 

`Parameters_PublicCompetition.R` for parameters  
`ProbDistr_thetaSOTA.R` for functions

                                            
### Section *4.2 Non-identical, independent classifiers* 

> nonidentical.R  

calls 

`Parameters_PublicCompetition.R` for parameters 

calls 

> indep_nonid_pmf_fun.R 

for functions  
`indep_nonid_pmf` - simulated pmf for non-identical, independent classifiers  
`nonid_cdf` - cdf for non-identical, independent classifiers  
`nonid_pmf` - pmf based on cdf
                                                  
### Section *4.3 Non-identical, dependent classifiers* 

> dependent_nonidentical.R  

calls `Parameters_PublicCompetition.R` for parameters 

> dep_nonid_pmf_fun.R 

for function  
`dep_nonid_pmf` - simulated pmf for dependent, nonidentical classifiers

calls `ProbDistr_thetaSOTA.R` for function `sim_ci` - simulated conf int



## Section *5 A kaggle challenge example* 

### Section *5.1 Strategies for estimating SOTA*

> AUROC.R 

for quick estimate of AUC_SOTA, ignoring the multiplicity effect

> SOTA_bootstrap.R 

for better estimate of AUC_SOTA

calls `dep_nonid_pmf_fun.R` for the function `dep_nonid_pmf` - simulated pmf
                                                
### Section *5.2 Simulations of AUCs* 

> Einar

auc_generator.R: creates AUC curves

melanoma.R: investigates the Melanoma 2020 kaggle competition
                                                
## Additional files:

SOTA_estimate.R: Creates toy example from the Melanoma 2020 kaggle competition to estimate theta_SOTA from observed performance

closed_form.R: Extra simulations

## Outdated files:

99Club_multiplicity.R: Contains text only
                         
                              
