# ninety-nine

Code for the manuscript 'What is the state of the art? Accounting for multiplicity in machine learning benchmark performance' 	
https://doi.org/10.48550/arXiv.2303.07272

Files needed to reproduce numbers/figures:

Parameters_PublicCompetition.R - all parameters used in the manuscript

## Section 'Two coin-flip examples': 

> CoinFlip.R  
> optional ProbDistr_thetaSOTA.R

## Section 'The probability distribution of p_SOTA':  

ProbDistr_thetaSOTA.R with functions.  

cdf       - cumulative distribution function  
pmf       - probability mass function. 
expect    - expectation. 
variance  - variance. 
sim_var - simulated variance. 
sim_mean - simulated mean. 
sim_ci - simulated confidence interval. 
                                                  
Section 'A simulated public competition example': PublicCompetition.R calls ProbDistr_thetaSOTA.R and Parameters_PublicCompetition.R
                                                  
                                                  
Section 'Dependent, identical classifiers': dependence.R calls ProbDistr_thetaSOTA.R and Parameters_PublicCompetition.R
                                            dep_id_pmf - simulated pmf for dependent, identical classifiers, fixed or random Y0
                                            
Section 'Non-identical, independent classifiers': nonidentical.R calls Parameters_PublicCompetition.R
                                                  indep_nonid_pmf_fun.R with function 
                                                  indep_nonid_pmf - simulated pmf for non-identical, independent classifiers
                                                  nonid_cdf - cdf for non-identical, independent classifiers
                                                  nonid_pmf - pmf based on cdf
                                                  
Section 'Non-identical, dependent classifiers': dependent_nonidentical.R calls Parameters_PublicCompetition.R and ProbDistr_thetaSOTA.R
                                                dep_nonid_pmf_fun.R with function
                                                dep_nonid_pmf - simulated pmf for dependent, nonidentical classifiers

Section 'A kaggle challenge example': AUROC.R for quick estimate of AUC_SOTA
                                      SOTA_bootstrap.R for better estimate of AUC_SOTA
                                                
Section 'Discussion': auc_generator.R: creates AUC curves
                      melanoma.R: investigates the Melanoma 2020 kaggle competition
                                                
Additional files:

SOTA_estimate.R: Creates toy example from the Melanoma 2020 kaggle competition to estimate theta_SOTA from observed performance
closed_form.R: Extra simulations

Outdated files:

99Club_multiplicity.R: Contains text only
                         
                              
