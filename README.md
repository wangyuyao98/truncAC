# Learning Treatment Effects under Covariate Dependent Left Truncation and Right Censoring

## Contents of the Repository

- `src/`: Code for implementing the estimators.

- `OSG/`: Code used to run simulations on OSG.
  - `ATE/`: Code used for ATE simulation.
    - `main_simu2_boot_OSG.R`: Main code for simulating data and computing the 'dr' estimator.
    - `main_cfdr_simu2_boot_OSG.R`: Main code for simulating data and computing the 'cf' estimator.
    - `organize_OSG_results_ATE.R`: Code for organizing the results collected from OSG.
    - `prepare_inputs.R`: Code for preparing seeds to run simulations on OSG.
    
  - `CATE/`: code for CATE simulation.
    - `main_HTEx_cf.R`: Main code for simulating data from setting 'HTEx' and compute the CATE estimator. 
    'HTE3', 'HTE5', 'HTE4' correspon to Scenarios (i), (ii), and (iii), respectively.
    - `organize_OSG_results_HTE.R`: Code for organizing the results collected from OSG.
    
    
    
    