# Learning Treatment Effects under Covariate Dependent Left Truncation and Right Censoring

## Contents of the Repository

- `src/`: Code for implementing the estimators.

- `src_R/`: Old code with computation all in R. The current `src/` folder has part of the computation in Cpp to speed up.

- `OSG/`: Code used to run simulations on OSG.
  - `ATE/`: Code used for ATE simulation.
    - `main_simu2_boot_OSG.R`: Main code for simulating data and computing the 'dr' estimator.
    - `main_cfdr_simu2_boot_OSG.R`: Main code for simulating data and computing the 'cf' estimator.
    - `organize_OSG_results_ATE.R`: Code for organizing the results collected from OSG.
    - `prepare_inputs.R`: Code for preparing seeds to run simulations on OSG.
    
  - `CATE/`: code for CATE simulation.
    - `main_HTEx_cf.R`: Main code for simulating data from setting 'HTEx' and compute the CATE estimator. 
    'HTE3', 'HTE5', 'HTE4' correspond to Scenarios (i), (ii), and (iii), respectively.
    - `organize_OSG_results_HTE.R`: Code for organizing the results collected from OSG.
    
  - `HAAS_analysis/`: code for HAAS data analysis.
    - `main_HAAS_educ2.R`: The main code for CATE analysis.
    - `main_HAAS_educ2_ATE_OSG.R`: The main code for estimating ATE on OSG.
    
    