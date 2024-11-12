# ENACT-BigCACTUS

## Overview
This workflow provides the code alongside detailed instructions to replicate the analysis described in the paper **Flight, L., Brennan, A., Wilson, I., & Chick, S. E. (2024). A Tutorial on Value-Based Adaptive Designs: Could a Value-Based Sequential 2-Arm Design Have Created More Health Economic Value for the Big CACTUS Trial? Value in Health. https://doi.org/10.1016/j.jval.2024.06.004**. This work is part of the wider ENACT project (https://www.sheffield.ac.uk/ctru/completed-trials/enact).

It could be generalised to fit the needs of your clinical trial. These instructions highlight how it was carried out for the Big CACTUS trial (https://www.sheffield.ac.uk/ctru/completed-trials/big-cactus). Anonymised and simulated data is provided, and sample sizes of 95 and 435 have been used.

## Information

Further details on the provided, created and derived data, as well as the programs and functions used, will be included in the [instructions](#instructions).

### Data provided:
- pilot_data_set.csv
- data_set.csv
- BigCACTUS_3arm_costs.csv

Note: This data has been appropriately anonymised in accordance with Information Governance guidelines for sharing purposes. It serves as an illustrative representation of the real data, preserving similar patterns and key correlations. 

### Data created/derived:
- CACTUS_spend_data_matlab.csv (1 total)
- bndlower.csv, bndupper.csv, muvec.csv, threspoint.csv, tvec.csv for Tmax 95 & 435 (10 total)
- Bootstraps CSVs (4 total)

### Programs:
- BC_creating_bootstrap_samples.R
- BC_running_analyses.R
- BigCACTUS_DelayDriver_Edit.m
- CACTUS_HE_and_cost_analysis.R 
- DoBoundaryPlotsBC.m
- DoSectionFivePlotsBC.m
- Modelling_pilot_data.R
- SetBC_basecase.m
- https://github.com/sechick/htadelay

### Functions:
- BC_bootstrap_analysis_function.R
- BC_model_param_function.R
- Big_CACTUS_R_model.R
- CACTUS_HE_model.R 
- Unadjusted_pilot_data_analysis.R

## Instructions

1. Download the required software packages:

     a) MATLAB (original analyses conducted in R2024a) - https://uk.mathworks.com/products/matlab.html
     
     b) RStudio and R (original analyses conducted in R 4.4.1) - https://rstudio-education.github.io/hopr/starting.html
     
     c) In addition, there will be software modules that will operate within each of these tools. The htadelay package has MATLAB code to generate stopping boundaries and make other general computations for value-adaptive trials. This ZIP has additional R code that will be used to manage clinical trial data and run bootstrap simulations and computations reported regarding the CACTUS trial

 2. Set up workspace:

      a) Download the folder from GitHub and place it in your working directory
      
      b) Ensure you know your working directory's file path (where the unzipped folder is located)
      
      c) When unzipped, the following folder organisation will be created. Some of the R code explicitly uses this structure (but you can edit the directory structure in the code if you wish to adapt it further)
```
Value-Based Sequential 2-Arm Design
├── Data
│    ├── Derived
│    │   ├── Bootstraps
│    │   ├──Tmax 95
│    │   └──Tmax 435
│    └── Raw
└── Programs
     ├── Functions
     └── htadelay-master
```

3. Modelling_pilot_data.R using RStudio:

    _The script uses pilot_data_set.csv, which contains information for the 34 (anonymised and simulated) participants from the CACTUS pilot trial, including data on treatment allocation, utilities, primary outcome, quality-adjusted life year, and costs._
    
    _It loads the dataset described, analyses the pilot data, runs the health economic model (including probabilistic sensitivity analysis), and performs the VOI analysis to calculate the fixed VOI sample size (n=435)._
    
    _It requires function Unadjusted_pilot_data_analysis.R, which analyses the pilot data to calculate the necessary parameters for the CACTUS pilot health economic model and function CACTUS_HE_model.R, which conducts the CACTUS pilot health economic model (with simplifications)._
    
    _Important outputs from this script are ‘prob_mean_INB’ and ‘prob_sd_INB’._

    a) Open the script
    
    b) Line 23 needs to be updated to your working directory; change xxxxxx to match your file path
    
    c) Run the script and make a note of your **‘prob_mean_INB’ and ‘prob_sd_INB’**

4. Download the htadelay MATLAB package:

    _This package calculates the value-based sequential design stopping boundaries._
   
    a) Follow the link (https://github.com/sechick/htadelay); click ‘code’; click ‘download zip’
   
    b) This should be saved in the ‘Programs’ folder within your working directory
   
    c) Move the provided MATLAB files (BigCACTUS_DelayDriver_Edit.m, DoBoundaryPlotsBC.m, DoSectionFivePlotsBC.m, SetBC_basecase.m) into the downloaded ‘htadelay-master’ folder

5. SetBC_basecase.m and BigCACTUS_DelayDriver_Edit.m using MATLAB:

    a) In BigCACTUS_DelayDriver_Edit.m, lines 9 and 60-64 need to be updated to your working directory; change xxxxxx to match your file path
   
    b) In SetBC_basecase.m (line 35) and BigCACTUS_DelayDriver_Edit.m (lines 60-64), change the value of xx in basic.TMax = xx and Tmax xx to 95
   
    c) In SetBC_basecase.m set the value of basic.mu0 in line 29 to match your ‘prob_mean_INB’
   
    d) Run BigCACTUS_DelayDriver_Edit.m (this uses the htadelay package, SetBC_basecase.m, and DoBoundaryPlotsBC.m)
   
    e) Repeat steps a-c for Tmax=435. You should now have **5 CSVs in each of the Tmax 95 and Tmax 435 folders** (located in Value-Based Sequential 2-Arm Design/Data/Derived)

6. CACTUS_HE_and_cost_analysis.R using RStudio:

    _This script uses data_set.csv, which contains information on the 270 (anonymised & simulated) participants from the Big CACTUS trial, including the imputed EQ-5D data. It also uses BigCACTUS_3arm_costs.csv, which includes the calculated fixed and variable costs of conducting the research for the Big CACTUS trial._
   
    _It runs the deterministic health economic model on all the data and runs code for sequential analysis of outcomes, as well as calculating costs from the grant application._
   
    _It requires function Big_CACTUS_R_model.R, which conducts the Big CACTUS health economic model (with some simplifications), and function BC_model_param_function.R, which calculates the health economic model parameters required to run the previous function._

    _The CSV CACTUS_spend_data_matlab.csv will be generated and exported (located in Value-Based Sequential 2-Arm Design/Data/Derived)._

    a) Open the script
   
    b) Line 23 needs to be updated to your working directory; change xxxxxx to match your file path
   
    c) Run the script

7. BC_creating_bootstrap_samples.R using RStudio:

    _This script creates the bootstrap samples using data_set.csv. Note that lines 290-323 (each scenario) need to be commented and uncommented out as required, as the bootstraps were run in batches of 2500 due to the time the code took to run._
   
    a) Open the script
   
    b) Line 24 needs to be updated to your working directory; change xxxxxx to match your file path
   
    c) Change the value of ‘fixed_mu_n’ in line 128 and ‘mu_n’ in line 265 to match your ‘prob_mean_INB’
   
    d) Run the script (a total of four times). You should now have **4 CSVs in the Bootstraps folder** (located in Value-Based Sequential 2-Arm Design/Data/Derived)

8. BC_running_analyses.R using RStudio:

    _This script merges the bootstrap samples (boot_95 and boot_435), analyses, and creates some useful output (output_95 and output_435)._
   
    _It uses data_set.csv, MATLAB-created CSVs, and the bootstrap CSVs._
   
    _It requires functions Big_CACTUS_R_model.R, BC_model_param_function.R, and BC_bootstrap_analysis_function.R, which reads the MATLAB output into R._
   
    _Important outputs from this script are ‘boot_95’ and ‘boot_435’._
   
    a) Open the script
   
    b) Line 23 needs to be updated to your working directory; change xxxxxx to match your file path
   
    c) Change the value of ‘priormean’ in lines 50 and 74 to match your ‘prob_mean_INB’
   
    d) Run the script
