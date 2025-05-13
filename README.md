# 580.430.SP25_G10_Valsartan
## Valsartan Pharmacokinetics and Pharmacodynamics Analysis

This repository contains MATLAB scripts, R scripts, and data files for analyzing the pharmacokinetics (PK) and pharmacodynamics (PD) of Valsartan. The provided code produces key figures and data outputs for understanding the drug's concentration-time profile, receptor binding dynamics, and blood pressure effects under various scenarios, including missed doses and dose-dependent variability.

### Folder Structure and Descriptions:

1. **PKPD driver case**
   - Contains MATLAB scripts and data for running the main pharmacokinetic and pharmacodynamic simulations.
   - **Files:**
     - `driver_case.m`: Contains multiple RunCase functions to generate figures and data outputs. The detailed RunCases are as follows:
       - **RunCase 1:** Takes `val_exp.csv` (empirical free valsartan plasma concentration vs. time) as input and produces **Figure 5**.
       - **RunCase 2:** Produces **Figure 7**.
       - **RunCase 3:** Produces **Figures 2d and 6**.
       - **RunCase 4:** Allows specification of missed dose(s) via `miss_dose_idx` variable. Produces **Figure 8**.
       - **RunCase 5:** Outputs 40 data files used for ShinyApp interactive visualization.
       - **RunCase 6:** Generates `AngII_ATR1_Complex_Timepoints.csv`, representing data for **Figure 3b**.
       - **RunCase 7:** Produces **Figure 9**.
       - **RunCase 8:** Takes `bp_exp.csv` (empirical blood pressure data for untreated patients) as input and produces **Figures 10a and 10c**.
     - `eqns_v2.m`: Setup of ode equations for the Valsartan PK/PD model.
     - `sim0_v2.m`: Main simulation script for running individual patient simulations.
     - `singlepatient_v2.m`: Testing simulation with all compartments plotted.
     - `singlepatient_v3.m`: Plotting the mass balance for free valsartan and produces **Figure 4**.
     - `bp_exp.csv`: Clinical blood pressure data for untreated patients.
     - `val_exp.csv`: Clinical free valsartan plasma concentration data.

2. **PopPK**
   - Contains scripts and data for population pharmacokinetics analysis, sensitivity analysis, and data visualization.
   - **Files:**
     - `poppk.m`: Population pharmacokinetics simulation script and produces **Figure 11**.
     - `eqns_v2.m`: System of equations for the population PK/PD model.
     - `sensitivity_with_output.m`: Takes `valsartan_random_pop_case1_case2.mat` and generates contour and scatter plots shown in **Figure 12**.
     - `Data_3Dplot.m`: Generates data (.csv files) for 3D interactive plots.
     - `beeswarm_sensitivity.r`: Generates beeswarm plots for sensitivity analysis for **Figure 13**.
     - `valsartan_generate_random_population.m`: Generates dataset for population simulations through mean +/- SD of weight and baseline receptor concentrations.
     - `valsartan_relative_sensitivity.mat`: Dataset for relative sensitivity analysis.

3. **sys_pharm_g10**
   - Contains the ShinyApp for interactive visualization of missed-dose scenarios and population sensitivity analysis.
   - **Files:**
     - `app.r`: Main R script for building and running the ShinyApp.

### Usage Instructions:
1. Ensure that MATLAB and R are installed with necessary toolboxes and libraries.
2. Navigate to the `PKPD driver case` folder and run `driver_case.m` to generate required data and figures.
3. Move to the `PopPK` folder to run population pharmacokinetics simulations and sensitivity analysis.

