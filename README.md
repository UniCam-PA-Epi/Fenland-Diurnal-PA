# Fenland-Diurnal-PA

<p align="center">
  <img src="Code/Legacy/FENST-logo-rgb-stacked-green-CROPSQUARE.png" alt="Logo for the Fenland study, showing a stylized windmill in Cambridge blue color.">
</p>

This repository provides the Stata code used to conduct the analysis for the manuscript "Fenland Diurnal PA and Metabolic Risk."  The research examines how total physical activity energy expenditure and variations in physical activity throughout the day relate to metabolic health risk factors, such as blood glucose, insulin levels, and blood pressure, in participants from the Fenland study population cohort.

## :open_book: Table of Contents
* [Usage](#hammer_and_wrench-usage)
* [Code File Descriptions](#bookmark_tabs-code-file-descriptions)
* [Data Availability](#bar_chart-data-availability)
* [Contributing](#raising_hand_woman-raising_hand_man-contributing)
* [Contact](#telephone_receiver-contact)
* [License](#balance_scale-license)

## :hammer_and_wrench: Usage

This repository contains Stata code to reproduce the analysis for the "Fenland Diurnal PA and Metabolic Risk" manuscript. The main analysis is conducted using the `0_Main.do` file. 

**To run the code:**

1. **Prerequisites:** Ensure you have Stata SE version 17.0 or above installed on your system and have access to the Fenland study dataset.
2. **Clone the repository:** `git clone https://github.com/UniCam-PA-Epi/Fenland-Diurnal-PA.git`
3. **Navigate to the project code directory:** `cd Fenland-Diurnal-PA/Code`
4. **Run the main do file:**  `stata do 0_Main.do "path to Fenland release"` (replace `"path to Fenland release"` with the actual file path to the Fenland study dataset)

**Important Notes:**

* The `0_Main.do` file calls other do-files within the repository, so ensure all files are kept in the same directory.
* The code expects the Fenland data to be in a specific format. Please refer to the manuscript or any accompanying documentation for details on data preparation.
* The code assumes that the input PAEE data has been pre-processed and represents average daily values (i.e., averaged over 24-hour periods).
* The code generates output files (e.g., tables, figures) in project code directory.

## :bookmark_tabs: Code File Descriptions

This section provides a brief overview of the functionality of each code file in the analysis sequence.
<br> <br> 

**[0_Main.do](Code/0_Main.do)**

`0_Main.do` is the primary script that executes the complete analysis pipeline for the "Fenland Diurnal PA and Metabolic Risk" manuscript. Users can reproduce the full analysis by running this single file, which automatically calls and executes all other code files in the correct sequence.
<br> <br> 

**[1_initialiseCovariates.do](Code/1_initialiseCovariates.do)**

`1_initialiseCovariates.do` performs initial data preparation for the main covariates used in the analysis, including:

* **Variable Renaming:** Simplifies variable names from the Fenland dataset to improve clarity and consistency within the analysis code. 
* **Categorical Variable Creation:** Generates categorical variables for sociodemographic factors (e.g., age groups, education level) and anthropometric measures (e.g., BMI categories).
<br> <br>

**[2_initialiseOutcomes.do](Code/2_initialiseOutcomes.do)**

`2_initialiseOutcomes.do` prepares the outcome variables for analysis:

* **Variable Renaming:** As in the previous file, this step renames outcome variables from the Fenland dataset for simplicity.
* **Outcome Variable Definitions:** Creates the outcome variables used in the analysis. These include:
    * **Metabolic Measures:** 
        * `glucose120`: 2-hour glucose (mmol/L)
        * `insulin`: Fasting insulin (pmol/L)
        * `leptin`: Leptin (ng/ml) 
        * `nefa`: Non-esterified fatty acids (μmol/L)
        * `adiponectin`: Adiponectin (μg/ml)
        * `ldl`: Low-density lipoprotein cholesterol (mmol/L)
        * `hdl`: High-density lipoprotein cholesterol (mmol/L)
    * **Cardiovascular Measures:**
        * `mbpsys`: Systolic blood pressure (mmHg)
        * `mbpdia`: Diastolic blood pressure (mmHg)
    * **Inflammatory Measures:**
        * `crp`: C-reactive protein (mg/L)     
<br> 

**[3_applyCosinorModel.do](Code/3_applyCosinorModel.do)**

`3_applyCosinorModel.do` applies a cosinor model to analyse diurnal patterns in individual-level physical activity energy expenditure (PAEE) collected over approximately 6 days of continuous monitoring.

* **Cosinor Model:** To quantify rhythmic patterns in PAEE over a 24-hour period, we use a multi-frequency cosinor model. This model is implemented using a gamma Generalized Linear Model (GLM) with a log link function. PAEE data are strictly positive and exhibit a right-skewed distribution. A gamma GLM is well-suited for this type of data, providing a better statistical fit than a Gaussian (normal) distribution. The main features of the cosinor model are:
  * **Mesor ($M$):** This relates to the average level of PAEE over the 24-hour cycle on the log scale.  
  * **Amplitude ($\alpha$):** This relates to the difference between the peak PAEE and the mesor on the log scale, representing the strength of the rhythm.
  * **Acrophase ($\phi$):** The time of day at which the peak of the rhythm occurs.
    
  **Important Note:** Due to the log link function, the mesor and amplitude values reflect the patterns on the logarithm of PAEE prior to exponentiation.

* **Multi-frequency Rhythms:** To capture both circadian (24-hour) and ultradian (shorter than 24-hour) rhythms, we jointly fit cosinor models to 24-hour, 12-hour, and 8-hour cycles. This allows us to identify potential rhythms with different periods. The combined cosinor model is represented as:
    * $f_{24}(t) = \alpha_{24} \cos (\frac{2 \pi}{24} (t - \phi_{24})$
    * $f_{12}(t) = \alpha_{12} \cos (\frac{2 \pi}{12} (t - \phi_{12})$
    * $f_{8}(t) = \alpha_{8} \cos (\frac{2 \pi}{8} (t - \phi_{8})$    
    * $f_{PAEE}(t) = \exp (f_{24}(t) + f_{12}(t) + f_{8}(t) + M)$
 
  where:
  *  $t$ represents clock time in hours, ranging from 0 (midnight) to 23.
  *  $f_{24}(t)$, $f_{12}(t)$, $f_{8}(t)$ are cosinor models for the 24-hour, 12-hour, and 8-hour cycles.
  *  $\alpha_{24}$, $\alpha_{12}$, and $\alpha_{8}$ are the estimated amplitudes for the 24-hour, 12-hour, and 8-hour cycles, respectively.
  *  $\phi_{24}$, $\phi_{12}$, and $\phi_{8}$ are the estimated acrophases for the respective rhythms.
  *  $f_{PAEE}(t)$ is the combined model representing superimposition and exponentiation of the cosinor models to produce a PAEE profile. The exponentiation, incorporating the estimated mesor value ($M$), ensures that the predicted PAEE values are always positive and allows for the modeling of varying average levels of PAEE.

* **Additional PAEE Metrics:** In addition to the cosinor parameters, the code also estimates:
  
  * **Total 24-hour PAEE:** The average total PAEE over a 24-hour period.
  * **Maximum Achieved PAEE:** The highest PAEE value within the 24-hour period.
  * **Hour of Maximum PAEE:** The hour of the day when PAEE is at its peak.
   
* **Stata Implementation:** The Stata command used to fit this multi-frequency cosinor model is:
  * `glm paee_hour sin24 cos24 sin12 cos12 sin8 cos8 [aw=pwear_hour], family(gamma) link(log)`. In this model:
    * `paee_hour` represents the PAEE level for a given hour.
    * `sin24`, `cos24`, `sin12`, `cos12`, `sin8`, and `cos8` are the sine and cosine transformations of the 24-hour clocktime values. These transformations convert the angular hour data into Cartesian coordinates, allowing the cyclical patterns to be modelled using a linear framework.
    * `pwear_hour` is the probability that an activity monitor was worn in a given hour, which is used as an analytic weight to account for potential variations in wear time.
      
  * After using the multi-frequency cosinor model to estimate regression coefficients ``_b[`sin`p']`` and ``_b[cos`p']`` for each `p` value in `24 12 8`, the following cosinor model feautres are estimated for each `p`:
    * **Amplitude ($\alpha$)** values for each rhythm ($\alpha_{24}$, $\alpha_{12}$, and $\alpha_{8}$) are estimated using the following Stata command:
      * ``nlcom sqrt(_b[sin`p']^2+_b[cos`p']^2)``.  
    * **Acrophase ($\phi$)** values for each rhythm ($\phi_{24}$, $\phi_{12}$, and $\phi_{8}$) are estimated using the following Stata command:
      * ``nlcom cond(_b[sin`p']<0,`p',0) + atan2(_b[sin`p'],_b[cos`p'])*`p'/(2*_pi)``.
    
* **Outlier Removal:** An IQR-based method to remove outliers in the PAEE data is applied before fitting the cosinor model. 
<br><br>

**[4_applyExclusions.do](Code/4_applyExclusions.do)**

`4_applyExclusions.do` applies exclusion criteria to the study sample for the main analysis. Participants are excluded based on the following criteria:

* **Cosinor Model Fit:** Removes participants for whom the cosinor model did not adequately fit the PAEE data. This is determined by examining the p-values associated with the sine and cosine terms in the 24-hour, 12-hour, and 8-hour cycles.
* **Insufficient Wear Time:** Excludes participants with insufficient wear time of the accelerometer device.
* **Missing Body Composition Data:** Removes participants without measurements for fat mass and fat-free mass.

**Specific Exclusion Rules:**

* Participants are excluded if the p-values for all sine and cosine terms in the 24-hour, 12-hour, and 8-hour cosinor models are greater than 0.05.
* Participants are excluded if the p-value for the total PAEE estimate is greater than 0.05.
* Participants with missing values for any of the estimated cosinor parameters (mesor, amplitudes, acrophases) are excluded.
* Participants with missing standard errors for the acrophases are excluded.
* Participants with less than 72 hours of consolidated wear time are excluded.
* Participants with missing fat mass or fat-free mass measurements are excluded.
<br> <br> 

**[5_clusterAnalysis.do](Code/5_clusterAnalysis.do)**

`5_clusterAnalysis.do` performs a k-means cluster analysis to group participants based on their diurnal PAEE patterns.

* **K-means Clustering:** Performs k-means clustering on the standardised cosinor parameters using the L2 distance metric. The starting cluster centers are randomly assigned using a fixed seed.
* **Standardisation:** Standardises the cosinor parameters (sine and cosine terms for 24-hour, 12-hour, and 8-hour cycles, and the mesor) to ensure that all variables contribute equally to the distance calculations in the k-means algorithm.
* **Optimal k Determination:** Determines the optimal number of clusters (k) using the elbow method with within-cluster sum of squares (WCSS) and linear splines regression.
* **Cluster Visualisation:** Generates a [panel figure](Code/Figures/kClusterPAEEProfiles.png) of the average PAEE profiles for each cluster.
<br> <br> 

**[6_descriptivesTables.do](Code/6_descriptivesTables.do)**

`6_descriptivesTables.do` generates the descriptive statistics tables for the manuscript by summarising study participant characteristics.

* **Descriptive Statistics:** Calculates median and interquartile range (IQR) for continuous variables and frequencies (counts and percentages) for categorical variables.
* **Variable Grouping:** Presents descriptive statistics for the total sample and stratified by subgroups (sex and cluster groups derived from the k-means cluster analysis).
* **Cosinor Parameter Transformation:** Transforms circular acrophase variables (which represent time of day) to their sine and cosine components for calculation of median and IQR.
* **Output:** Exports the descriptive statistics to an [Excel file](Code/Tables/6_descriptiveTables.xlsx) with separate sheets for each subgroup analysis.
<br> <br> 

**[7_violinPlots.do](Code/7_violinPlots.do)**

`7_violinPlots.do` creates violin plots to visualize the distributions of outcome variables across tertiles of total PAEE, stratified by sex.

* **Violin Plots:** Uses the `violinplot` command (from the `violinplot` package) to generate violin plots, which display the density distribution of the outcome variables.
* **Data Transformation:** Applies a natural logarithmic transformation to the outcome variables to improve visualisation and interpretation.
* **PAEE Tertiles:** Divides the participants into tertiles based on their total PAEE (`totalPAEE_hat`), separately for males and females.
* **Stratification:** Generates separate violin plots for each sex, allowing for comparison of outcome distributions across PAEE tertiles within each sex.
* **Output:** Saves the individual violin plots as graph files and combines them into a [single figure](Code/Figures/violinPlots.png).
<br> <br> 

**[8_cosinorFeatureAnalysis.do](Code/8_cosinorFeatureAnalysis.do)**

`8_cosinorFeatureAnalysis.do` investigates the associations between cardiometabolic risk factors and features of the multi-component cosinor model, using a nested GLM approach with predictions at shifted PAEE profiles.

* **Cosinor Component Decomposition:** Decomposes the circular acrophase variables into sine and cosine components for appropriate inclusion in linear models.
* **Nested GLMs:** Constructs and fits nested generalized linear models (GLMs) with different sets of covariates to assess the contribution of various predictors to the outcome variables. This includes baseline characteristics, body composition measures, and interactions between cosinor features (mesor, amplitude, acrophase) and other covariates.
* **Model Evaluation:** Uses likelihood ratio tests to evaluate the significance of adding different blocks of predictors to the model.
* **Predicted Outcome Curves:** Generates predicted outcome values for time-shifted PAEE profiles by systematically adjusting the acrophase of each cosinor component (24h, 12h, 8h). This allows for visualization of how shifting the timing of PAEE affects the predicted outcomes.
* **Visualization:** Creates plots of the predicted outcome curves, stratified by sex, to illustrate the impact of shifting the PAEE profile on each outcome.
* **Model Diagnostics:** Saves the results of likelihood ratio testing to an [Excel file](Code/Tables/8_cosinorFeatureAnalysis.xlsx).
* **Outcome-Specific Models:** Adapts the GLM family (Gaussian or inverse-Gaussian) and link function (log or identity) based on the characteristics of each outcome variable.
<br> <br> 

**[9_totalPAEEAnalysis.do](Code/9_totalPAEEAnalysis.do)**

`9_totalPAEEAnalysis.do` examines the associations between total PAEE and cardiometabolic risk factors using GLMs and generates predicted outcome curves across the range of total PAEE values.

* **Data Transformation:** Transforms total PAEE from J/min/kg to kJ/hour/kg for easier interpretation.
* **GLMs:** Fits generalized linear models (GLMs) to assess the relationship between total PAEE and each outcome variable. The models include total PAEE, its quadratic term, and interactions with sex, along with baseline characteristics and body composition measures as covariates.
* **Model Evaluation:** Performs likelihood ratio tests to compare the full model (including total PAEE and its interactions) to a base model with only baseline covariates.
* **Predicted Outcome Curves:**  Generates predicted outcome values across a range of total PAEE values (from 20 to 100 kJ/hour/kg) to visualize the relationship between total PAEE and each outcome.
* **Visualization:**  Creates plots of the predicted outcome curves, stratified by sex, to illustrate the impact of different total PAEE levels on the outcomes.
* **Model Diagnostics:** Saves the results of likelihood ratio testing to an [Excel file](Code/Tables/9_totalPAEEAnalysis.xlsx).
* **Outcome-Specific Models:** Adapts the GLM family (Gaussian or inverse-Gaussian) and link function (log or identity) based on the characteristics of each outcome variable.
<br> <br>

**[10_panelPlots.do](Code/10_panelPlots.do)**

`10_panelPlots.do` creates panel plots by combining the individual outcome plots generated in the `8_cosinorFeatureAnalysis.do` and `9_totalPAEEAnalysis.do` scripts.

* **Plot Combination:** Combines the individual graphs produced in the previous analysis steps (`cosinorFeatureAnalysis` and `totalPAEEAnalysis`) for each outcome variable, cluster group (including the pooled sample), and model adjustment level.
* **Panel Layout:** Arranges the combined plots into a 2x5 grid (2 rows, 5 columns), creating a panel figure for each combination of cluster group, analysis type, and model level.
* **Output:** Exports the panel plots as PNG image files (`Code/Figures/<group>_<analysis>_<model>.png`).
<br> <br>

## :bar_chart: Data Availability

The Fenland study dataset is not publicly available due to participant confidentiality and data sharing agreements. However, researchers can apply for access to the data through the Medical Research Council (MRC) Epidemiology Unit at the University of Cambridge. 

For more information on the data access process and to submit a data request, please visit the MRC Epidemiology Unit's data sharing webpage: [http://www.mrc‐epid.cam.ac.uk/research/data‐sharing/](http://www.mrc‐epid.cam.ac.uk/research/data‐sharing/)

## :raising_hand_woman: :raising_hand_man: Contributing

This repository primarily serves as a resource for reproducing the analysis presented in the "Fenland Diurnal PA and Metabolic Risk" manuscript. However, we welcome feedback and suggestions for improvements! If you encounter any issues with the code or have suggestions for enhancements, please open an issue on the GitHub repository.

## :telephone_receiver: Contact

If you have any questions, or if you are interested in collaborating on future research related to this project, please contact:

  * **Dr Soren Brage, PhD**
    * **Affiliation:** MRC Epidemiology Unit, University of Cambridge, UK
    * **Email:** soren.brage@mrc-epid.cam.ac.uk 
    * **ORCID:** [https://orcid.org/0000-0002-1265-7355](https://orcid.org/0000-0002-1265-7355)
      
  * **Dr Philip Lewis, PhD**
    * **Affiliation:** Institute & Policlinic for Occupational Medicine, University Hospital of Cologne, Germany 
    * **Email:** philip.lewis@uk-koeln.de
    * **ORCID:** [https://orcid.org/0000-0002-5881-6032](https://orcid.org/0000-0002-5881-6032)
    
  * **Dr Tomas Gonzales, PhD**
    * **Affiliation:** MRC Epidemiology Unit, University of Cambridge, UK
    * **Email:** tomas.gonzales@mrc-epid.cam.ac.uk 
    * **ORCID:** [https://orcid.org/0000-0003-0085-8771](https://orcid.org/0000-0003-0085-8771)

## :balance_scale: License

The code contained in this repository is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
