# Fenland-Diurnal-PA

This repository provides the Stata code used to conduct the analysis for the manuscript "Fenland Diurnal PA and Metabolic Risk."  The research examines how total physical activity energy expenditure and variations in physical activity throughout the day relate to metabolic health outcomes, such as blood glucose, insulin levels, and blood pressure, in participants from the Fenland cohort.

## Table of Contents
* [Usage](#usage)
* [Code File Descriptions](#code-file-descriptions)
* [Data Availability](#data-availability)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)

## Usage

This repository contains Stata code to reproduce the analysis for the "Fenland Diurnal PA and Metabolic Risk" manuscript. The main analysis is conducted using the `0_Main.do` file. 

**To run the code:**

1. **Prerequisites:** Ensure you have Stata (at least version 17.0) installed on your system and have access to the Fenland study dataset.
2. **Clone the repository:** `git clone https://github.com/UniCam-PA-Epi/Fenland-Diurnal-PA.git`
3. **Navigate to the project code directory:** `cd Fenland-Diurnal-PA/Code`
4. **Run the main do file:**  `stata do 0_Main.do "path to Fenland release"` (replace `"path to Fenland release"` with the actual file path to the Fenland study dataset)

**Important Notes:**

* The `0_Main.do` file calls other do-files within the repository, so ensure all files are kept in the same directory.
* The code expects the Fenland data to be in a specific format. Please refer to the manuscript or any accompanying documentation for details on data preparation.
* PAEE data should have already been pre-processed and averaged to 24-hour periods.
* The code generates output files (e.g., tables, figures) in project code directory.

## Code File Descriptions

This section provides a brief overview of the functionality of each code file in the analysis sequence.
<br> <br> 

**[0_Main.do](Code/0_Main.do)**

`0_Main.do` executes the complete analysis pipeline for the "Fenland Diurnal PA and Metabolic Risk" manuscript.
<br> <br> 

**[1_initialiseCovariates.do](Code/1_initialiseCovariates.do)**

`1_initialiseCovariates.do` performs initial data preparation for the main covariates used in the analysis, including:

* **Variable Renaming:**  Simplifies variable names from the Fenland dataset to improve clarity and consistency within the analysis code. 
* **Categorical Variable Creation:**  Generates categorical variables for sociodemographic factors (e.g., age groups, education level) and anthropometric measures (e.g., BMI categories).
<br> <br>

**[2_initialiseOutcomes.do](Code/2_initialiseOutcomes.do)**

`2_initialiseOutcomes.do` prepares the outcome variables for analysis:

* **Variable Renaming:**  As in the previous file, this step renames outcome variables from the Fenland dataset for simplicity.
* **Outcome Variable Definitions:**  Creates the outcome variables used in the analysis. These include:
    * **Metabolic Measures:** 
        * `glucose120`: 2-hour glucose 
        * `insulin`: Fasting insulin 
        * `leptin`: Blood leptin 
        * `nefa`: Non-esterified fatty acids 
        * `adiponectin`: Blood adiponectin 
        * `ldl`: Low-density lipoprotein cholesterol 
        * `hdl`: High-density lipoprotein cholesterol 
    * **Cardiovascular Measures:**
        * `mbpsys`: Systolic blood pressure 
        * `mbpdia`: Diastolic blood pressure 
        * `crp`: C-reactive protein     
<br> 

**[3_applyCosinorModel.do](Code/3_applyCosinorModel.do)**

`3_applyCosinorModel.do` applies a cosinor model to analyse diurnal patterns in physical activity energy expenditure (PAEE) collected over approximately 6 days of continuous monitoring.

* **Cosinor Analysis:**  Utilises multiple superimposed cosinor models to analyze the rhythmic patterns of PAEE throughout a 24-hour period. This involves fitting a gamma Generalized Linear Model (GLM) with a log link function to estimate parameters of the cyclical pattern, including:
    * **Mesor:** The average level of PAEE over the 24-hour cycle.
    * **Amplitude:** The difference between the peak PAEE and the mesor, representing the strength of the diurnal rhythm.
    * **Acrophase:** The time of day at which the peak PAEE occurs.
* **Multi-frequency Rhythms:** The analysis estimates PAEE rhythms at multiple frequencies by jointly fitting cosinor models to 24-hour, 12-hour, and 8-hour cycles. This allows for the quantification ultradian rhythms in addition to the main 24-hour circadian rhythm.
* **PAEE Metrics:**  In addition to the cosinor parameters, the code estimates:
    * **Maximum Achieved PAEE:** The highest PAEE value within the 24-hour period.
    * **Hour of Maximum PAEE:** The hour of the day when PAEE is highest.
    * **Total 24-hour PAEE:** The average total PAEE over a 24-hour period.
* **Outlier Removal:** An IQR-based method to remove outliers in the PAEE data is applied before fitting the cosinor model. 
<br> <br>

**[4_applyExclusions.do](Code/4_applyExclusions.do)**

`4_applyExclusions.do` applies exclusion criteria to the study sample for the main analysis. Participants are excluded based on the following criteria:

* **Cosinor Model Fit:**  Removes participants for whom the cosinor model did not adequately fit the PAEE data. This is determined by examining the p-values associated with the sine and cosine terms in the 24-hour, 12-hour, and 8-hour cycles.
* **Insufficient Wear Time:**  Excludes participants with insufficient wear time of the accelerometer device.
* **Missing Body Composition Data:**  Removes participants without measurements for fat mass and fat-free mass.

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
* **Optimal k Determination:**  Determines the optimal number of clusters (k) using the elbow method with within-cluster sum of squares (WCSS) and linear splines regression.
* **Cluster Visualisation:** Generates plots of the average PAEE profiles for each cluster.
<br> <br> 

**[6_descriptivesTables.do](Code/6_descriptivesTables.do)**

**[7_violinPlots.do](Code/7_violinPlots.do)**

**[8_cosinorFeatureAnalysis.do](Code/8_cosinorFeatureAnalysis.do)**

**[9_totalPAEEAnalysis.do](Code/9_totalPAEEAnalysis.do)**

## Data Availability

The Fenland study dataset analysed in this research is not publicly available due to participant confidentiality and data sharing agreements. However, researchers can apply for access to the data through the Medical Research Council (MRC) Epidemiology Unit at the University of Cambridge. 

For more information on the data access process and to submit a data request, please visit the MRC Epidemiology Unit's data sharing webpage: [http://www.mrc‐epid.cam.ac.uk/research/data‐sharing/](http://www.mrc‐epid.cam.ac.uk/research/data‐sharing/)

## Contributing

This repository primarily serves as a resource for reproducing the analysis presented in the "Fenland Diurnal PA and Metabolic Risk" manuscript. However, we welcome feedback and suggestions for improvements! If you encounter any issues with the code or have suggestions for enhancements, please open an issue on the GitHub repository.

If you are interested in collaborating on future research related to this project, please contact **Soren Brage** at soren.brage@mrc-epid.cam.ac.uk or **Tomas Gonzales** at tomas.gonzales@mrc-epid.cam.ac.uk.

## Contact

**Tomas  Gonzales, MRC Epidemiology Unit, University of Cambridge**

* **Email:** tomas.gonzales@mrc-epid.cam.ac.uk 
* **ORCID:** [https://orcid.org/0000-0003-0085-8771](https://orcid.org/0000-0003-0085-8771)

## License

The code contained in this repository is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
