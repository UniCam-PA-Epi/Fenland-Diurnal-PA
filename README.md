# Fenland-Diurnal-PA
This repository houses code for the "Fenland Diurnal PA and Metabolic Risk" manuscript.

## Table of Contents
* [Usage](#usage)
* [Code File Descriptions](#code-file-descriptions)
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

**1. `1_initialiseCovariates.do`**

This do-file performs initial data preparation for the main covariates used in the analysis, including:

* **Variable Renaming:**  Simplifies variable names from the Fenland dataset to improve clarity and consistency within the analysis code. 
* **Categorical Variable Creation:**  Generates categorical variables for sociodemographic factors (e.g., age groups, education level) and anthropometric measures (e.g., BMI categories).

**2. `2_initialiseOutcomes.do`**

This do-file prepares the outcome variables for analysis:

* **Variable Renaming:**  As in the previous file, this step renames outcome variables from the Fenland dataset for simplicity.
* **Outcome Variable Definitions:**  Creates the key outcome variables used in the analysis. These include:
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

**3. `3_applyCosinorModel.do`**

This do-file applies a cosinor model to analyse diurnal patterns in physical activity energy expenditure (PAEE) collected over approximately 6 days of continuous monitoring.

* **Cosinor Analysis:**  Utilises multiple superimposed cosinor models to analyze the rhythmic patterns of PAEE throughout a 24-hour period. This involves fitting a gamma Generalized Linear Model (GLM) with a log link function to estimate parameters of the cyclical pattern, including:
    * **Mesor:** The average level of PAEE over the 24-hour cycle.
    * **Amplitude:** The difference between the peak PAEE and the mesor, representing the strength of the diurnal rhythm.
    * **Acrophase:** The time of day at which the peak PAEE occurs.
* **Multi-frequency Rhythms:** The analysis explores PAEE rhythms at multiple frequencies by jointly fitting cosinor models to 24-hour, 12-hour, and 8-hour cycles. This allows for the identification of potential ultradian rhythms in addition to the main circadian rhythm.
* **PAEE Metrics:**  In addition to the cosinor parameters, the code estimates:
    * **Maximum Achieved PAEE:** The highest PAEE value within the 24-hour period.
    * **Hour of Maximum PAEE:** The hour of the day when PAEE is highest.
    * **Total 24-hour PAEE:** The average total PAEE over a 24-hour period.
* **Outlier Removal:** An IQR-based method to remove outliers in the PAEE data is applied before fitting the cosinor model. 


* ## License

The code contained in this repository is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
