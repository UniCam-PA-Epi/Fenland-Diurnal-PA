//Fenland code



***************PA data

// PAEE: p_paee_consolidated

// Net value of MET per minute: stdMET_highIC_Branch 



******Exclude if <72h or quads of 8h

/*
- Participants were excluded if they had worn their sensors for < 72 h overall, or had not worn their sensors for a combined total 
of at least 8 h throughout each quadrant of the day. Day quadrants were defined as time blocks between 3 am to 9 am, 9 am to 3 pm, 
3 pm to 9 pm, and 9 pm to 3 am, thus expecting 36 h in each of these for a perfectly balanced and fully compliant 6-day wear. 
*/

local pwear 72
local pwear_quad 8

gen include_sensor = .
replace include_sensor = 1 if p_pwear_consolidated >= `pwear' & p_pwear_consolidated != .
replace include_sensor = 1 if p_pwear_morning >= `pwear_quad' & p_pwear_noon >= `pwear_quad' & p_pwear_afternoon >= `pwear_quad' & p_pwear_night >= `pwear_quad' & q_int_stdmet_highic_branch_cat0 != .

codebook  include_sensor 

drop if include_sensor == .


**********************Data Cleaning P1

***drop pregnant women

clonevar preg = gq_preg_der
		recode preg (-7 = .) 
		recode preg (-1 = .)
	label var preg "Pregnant; 2=Yes, 3=No"
	label define preg_lbl1 2 "Yes" 3 "No" 
	label values preg preg_lbl1    
	tab preg

drop if preg == 2

drop gq_preg_der



**-------------------------
**COVARIATES
**-------------------------



***SEX

gen sex = 0 if e_sex =="F" 
	replace sex = 1 if e_sex == "M" 
	label var sex "0 Women, 1 Men"
	label define SEX_lbl 0 "Women" 1 "Men"
	label values sex SEX_lbl
	tab sex
drop e_sex


**AGE
clonevar age_p1 = g_ageattest_attended  
	label variable age_p1 "Age, years"
drop age g_ageattest_attended
	
	**Age categories
	egen age_cat_p1 = cut(age_p1), at(20 40 45 50 55 70)
	recode age_cat_p1 (20 = 1) (40 = 2) (45 = 3) (50 = 4) (55 = 5) (70 = 6)
	lab define age2 1 "< 40" 2 "40-45" 3 "45-50" 4 "50-55" 5 "> 55"
	lab val age_cat_p1 age2
	label var age_cat_p1 "AGE"

**Weight
clonevar weight_p1 = e_weight 
	replace weight_p1 = . if e_weight == -9
drop e_weight

	
**Waist circumference 
clonevar waist_p1 = e_average_waist 
	replace waist_p1 = . if e_average_waist == -1
drop e_average_waist


**Hip circumference 
clonevar hip_p1 = e_average_hip
	replace hip_p1 = . if e_average_hip == -1
drop e_average_hip

	
**Waist:Hip ratio
gen waisthipratio_p1 = waist_p1/hip_p1  
 							 

*BMI*
gen bmi_p1 = .
	replace bmi_p1 = e_bmi if e_bmi!=.			
	replace bmi_p1 = . if e_bmi == -1 	
		label var bmi_p1 "BMI (kg/m2: E_Weight / E_Height Squared)"
drop e_bmi

	//CATEGORIES
	gen bmi_cat_p1 = .
	replace bmi_cat_p1 = 3 if bmi_p1 >= 30 & bmi_p1 !=. 
	replace bmi_cat_p1 = 2 if bmi_p1 <30 & bmi_p1 !=. 
	replace bmi_cat_p1 = 1 if bmi_p1 <25 & bmi_p1 !=. 
	sum bmi_cat_p1, detail 
	label define cat_BMI_label 1 "< 25" 2 "25 - 30" 3 "> 30" 
	label values bmi_cat_p1 cat_BMI_label 
	label var bmi_cat_p1 "BMI"

**Worktype
gen worktype_p1 = 0
	replace worktype_p1 = 1 if worktype_clean == 1
	replace worktype_p1 = 2 if worktype_clean == 2
	replace worktype_p1 = 3 if worktype_clean == 3
	replace worktype_p1 = 3 if worktype_clean == 4
	replace worktype_p1 = 6 if worktype_clean == -1 | worktype_clean == . 
	replace worktype_p1 = 4 if worktype_clean == -1 & gq_job_retired_der == 2
	replace worktype_p1 = 5 if worktype_clean == -1 & gq_job_unempl_der == 2 

	label define WORKTYPE_label 1 "Sedentary" 2 "Standing" 3 "Manual work" 4 "Retired" 5 "Unemployed" 6 "Unknown"
	label values worktype_p1 WORKTYPE_label 
	label var worktype_p1 "WORK TYPE"
	tab worktype_p1



**Marital status -- Soren's suggested method 
gen marital_p1 = gq_marit_der
	recode marital_p1 (5=3) (4=3) (-8=6) (-1=6) (.=6)
	label var marital_p1 "What is your marital status?"
	label define marital_label 1 "Single" 2 "Married/living as married" 3 "Widowed/separated/divorced" 6 "Unknown" 
	label values marital_p1 marital_label
	tab marital_p1
//drop gq_marit gq_marit_der	
	
		*TRUE-MISSING METHOD
		clonevar marital2_p1 = marital_p1
			recode marital2_p1 (6 = .)

**Education
	gen education_p1 = 0
	replace education_p1 = 1 if gq_edu_none_der ==2 | gq_edu_slc_der ==2 | gq_edu_cse_der ==2 | gq_edu_gcse_der ==2
	replace education_p1 = 2 if gq_edu_matricul_der ==2 | gq_edu_alev_der ==2 | gq_edu_techcg_der ==2 | gq_edu_apprentice_der ==2 | gq_edu_hnc_der ==2 | gq_edu_hnd_nvq_der ==2 | gq_edu_secret_der ==2 | gq_edu_trade_der ==2
	replace education_p1 = 3 if gq_edu_degree_der ==2 
	replace education_p1 = 1 if education_p1 == 0 
	label var education_p1 "EDUCATION"
	label define EDUCATION_label 1 "Compulsory" 2 "Further" 3 "Higher"
	label values education_p1 EDUCATION_label 
	tab education_p1
//drop gq_edu_none_der gq_edu_slc_der gq_edu_cse_der gq_edu_gcse_der gq_edu_matricul_der gq_edu_alev_der gq_edu_techcg_der gq_edu_apprentice_der gq_edu_hnc_der gq_edu_hnd_nvq_der gq_edu_secret_der gq_edu_trade_der gq_edu_degree_der
	

**Income -- Soren's suggested method 
clonevar income_p1 = gq_income_der
	recode income_p1 (-7 = -1)
	recode income_p1 (. = -1)
	recode income_p1 (-1 = 4)
	label var income_p1 "INCOME"
	label define income_p1_label 1 "< £20000" 2 "£20000 - £40000" 3 ">£40000" 4 "NA"
	label values income_p1 income_p1_label  
	tab income_p1
//drop gq_income_der

		*TRUE-MISSING METHOD
		clonevar income2_p1 = income_p1
			recode income2_p1 (4 = .)

**Smoking -- Soren's suggested method 
gen smoking_p1 = g_smoke 
	recode smoking_p1 (3 = 5)
	recode smoking_p1 (9 = 5)
	recode smoking_p1 (1 = 5)
	recode smoking_p1 (. = 5)
	label var smoking_p1 "SMOKING"
	label define smoking_p1_label 0 "Never smoked" 2 "Ex smoker" 4 "Current smoker" 5 "Unknown"
	label values smoking_p1 smoking_p1_label
	tab smoking_p1
//drop g_smoke

		*TRUE-MISSING METHOD
		clonevar smoking2_p1 = smoking_p1
			recode smoking2_p1 (5 = .)	

			
**Alcohol  - Standard drinks/wk   (10g = 1 standard drink)
gen alcohol_unitwk_p1 = (alcohol / 8)*7
label var alcohol_unitwk_p1 "ALCOHOL"



**Ethnic origin 
clonevar ethnic = gq_eth_5c_der
	recode ethnic (-8 = -1) 
	recode ethnic (. = -1) 
	recode ethnic (-1 = 6)
	label var ethnic "1= White; 2= South Asian; 3= Black; 4= East Asian; 5= Others; 6= NA"
	label define ETHNIC_label 1 "White" 2 "South Asian" 3 "Black" 4 "East Asian" 5 "Others" 6 "NA" 
	label values ethnic ETHNIC_label
	tab ethnic 
*drop GQ_eth gq_eth_DER gq_eth_5c_DER

		*TRUE-MISSING METHOD
		clonevar ethnic2 = ethnic 
			recode ethnic2 (6 = .)

**DIET
clonevar diet_p1 = mdspyr_v2


**Fam history of DM

destring fh04_motherdiabetic, gen(MotherDM) force
	recode MotherDM (-1 = .) (4 = .)
	label define motherDM 2 "Yes" 3 "No"
	label values MotherDM motherDM
	
destring fh06_fatherdiabetic, gen(FatherDM) force
	recode FatherDM (-1 = .) (4 = .)
	label define fatherDM 2 "Yes" 3 "No"
	label values FatherDM fatherDM	



**-------------------------
**OUTCOMES 
**-------------------------



// systolic blood pressure  ( BPSys1, BPSys2, BPSys3) - BPSysAVG //average of all available values - BUT, some only have 1 value -- NEEDED TO FIX TO AT LEAST x 2
egen SBP_p1 = rowmean(bpsys?) if (missing(bpsys1) + missing(bpsys2) + missing(bpsys3)) <= 1
replace SBP_p1 =round(SBP_p1, 0.1)

// diastolic blood pressure  ( BPDia1, BPDia2, BPDia3) - BPDiaAVG //average of all available values (7 missing) - BUT, some only have 1 value -- NEEDED TO FIX TO AT LEAST x 2
egen DBP_p1 = rowmean(bpdia?) if (missing(bpdia1) + missing(bpdia2) + missing(bpdia3)) <= 1
replace DBP_p1 =round(DBP_p1, 0.1)


*** Blood biochemistry markers ***

// Fasting and PP glucose 
gen glucose_0_p1 = .
	replace glucose_0_p1 = glucose0 if glucose0!=.				
	replace glucose_0_p1 = repeatglucose0 if glucose0==.
		label var glucose_0_p1 "Serum fasting Glucose (mmol/L)"

gen glucose_120_p1 = .
	replace glucose_120_p1 = glucose120 if glucose120!=.		
	replace glucose_120_p1 = repeatglucose120 if glucose120==.
		label var glucose_120_p1 "2 hour Glucose measurement (mmol/L)"
		
		
// Triglycerides
gen TG_p1 = .
	replace TG_p1 = triglyceride0 if triglyceride0!=.			
		label var TG_p1 "Serum fasting Triglyceride (mmol/L)"	

// Total Cholesterol
gen Chol_p1 = .
	replace Chol_p1 = chol0 if chol0!=.						
		label var Chol_p1 "Serum fasting Cholesterol (mmol/L)"	

// Fasting HDL-cholesterol 
gen HDL_p1 = .
	replace HDL_p1 = hdl0 if hdl0!=.								
		label var HDL_p1 "Serum fasting HDL (mmol/L)"	

// Fasting LDL-cholesterol 
gen LDL_p1 = .
	replace LDL_p1 = ldl0 if ldl0!=.								
		label var LDL_p1 "Serum fasting LDL (mmol/L)"	

// TOTAL CHOL / HDL ratio  (lower better)
gen CholHDL_p1 = .
	replace CholHDL_p1 = chol_hdlcholratio0 if chol_hdlcholratio0!=.		
		label var CholHDL_p1 "Chol to HDL ratio (mmol/L)"

// LDL / HDL ratio (lower better)
gen LDLHDL_p1 = LDL_p1 / HDL_p1
	label var LDLHDL_p1 "LDL to HDL ratio (mmol/L)"	
	
	
// HbA1c 
		** Hba1C0_PC  // 4k missing due to different DCCT vs IFCC output
		** Hba1C0_MM  // 7k missing due to different DCCT vs IFCC output

** Convert % to MM/MOL
	gen hba1c0_pc_convert = (hba1c0_mm / 10.929) + 2.15
	replace hba1c0_pc_convert=round(hba1c0_pc_convert, 0.1)
	
//Using this
** Convert MM/MOL to %
	gen hba1c0_mm_convert = (hba1c0_pc - 2.15) * 10.929
	replace hba1c0_mm_convert=round(hba1c0_mm_convert, 1)
	
**COMBINE
	clonevar HBA1C_PERC = hba1c0_pc
	replace HBA1C_PERC = hba1c0_pc_convert if HBA1C_PERC==.

//and this
	clonevar HbA1c_p1 = hba1c0_mm
	replace HbA1c_p1 = hba1c0_mm_convert if HbA1c_p1 == .
	




**********************Data Cleaning P1P2

***drop pregnant women

clonevar preg_p1 = gq_preg_der
		recode preg_p1 (-7 = .) 
		recode preg_p1 (-1 = .) 
	label var preg_p1 "Pregnant; 2=Yes, 3=No"
	label define preg_lbl 2 "Yes" 3 "No" 
	label values preg_p1 preg_lbl    
	tab preg_p1


gen preg_p2 = .
replace preg_p2 = 2 if gq_a3kpregnant_p2 == "2"
replace preg_p2 = 3 if gq_a3kpregnant_p2 == "3"
label var preg_p2 "Pregnant; 2=Yes, 3=No"
label define preg_lbl2 2 "Yes" 3 "No" 
label values preg_p2 preg_lbl2
tab preg_p2

drop if preg_p1 == 2
drop if preg_p2 == 2

drop gq_preg_der gq_a3kpregnant_p2



**-------------------------
**COVARIATES
**-------------------------

***AGE
rename p_age_p2 age_p2


**Weight
clonevar weight_p2 = cl_weight_p2 
drop cl_weight_p2 cl_weight_code_p2 p_weight_p2 p_weighted_obs_p2

**Waist circumference 
clonevar waist_p2 = e_average_waist_p2 
drop e_waist1_p2 e_waist2_p2 e_waist3_p2 e_average_waist_p2 


**Hip circumference 
clonevar hip_p2 = e_average_hip_p2
	replace hip_p2 = . if e_average_hip_p2 == -1
drop e_average_hip_p2 e_hip1_p2 e_hip2_p2 e_hip3_p2


**Waist:Hip ratio
gen waisthipratio_p2 = waist_p2/hip_p2


*BMI*
gen bmi_p2 = .
	replace bmi_p2 = cl_bmi_p2 if cl_bmi_p2 != .			
	replace bmi_p2 = . if cl_bmi_p2 == -1 	
		label var bmi_p2 "BMI (kg/m2: E_Weight / E_Height Squared)"
drop cl_bmi_p2

	//CATEGORIES
	gen bmi_cat_p2 = .
	replace bmi_cat_p2 = 3 if bmi_p2 >= 30 & bmi_p2 !=. 
	replace bmi_cat_p2 = 2 if bmi_p2 <30 & bmi_p2 !=. 
	replace bmi_cat_p2 = 1 if bmi_p2 <25 & bmi_p2 !=. 
	//label define cat_BMI_label 1 "< 25" 2 "25 - 30" 3 "> 30" 
	label values bmi_cat_p2 cat_BMI_label 



**Worktype
gen worktype_p2 = .
	replace worktype_p2 = 1 if worktype_clean_p2 == "1"
	replace worktype_p2 = 2 if worktype_clean_p2 == "2"
	replace worktype_p2 = 3 if worktype_clean_p2 == "3"
	replace worktype_p2 = 3 if worktype_clean_p2 == "4"
	replace worktype_p2 = 4 if worktype_p2 == . & gq_b1dretired_p2 == "2"
	replace worktype_p2 = 5 if worktype_p2 == . & gq_b1funemployed_p2 == "2" 

	label define WORKTYPE_label2 1 "Sedentary" 2 "Standing" 3 "Manual work" 4 "Retired" 5 "Unemployed" 
	label values worktype_p2 WORKTYPE_label2
	tab worktype_p2


**Marital status
			
gen marital_p2 = .
	replace marital_p2 = 1 if gq_b13maritalstatus_p2 == "1"
	replace marital_p2 = 2 if gq_b13maritalstatus_p2 == "2"
	replace marital_p2 = 3 if gq_b13maritalstatus_p2 == "3"
	replace marital_p2 = 3 if gq_b13maritalstatus_p2 == "4"
	replace marital_p2 = 3 if gq_b13maritalstatus_p2 == "5"
	label var marital_p2 "What is your marital status?"
	label define marital_label2 1 "Single" 2 "Married/living as married" 3 "Widowed/separated/divorced" 
	label values marital_p2 marital_label2
	tab marital_p2


**Education

gen education_p2 = 0
	replace education_p2 = 1 if gq_b6nnone_p2 =="2" | gq_b6aslc_p2 ==2 | gq_b6bcse_p2 ==2 | gq_b6cgcse_p2 ==2
	replace education_p2 = 2 if gq_b6dmatriculation_p2 ==2 | gq_b6ealevels_p2 ==2 | gq_b6fc_g_p2 ==2 | gq_b6happrenticeship_p2 ==2 | gq_b6jhnc_p2 ==2 | gq_b6ghnd_nvq_p2 ==2 | gq_b6isecretarial_p2 ==2 | gq_b6ltradecertificate_p2 ==2
	replace education_p2 = 3 if gq_b6kdegree_p2 =="2"
	replace education_p2 = 1 if education_p2 == 0 
	label var education_p2 "1=compulsory 2=further 3=higher"
	label define EDUCATION_label2 1 "Basic" 2 "Further" 3 "Higher"
	label values education_p2 EDUCATION_label2 
	tab education_p2

**Income
gen income_p2 = .
	replace income_p2 = 1 if gq_b8householdincome_p2 == "1"
	replace income_p2 = 2 if gq_b8householdincome_p2 == "2"
	replace income_p2 = 3 if gq_b8householdincome_p2 == "3"
	replace income_p2 = 4 if gq_b8householdincome_p2 == "4"
	replace income_p2 = 5 if gq_b8householdincome_p2 == "5"
	replace income_p2 = 6 if income_p2 == .
	label var income_p2 "Please can you indicate what your household income is?"
	label define income_p2_label 1 "< £20000" 2 "£20000 - £40000" 3 "£40000 - £60000" 4 "£60000 - £80000" 5 "£80000 and above" 6 "NA"
	label values income_p2 income_p2_label  
	tab income_p2
*drop gq_b8householdincome_p2 

		*TRUE-MISSING METHOD
		clonevar income2_p2 = income_p2
			recode income2_p2 (6 = .)


**Smoking -- Soren's suggested method 
clonevar smoking_p2 = g_smoke_p2 
	recode smoking_p2 (3 = 5)
	recode smoking_p2 (9 = 5)
	recode smoking_p2 (1 = 5)
	recode smoking_p2 (. = 5)
	label var smoking_p2 "Current Smoking Consumption History / Status"
	label define smoking_p2_label 0 "Never smoked" 2 "Ex smoker" 4 "Current smoker" 5 "Unknown"
	label values smoking_p2 smoking_p2_label
	tab smoking_p2
//drop g_smoke_p2

		*TRUE-MISSING METHOD
		clonevar smoking2_p2 = smoking_p2
			recode smoking2_p2 (5 = .)	


**Alcohol  - Standard drinks/wk   (10g = 1 standard drink)
gen alcohol_unitwk_p2 = (alcohol_v2_p2 / 8)*7


**DIET
clonevar diet_p2 = mdspyr_v2_p2



**-------------------------
**OUTCOMES 
**-------------------------



// systolic blood pressure  ( BPSys1, BPSys2, BPSys3) - BPSysAVG //average of all available values - BUT, some only have 1 value -- NEEDED TO FIX TO AT LEAST x 2
egen SBP_p2 = rowmean(bpsys?_p2) if (missing(bpsys1_p2) + missing(bpsys2_p2) + missing(bpsys3_p2)) <= 1
replace SBP_p2 = round(SBP_p2, 0.1)

// diastolic blood pressure  ( BPDia1, BPDia2, BPDia3) - BPDiaAVG //average of all available values (7 missing) - BUT, some only have 1 value -- NEEDED TO FIX TO AT LEAST x 2
egen DBP_p2 = rowmean(bpdia?_p2) if (missing(bpdia1_p2) + missing(bpdia2_p2) + missing(bpdia3_p2)) <= 1
replace DBP_p2 = round(DBP_p2, 0.1)



*** Blood biochemistry markers ***

// Fasting and PP glucose 
gen glucose_0_p2 = .
	replace glucose_0_p2 = glucose0_p2 if glucose0_p2!=.				
	replace glucose_0_p2 = repeatglucose0_p2 if glucose0_p2==.
		label var glucose_0_p2 "Serum fasting Glucose (mmol/L)"

gen glucose_120_p2 = .
	replace glucose_120_p2 = glucose120_p2 if glucose120_p2!=.		
	replace glucose_120_p2 = repeatglucose120_p2 if glucose120_p2==.
		label var glucose_120_p2 "2 hour Glucose measurement (mmol/L)"
	
		
// Triglycerides
gen TG_p2 = .
	replace TG_p2 = triglyceride0_p2 if triglyceride0_p2!=.
		label var TG_p2 "Serum fasting Triglyceride (mmol/L)"	

// Total Cholesterol
gen Chol_p2 = .
	replace Chol_p2 = chol0_p2 if chol0_p2!=.						
		label var Chol_p2 "Serum fasting Cholesterol (mmol/L)"	

// Fasting HDL-cholesterol 
gen HDL_p2 = .
	replace HDL_p2 = hdl0_p2 if hdl0_p2!=.								
		label var HDL_p2 "Serum fasting HDL (mmol/L)"	

// Fasting LDL-cholesterol 
gen LDL_p2 = .
	replace LDL_p2 = ldl0_p2 if ldl0_p2!=.								
		label var LDL_p2 "Serum fasting LDL (mmol/L)"	

// TOTAL CHOL / HDL ratio  (lower better)
gen CholHDL_p2 = .
	replace CholHDL_p2 = chol_hdlcholratio0_p2 if chol_hdlcholratio0_p2!=.
		label var CholHDL_p2 "Chol to HDL ratio (mmol/L)"

// LDL / HDL ratio (lower better)
gen LDLHDL_p2 = LDL_p2 / HDL_p2
	label var LDLHDL_p2 "LDL to HDL ratio (mmol/L)"		

	
// HbA1c 
gen HbA1c_p2 = hba1c0_mm_p2
	replace HbA1c_p2 = round(HbA1c_p2, 0.1)
*drop hba1c0_mm_p2 hba1c0_pc_p2





**************CMRS


gen BP_p1 = (SBP_p1 + DBP_p1) / 2
gen BP_p2 = (SBP_p2 + DBP_p2) / 2


*************Paddy's code for CMRS


//normalise when the data is skewed, not all of them are, so double check and fix the rest of the code accordingly 

local variables SBP_p1 DBP_p1 BP_p1 glucose_0_p1 TG_p1 HDL_p1 waist_p1 SBP_p2 DBP_p2 BP_p2 glucose_0_p2 TG_p2 HDL_p2 waist_p2 CholHDL_p1 CholHDL_p2 HbA1c_p1 HbA1c_p2

foreach var of local variables {
    
	gen log10_`var' = log10(`var')
}


//standardisation
foreach v of varlist log10_SBP_p1 log10_DBP_p1 log10_BP_p1 log10_glucose_0_p1 log10_TG_p1 log10_HDL_p1 log10_waist_p1 log10_SBP_p2 log10_DBP_p2 log10_BP_p2 log10_glucose_0_p2 log10_TG_p2 log10_HDL_p2 log10_waist_p2 log10_CholHDL_p1 log10_CholHDL_p2 log10_HbA1c_p1 log10_HbA1c_p2 { 
	
	egen z_`v'_W = std(`v') if sex == 0
	}

foreach v of varlist log10_SBP_p1 log10_DBP_p1 log10_BP_p1 log10_glucose_0_p1 log10_TG_p1 log10_HDL_p1 log10_waist_p1 log10_SBP_p2 log10_DBP_p2 log10_BP_p2 log10_glucose_0_p2 log10_TG_p2 log10_HDL_p2 log10_waist_p2 log10_CholHDL_p1 log10_CholHDL_p2 log10_HbA1c_p1 log10_HbA1c_p2 {
	
    egen z_`v'_M = std(`v') if sex == 1
}

foreach v of varlist log10_SBP_p1 log10_DBP_p1 log10_BP_p1 log10_glucose_0_p1 log10_TG_p1 log10_HDL_p1 log10_waist_p1 log10_SBP_p2 log10_DBP_p2 log10_BP_p2 log10_glucose_0_p2 log10_TG_p2 log10_HDL_p2 log10_waist_p2 log10_CholHDL_p1 log10_CholHDL_p2 log10_HbA1c_p1 log10_HbA1c_p2 {
	
    gen z_`v' = .
}

foreach v of varlist log10_SBP_p1 log10_DBP_p1 log10_BP_p1 log10_glucose_0_p1 log10_TG_p1 log10_HDL_p1 log10_waist_p1 log10_SBP_p2 log10_DBP_p2 log10_BP_p2 log10_glucose_0_p2 log10_TG_p2 log10_HDL_p2 log10_waist_p2 log10_CholHDL_p1 log10_CholHDL_p2 log10_HbA1c_p1 log10_HbA1c_p2 {
	
    replace z_`v' = z_`v'_W if sex == 0
}

foreach v of varlist log10_SBP_p1 log10_DBP_p1 log10_BP_p1 log10_glucose_0_p1 log10_TG_p1 log10_HDL_p1 log10_waist_p1 log10_SBP_p2 log10_DBP_p2 log10_BP_p2 log10_glucose_0_p2 log10_TG_p2 log10_HDL_p2 log10_waist_p2 log10_CholHDL_p1 log10_CholHDL_p2 log10_HbA1c_p1 log10_HbA1c_p2 {
	
    replace z_`v' = z_`v'_M if sex == 1
}


//final calculation

********USING HbA1c (better not according to Soren)

gen CMRS_p1_M = (z_log10_BP_p1_M + z_log10_TG_p1_M + z_log10_CholHDL_p1_M + z_log10_HbA1c_p1_M + z_log10_waist_p1_M)
gen CMRS_p1_W = (z_log10_BP_p1_W + z_log10_TG_p1_W + z_log10_CholHDL_p1_W + z_log10_HbA1c_p1_W + z_log10_waist_p1_W)
gen CMRS_p1 = (z_log10_BP_p1 + z_log10_TG_p1 + z_log10_CholHDL_p1 + z_log10_HbA1c_p1 + z_log10_waist_p1)

gen CMRS_p2_M = (z_log10_BP_p2_M + z_log10_TG_p2_M + z_log10_CholHDL_p2_M + z_log10_HbA1c_p2_M + z_log10_waist_p2_M)
gen CMRS_p2_W = (z_log10_BP_p2_W + z_log10_TG_p2_W + z_log10_CholHDL_p2_W + z_log10_HbA1c_p2_W + z_log10_waist_p2_W)
gen CMRS_p2 = (z_log10_BP_p2 + z_log10_TG_p2 + z_log10_CholHDL_p2 + z_log10_HbA1c_p2 + z_log10_waist_p2)


********USING 2h Glucose (to use)

gen CMRS_p1_M = (z_BP_p1_M + z_log10_TG_p1_M + z_log10_CholHDL_p1_M + z_log10_glucose_120_p1_M + z_log10_waist_p1_M)
gen CMRS_p1_W = (z_BP_p1_W + z_log10_TG_p1_W + z_log10_CholHDL_p1_W + z_log10_glucose_120_p1_W + z_log10_waist_p1_W)
gen CMRS_p1 = (z_BP_p1 + z_log10_TG_p1 + z_log10_CholHDL_p1 + z_log10_glucose_120_p1 + z_log10_waist_p1)

gen CMRS_p2_M = (z_BP_p2_M + z_log10_TG_p2_M + z_log10_CholHDL_p2_M + z_log10_glucose_120_p2_M + z_log10_waist_p2_M)
gen CMRS_p2_W = (z_BP_p2_W + z_log10_TG_p2_W + z_log10_CholHDL_p2_W + z_log10_glucose_120_p2_W + z_log10_waist_p2_W)
gen CMRS_p2 = (z_BP_p2 + z_log10_TG_p2 + z_log10_CholHDL_p2 + z_log10_glucose_120_p2 + z_log10_waist_p2)

