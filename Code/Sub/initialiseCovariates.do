********************************************************************************
* STEP #3: Organise Covariates
*sex, age, ethnic, testsite, season, 
*marital_status, income, edeucation, work_status,
*smoking, diet, alcohol,cardiometabol_med
********************************************************************************
//Sex
label define sexlab 0 "Women" 1 "Men"
label values sex sexlab
tab sex, missing // no missing

//Age
gen age = AgeAtTest_DM_Attended // no missing
	  
//Ethnicity
gen ethnic = gq_eth_DER 
replace ethnic = 0 if ethnic==1 | ethnic==2 | ethnic==3
replace ethnic = 1 if ethnic !=0
label define ethlab 0 "white" 1 "non-white"
label values ethnic ethlab //n=587 missing (categorised as non-white)

//Test Site	  
encode TestSite, gen(testsite) // no missing data

//Season: centred around solstices and equinox
gen date = date(P_startdate, "DMY", 2022)
format date %td
gen doy=doy(date)
gen season = 0 if doy>=80-45 & doy<=80+46 // approx autumn equinox
replace season = 1 if doy>=172-45 & doy<=172+47 // approx summer solstice
replace season = 2 if doy>=266-46 & doy<=266+45 // approx autumn equinox
replace season = 3 if doy>=355-45 | doy<=80-46 // approx winter solstice
drop doy
label define seasonlab 0 "spring" 1 "summer" 2 "autumn" 3 "winter"
label values season seasonlab // no missing data

//Marital Status
gen marital_status = gq_marit_DER
label define marlab 1 "single" 2 "married/living as married" 3 "widowed" ///
4 "separated" 5 "divorced" 0 "unknown"
label values marital_status marlab
replace marital_status=0 if marital_status==-1 | marital_status==-8 | ///
marital_status==. 
// n=2642 missing (labelled as unknown)

//Income
gen income = gq_income_DER
label define incomelab 1 "<20k" 2 "20-40k" 3 ">40k" 0 "unknown"
label values income incomelab
replace income=0 if income==. | income==-7 | income==-1 
// n=318 missing (labelled as unknown)

//Education
rename EDUCATION education
label define edulab 0 "none" 1 "compulsory" 2 "further" 3 "higher"
label values education edulab
tab education, missing	// no missing data

//Work Status
encode Worktype, gen(worktype2) // note: -1=26, 1=27, 2=28, 3=29, 4=30
gen worktype=1 if worktype2== 27 | worktype2== 1 | worktype2== 2 | ///
worktype2== 5 | worktype2== 6 | worktype2== 8 | worktype2== 9 | ///
worktype2== 10 | worktype2== 17
replace worktype=2 if worktype2== 28 | worktype2== 7| worktype2== 11 | ///
worktype2== 13 | worktype2== 14 | worktype2== 18| worktype2== 19 | ///
worktype2== 20 | worktype2== 22
replace worktype=3 if worktype2== 29 | worktype2== 30| worktype2== 12 | ///
worktype2== 15| worktype2== 16 | worktype2== 21 | worktype2== 23 | ///
worktype2== 24 | worktype2== 25 
replace worktype=5 if worktype2== 26 | worktype2== . | worktype2== 3 | ///
worktype2== 4 
replace worktype= 4 if gq_job_retired_DER==2 & worktype==5
replace worktype= 4 if gq_job_unempl_DER==2 & worktype==5
replace worktype= 4 if gq_job_sick_DER==2 & worktype==5
gen work_status = worktype
label define worklab 1 "Sedentary" 2 "Standing" 3 "Manual" 4 "Not working" ///
5 "Unknown" 
label values work_status worklab
// n=637 missing (labelled as unknown)

//Smoking

// Changed to categorical var

replace SMOKE = 3 if SMOKE == 9 | SMOKE == .

tab SMOKE

#delimit ;

label 	define smoke_labels 
        0 "Never smoked"
        1 "Ex smoker"
        2 "Current smoker"
        3 "Missing/Unknown"
        ,
        replace
        ;
                
#delimit cr

label values SMOKE smoke_labels

rename SMOKE smoke

/*
gen smoking = gq_sm_pckyrs_DER
replace smoking=. if smoking==-8 | smoking==-7 | smoking==. //n=877 missing
*/



//DASH
gen diet = DASH_AccordanceScore //n=7 missing

//Alcohol
gen alcohol=gq_alc_freq_DER
replace alcohol=. if alcohol==-7 | alcohol==-1 //n=195 missing  

//Medications
gen cardiometabol_med = gq_med_statin_DER
replace cardiometabol_med=0 if cardiometabol_med==-1
replace cardiometabol_med=1 if cardiometabol_med==2
replace cardiometabol_med=1 if gq_med_cardio_DER==2 | gq_med_cardio_DER==-3 | gq_med_cardio_DER==-5
replace cardiometabol_med=1 if gq_med_acei_DER==2
replace cardiometabol_med=1 if gq_med_arb_DER==2
replace cardiometabol_med=1 if gq_med_dm_DER==2
replace cardiometabol_med=1 if gq_med_lipid_DER==2
replace cardiometabol_med=1 if gq_med_betab_DER==2 | gq_med_betab_DER==-3 | gq_med_betab_DER==-5   
replace cardiometabol_med=1 if gq_med_diuret_DER==2
replace cardiometabol_med=1 if gq_med_nitr_DER==2
replace cardiometabol_med=1 if BNF_BetaBlocker==1
replace cardiometabol_med=1 if BNF_LipidLowering==1
replace cardiometabol_med=1 if BNF_OralAntidiabetic==1
mdesc cardiometabol_med //n=7 missing