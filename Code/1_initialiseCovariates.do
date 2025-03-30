version 17.0

***************************
** Initialize covariates **
***************************

*********
** Sex **
*********

label define sexlab 0 "Women" 1 "Men"
label values sex sexlab

*********
** Age **
*********

gen age = AgeAtTest_DM_Attended

rename AGE_GROUP agecat3
label define agelab3 1 "<45" 2 ">=45 & <55" 3 ">=55", replace
label values agecat3 agelab3

gen agecat2 = .
replace agecat2 = 1 if age <  50 & age != .
replace agecat2 = 2 if age >= 50 & age != .
label define agelab2 1 "<50" 2 ">=50", replace
label values agecat2 agelab2

*********
** BMI **
*********

rename BMI bmi
rename BMI_GROUP bmicat3
label define bmilab3 1 "Normal" 2 "Overweight" 3 "Obese", replace
label values bmicat3 bmilab3

gen bmicat2 = .
replace bmicat2 = 1 if bmicat3 == 1
replace bmicat2 = 2 if bmicat3 == 2 | bmicat3 == 3
label define bmilab2 1 "Normal" 2 "Overweight/Obese", replace
label values bmicat2 bmilab2

*********
** RHR **
*********

rename RestingHeartRate rhr
egen rhr2 = rowmedian(BPPR1 BPPR2 BPPR3)
replace rhr = rhr2 if rhr == . & rhr2 != .
drop rhr2



***************         
** Ethnicity **
***************

gen ethnic3 = gq_eth_DER 


replace ethnic3 = 0 if ethnic3==1 | ethnic3==2 | ethnic3==3
replace ethnic3 = 1 if ethnic3 !=0 & ethnic3 !=.
replace ethnic3 = 2 if ethnic3 == .

label define ethlab3 0 "White" 1 "Non-white" 2 "Missing/Unknown", replace
label values ethnic3 ethlab3

gen ethnic2 = .
replace ethnic2 = 0 if ethnic3 == 0
replace ethnic2 = 1 if ethnic3 == 1 | ethnic3 == 2
label define ethlab2 0 "White" 1 "Non-white/Missing/Unknown", replace
label values ethnic2 ethlab2

***************    
** Test Site **
*************** 

encode TestSite, gen(testsite)

**************************************************
** Season: centred around solstices and equinox **
**************************************************

gen date = date(P_startdate, "DMY", 2022)
format date %td
gen doy=doy(date)
gen season = 0 if doy>=80-45 & doy<=80+46       // approx autumn equinox
replace season = 1 if doy>=172-45 & doy<=172+47 // approx summer solstice
replace season = 2 if doy>=266-46 & doy<=266+45 // approx autumn equinox
replace season = 3 if doy>=355-45 | doy<=80-46  // approx winter solstice
drop doy
label define seasonlab 0 "Spring" 1 "Summer" 2 "Autumn" 3 "Winter"
label values season seasonlab


********************
** Marital Status **
********************

/*
gen marital_status = gq_marit_DER
replace marital_status=6 if marital_status==-1 | marital_status==-8 | marital_status==. 
label define marlab 1 "Single" 2 "Married/living as married" 3 "Widowed" 4 "Separated" 5 "Divorced" 6 "Missing/Unknown"
label values marital_status marlab
*/

gen     marital_status4 = .
replace marital_status4 = 1 if gq_marit_DER == 1
replace marital_status4 = 2 if gq_marit_DER == 2
replace marital_status4 = 3 if gq_marit_DER == 3 | gq_marit_DER == 4 | gq_marit_DER == 5
replace marital_status4 = 4 if gq_marit_DER ==-1 | gq_marit_DER ==-8 | gq_marit_DER == .
label define marlab4 1 "Single" 2 "Married/living as married" 3 "Widowed/Separated/Divorced" 4 "Missing/Unknown", replace
label values marital_status4 marlab4

gen marital_status3 = . 
replace marital_status3 = 1 if marital_status4 == 2
replace marital_status3 = 2 if marital_status4 == 1 | marital_status4 == 3
replace marital_status3 = 3 if marital_status4 == 4
label define marlab3 1 "Married/living as married" 2 "Single/Widowed/Separated/Divorce" 3 "Missing/Unknown", replace
label values marital_status3 marlab3

************
** Income **
************

gen income4 = gq_income_DER
replace income4=4 if income==. | income==-7 | income==-1 
label define incomelab4 1 "<20k" 2 "20-40k" 3 ">40k" 4 "Missing/Unknown", replace
label values income4 incomelab4

gen income3 = .
replace income3 = 1 if income4 == 3
replace income3 = 2 if income4 == 1 | income4 == 2
replace income3 = 3 if income4 == 4
label define incomelab3 1 ">40k" 2 "<=40k" 3 "Missing/Unknown", replace
label values income3 incomelab3

***************
** Education **
***************

rename EDUCATION education4
replace education4 = 4 if education4 == .
label define edulab4 0 "None" 1 "Compulsory" 2 "Further" 3 "Higher" 4 "Missing/Unknown", replace 
label values education4 edulab4

gen education3 = .
replace education3 = 0 if education4 == 3
replace education3 = 1 if education4 == 2
replace education3 = 2 if education4 == 1 | education4 == 0
replace education3 = 3 if education4 == 4
label define edulab3 0 "Higher" 1 "Further" 2 "None/Compulsory" 3 "Missing/Unknown", replace 
label values education3 edulab3


*****************
** Work Status **
*****************


gen work_status6 = 0
replace work_status6 = 1 if Worktype == "1" 
replace work_status6 = 2 if Worktype == "2" 
replace work_status6 = 3 if Worktype == "3" 
replace work_status6 = 3 if Worktype == "4"
replace work_status6 = 6 if Worktype == "-1" | Worktype == "(-10) N/A" | Worktype == "" | Worktype == "(-10) None" 
replace work_status6 = 4 if Worktype == "-1" & gq_job_retired_DER == 2
replace work_status6 = 5 if Worktype == "-1" & gq_job_unempl_DER == 2 

*All mixed recoded to most sendentary selected

#delimit ;

replace work_status6 = 1 if  Worktype == "(-5) 1 AND 2" 												| 
                            Worktype == "(-5) 1 and 2" 												| 
                            Worktype == "(-5) 1 and 2 and 3" 										| 
                            Worktype == "(-5) 1 AND 3" 												| 
                            Worktype == "(-5) 1 and 3" 												| 
                            Worktype == "(-5) 1 AND 2 AND 3" 										| 
                            Worktype == "(-5) 1 AND 3 AND 4" 										| 
                            Worktype == "(-5) 1 AND 2 AND 3 AND 4" 									| 
                            Worktype == "(-10) TICK IN BOX 1 '8 MONTHS'. TICK IN BOX 4 '4 MONTHS'" 	| 
                            Worktype == "(-5) 1 AND 2 (-10) I WORK FROM HOME" 						| 
                            Worktype == "(-10) 1 AND 2 HALFWAY" 									| 
                            Worktype == "(-10) 1 AND 2 50/50" 										| 
                            Worktype == "(-5) 1 AND 3 (-10) 50/50 SPLIT VOLUNTEER CAN'T DECIDE" 	| 
                            Worktype == "(-10) job 2 box 1 and job 1 box 3" 						| 
                            Worktype == "(-5) 1 AND 2 (-10) A BIT OF BOTH" 							| 
                            Worktype == "(-10) TICK 1 45 AND TICK 2 18" 							| 
                            Worktype == "(-5) 1 AND 4 (-10) 50%/50%"								;

                                
replace work_status6 = 2 if  Worktype == "(-5) 2 AND 3" 	| 
                            Worktype == "(-5) 2 and 3" 	| 
                            Worktype == "(-5) 2 and 4" 	| 
                            Worktype == "(-5) 2 AND 4"	;
                        
replace work_status6 = 3 if  Worktype == "(-5) 3 AND 4" 	| 
                            Worktype == "(-5) 3 and 4"	;

replace work_status6 = 6 if  work_status == 0 |
                            work_status == . ;

#delimit cr


label define work_label6 1 "Sedentary" 2 "Standing" 3 "Manual work" 4 "Retired" 5 "Unemployed" 6 "Missing/Unknown", replace
label values work_status6 work_label6 

gen work_status5 = .
replace work_status5 = 1 if work_status6 == 1
replace work_status5 = 2 if work_status6 == 2
replace work_status5 = 3 if work_status6 == 3
replace work_status5 = 4 if work_status6 == 4 | work_status6 == 5
replace work_status5 = 5 if work_status6 == 6

label define work_label5 1 "Sedentary" 2 "Standing" 3 "Manual work" 4 "Retired/Unemployed" 5 "Missing/Unknown", replace
label values work_status5 work_label5 

*************
** Smoking **
*************

rename SMOKE smoke4
replace smoke4 = 3 if smoke4 == 9 | smoke4 == .
label define smoke_labels4 0 "Never smoked" 1 "Ex smoker" 2 "Current smoker" 3 "Missing/Unknown", replace
label values smoke4 smoke_labels4

gen     smoke3 = .
replace smoke3 = 0 if smoke4 == 0
replace smoke3 = 1 if smoke4 == 1 | smoke4 == 2
replace smoke3 = 2 if smoke4 == 3
label define smoke_labels3 0 "Never smoker" 1 "Ever smoker" 2 "Missing/Unknown", replace
label values smoke3 smoke_labels3

**********
** Diet **
**********

xtile diet = DASH_AccordanceScore, nquantiles(2)
replace diet = 1 if diet == .
//label define diet_labels 1 "1st quintile of DASH" 2 "2nd quintile of DASH" 3 "3rd quintile of DASH" 4 "4th quintile of DASH" 5 "5th quintile of DASH" 6 "Missing/Unknown"
label define diet_labels 1 "Below median of DASH/Missing/Unknown" 2 "Above median of DASH"
label values diet diet_labels

*************
** Alcohol **
*************

recode gq_alc_freqc_DER (1/2=1 "<1/wk") (3/4=2 "1-4/wk") (5=3 "Almost-daily") (-7/.=4 "Unknown") , gen(alcohol4)
label define alcohol_labels4 1 "<1 per week" 2 "1-4 per week" 3 "Almost daily" 4 "Missing/Unknown", replace 
label values alcohol4 alcohol_labels4

gen alcohol3 = .
replace alcohol3 = 1 if alcohol4 == 1
replace alcohol3 = 2 if alcohol4 == 2 | alcohol4 == 3
replace alcohol3 = 3 if alcohol4 == 4
label define alcohol_labels3 1 "<1 per week" 2 ">=1 per week" 3 "Missing/Unknown", replace 
label values alcohol3 alcohol_labels3

*********************************
** Cardiometobolic medications **
*********************************

gen     cardiometabol_med = gq_med_statin_DER
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
replace cardiometabol_med=2 if cardiometabol_med == .
label define cardiometabol_labels 0 "No cardiometabol med" 1 "Taking cardiometabol med" 2 "Missing/Unknown" 
label values cardiometabol_med cardiometabol_labels

**********************************************************************************************
** Time of testing clinical assessments started, approximately (using rest test start time) **
**********************************************************************************************

gen testTime = mm(clock(RestStart, "hm")) + 60*hh(clock(RestStart, "hm"))
replace testTime = P_startminofday if testTime == .
replace testTime = mm(clock(TimeOf120BloodSample, "hm")) + 60*hh(clock(TimeOf120BloodSample, "hm")) if testTime < 360

gen testTime_sin = sin(testTime*2*_pi/1440)
gen testTime_cos = cos(testTime*2*_pi/1440)

drop testTime

// also adding testDay

gen testDay     = doy(date(P_startdate, "DMY", 2022))
gen testDay_sin = sin(testDay*2*_pi/365.25)
gen testDay_cos = cos(testDay*2*_pi/365.25)

drop testDay





