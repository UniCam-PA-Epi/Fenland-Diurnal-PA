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

rename AGE_GROUP agecat
label define agelab 1 "<45" 2 ">=45 & <55" 3 ">=55"
label values agecat agelab

*********
** BMI **
*********

rename BMI bmi
rename BMI_GROUP bmicat
label define bmilab 1 "Normal" 2 "Overweight" 3 "Obese"
label values bmicat bmilab

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

gen ethnic = gq_eth_DER 


replace ethnic = 0 if ethnic==1 | ethnic==2 | ethnic==3
replace ethnic = 1 if ethnic !=0 & ethnic !=.
replace ethnic = 2 if ethnic == .

label define ethlab 0 "White" 1 "Non-white" 2 "Missing/Unknown"
label values ethnic ethlab

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

gen marital_status = gq_marit_DER
replace marital_status=6 if marital_status==-1 | marital_status==-8 | marital_status==. 
label define marlab 1 "Single" 2 "Married/living as married" 3 "Widowed" 4 "Separated" 5 "Divorced" 6 "Missing/Unknown"
label values marital_status marlab


************
** Income **
************

gen income = gq_income_DER
replace income=4 if income==. | income==-7 | income==-1 
label define incomelab 1 "<20k" 2 "20-40k" 3 ">40k" 4 "Missing/Unknown"
label values income incomelab


***************
** Education **
***************

rename EDUCATION education
replace education = 4 if education == .
label define edulab 0 "None" 1 "Compulsory" 2 "Further" 3 "Higher" 4 "Missing/Unknown"
label values education edulab

*****************
** Work Status **
*****************

gen work_status = 0
replace work_status = 1 if Worktype == "1" 
replace work_status = 2 if Worktype == "2" 
replace work_status = 3 if Worktype == "3" 
replace work_status = 3 if Worktype == "4"
replace work_status = 6 if Worktype == "-1" | Worktype == "(-10) N/A" | Worktype == "" | Worktype == "(-10) None" 
replace work_status = 4 if Worktype == "-1" & gq_job_retired_DER == 2
replace work_status = 5 if Worktype == "-1" & gq_job_unempl_DER == 2 

*All mixed recoded to most sendentary selected

#delimit ;

replace work_status = 1 if  Worktype == "(-5) 1 AND 2" 												| 
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

                                
replace work_status = 2 if  Worktype == "(-5) 2 AND 3" 	| 
                            Worktype == "(-5) 2 and 3" 	| 
                            Worktype == "(-5) 2 and 4" 	| 
                            Worktype == "(-5) 2 AND 4"	;
                        
replace work_status = 3 if  Worktype == "(-5) 3 AND 4" 	| 
                            Worktype == "(-5) 3 and 4"	;

replace work_status = 6 if  work_status == 0 |
                            work_status == . ;

#delimit cr

label define work_label 1 "Sedentary" 2 "Standing" 3 "Manual work" 4 "Retired" 5 "Unemployed" 6 "Missing/Unknown"
label values work_status work_label 
tab work_status

*************
** Smoking **
*************

rename SMOKE smoke
replace smoke = 3 if smoke == 9 | smoke == .
label define smoke_labels 0 "Never smoked" 1 "Ex smoker" 2 "Current smoker" 3 "Missing/Unknown"
label values smoke smoke_labels

**********
** Diet **
**********

xtile diet = DASH_AccordanceScore, nquantiles(5)
replace diet = 6 if diet == .
label define diet_labels 1 "1st quintile of DASH" 2 "2nd quintile of DASH" 3 "3rd quintile of DASH" 4 "4th quintile of DASH" 5 "5th quintile of DASH" 6 "Missing/Unknown"
label values diet diet_labels

*************
** Alcohol **
*************

recode gq_alc_freqc_DER (1/2=1 "<1/wk") (3/4=2 "1-4/wk") (5=3 "Almost-daily") (-7/.=4 "Unknown") , gen(alcohol)
label define alcohol_labels 1 "<1 per week" 2 "1-4 per week" 3 "Almost daily" 4 "Missing/Unknown" 
label values alcohol alcohol_labels

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





