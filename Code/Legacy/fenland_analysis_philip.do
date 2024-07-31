
cd "C:\Users\Mona\Desktop\Cambridge"
use "fenland_raw"
save "fenland_working", replace
/*
CONTENTS (CTRL+F & one of the following words to find sections):
STEP #1: Organise PAEE
STEP #2: Organise Outomes
STEP #3: Organise Covariates
STEP #4: Descriprive (Table 1)
STEP #5: Joint Associations (Table 2)
STEP #6: Joint Associations - Adiposity Adjusted (Table 3)
STEP #7: Volume Associations & Quintile Graphs
STEP #8: CoDA & Graphs (Figure 3)
STEP #9: CoDA & Graphs - Adiposity Adjusted (Suppl)
*/

********************************************************************************
* STEP #1: Organise PAEE
********************************************************************************
clear
use "fenland_working"

desc
// Obs: n=12,435 

//Remove individiuals with missing hours data
egen nmis=rmiss(paee_hour*)
drop if nmis!=0 // n=12,278 (n=157 deleted)
drop nmis 

//Remove individuals with less than 4 days wear time in any time window 
foreach var of varlist pwear_hour1-pwear_hour24 {
gen low`var'=1 if `var'<=3.999999
replace low`var'=0 if low`var'!=1
}
egen lowpwearT = rsum(lowpwear*)
drop if lowpwearT!=0 // n=11,766 (n=512 deleted)	  
drop lowpwear*

//Generate total PAEE, time-block PAEE, and percentage PAEE
egen paeeT = rsum(paee_hour*) // total
gen paeeT100 = paeeT/100
gen paeeTt = (paeeT*60)/1000

egen early=rsum(paee_hour5-paee_hour8) // time-block 
egen earlymid=rsum(paee_hour9-paee_hour12)
egen mid=rsum(paee_hour13-paee_hour16)
egen latemid=rsum(paee_hour17-paee_hour20)
egen late=rsum(paee_hour21-paee_hour24)
egen night=rsum(paee_hour1-paee_hour4)
gen earlyp=(early/paeeT)*100 // time-block percentage of total
gen earlymidp=(earlymid/paeeT)*100
gen midp=(mid/paeeT)*100
gen latemidp=(latemid/paeeT)*100
gen latep=(late/paeeT)*100
gen nightp=(night/paeeT)*100
gen earlyp5=earlyp/5 //time-block percentage of total divided by 5
gen earlymidp5=earlymidp/5
gen midp5=midp/5
gen latemidp5=latemidp/5
gen latep5=latep/5
gen nightp5=nightp/5

//Amplitude
gen amplitude = ((P_MORNING_adj_coefficient^2) + (P_MIDNIGHT_adj_coefficient^2))^0.5 
gen log_rel_amplitude = ln(amplitude/paeeTt)

//Acrophase
gen acrophase_rad = atan2(P_MORNING_adj_coefficient, P_MIDNIGHT_adj_coefficient)
drop if acrophase_rad==.
replace acrophase_rad = acrophase_rad + 2*_pi if acrophase_rad<0
gen sin_acro = sin(acrophase_rad)
gen cos_acro = cos(acrophase_rad)
gen acrophase_deg = (acrophase_rad)*57.29
gen acrophase_hours = (acrophase_deg)/15
sum acrophase_hours, detail
gen acrophase = acrophase_hours
gen acrophase_hour = (24*acrophase_rad)/(2*_pi)
********************************************************************************
* STEP #2: Organise Outomes
********************************************************************************
//Bodyfat
rename BodyFatPercent_Cons bodyfat
replace bodyfat=. if bodyfat==-1 //n=9 missing

//Insulin & log_Insulin
rename G_Insulin insulin
tab insulin, missing
//one extremely high value (1370, next is 747). n=2,323 (19.74%) missing
gen log_insulin=log(insulin)

//CRP & log_CRP
rename G_hs_CRP crp
replace crp=. if G_hs_CRP_Threshold==2
mdesc crp // n=2,546 missing (21.64%)
gen log_crp=log(crp)

//Leptin & log_Leptin
rename G_Leptin leptin
replace leptin=. if G_Leptin_Threshold!=1
mdesc leptin // n=2,348 missing (19.96%)
gen log_leptin=log(leptin)

//NEFA & log_NEFA
rename G_NEFA nefa
replace nefa=. if G_NEFA_Threshold!=1
mdesc nefa // n=2,345 missing (19.93%)
gen log_nefa=log(nefa)

//Adiponectin & log_Adiponectin
rename G_Adiponectin adiponectin
replace adiponectin=. if G_Adiponectin_Threshold!=1
mdesc adiponectin // n=2,340 missing (19.89%)
gen log_adiponectin = log(adiponectin)

//LDL
rename LDL0 ldl
mdesc ldl // n=158 missing (1,34%)

//Blood Pressures: generate averages for persons with 3 measures
mdesc BPDia1 
sum BPDia1 if BPDia2==. 
sum BPDia1 if BPDia3==. 
// n=5 with 1, n=70 with 2, n=3 with 3 missing measures
mdesc BPSys1 
sum BPSys1 if BPSys2==. 
sum BPSys1 if BPSys3==. 
// n=5 with 1, n=99 with 2, n=2 with 3 missing measures
gen mbpsys = (BPSys1+BPSys2+BPSys3)/3
gen mbpdia = (BPDia1+BPDia2+BPDia3)/3 
mdesc mbpdia // n=78 missing (0.66%)
mdesc mbpsys //  n=106 missing (0.9%)
sum mbpsys if mbpsys >=140 // ~1500 could be classified as having hypertension

//Glucose: convert to binary variable with levels ≥ 5.6 mmol/L defined 
//as “high glucose level” and < 5.6 mmol/L defined as “low glucose level” 
gen glucose0 = Glucose0
replace glucose0 = RepeatGlucose0 if RepeatGlucose0<=Glucose0
sum glucose0 if glucose0>=5.6 // i.e., high glucose (n=705)
gen glucose_cat = 1 if glucose0>=5.6
replace glucose_cat = 0 if glucose_cat!=1 & glucose0!=.
replace glucose_cat = . if glucose0==.
mdesc glucose_cat // 92 missing
mdesc glucose0 // 92 missing

********************************************************************************
* STEP #3: Organise Covariates
*sex, age, ethnic, testsite, season, 
*marital_status, income, edeucation, work_status,
*smoking, diet, alcohol,cardiometabol_med
********************************************************************************
//Sex
label define sexlab 0 "female" 1 "male"
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
gen smoking = gq_sm_pckyrs_DER
replace smoking=. if smoking==-8 | smoking==-7 | smoking==. //n=877 missing

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

save "working1.dta", replace
********************************************************************************
* STEP #4: Descriprive (Table 1)
********************************************************************************
clear
use "working1.dta"

tab sex
bysort sex: sum age
bysort sex: tab sex ethnic, row
bysort sex: tab sex education, row
bysort sex: tab sex income, row
bysort sex: tab sex work_s, row
bysort sex: tab sex marital_s, row
bysort sex: tab sex season, row
bysort sex: sum smoking
bysort sex: sum diet
bysort sex: sum alcohol
bysort sex: tab sex cardiometabol_med, row
bysort sex: tab sex testsite, row
bysort sex: sum paeeTt
bysort sex: sum bodyfat
bysort sex: sum insulin
bysort sex: sum leptin
bysort sex: sum nefa
bysort sex: sum ldl
bysort sex: sum adiponectin
bysort sex: sum crp
bysort sex: sum mbpdia
bysort sex: sum mbpsys
bysort sex: tab sex glucose_cat,row

********************************************************************************
* STEP #5: Joint Associations (Table 2)
*Ignore Acrophase Estimates & Acrophase Graphs               						
********************************************************************************
clear
use "working1.dta"
*************************************************
*ADIPOSITY
*************************************************
*Adiposty Female
regress bodyfat paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fadipo_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fadipo_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fadipo_acro_ul = fadipo_acro_risk + (1.96*fadipo_acro_se)
gen fadipo_acro_ll = fadipo_acro_risk - (1.96*fadipo_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

*Adiposty Male
regress bodyfat paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen madipo_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen madipo_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen madipo_acro_ul = madipo_acro_risk + (1.96*madipo_acro_se)
gen madipo_acro_ll = madipo_acro_risk - (1.96*madipo_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

quietly: twoway (connect fadipo_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fadipo_acro_ul fadipo_acro_ll acrophase, color(red%10)) ///
		(connect madipo_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap madipo_acro_ul madipo_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (%)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) ylabel(-3(1)3) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:A.     Adiposity}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_adipo", replace

		
*************************************************
*Insulin
*************************************************
*Insulin Female
regress log_insulin paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen finsulin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen finsulin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen finsulin_acro_risk = (exp(finsulin_acro_risk_log)-1)*100
gen finsulin_acro_ul = (exp(finsulin_acro_risk_log)+(1.96*finsulin_acro_se)-1)*100
gen finsulin_acro_ll = (exp(finsulin_acro_risk_log) - (1.96*finsulin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100


lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*Insulin Male
regress log_insulin paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen minsulin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen minsulin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen minsulin_acro_risk = (exp(minsulin_acro_risk_log)-1)*100
gen minsulin_acro_ul = (exp(minsulin_acro_risk_log)+(1.96*minsulin_acro_se)-1)*100
gen minsulin_acro_ll = (exp(minsulin_acro_risk_log) - (1.96*minsulin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


quietly: twoway (connect finsulin_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap finsulin_acro_ul finsulin_acro_ll acrophase, color(red%10)) ///
		(connect minsulin_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap minsulin_acro_ul minsulin_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (%)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:B.     Insulin}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_insulin", replace

*************************************************
*Diastolic BP
*************************************************
*Diastolic BP Female
regress mbpdia paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fdia_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fdia_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fdia_acro_ul = fdia_acro_risk + (1.96*fdia_acro_se)
gen fdia_acro_ll = fdia_acro_risk - (1.96*fdia_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase


*Diastolic BP Male
regress mbpdia paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen mdia_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mdia_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mdia_acro_ul = mdia_acro_risk + (1.96*mdia_acro_se)
gen mdia_acro_ll = mdia_acro_risk - (1.96*mdia_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

quietly: twoway (connect fdia_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fdia_acro_ul fdia_acro_ll acrophase, color(red%10)) ///
		(connect mdia_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap mdia_acro_ul mdia_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (mmHg)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:C.     Diastolic BP}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_mbpdia", replace

*************************************************
*Glucose
*************************************************
*Glucose Female
logistic glucose_cat paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fgluc_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fgluc_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fgluc_acro_ul = fgluc_acro_risk + (1.96*fgluc_acro_se)
gen fgluc_acro_ll = fgluc_acro_risk - (1.96*fgluc_acro_se)

di exp(_b[log_rel_amplitude])*0.01 // amplitude
di exp((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))*0.01 // amplitude
di exp((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))*0.01 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase


*Glucose Male
logistic glucose_cat paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen mgluc_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mgluc_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mgluc_acro_ul = mgluc_acro_risk + (1.96*mgluc_acro_se)
gen mgluc_acro_ll = mgluc_acro_risk - (1.96*mgluc_acro_se)

di exp(_b[log_rel_amplitude])*0.01 // amplitude
di exp((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))*0.01 // amplitude
di exp((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))*0.01 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

quietly: twoway (connect fgluc_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fgluc_acro_ul fgluc_acro_ll acrophase, color(red%10)) ///
		(connect mgluc_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap mgluc_acro_ul mgluc_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (OR)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:D.     High Glucose}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_gluc", replace


grc1leg2 "C_acro_adipo" "C_acro_insulin" "C_acro_mbpdia" "C_acro_gluc", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(2) legendfrom("C_acro_adipo") position(3) lcols(1)
graph save "C_acro_combined_paper1", replace
graph export C_acro_combined_paper1.tif, replace


*************************************************
*LDL
*************************************************
*LDL Female
regress ldl paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fldl_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fldl_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fldl_acro_ul = fldl_acro_risk + (1.96*fldl_acro_se)
gen fldl_acro_ll = fldl_acro_risk - (1.96*fldl_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

*LDL Male
regress ldl paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen mldl_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mldl_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mldl_acro_ul = mldl_acro_risk + (1.96*mldl_acro_se)
gen mldl_acro_ll = mldl_acro_risk - (1.96*mldl_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

quietly: twoway (connect fldl_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fldl_acro_ul fldl_acro_ll acrophase, color(red%10)) ///
		(connect mldl_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap mldl_acro_ul mldl_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (%)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:A.     LDL}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_ldl", replace


*************************************************
*NEFA
*************************************************
*NEFA Female
regress log_nefa paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fnefa_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fnefa_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fnefa_acro_risk = (exp(fnefa_acro_risk_log)-1)*100
gen fnefa_acro_ul = (exp(fnefa_acro_risk_log)+(1.96*fnefa_acro_se)-1)*100
gen fnefa_acro_ll = (exp(fnefa_acro_risk_log) - (1.96*fnefa_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


*NEFA Male
regress log_nefa paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen mnefa_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mnefa_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mnefa_acro_risk = (exp(mnefa_acro_risk_log)-1)*100
gen mnefa_acro_ul = (exp(mnefa_acro_risk_log)+(1.96*mnefa_acro_se)-1)*100
gen mnefa_acro_ll = (exp(mnefa_acro_risk_log) - (1.96*mnefa_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


quietly: twoway (connect fnefa_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fnefa_acro_ul fnefa_acro_ll acrophase, color(red%10)) ///
		(connect mnefa_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap mnefa_acro_ul mnefa_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (%)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:B.     NEFAs}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_nefa", replace

*************************************************
*Leptin
*************************************************
*Leptin Female
regress log_leptin paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fleptin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fleptin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fleptin_acro_risk = (exp(fleptin_acro_risk_log)-1)*100
gen fleptin_acro_ul = (exp(fleptin_acro_risk_log)+(1.96*fleptin_acro_se)-1)*100
gen fleptin_acro_ll = (exp(fleptin_acro_risk_log) - (1.96*fleptin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


*Leptin Male
regress log_leptin paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen mleptin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mleptin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mleptin_acro_risk = (exp(mleptin_acro_risk_log)-1)*100
gen mleptin_acro_ul = (exp(mleptin_acro_risk_log)+(1.96*mleptin_acro_se)-1)*100
gen mleptin_acro_ll = (exp(mleptin_acro_risk_log) - (1.96*mleptin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


quietly: twoway (connect fleptin_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fleptin_acro_ul fleptin_acro_ll acrophase, color(red%10)) ///
		(connect mleptin_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap mleptin_acro_ul mleptin_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (%)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:C.     Leptin}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_leptin", replace

*************************************************
*Systolic BP
*************************************************
*Systolic BP Female
regress mbpsys paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fsys_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fsys_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fsys_acro_ul = fsys_acro_risk + (1.96*fsys_acro_se)
gen fsys_acro_ll = fsys_acro_risk - (1.96*fsys_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

*Systolic BP Male
regress mbpsys paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen msys_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen msys_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen msys_acro_ul = msys_acro_risk + (1.96*msys_acro_se)
gen msys_acro_ll = msys_acro_risk - (1.96*msys_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase
quietly: twoway (connect fleptin_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fleptin_acro_ul fleptin_acro_ll acrophase, color(red%10)) ///
		(connect mleptin_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap mleptin_acro_ul mleptin_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (mmHg)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:D.     Systolic BP}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_mbpsys", replace


*************************************************
*Adiponectin
*************************************************
*adiponectin Female
regress log_adiponectin paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fadiponectin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fadiponectin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fadiponectin_acro_risk = (exp(fadiponectin_acro_risk_log)-1)*100
gen fadiponectin_acro_ul = (exp(fadiponectin_acro_risk_log)+(1.96*fadiponectin_acro_se)-1)*100
gen fadiponectin_acro_ll = (exp(fadiponectin_acro_risk_log) - (1.96*fadiponectin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


*adiponectin Male
regress log_adiponectin paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen madiponectin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen madiponectin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen madiponectin_acro_risk = (exp(madiponectin_acro_risk_log)-1)*100
gen madiponectin_acro_ul = (exp(madiponectin_acro_risk_log)+(1.96*madiponectin_acro_se)-1)*100
gen madiponectin_acro_ll = (exp(madiponectin_acro_risk_log) - (1.96*madiponectin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


quietly: twoway (connect fadiponectin_acro_risk acrophase, sort lcolor(red) msymbol(Oh) mcolor(red%30)) ///
		(rcap fadiponectin_acro_ul fadiponectin_acro_ll acrophase, color(red%10)) ///
		(connect madiponectin_acro_risk acrophase, sort lcolor(blue) msymbol(Dh) mcolor(blue%30)) ///
		(rcap madiponectin_acro_ul madiponectin_acro_ll acrophase, color(blue%10)), ///
		ytitle("Difference (%)", color(black) margin(medium)) ///
		xtitle("acrophase (clock time)", color(black) margin(medium)) ///
		xlabel(0(4)24) yline(0, lcolor(black)) ///
		legend(order(1 "female" 3 "male") cols(1) pos(3)) /// 
		note("{bf:E.     Adiponectin}",  position(11) ring(50) size(large)) ///
		graphregion(fcolor(white) lcolor(white)) bgcolor(white) ///
		plotregion(fcolor(white))
graph save "C_acro_adiponectin", replace


*************************************************
*Crp
*************************************************
*crp Female
regress log_crp paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fcrp_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fcrp_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fcrp_acro_risk = (exp(fcrp_acro_risk_log)-1)*100
gen fcrp_acro_ul = (exp(fcrp_acro_risk_log)+(1.96*fcrp_acro_se)-1)*100
gen fcrp_acro_ll = (exp(fcrp_acro_risk_log) - (1.96*fcrp_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*crp Male
regress log_crp paeeTt log_rel_amplitude sin_acro cos_acro age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen mcrp_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mcrp_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mcrp_acro_risk = (exp(mcrp_acro_risk_log)-1)*100
gen mcrp_acro_ul = (exp(mcrp_acro_risk_log)+(1.96*mcrp_acro_se)-1)*100
gen mcrp_acro_ll = (exp(mcrp_acro_risk_log) - (1.96*mcrp_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


********************************************************************************
* STEP #6: Joint Associations - Adiposity Adjusted (Table 3)
********************************************************************************
clear
use "working1.dta"
*************************************************
*Insulin
*************************************************
*Insulin Female
regress log_insulin paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen finsulin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen finsulin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen finsulin_acro_risk = (exp(finsulin_acro_risk_log)-1)*100
gen finsulin_acro_ul = (exp(finsulin_acro_risk_log)+(1.96*finsulin_acro_se)-1)*100
gen finsulin_acro_ll = (exp(finsulin_acro_risk_log) - (1.96*finsulin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100


lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*Insulin Male
regress log_insulin paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen minsulin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen minsulin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen minsulin_acro_risk = (exp(minsulin_acro_risk_log)-1)*100
gen minsulin_acro_ul = (exp(minsulin_acro_risk_log)+(1.96*minsulin_acro_se)-1)*100
gen minsulin_acro_ll = (exp(minsulin_acro_risk_log) - (1.96*minsulin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*************************************************
*Diastolic BP
*************************************************
*Diastolic BP Female
regress mbpdia paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fdia_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fdia_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fdia_acro_ul = fdia_acro_risk + (1.96*fdia_acro_se)
gen fdia_acro_ll = fdia_acro_risk - (1.96*fdia_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase


*Diastolic BP Male
regress mbpdia paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen mdia_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mdia_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mdia_acro_ul = mdia_acro_risk + (1.96*mdia_acro_se)
gen mdia_acro_ll = mdia_acro_risk - (1.96*mdia_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

*************************************************
*Glucose
*************************************************
*Glucose Female
logistic glucose_cat paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fgluc_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fgluc_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fgluc_acro_ul = fgluc_acro_risk + (1.96*fgluc_acro_se)
gen fgluc_acro_ll = fgluc_acro_risk - (1.96*fgluc_acro_se)

di exp(_b[log_rel_amplitude])*0.01 // amplitude
di exp((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))*0.01 // amplitude
di exp((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))*0.01 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase


*Glucose Male
logistic glucose_cat paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen mgluc_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mgluc_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mgluc_acro_ul = mgluc_acro_risk + (1.96*mgluc_acro_se)
gen mgluc_acro_ll = mgluc_acro_risk - (1.96*mgluc_acro_se)

di exp(_b[log_rel_amplitude])*0.01 // amplitude
di exp((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))*0.01 // amplitude
di exp((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))*0.01 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase


*************************************************
*LDL
*************************************************
*LDL Female
regress ldl paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fldl_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fldl_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fldl_acro_ul = fldl_acro_risk + (1.96*fldl_acro_se)
gen fldl_acro_ll = fldl_acro_risk - (1.96*fldl_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

*LDL Male
regress ldl paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen mldl_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mldl_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mldl_acro_ul = mldl_acro_risk + (1.96*mldl_acro_se)
gen mldl_acro_ll = mldl_acro_risk - (1.96*mldl_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase


*************************************************
*NEFA
*************************************************
*NEFA Female
regress log_nefa paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fnefa_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fnefa_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fnefa_acro_risk = (exp(fnefa_acro_risk_log)-1)*100
gen fnefa_acro_ul = (exp(fnefa_acro_risk_log)+(1.96*fnefa_acro_se)-1)*100
gen fnefa_acro_ll = (exp(fnefa_acro_risk_log) - (1.96*fnefa_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


*NEFA Male
regress log_nefa paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen mnefa_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mnefa_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mnefa_acro_risk = (exp(mnefa_acro_risk_log)-1)*100
gen mnefa_acro_ul = (exp(mnefa_acro_risk_log)+(1.96*mnefa_acro_se)-1)*100
gen mnefa_acro_ll = (exp(mnefa_acro_risk_log) - (1.96*mnefa_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*************************************************
*Leptin
*************************************************
*Leptin Female
regress log_leptin paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fleptin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fleptin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fleptin_acro_risk = (exp(fleptin_acro_risk_log)-1)*100
gen fleptin_acro_ul = (exp(fleptin_acro_risk_log)+(1.96*fleptin_acro_se)-1)*100
gen fleptin_acro_ll = (exp(fleptin_acro_risk_log) - (1.96*fleptin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


*Leptin Male
regress log_leptin paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen mleptin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mleptin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mleptin_acro_risk = (exp(mleptin_acro_risk_log)-1)*100
gen mleptin_acro_ul = (exp(mleptin_acro_risk_log)+(1.96*mleptin_acro_se)-1)*100
gen mleptin_acro_ll = (exp(mleptin_acro_risk_log) - (1.96*mleptin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*************************************************
*Systolic BP
*************************************************
*Systolic BP Female
regress mbpsys paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen fsys_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fsys_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fsys_acro_ul = fsys_acro_risk + (1.96*fsys_acro_se)
gen fsys_acro_ll = fsys_acro_risk - (1.96*fsys_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

*Systolic BP Male
regress mbpsys paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di _b[paeeTt] // volume
di _b[paeeTt]+(1.96*_se[paeeTt]) // volume
di _b[paeeTt]-(1.96*_se[paeeTt]) // volume

gen msys_acro_risk = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen msys_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen msys_acro_ul = msys_acro_risk + (1.96*msys_acro_se)
gen msys_acro_ll = msys_acro_risk - (1.96*msys_acro_se)

di _b[log_rel_amplitude]/100
di ((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))/100 // amplitude
di ((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))/100 // amplitude

lincom _b[sin_acro] + _b[cos_acro]
di r(estimate) // acrophase
di r(estimate) + (1.96*r(se)) // acrophase
di r(estimate) - (1.96*r(se)) // acrophase

*************************************************
*Adiponectin
*************************************************
*adiponectin Female
regress log_adiponectin paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fadiponectin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fadiponectin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fadiponectin_acro_risk = (exp(fadiponectin_acro_risk_log)-1)*100
gen fadiponectin_acro_ul = (exp(fadiponectin_acro_risk_log)+(1.96*fadiponectin_acro_se)-1)*100
gen fadiponectin_acro_ll = (exp(fadiponectin_acro_risk_log) - (1.96*fadiponectin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


*adiponectin Male
regress log_adiponectin paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen madiponectin_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen madiponectin_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen madiponectin_acro_risk = (exp(madiponectin_acro_risk_log)-1)*100
gen madiponectin_acro_ul = (exp(madiponectin_acro_risk_log)+(1.96*madiponectin_acro_se)-1)*100
gen madiponectin_acro_ll = (exp(madiponectin_acro_risk_log) - (1.96*madiponectin_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*************************************************
*Crp
*************************************************
*crp Female
regress log_crp paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 0

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen fcrp_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen fcrp_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen fcrp_acro_risk = (exp(fcrp_acro_risk_log)-1)*100
gen fcrp_acro_ul = (exp(fcrp_acro_risk_log)+(1.96*fcrp_acro_se)-1)*100
gen fcrp_acro_ll = (exp(fcrp_acro_risk_log) - (1.96*fcrp_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase

*crp Male
regress log_crp paeeTt log_rel_amplitude sin_acro cos_acro bodyfat age i.ethnic ///
i.testsite i.season i.education i.income i.work_s smoking diet alcohol ///
i.cardiometabol_med if sex == 1

di (exp(_b[paeeTt])-1)*100 // volume
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100 // volume
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 // volume

gen mcrp_acro_risk_log = (_b[sin_acro]*sin_acro) + (_b[cos_acro]*cos_acro)
gen mcrp_acro_se = sqrt((_se[sin_acro]*sin_acro)^2 + (_se[cos_acro]*cos_acro)^2)
gen mcrp_acro_risk = (exp(mcrp_acro_risk_log)-1)*100
gen mcrp_acro_ul = (exp(mcrp_acro_risk_log)+(1.96*mcrp_acro_se)-1)*100
gen mcrp_acro_ll = (exp(mcrp_acro_risk_log) - (1.96*mcrp_acro_se)-1)*100

di ((1.01^(_b[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])+(1.96*_se[log_rel_amplitude]))-1)*100
di (1.01^((_b[log_rel_amplitude])-(1.96*_se[log_rel_amplitude]))-1)*100

lincom _b[sin_acro] + _b[cos_acro]
di (exp(r(estimate))-1)*100 // acrophase
di (exp(r(estimate) + (1.96*r(se)))-1)*100 // acrophase
di (exp(r(estimate) - (1.96*r(se)))-1)*100 // acrophase


********************************************************************************
* STEP #7: Volume Associations & Quintile Graphs
********************************************************************************
************************
*Total PAEE & Bodyfat  *
************************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress bodyfat paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di _b[paeeTt]
di _b[paeeTt]+(1.96*_se[paeeTt])
di _b[paeeTt]-(1.96*_se[paeeTt])
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=bodyfat (semean) se=bodyfat, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (%)", color(black) margin(medium)) ///
ylabel(25(5)40, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Adiposity", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:A.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_bodyfat", replace
*Male Graph
clear

use "working1.dta"
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=bodyfat (semean) se=bodyfat, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (%)", color(black) margin(medium)) ///
ylabel(25(5)40, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Adiposity", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:A.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_bodyfat", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=bodyfat (semean) se=bodyfat, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (%)", color(black) margin(medium)) ///
ylabel(20(5)40, nogrid) xlabel(0(20)100) ///
legend(order(4 "female" 6 "male")) xtitle("") title("Adiposity", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:A.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_bodyfat", replace
*******************
*Total PAEE & LDL *
*******************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress ldl paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di _b[paeeTt]
di _b[paeeTt]+(1.96*_se[paeeTt])
di _b[paeeTt]-(1.96*_se[paeeTt])
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=ldl (semean) se=ldl, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mmol/L)", color(black) margin(medium)) ///
ylabel(3.0(0.1)3.6, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("LDL", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:A.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_ldl", replace
*Male Graph
clear

use "working1.dta"
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=ldl (semean) se=ldl, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mmol/L)", color(black) margin(medium)) ///
ylabel(3(0.1)3.6, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("LDL", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:B.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_ldl", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=ldl (semean) se=ldl, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (mmol/L)", color(black) margin(medium)) ///
ylabel(3(0.1)3.6, nogrid) xlabel(0(20)100) ///
legend(order(4 "female" 6 "male")) xtitle("") title("LDL", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:A.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_ldl", replace
*********************
*Total PAEE & NEFAs *
*********************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress log_nefa paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di (exp(_b[paeeTt])-1)*100
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=nefa (semean) se=nefa, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (μmol/L)", color(black) margin(medium)) ///
ylabel(275(50)425, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("NEFAs", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:C.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_nefa", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=nefa (semean) se=nefa, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (μmol/L)", color(black) margin(medium)) ///
ylabel(275(50)425, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("NEFAs", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:C.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_nefa", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=nefa (semean) se=nefa, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (μmol/L)", color(black) margin(medium)) ///
ylabel(275(50)425, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("NEFAs", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:B.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_nefa", replace
***********************
*Total PAEE & Insulin *
***********************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress log_insulin paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di (exp(_b[paeeTt])-1)*100 
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=insulin (semean) se=insulin, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (pmol/L)", color(black) margin(medium)) ///
ylabel(25(10)75, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Insulin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:D.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_insulin", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=insulin (semean) se=insulin, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (pmol/L)", color(black) margin(medium)) ///
ylabel(25(10)75, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Insulin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:D.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_insulin", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=insulin (semean) se=insulin, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (pmol/L)", color(black) margin(medium)) ///
ylabel(25(10)75, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Insulin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:B.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_insulin", replace
***************************
*Total PAEE & Adiponectin *
***************************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress log_adiponectin paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di (exp(_b[paeeTt])-1)*100
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=adiponectin (semean) se=adiponectin, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (μg/ml)", color(black) margin(medium)) ///
ylabel(4.5(1)9.5, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Adiponectin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:E.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_adiponectin", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=adiponectin (semean) se=adiponectin, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (μg/ml)", color(black) margin(medium)) ///
ylabel(4.5(1)9.5, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Adiponectin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:E.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_adiponectin", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=adiponectin (semean) se=adiponectin, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (μg/ml)", color(black) margin(medium)) ///
ylabel(4.5(1)9.5, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Adiponectin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:E.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_adiponectin", replace
**********************
*Total PAEE & Leptin *
**********************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress log_leptin paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di (exp(_b[paeeTt])-1)*100
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100 
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=leptin (semean) se=leptin, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (ng/mL)", color(black) margin(medium)) ///
ylabel(0(10)40, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Leptin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:F.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_leptin", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=leptin (semean) se=leptin, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (ng/mL)", color(black) margin(medium)) ///
ylabel(0(10)40, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Leptin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:F.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_leptin", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=leptin (semean) se=leptin, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (ng/mL)", color(black) margin(medium)) ///
ylabel(0(10)40, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Leptin", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:C.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_leptin", replace
*******************
*Total PAEE & CRP *
*******************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress log_crp paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di (exp(_b[paeeTt])-1)*100 
di (exp((_b[paeeTt])+(1.96*_se[paeeTt]))-1)*100
di (exp((_b[paeeTt])-(1.96*_se[paeeTt]))-1)*100
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=crp (semean) se=crp, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mg/ml)", color(black) margin(medium)) ///
ylabel(0(1)6, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("CRP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:G.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_crp", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=crp (semean) se=crp, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mg/ml)", color(black) margin(medium)) ///
ylabel(0(1)6, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("CRP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:G.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_crp", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=crp (semean) se=crp, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (mg/ml)", color(black) margin(medium)) ///
ylabel(0(1)6, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("CRP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:F.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_crp", replace
**********************
*Total PAEE & MBPDIA *
**********************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress mbpdia paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di _b[paeeTt]
di _b[paeeTt]+(1.96*_se[paeeTt])
di _b[paeeTt]-(1.96*_se[paeeTt])
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=mbpdia (semean) se=mbpdia, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mm/Hg)", color(black) margin(medium)) ///
ylabel(68(2)80, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Diastolic BP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:H.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_mbpdia", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=mbpdia (semean) se=mbpdia, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mm/Hg)", color(black) margin(medium)) ///
ylabel(68(2)80, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Diastolic BP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:H.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_mbpdia", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=mbpdia (semean) se=mbpdia, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (mm/Hg)", color(black) margin(medium)) ///
ylabel(68(2)80, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Diastolic BP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:C.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_mbpdia", replace
**********************
*Total PAEE & MBPSYS *
**********************
*Table
clear

use "working1.dta"
forvalues x = 0/1 {
regress mbpsys paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med ///
if sex == `x'
di _b[paeeTt]
di _b[paeeTt]+(1.96*_se[paeeTt])
di _b[paeeTt]-(1.96*_se[paeeTt])
}
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=mbpsys (semean) se=mbpsys, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mm/Hg)", color(black) margin(medium)) ///
ylabel(110(5)135, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Systolic BP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:I.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_mbpsys", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=mbpsys (semean) se=mbpsys, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mm/Hg)", color(black) margin(medium)) ///
ylabel(110(5)135, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Systolic BP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:I.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_mbpsys", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=mbpsys (semean) se=mbpsys, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (mm/Hg)", color(black) margin(medium)) ///
ylabel(110(5)135, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Systolic BP", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:D.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_mbpsys", replace
***********************
*Total PAEE & GLUCOSE *
***********************
*Table
clear

use "working1.dta"
bysort sex: logistic glucose_cat paeeTt age i.ethnic i.testsite i.season ///
i.education i.income i.work_s smoking diet alcohol i.cardiometabol_med
*Female Graph
use "working1.dta", clear
drop if sex==1
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=glucose0 (semean) se=glucose0, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mmol/L)", color(black) margin(medium)) ///
ylabel(4.5(0.2)5.3, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Glucose", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:D.}",  position(11) ring(50) size(large))
graph save "F_totalPAEE_glucose", replace
*Male Graph
use "working1.dta", clear
drop if sex==0
xtile paeeTx = paeeTt, nq(5)
gen quint_meanPAEE = .
forvalues x = 1/5 {
egen  quint_meanPAEE`x' = mean(paeeTt) if paeeTx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEE = quint_meanPAEE`x' if paeeTx==`x'
}
collapse (mean) y=glucose0 (semean) se=glucose0, by(quint_meanPAEE)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEE, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEE, color(black)) ///
(connect y quint_meanPAEE, color(black) msymbol(Oh)), ///
ytitle("mean (mmol/L)", color(black) margin(medium)) ///
ylabel(4.5(0.2)5.3, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Glucose", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:D.}",  position(11) ring(50) size(large))
graph save "M_totalPAEE_glucose", replace
*Combined Graph
use "working1.dta", clear
xtile paeeTfx = paeeTt if sex==0, nq(5)
gen quint_meanPAEEf = .
xtile paeeTmx = paeeTt if sex==1, nq(5)
gen quint_meanPAEEm = .
forvalues x = 1/5 {
egen  quint_meanPAEEf`x' = mean(paeeTt) if paeeTfx==`x'
egen  quint_meanPAEEm`x' = mean(paeeTt) if paeeTmx==`x'
}
forvalues x = 1/5 {
replace quint_meanPAEEf = quint_meanPAEEf`x' if paeeTfx==`x'
replace quint_meanPAEEm = quint_meanPAEEm`x' if paeeTmx==`x'
}
collapse (mean) y=glucose0 (semean) se=glucose0, by(quint_meanPAEEf quint_meanPAEEm)
gen yu= y+1.96*se
gen yl= y-1.96*se
twoway (scatter y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(scatter y quint_meanPAEEm, color(black) msymbol(D)) ///
(rcap yu yl quint_meanPAEEf, color(black)) ///
(connect y quint_meanPAEEf, color(black) msymbol(Oh)) ///
(rcap yu yl quint_meanPAEEm, color(black)) ///
(connect y quint_meanPAEEm, color(black) lpatter(dash) msymbol(D)), ///
ytitle("mean (mmol/L)", color(black) margin(medium)) ///
ylabel(4.5(0.2)5.3, nogrid) xlabel(0(20)100) ///
legend(off) xtitle("") title("Glucose", color(black) margin(medium)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white)) ///
note("{bf:D.}",  position(11) ring(50) size(large))
graph save "C_totalPAEE_glucose", replace

*****************
*Combine graphs *
*****************
clear

use "working1.dta"
graph combine "F_totalPAEE_bodyfat" "F_totalPAEE_ldl" "F_totalPAEE_nefa" ///
"F_totalPAEE_insulin" "F_totalPAEE_adiponectin" "F_totalPAEE_leptin" ///
"F_totalPAEE_crp" "F_totalPAEE_mbpdia" "F_totalPAEE_mbpsys" ///
"F_totalPAEE_glucose", graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(5) b1("mean PAEE (kJ/day/kg) of quintiles", color(black) margin(medium)) ///
xcom
graph save "F_totalPAEE_combined",replace
graph export F_totalPAEE_combined.tif, replace

graph combine "M_totalPAEE_bodyfat" "M_totalPAEE_ldl" "M_totalPAEE_nefa" ///
"M_totalPAEE_insulin" "M_totalPAEE_adiponectin" "M_totalPAEE_leptin" ///
"M_totalPAEE_crp" "M_totalPAEE_mbpdia" "M_totalPAEE_mbpsys" ///
"M_totalPAEE_glucose", graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(5) b1("mean PAEE (kJ/day/kg) of quintiles", color(black) margin(medium)) ///
xcom
graph save "M_totalPAEE_combined", replace
graph export M_totalPAEE_combined.tif, replace

graph combine "C_totalPAEE_bodyfat" "C_totalPAEE_ldl" "C_totalPAEE_nefa" ///
"C_totalPAEE_insulin" "C_totalPAEE_adiponectin" "C_totalPAEE_leptin" ///
"C_totalPAEE_crp" "C_totalPAEE_mbpdia" "C_totalPAEE_mbpsys" ///
"C_totalPAEE_glucose", graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(5) b1("mean PAEE (kJ/day/kg) of quintiles" "(female = open circle & full line; male = filled diamond & dashed line)", color(black) margin(medium)) ///
xcom
graph save "C_totalPAEE_combined", replace
graph export C_totalPAEE_combined.tif, replace

grc1leg2 "C_totalPAEE_bodyfat" "C_totalPAEE_insulin" ///
"C_totalPAEE_mbpdia" "C_totalPAEE_glucose", graphregion(fcolor(white) ///
lcolor(white)) plotregion(fcolor(white)) cols(2) ///
b1("PAEE volume (kJ/day/kg) of quintiles", color(black) margin(med)) ///
legendfrom("C_totalPAEE_bodyfat") position(3) lcols(1)
graph save "C_totalPAEE_combined_paper", replace
graph export C_totalPAEE_combined_paper.tif, replace

clear

use "working1.dta"
grc1leg2 "C_totalPAEE_ldl" "C_totalPAEE_nefa" ///
"C_totalPAEE_leptin" "C_totalPAEE_mbpsys" "C_totalPAEE_adiponectin" ///
"C_totalPAEE_crp", graphregion(fcolor(white) ///
lcolor(white)) plotregion(fcolor(white)) cols(3) ///
b1("PAEE volume (kJ/day/kg) of quintiles", color(black) margin(med)) ///
legendfrom("C_totalPAEE_ldl") position(3) lcols(1)
graph save "C_totalPAEE_combined_supp1", replace
graph export C_totalPAEE_combined_supp1.tif, replace


********************************************************************************
* STEP #8: CoDA & Graphs (Figure 3)
********************************************************************************
clear

use "working1.dta"

******************
*CoDA & Bodyfat  *
******************
* Bodyfat from 4am-8am
bysort sex: regress bodyfat earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress bodyfat earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress bodyfat earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* Bodyfat from 8am-12noon
bysort sex: regress bodyfat earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress bodyfat earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress bodyfat earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* Bodyfat from 12noon-4pm
bysort sex: regress bodyfat earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress bodyfat earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress bodyfat earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* Bodyfat from 4pm-8pm
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* Bodyfat from 8pm-12midnight
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* Bodyfat from 12midnight-4am
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress bodyfat earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

**************
*CoDA & LDL  *
**************
* ldl from 4am-8am
bysort sex: regress ldl earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress ldl earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress ldl earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 8am-12noon
bysort sex: regress ldl earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress ldl earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress ldl earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 12noon-4pm
bysort sex: regress ldl earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress ldl earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress ldl earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 4pm-8pm
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 8pm-12midnight
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 12midnight-4am
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

*****************
*CoDA & Leptin  *
*****************
* log_leptin from 4am-8am
bysort sex: regress log_leptin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_leptin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_leptin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 8am-12noon
bysort sex: regress log_leptin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_leptin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_leptin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 12noon-4pm
bysort sex: regress log_leptin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_leptin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_leptin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 4pm-8pm
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 8pm-12midnight
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 12midnight-4am
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

**********************
*CoDA & Adiponectin  *
**********************
* log_adiponectin from 4am-8am
bysort sex: regress log_adiponectin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_adiponectin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_adiponectin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 8am-12noon
bysort sex: regress log_adiponectin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_adiponectin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_adiponectin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 12noon-4pm
bysort sex: regress log_adiponectin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_adiponectin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_adiponectin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 4pm-8pm
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 8pm-12midnight
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 12midnight-4am
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

****************
*CoDA & NEFAs  *
****************
* log_nefa from 4am-8am
bysort sex: regress log_nefa earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_nefa earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_nefa earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 8am-12noon
bysort sex: regress log_nefa earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_nefa earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_nefa earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 12noon-4pm
bysort sex: regress log_nefa earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_nefa earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_nefa earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 4pm-8pm
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 8pm-12midnight
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 12midnight-4am
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

*****************
*CoDA & Sys BP  *
*****************
* mbpsys from 4am-8am
bysort sex: regress mbpsys earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpsys earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpsys earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 8am-12noon
bysort sex: regress mbpsys earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpsys earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpsys earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 12noon-4pm
bysort sex: regress mbpsys earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpsys earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpsys earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 4pm-8pm
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 8pm-12midnight
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 12midnight-4am
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

*****************
*CoDA & Dia BP  *
*****************
* mbpdia from 4am-8am
bysort sex: regress mbpdia earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpdia earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpdia earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 8am-12noon
bysort sex: regress mbpdia earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpdia earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpdia earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 12noon-4pm
bysort sex: regress mbpdia earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpdia earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpdia earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 4pm-8pm
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 8pm-12midnight
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 12midnight-4am
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

**************
*CoDA & CRP  *
**************
* log_crp from 4am-8am
bysort sex: regress log_crp earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_crp earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_crp earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 8am-12noon
bysort sex: regress log_crp earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_crp earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_crp earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 12noon-4pm
bysort sex: regress log_crp earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_crp earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_crp earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 4pm-8pm
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 8pm-12midnight
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 12midnight-4am
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

******************
*CoDA & Insulin  *
******************
* log_insulin from 4am-8am
bysort sex: regress log_insulin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_insulin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_insulin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 8am-12noon
bysort sex: regress log_insulin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_insulin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_insulin earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 12noon-4pm
bysort sex: regress log_insulin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_insulin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_insulin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 4pm-8pm
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 8pm-12midnight
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 12midnight-4am
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med

******************
*CoDA & Glucose  *
******************
* glucose_cat from 4am-8am
bysort sex: logistic glucose_cat earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: logistic glucose_cat earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: logistic glucose_cat earlymidp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 8am-12noon
bysort sex: logistic glucose_cat earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: logistic glucose_cat earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: logistic glucose_cat earlyp5 midp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 12noon-4pm
bysort sex: logistic glucose_cat earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: logistic glucose_cat earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: logistic glucose_cat earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 4pm-8pm
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latep5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 8pm-12midnight
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season 
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 12midnight-4am
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season 
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeT age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med


**********************************
**# CoDA Graphs from excel (betas copied from analyses above*  
*		 | Female....|....Male   *
*Body fat| Sheet 1   |  Sheet 2  *
*LDL     | Sheet 3   |  Sheet 4  *
*Leptin  | Sheet 5   |  Sheet 6  *
*Adipon  | Sheet 7   |  Sheet 8  *
*NEFA    | Sheet 9   |  Sheet 10 *
*Sys BP  | Sheet 11  |  Sheet 12 *
*Dia BP  | Sheet 13  |  Sheet 14 *
*CRP     | Sheet 15  |  Sheet 16 *
*Insulin | Sheet 17  |  Sheet 18 *
*Glucose | Sheet 19  |  Sheet 20 *
**********************************

*Bodyfat Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet1") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -1.5 0 -1.5 2 1.5 2 1.5 0 -1.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -1.5 4 -1.5 6 1.5 6 1.5 4 -1.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -1.5 8 -1.5 10 1.5 10 1.5 8 -1.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-1.5(0.5)1.5, angle(0) nogrid norescale) yscale(titlegap(*-15)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin))  note("{bf:A.     Adiposity}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_bodyfat", replace
graph export f_coda_bodyfat.tif, replace

*Bodyfat Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet2") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -1.5 0 -1.5 2 1.5 2 1.5 0 -1.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -1.5 4 -1.5 6 1.5 6 1.5 4 -1.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -1.5 8 -1.5 10 1.5 10 1.5 8 -1.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-1.5(0.5)1.5, angle(0) nogrid norescale) yscale(titlegap(*-15)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "5% PAEE from 4am-8am" 4 "5% PAEE from 8am-12pm" 6 "5% PAEE from 12pm-4pm" 9 "5% PAEE from 4pm-8pm" 11 "5% PAEE from 8pm-12am" 14  "5% PAEE from 12am-4am") cols(1) pos(3))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:B.     Adiposity}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_bodyfat", replace
graph export m_coda_bodyfat.tif, replace

*LDL Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet3") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -0.15 0 -0.15 2 0.15 2 0.15 0 -0.15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -0.15 4 -0.15 6 0.15 6 0.15 4 -0.15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -0.15 8 -0.15 10 0.15 10 0.15 8 -0.15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmol/L ± 95% CI)", margin(medium) color(black)) ylabel(-0.15(0.05)0.15, angle(0) nogrid norescale) yscale(titlegap(*-35)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:A.     LDL}",  position(11) ring(50) size(large))   ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_ldl", replace
graph export f_coda_ldl.tif, replace

*LDL Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet4") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -0.15 0 -0.15 2 0.15 2 0.15 0 -0.15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -0.15 4 -0.15 6 0.15 6 0.15 4 -0.15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -0.15 8 -0.15 10 0.15 10 0.15 8 -0.15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmol/L ± 95% CI)", margin(medium) color(black)) ylabel(-0.15(0.05)0.15, angle(0) nogrid norescale) yscale(titlegap(*-35)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin))  note("{bf:B.     LDL}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_ldl", replace
graph export m_coda_ldl.tif, replace

*Leptin Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet5") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -15 0 -15 2 15 2 15 0 -15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -15 4 -15 6 15 6 15 4 -15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -15 8 -15 10 15 10 15 8 -15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-15(5)15, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:E.     Leptin}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_leptin", replace
graph export f_coda_leptin.tif, replace

*Leptin Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet6") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -15 0 -15 2 15 2 15 0 -15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -15 4 -15 6 15 6 15 4 -15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -15 8 -15 10 15 10 15 8 -15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-15(5)15, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:F.     Leptin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_leptin", replace
graph export m_coda_leptin.tif, replace

*Adiponectin Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet7") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -7.5 0 -7.5 2 7.5 2 7.5 0 -7.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -7.5 4 -7.5 6 7.5 6 7.5 4 -7.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -7.5 8 -7.5 10 7.5 10 7.5 8 -7.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-7.5(2.5)7.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:C.     Adiponectin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adiponectin", replace
graph export f_coda_adiponectin.tif, replace

*Adiponectin Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet8") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -7.5 0 -7.5 2 7.5 2 7.5 0 -7.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -7.5 4 -7.5 6 7.5 6 7.5 4 -7.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -7.5 8 -7.5 10 7.5 10 7.5 8 -7.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-7.5(2.5)7.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:D.     Adiponectin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adiponectin", replace
graph export m_coda_adiponectin.tif, replace

*NEFA Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet9") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:C.     NEFAs}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_nefa", replace
graph export f_coda_nefa.tif, replace

*NEFA Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet10") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:D.     NEFAs}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_nefa", replace
graph export m_coda_nefa.tif, replace

*SysBP Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet11") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -4 0 -4 2 4 2 4 0 -4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -4 4 -4 6 4 6 4 4 -4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -4 8 -4 10 4 10 4 8 -4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-4(1)4, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:A.     Systolic BP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_sysbp", replace
graph export f_coda_sysbp.tif, replace

*SysBP Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet12") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -4 0 -4 2 4 2 4 0 -4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -4 4 -4 6 4 6 4 4 -4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -4 8 -4 10 4 10 4 8 -4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-4(1)4, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:B.     Systolic BP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_sysbp", replace
graph export m_coda_sysbp.tif, replace

*DiaBP Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet13") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -2.5 0 -2.5 2 2.5 2 2.5 0 -2.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -2.5 4 -2.5 6 2.5 6 2.5 4 -2.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -2.5 8 -2.5 10 2.5 10 2.5 8 -2.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-2.5(0.5)2.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:C.     Diastolic BP}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_diabp", replace
graph export f_coda_diabp.tif, replace

*DiaBP Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet14") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -2.5 0 -2.5 2 2.5 2 2.5 0 -2.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -2.5 4 -2.5 6 2.5 6 2.5 4 -2.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -2.5 8 -2.5 10 2.5 10 2.5 8 -2.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-2.5(0.5)2.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:D.     Diastolic BP}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_diabp", replace
graph export m_coda_diabp.tif, replace

*CRP Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet15") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -20 0 -20 2 20 2 20 0 -20 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -20 4 -20 6 20 6 20 4 -20 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -20 8 -20 10 20 10 20 8 -20 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-20(5)20, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:E.     CRP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_crp", replace
graph export f_coda_crp.tif, replace

*CRP Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet16") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -20 0 -20 2 20 2 20 0 -20 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -20 4 -20 6 20 6 20 4 -20 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -20 8 -20 10 20 10 20 8 -20 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-20(5)20, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:F.     CRP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_crp", replace
graph export m_coda_crp.tif, replace

*Insulin Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet17") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:E.     Insulin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_insulin", replace
graph export f_coda_insulin.tif, replace

*Insulin Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet18") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:F.     Insulin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_insulin", replace
graph export m_coda_insulin.tif, replace

*Glucose Female
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet19") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri 0.4 0 0.4 2 2.2 2 2.2 0 0.4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri 0.4 4 0.4 6 2.2 6 2.2 4 0.4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri 0.4 8 0.4 10 2.2 10 2.2 8 0.4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (OR ± 95% CI)", margin(medium) color(black)) ylabel(0.4(0.2)2.2, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(1, lcolor(black) lwidth(thin)) note("{bf:G.     High Glucose}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_glucose", replace
graph export f_coda_glucose.tif, replace

*Glucose Male
clear

import excel using "Coda_data for graphs_27.10.2023.xlsx", sheet("Sheet20") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri 0.4 0 0.4 2 2.2 2 2.2 0 0.4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri 0.4 4 0.4 6 2.2 6 2.2 4 0.4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri 0.4 8 0.4 10 2.2 10 2.2 8 0.4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (OR ± 95% CI)", margin(medium) color(black)) ylabel(0.4(0.2)2.2, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(1, lcolor(black) lwidth(thin)) note("{bf:H.     High Glucose}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_glucose", replace
graph export m_coda_glucose.tif, replace


*Graph for Paper
clear

grc1leg2 "F_coda_bodyfat" "F_coda_insulin" "F_coda_diabp" "F_coda_glucose", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Female, color(black) size(medsmall)) ///
legendfrom("F_coda_bodyfat") pos(3) lcols(1) name(gr1, replace)
grc1leg2 "M_coda_bodyfat" "M_coda_insulin" "M_coda_diabp" "M_coda_glucose", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Male, color(black) size(medsmall)) ///
legendfrom("M_coda_bodyfat") pos(3) lcols(1) name(gr2, replace)
grc1leg2  gr1 gr2, col(2) graphregion(fcolor(white) lcolor(white)) ///
plotregion(fcolor(white)) legendfrom(gr1) pos(3) lcols(1) ///
b1("X-axes: Time blocks into which 5% PAEE is redistributed", color(black) margin(small) size(vsmall)) 
graph save "C_coda_combined_paper", replace
graph export C_coda_combined_paper.tif, replace

*Graphs for Supplement
clear

grc1leg2 "F_coda_ldl" "F_coda_nefa" "F_coda_leptin", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Female, color(black) size(medsmall)) ///
legendfrom("F_coda_ldl") pos(3) lcols(1) name(gr1, replace)
grc1leg2 "M_coda_ldl" "M_coda_nefa" "M_coda_leptin", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Male, color(black) size(medsmall)) ///
legendfrom("M_coda_ldl") pos(3) lcols(1) name(gr2, replace)
grc1leg2  gr1 gr2, col(2) graphregion(fcolor(white) lcolor(white)) ///
plotregion(fcolor(white)) legendfrom(gr1) pos(3) lcols(1) ///
b1("X-axes: Time blocks into which 5% PAEE is redistributed", color(black) margin(small) size(vsmall)) 
graph save "C_coda_combined_supp1", replace
graph export C_coda_combined_supp1.tif, replace

clear

grc1leg2 "F_coda_sysbp" "F_coda_adiponectin" "F_coda_crp", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Female, color(black) size(medsmall)) ///
legendfrom("F_coda_sysbp") pos(3) lcols(1) name(gr1, replace)
grc1leg2 "M_coda_sysbp" "M_coda_adiponectin" "M_coda_crp", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Male, color(black) size(medsmall)) ///
legendfrom("M_coda_sysbp") pos(3) lcols(1) name(gr2, replace)
grc1leg2  gr1 gr2, col(2) graphregion(fcolor(white) lcolor(white)) ///
plotregion(fcolor(white)) legendfrom(gr1) pos(3) lcols(1) ///
b1("X-axes: Time blocks into which 5% PAEE is redistributed", color(black) margin(small) size(vsmall)) 
graph save "C_coda_combined_supp2", replace
graph export C_coda_combined_supp2.tif, replace

********************************************************************************
* STEP #9: CoDA & Graphs - Adiposity Adjusted (Suppl)
********************************************************************************
clear
use "working1.dta"
**************
*CoDA & LDL  *
**************
* ldl from 4am-8am
bysort sex: regress ldl earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 8am-12noon
bysort sex: regress ldl earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 12noon-4pm
bysort sex: regress ldl earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 4pm-8pm
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 8pm-12midnight
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* ldl from 12midnight-4am
bysort sex: regress ldl earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
*****************
*CoDA & Leptin  *
*****************
* log_leptin from 4am-8am
bysort sex: regress log_leptin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 8am-12noon
bysort sex: regress log_leptin earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 12noon-4pm
bysort sex: regress log_leptin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 4pm-8pm
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 8pm-12midnight
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_leptin from 12midnight-4am
bysort sex: regress log_leptin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
**********************
*CoDA & Adiponectin  *
**********************
* log_adiponectin from 4am-8am
bysort sex: regress log_adiponectin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 8am-12noon
bysort sex: regress log_adiponectin earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 12noon-4pm
bysort sex: regress log_adiponectin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 4pm-8pm
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 8pm-12midnight
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_adiponectin from 12midnight-4am
bysort sex: regress log_adiponectin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
****************
*CoDA & NEFAs  *
****************
* log_nefa from 4am-8am
bysort sex: regress log_nefa earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 8am-12noon
bysort sex: regress log_nefa earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 12noon-4pm
bysort sex: regress log_nefa earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 4pm-8pm
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 8pm-12midnight
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_nefa from 12midnight-4am
bysort sex: regress log_nefa earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
*****************
*CoDA & Sys BP  *
*****************
* mbpsys from 4am-8am
bysort sex: regress mbpsys earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 8am-12noon
bysort sex: regress mbpsys earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 12noon-4pm
bysort sex: regress mbpsys earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 4pm-8pm
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 8pm-12midnight
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpsys from 12midnight-4am
bysort sex: regress mbpsys earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
*****************
*CoDA & Dia BP  *
*****************
* mbpdia from 4am-8am
bysort sex: regress mbpdia earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 8am-12noon
bysort sex: regress mbpdia earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 12noon-4pm
bysort sex: regress mbpdia earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 4pm-8pm
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 8pm-12midnight
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* mbpdia from 12midnight-4am
bysort sex: regress mbpdia earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
**************
*CoDA & CRP  *
**************
* log_crp from 4am-8am
bysort sex: regress log_crp earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 8am-12noon
bysort sex: regress log_crp earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 12noon-4pm
bysort sex: regress log_crp earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 4pm-8pm
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 8pm-12midnight
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_crp from 12midnight-4am
bysort sex: regress log_crp earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
******************
*CoDA & Insulin  *
******************
* log_insulin from 4am-8am
bysort sex: regress log_insulin earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 8am-12noon
bysort sex: regress log_insulin earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 12noon-4pm
bysort sex: regress log_insulin earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 4pm-8pm
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 8pm-12midnight
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* log_insulin from 12midnight-4am
bysort sex: regress log_insulin earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
******************
*CoDA & Glucose  *
******************
* glucose_cat from 4am-8am
bysort sex: logistic glucose_cat earlymidp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 8am-12noon
bysort sex: logistic glucose_cat earlyp5 midp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 12noon-4pm
bysort sex: logistic glucose_cat earlyp5 earlymidp5 latemidp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 4pm-8pm
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latep5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 8pm-12midnight
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 nightp5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med
* glucose_cat from 12midnight-4am
bysort sex: logistic glucose_cat earlyp5 earlymidp5 midp5 latemidp5 latep5 paeeTt bodyfat age i.ethnic i.testsite i.season i.education i.income i.worktype smoking diet alcohol i.cardiometabol_med


**********************************
**# CoDA Graphs (bodyfat adjusted*  
*		 | Female....|....Male   *
*LDL     | Sheet 3   |  Sheet 4  *
*Leptin  | Sheet 5   |  Sheet 6  *
*Adipon  | Sheet 7   |  Sheet 8  *
*NEFA    | Sheet 9   |  Sheet 10 *
*Sys BP  | Sheet 11  |  Sheet 12 *
*Dia BP  | Sheet 13  |  Sheet 14 *
*CRP     | Sheet 15  |  Sheet 16 *
*Insulin | Sheet 17  |  Sheet 18 *
*Glucose | Sheet 19  |  Sheet 20 *
**********************************

*LDL Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet3") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -0.15 0 -0.15 2 0.15 2 0.15 0 -0.15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -0.15 4 -0.15 6 0.15 6 0.15 4 -0.15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -0.15 8 -0.15 10 0.15 10 0.15 8 -0.15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmol/L ± 95% CI)", margin(medium) color(black)) ylabel(-0.15(0.05)0.15, angle(0) nogrid norescale) yscale(titlegap(*-35)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:A.     LDL}",  position(11) ring(50) size(large))   ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_ldl", replace
graph export F_coda_adipo_ldl.tif, replace

*LDL Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet4") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -0.15 0 -0.15 2 0.15 2 0.15 0 -0.15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -0.15 4 -0.15 6 0.15 6 0.15 4 -0.15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -0.15 8 -0.15 10 0.15 10 0.15 8 -0.15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmol/L ± 95% CI)", margin(medium) color(black)) ylabel(-0.15(0.05)0.15, angle(0) nogrid norescale) yscale(titlegap(*-35)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin))  note("{bf:B.     LDL}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_ldl", replace
graph export M_coda_adipo_ldl.tif, replace

*Leptin Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet5") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -15 0 -15 2 15 2 15 0 -15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -15 4 -15 6 15 6 15 4 -15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -15 8 -15 10 15 10 15 8 -15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-15(5)15, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:E.     Leptin}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_leptin", replace
graph export F_coda_adipo_leptin.tif, replace

*Leptin Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet6") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -15 0 -15 2 15 2 15 0 -15 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -15 4 -15 6 15 6 15 4 -15 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -15 8 -15 10 15 10 15 8 -15 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-15(5)15, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:F.     Leptin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_leptin", replace
graph export M_coda_adipo_leptin.tif, replace

*Adiponectin Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet7") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -7.5 0 -7.5 2 7.5 2 7.5 0 -7.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -7.5 4 -7.5 6 7.5 6 7.5 4 -7.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -7.5 8 -7.5 10 7.5 10 7.5 8 -7.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-7.5(2.5)7.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:C.     Adiponectin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_adiponectin", replace
graph export F_coda_adipo_adiponectin.tif, replace

*Adiponectin Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet8") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -7.5 0 -7.5 2 7.5 2 7.5 0 -7.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -7.5 4 -7.5 6 7.5 6 7.5 4 -7.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -7.5 8 -7.5 10 7.5 10 7.5 8 -7.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-7.5(2.5)7.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:D.     Adiponectin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_adiponectin", replace
graph export M_coda_adipo_adiponectin.tif, replace

*NEFA Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet9") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin))  note("{bf:C.     NEFAs}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_nefa", replace
graph export F_coda_adipo_nefa.tif, replace

*NEFA Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet10") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:D.     NEFAs}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_nefa", replace
graph export M_coda_adipo_nefa.tif, replace

*SysBP Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet11") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -4 0 -4 2 4 2 4 0 -4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -4 4 -4 6 4 6 4 4 -4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -4 8 -4 10 4 10 4 8 -4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-4(1)4, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:A.     Systolic BP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_sysbp", replace
graph export F_coda_adipo_sysbp.tif, replace

*SysBP Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet12") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -4 0 -4 2 4 2 4 0 -4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -4 4 -4 6 4 6 4 4 -4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -4 8 -4 10 4 10 4 8 -4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-4(1)4, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:B.     Systolic BP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_sysbp", replace
graph export M_coda_adipo_sysbp.tif, replace

*DiaBP Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet13") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -2.5 0 -2.5 2 2.5 2 2.5 0 -2.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -2.5 4 -2.5 6 2.5 6 2.5 4 -2.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -2.5 8 -2.5 10 2.5 10 2.5 8 -2.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-2.5(0.5)2.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:C.     Diastolic BP}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_diabp", replace
graph export F_coda_adipo_diabp.tif, replace

*DiaBP Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet14") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -2.5 0 -2.5 2 2.5 2 2.5 0 -2.5 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -2.5 4 -2.5 6 2.5 6 2.5 4 -2.5 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -2.5 8 -2.5 10 2.5 10 2.5 8 -2.5 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (mmHg ± 95% CI)", margin(medium) color(black)) ylabel(-2.5(0.5)2.5, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:D.     Diastolic BP}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_diabp", replace
graph export M_coda_adipo_diabp.tif, replace

*CRP Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet15") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -20 0 -20 2 20 2 20 0 -20 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -20 4 -20 6 20 6 20 4 -20 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -20 8 -20 10 20 10 20 8 -20 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-20(5)20, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:E.     CRP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_crp", replace
graph export F_coda_adipo_crp.tif, replace

*CRP Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet16") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -20 0 -20 2 20 2 20 0 -20 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -20 4 -20 6 20 6 20 4 -20 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -20 8 -20 10 20 10 20 8 -20 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-20(5)20, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(0, lcolor(black) lwidth(thin)) note("{bf:F.     CRP}",  position(11) ring(50) size(large))  ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_crp", replace
graph export M_coda_adipo_crp.tif, replace

*Insulin Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet17") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin))  note("{bf:A.     Insulin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_insulin", replace
graph export F_coda_adipo_insulin.tif, replace

*Insulin Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet18") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri -10 0 -10 2 10 2 10 0 -10 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri -10 4 -10 6 10 6 10 4 -10 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri -10 8 -10 10 10 10 10 8 -10 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (% ± 95% CI)", margin(medium) color(black)) ylabel(-10(2.5)10, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(order(1 "4am-8am" 4 "8am-12pm" 6 "12pm-4pm" 9 "4pm-8pm" 11 "8pm-12am" 14  "12am-4am") cols(1) pos(3) ///
title("Markers indicate the time block" "from where 5% PAEE is taken:", color(black) size(tiny)) size(tiny))  ///
yline(0, lcolor(black) lwidth(thin)) note("{bf:B.     Insulin}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_insulin", replace
graph export M_coda_adipo_insulin.tif, replace

*Glucose Female
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet19") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri 0.4 0 0.4 2 2.2 2 2.2 0 0.4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri 0.4 4 0.4 6 2.2 6 2.2 4 0.4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri 0.4 8 0.4 10 2.2 10 2.2 8 0.4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (OR ± 95% CI)", margin(medium) color(black)) ylabel(0.4(0.2)2.2, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(1, lcolor(black) lwidth(thin)) note("{bf:E.     High Glucose}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "F_coda_adipo_glucose", replace
graph export F_coda_adipo_glucose.tif, replace

*Glucose Male
clear

import excel using "Coda_AipoAdj_data for graphs_21.11.2023.xlsx", sheet("Sheet20") firstrow clear
replace time=time+5 if time==6 //11
replace time=time+4 if time==5 //9
replace time=time+3 if time==4 //7
replace time=time+2 if time==3 //5
replace time=time+1 if time==2 //3
replace time=time-0.60 if sub==1
replace time=time-0.40 if sub==2
replace time=time-0.20 if sub==3
replace time=time+0.0 if sub==4
replace time=time+0.20 if sub==5
replace time=time+0.40 if sub==6
twoway (scatter coeff time if sub==1, color(black) msymbol(O)) (rcap lci uci time if sub==1, color(black) lwidth(thin)) ///
(scatteri 0.4 0 0.4 2 2.2 2 2.2 0 0.4 0, recast(area) color(black%10)) ///
(scatter coeff time if sub==2, color(red) msymbol(X)) (rcap lci uci time if sub==2, color(red) lwidth(thin)) ///
(scatter coeff time if sub==3, color(blue) msymbol(T)) (rcap lci uci time if sub==3, color(blue) lwidth(thin)) ///
(scatteri 0.4 4 0.4 6 2.2 6 2.2 4 0.4 4, recast(area) color(black%10)) ///
(scatter coeff time if sub==4, color(orange) msymbol(D)) (rcap lci uci time if sub==4, color(orange) lwidth(thin)) ///
(scatter coeff time if sub==5, color(green) msymbol(V)) (rcap lci uci time if sub==5, color(green) lwidth(thin)) ///
(scatteri 0.4 8 0.4 10 2.2 10 2.2 8 0.4 8, recast(area) color(black%10)) ///
(scatter coeff time if sub==6, color(purple) msymbol(S)) (rcap lci uci time if sub==6, color(purple) lwidth(thin)), ///
ytitle("Difference (OR ± 95% CI)", margin(medium) color(black)) ylabel(0.4(0.2)2.2, angle(0) nogrid norescale) yscale(titlegap(*-10)) ///
xtitle("", margin(medium) color(black)) xlabel(1 "4am-8am" 3 "8am-12pm" 5 "12pm-4pm" 7 "4pm-8pm" 9 "8pm-12am" 11 "12am-4am", angle(0)) ///
legend(off) yline(1, lcolor(black) lwidth(thin)) note("{bf:F.     High Glucose}",  position(11) ring(50) size(large)) ///
graphregion(fcolor(white) lcolor(white)) bgcolor(white) plotregion(fcolor(white) margin(b=0 t=0 l=0 4=0))
graph save "M_coda_adipo_glucose", replace
graph export M_coda_adipo_glucose.tif, replace


*Graph for Paper
clear

grc1leg2 "F_coda_adipo_insulin" "F_coda_adipo_diabp" "F_coda_adipo_glucose", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Female, color(black) size(medsmall)) ///
legendfrom("F_coda_adipo_insulin") pos(3) lcols(1) name(gr1, replace)
grc1leg2 "M_coda_adipo_insulin" "M_coda_adipo_diabp" "M_coda_adipo_glucose", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Male, color(black) size(medsmall)) ///
legendfrom("M_coda_adipo_insulin") pos(3) lcols(1) name(gr2, replace)
grc1leg2  gr1 gr2, col(2) graphregion(fcolor(white) lcolor(white)) ///
plotregion(fcolor(white)) legendfrom(gr1) pos(3) lcols(1) ///
b1("X-axes: Time blocks into which 5% PAEE is redistributed", color(black) margin(small) size(vsmall)) 
graph save "C_coda_adipo_combined_paper", replace
graph export C_coda_adipo_combined_paper.tif, replace

*Graphs for Supplement
clear

grc1leg2 "F_coda_adipo_ldl" "F_coda_adipo_nefa" "F_coda_adipo_leptin", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Female, color(black) size(medsmall)) ///
legendfrom("F_coda_adipo_ldl") pos(3) lcols(1) name(gr1, replace)
grc1leg2 "M_coda_adipo_ldl" "M_coda_adipo_nefa" "M_coda_adipo_leptin", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Male, color(black) size(medsmall)) ///
legendfrom("M_coda_adipo_ldl") pos(3) lcols(1) name(gr2, replace)
grc1leg2  gr1 gr2, col(2) graphregion(fcolor(white) lcolor(white)) ///
plotregion(fcolor(white)) legendfrom(gr1) pos(3) lcols(1) ///
b1("X-axes: Time blocks into which 5% PAEE is redistributed", color(black) margin(small) size(vsmall)) 
graph save "C_coda_adipo_combined_supp1", replace
graph export C_coda_adipo_combined_supp1.tif, replace

clear

grc1leg2 "F_coda_adipo_sysbp" "F_coda_adipo_adiponectin" "F_coda_adipo_crp", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Female, color(black) size(medsmall)) ///
legendfrom("F_coda_adipo_sysbp") pos(3) lcols(1) name(gr1, replace)
grc1leg2 "M_coda_adipo_sysbp" "M_coda_adipo_adiponectin" "M_coda_adipo_crp", ///
graphregion(fcolor(white) lcolor(white)) plotregion(fcolor(white)) ///
cols(1) iscale(0.5) title(Male, color(black) size(medsmall)) ///
legendfrom("M_coda_adipo_sysbp") pos(3) lcols(1) name(gr2, replace)
grc1leg2  gr1 gr2, col(2) graphregion(fcolor(white) lcolor(white)) ///
plotregion(fcolor(white)) legendfrom(gr1) pos(3) lcols(1) ///
b1("X-axes: Time blocks into which 5% PAEE is redistributed", color(black) margin(small) size(vsmall)) 
graph save "C_coda_adipo_combined_supp2", replace
graph export C_coda_adipo_combined_supp2.tif, replace
*/

