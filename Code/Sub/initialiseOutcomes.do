********************************************************************************
* STEP #2: Organise Outomes
********************************************************************************
//Bodyfat
rename BodyFatPercent_Cons bodyfat
replace bodyfat=. if bodyfat==-1 //n=9 missing

// Fat mass and fat free mass

rename FM fatMass
rename FFM fatFreeMass

gen log_fatMass = ln(fatMass)
gen log_fatFreeMass = ln(fatFreeMass)


//Insulin & log_Insulin
rename G_Insulin insulin
tab insulin, missing
//one extremely high value (1370, next is 747). n=2,323 (19.74%) missing
gen log_insulin=ln(insulin)

//CRP & log_CRP
rename G_hs_CRP crp
replace crp=. if G_hs_CRP_Threshold==2
mdesc crp // n=2,546 missing (21.64%)
gen log_crp=ln(crp)

//Leptin & log_Leptin
rename G_Leptin leptin
replace leptin=. if G_Leptin_Threshold!=1
mdesc leptin // n=2,348 missing (19.96%)
gen log_leptin=ln(leptin)

//NEFA & log_NEFA
rename G_NEFA nefa
replace nefa=. if G_NEFA_Threshold!=1
mdesc nefa // n=2,345 missing (19.93%)
gen log_nefa=ln(nefa)

//Adiponectin & log_Adiponectin
rename G_Adiponectin adiponectin
replace adiponectin=. if G_Adiponectin_Threshold!=1
mdesc adiponectin // n=2,340 missing (19.89%)
gen log_adiponectin = ln(adiponectin)

//LDL
rename LDL0 ldl
mdesc ldl // n=158 missing (1,34%)

rename HDL0 hdl
rename Chol0 cholesterol

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

// new glucose //

gen glucose0 = cond(RepeatGlucose0<=Glucose0,RepeatGlucose0,Glucose0)
gen log_glucose0 = ln(glucose0)

/*
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
*/