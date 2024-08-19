version 17.0

**********************************
** Initialise outcome variables **
**********************************




********************************
** Fat mass and fat free mass **
********************************

rename FM fatMass
rename FFM fatFreeMass

*************
** Insulin **
*************

rename G_Insulin insulin

*********
** CRP **
*********

rename G_hs_CRP crp
replace crp = . if G_hs_CRP_Threshold==2


************
** Leptin **
************

rename G_Leptin leptin
replace leptin = . if G_Leptin_Threshold!=1

**********
** NEFA **
**********

rename G_NEFA nefa
replace nefa = . if G_NEFA_Threshold!=1

*****************
** Adiponectin **
*****************

rename G_Adiponectin adiponectin
replace adiponectin = . if G_Adiponectin_Threshold!=1

******************
** Blood lipids **
******************

rename LDL0     ldl
rename HDL0     hdl
rename Chol0    cholesterol

********************************************
** Blood Pressures: generate median value **
********************************************

egen mbpsys = rowmedian(BPSys1 BPSys2 BPSys3)
egen mbpdia = rowmedian(BPDia1 BPDia2 BPDia3) 

************* 
** Glucose **
************* 

gen glucose0 = cond(RepeatGlucose0<=Glucose0,RepeatGlucose0,Glucose0)
gen glucose120 = cond(RepeatGlucose120<=Glucose120,RepeatGlucose120,Glucose120)


