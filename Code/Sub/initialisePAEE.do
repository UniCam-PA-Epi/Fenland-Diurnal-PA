********************************************************************************
* STEP #1: Organise PAEE
********************************************************************************

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