********************************************************************************
* STEP #1: Organise PAEE
********************************************************************************

egen paeeT = rsum(paee_hour*)

keep ID paeeT paee_hour*

reshape long paee_hour, i(ID) j(hour)
asdf
foreach i in 24{

    gen sin`i' = sin(hour*2*_pi/`i')
    gen cos`i' = cos(hour*2*_pi/`i')
}

gen sin24_hat       = .
gen sin24_se        = .
gen sin24_lb        = .
gen sin24_ub        = .
gen sin24_p         = .

gen cos24_hat       = .
gen cos24_se        = .
gen cos24_lb        = .
gen cos24_ub        = .
gen cos24_p         = .

gen acrophase_hat   = .
gen acrophase_se    = .
gen acrophase_lb    = .
gen acrophase_ub    = .
gen acrophase_p     = .

gen amplitude_hat   = .
gen amplitude_se    = .
gen amplitude_lb    = .
gen amplitude_ub    = .
gen amplitude_p     = .

gen max_hat          = .    
gen max_se           = .         
gen max_lb           = .       
gen max_ub           = .   
gen max_p            = .

gen min_hat          = .    
gen min_se           = .         
gen min_lb           = .       
gen min_ub           = .   
gen min_p            = .

levelsof ID, local(IDList)
foreach curID of local IDList{

    glm paee_hour sin24 cos24 if ID == "`curID'", family(gaussian) link(log)
    estimates store m1

    predict x if ID == "`curID'"

    margin, dydx(sin24 cos24) post

    replace sin24_hat       = _b[sin24]                                     if ID == "`curID'"
    replace sin24_se        = _se[sin24]                                    if ID == "`curID'"
    replace sin24_lb        = _b[sin24] + invnormal(0.025)*_se[sin24]       if ID == "`curID'"
    replace sin24_ub        = _b[sin24] + invnormal(0.975)*_se[sin24]       if ID == "`curID'"
    replace sin24_p         = 2*normal(-abs(_b[sin24]/_se[sin24]))          if ID == "`curID'"

    replace cos24_hat       = _b[cos24]                                     if ID == "`curID'"
    replace cos24_se        = _se[cos24]                                    if ID == "`curID'"
    replace cos24_lb        = _b[cos24] + invnormal(0.025)*_se[cos24]       if ID == "`curID'"
    replace cos24_ub        = _b[cos24] + invnormal(0.975)*_se[cos24]       if ID == "`curID'"
    replace cos24_p         = 2*normal(-abs(_b[cos24]/_se[cos24]))          if ID == "`curID'"

    nlcom cond(_b[sin24]<0,24,0)+atan2(_b[sin24],_b[cos24])*24/(2*_pi) 
    replace acrophase_hat   = r(b)[1,1]                                         if ID == "`curID'"
    replace acrophase_se    = sqrt(r(V)[1,1])                                   if ID == "`curID'"
    replace acrophase_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace acrophase_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace acrophase_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         if ID == "`curID'"
    local curAcro = r(b)[1,1]
    
    estimates restore m1
    margin, at( sin24=(`=sin(`curAcro'*2*_pi/24-_pi)',`=sin(`curAcro'*2*_pi/24)')   ///
                cos24=(`=cos(`curAcro'*2*_pi/24-_pi)',`=cos(`curAcro'*2*_pi/24)')   ///
              ) post
    
    nlcom  _b[4._at]
    replace max_hat         = r(b)[1,1]                                         if ID == "`curID'"
    replace max_se          = sqrt(r(V)[1,1])                                   if ID == "`curID'"
    replace max_lb          = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace max_ub          = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace max_p           = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         if ID == "`curID'"

    nlcom  _b[1._at]
    replace min_hat         = r(b)[1,1]                                         if ID == "`curID'"
    replace min_se          = sqrt(r(V)[1,1])                                   if ID == "`curID'"
    replace min_lb          = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace min_ub          = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace min_p           = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         if ID == "`curID'"

    nlcom (_b[4._at]-_b[1._at])/2
    replace amplitude_hat   = r(b)[1,1]                                         if ID == "`curID'"
    replace amplitude_se    = sqrt(r(V)[1,1])                                   if ID == "`curID'"
    replace amplitude_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace amplitude_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      if ID == "`curID'"
    replace amplitude_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         if ID == "`curID'"

   // keep if ID == "R16601063"

    asdf
}














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