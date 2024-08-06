********************************************************************************
** Redo Cosinor Modelling
********************************************************************************


******************************************************
** Initialise dataset + outcome modelling variables **
******************************************************

keep ID paee_hour*
reshape long paee_hour, i(ID) j(hour)
drop if paee_hour == .

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

gen amplitude_hat   = .
gen amplitude_se    = .
gen amplitude_lb    = .
gen amplitude_ub    = .
gen amplitude_p     = .

gen paeet_hat       = .
gen paeet_se        = .
gen paeet_lb        = .
gen paeet_ub        = .
gen paeet_p         = .

gen qfactor_hat      = .
gen qfactor_se       = .         
gen qfactor_lb       = .       
gen qfactor_ub       = .   
gen qfactor_p        = .

gen upperOutlier    = .
gen lowerOutlier    = .



levelsof ID, local(IDList)

/*
foreach curID of local IDList{

    su paee_hour if ID == "`curID'", detail
    replace upperOutlier = cond((paee_hour - r(p75)) / (r(p75) - r(p25)) >=  2 ,1,0) if ID == "`curID'",
    replace lowerOutlier = cond((paee_hour - r(p25)) / (r(p75) - r(p25)) <= -2 ,1,0) if ID == "`curID'",

}

drop if upperOutlier == 1
drop if lowerOutlier == 1
drop upperOutlier lowerOutlier
*/

foreach curID of local IDList{

    count if ID == "`curID'"

    if r(N) >= 18{

        preserve
        keep if ID == "`curID'"
        
        ************************************************************************************************
        ** Apply cosinor GLM model, using log-link function to ensure paee values are always positive **
        ************************************************************************************************

        glm paee_hour sin24 cos24, family(gaussian) link(log)
        estimates store m1

        *******************************************************
        ** Estimate real-values for sin and cos coefficients **
        *******************************************************

        margin, dydx(sin24 cos24)

        replace sin24_hat       = _b[sin24]                                         
        replace sin24_se        = _se[sin24]                                        
        replace sin24_lb        = _b[sin24] + invnormal(0.025)*_se[sin24]           
        replace sin24_ub        = _b[sin24] + invnormal(0.975)*_se[sin24]           
        replace sin24_p         = 2*normal(-abs(_b[sin24]/_se[sin24]))              

        replace cos24_hat       = _b[cos24]                                         
        replace cos24_se        = _se[cos24]                                        
        replace cos24_lb        = _b[cos24] + invnormal(0.025)*_se[cos24]           
        replace cos24_ub        = _b[cos24] + invnormal(0.975)*_se[cos24]           
        replace cos24_p         = 2*normal(-abs(_b[cos24]/_se[cos24]))              

        ************************
        ** Estimate acrophase **
        ************************

        nlcom cond(_b[sin24]<0,24,0)+atan2(_b[sin24],_b[cos24])*24/(2*_pi) 
        replace acrophase_hat   = r(b)[1,1]                                         
        replace acrophase_se    = sqrt(r(V)[1,1])                                   
        replace acrophase_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
        replace acrophase_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
        replace acrophase_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         
        local curAcro = r(b)[1,1]

        *************************************************************************************************************
        ** Get predicted values at max and min hour (at 25 and 26), as well as hours 1 through 24 (at 1 through 24) **
        *************************************************************************************************************

        estimates restore m1
        margin, at(sin24=(`=sin(1*2*_pi/24)')               cos24=(`=cos(1*2*_pi/24)'))             ///
                at(sin24=(`=sin(2*2*_pi/24)')               cos24=(`=cos(2*2*_pi/24)'))             ///
                at(sin24=(`=sin(3*2*_pi/24)')               cos24=(`=cos(3*2*_pi/24)'))             ///
                at(sin24=(`=sin(4*2*_pi/24)')               cos24=(`=cos(4*2*_pi/24)'))             /// 
                at(sin24=(`=sin(5*2*_pi/24)')               cos24=(`=cos(5*2*_pi/24)'))             ///
                at(sin24=(`=sin(6*2*_pi/24)')               cos24=(`=cos(6*2*_pi/24)'))             ///
                at(sin24=(`=sin(7*2*_pi/24)')               cos24=(`=cos(7*2*_pi/24)'))             ///
                at(sin24=(`=sin(8*2*_pi/24)')               cos24=(`=cos(8*2*_pi/24)'))             ///
                at(sin24=(`=sin(9*2*_pi/24)')               cos24=(`=cos(9*2*_pi/24)'))             ///
                at(sin24=(`=sin(10*2*_pi/24)')              cos24=(`=cos(10*2*_pi/24)'))            ///
                at(sin24=(`=sin(11*2*_pi/24)')              cos24=(`=cos(11*2*_pi/24)'))            ///
                at(sin24=(`=sin(12*2*_pi/24)')              cos24=(`=cos(12*2*_pi/24)'))            /// 
                at(sin24=(`=sin(13*2*_pi/24)')              cos24=(`=cos(13*2*_pi/24)'))            ///
                at(sin24=(`=sin(14*2*_pi/24)')              cos24=(`=cos(14*2*_pi/24)'))            ///
                at(sin24=(`=sin(15*2*_pi/24)')              cos24=(`=cos(15*2*_pi/24)'))            ///
                at(sin24=(`=sin(16*2*_pi/24)')              cos24=(`=cos(16*2*_pi/24)'))            ///
                at(sin24=(`=sin(17*2*_pi/24)')              cos24=(`=cos(17*2*_pi/24)'))            ///
                at(sin24=(`=sin(18*2*_pi/24)')              cos24=(`=cos(18*2*_pi/24)'))            ///
                at(sin24=(`=sin(19*2*_pi/24)')              cos24=(`=cos(19*2*_pi/24)'))            /// 
                at(sin24=(`=sin(20*2*_pi/24)')              cos24=(`=cos(20*2*_pi/24)'))            ///
                at(sin24=(`=sin(21*2*_pi/24)')              cos24=(`=cos(21*2*_pi/24)'))            ///
                at(sin24=(`=sin(22*2*_pi/24)')              cos24=(`=cos(22*2*_pi/24)'))            ///
                at(sin24=(`=sin(23*2*_pi/24)')              cos24=(`=cos(23*2*_pi/24)'))            ///
                at(sin24=(`=sin(24*2*_pi/24)')              cos24=(`=cos(24*2*_pi/24)'))            ///
                at(sin24=(`=sin(`curAcro'*2*_pi/24)')       cos24=(`=cos(`curAcro'*2*_pi/24)'))     ///
                at(sin24=(`=sin(`curAcro'*2*_pi/24-_pi)')   cos24=(`=cos(`curAcro'*2*_pi/24-_pi)')) ///
                post    

        ************************
        ** Estimate max value **
        ************************

        nlcom  _b[25._at]
        replace max_hat         = r(b)[1,1]                                         
        replace max_se          = sqrt(r(V)[1,1])                                   
        replace max_lb          = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
        replace max_ub          = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
        replace max_p           = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         

        ************************
        ** Estimate min value **
        ************************

        nlcom  _b[26._at]
        replace min_hat         = r(b)[1,1]                                         
        replace min_se          = sqrt(r(V)[1,1])                                   
        replace min_lb          = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
        replace min_ub          = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
        replace min_p           = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         

        ************************
        ** Estimate amplitude **
        ************************

        nlcom (_b[25._at]-_b[26._at])/2
        replace amplitude_hat   = r(b)[1,1]                                         
        replace amplitude_se    = sqrt(r(V)[1,1])                                   
        replace amplitude_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
        replace amplitude_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
        replace amplitude_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         

        *************************
        ** Estimate total paee **
        *************************

        nlcom   _b[1._at]+_b[2._at]+_b[3._at]       +   ///
                _b[4._at]+_b[5._at]+_b[6._at]       +   ///
                _b[7._at]+_b[8._at]+_b[9._at]       +   ///
                _b[10._at]+_b[11._at]+_b[12._at]    +   ///
                _b[13._at]+_b[14._at]+_b[15._at]    +   ///
                _b[16._at]+_b[17._at]+_b[18._at]    +   ///
                _b[19._at]+_b[20._at]+_b[21._at]    +   ///
                _b[22._at]+_b[23._at]+_b[24._at]        

        replace paeet_hat       = r(b)[1,1]                                         
        replace paeet_se        = sqrt(r(V)[1,1])                                   
        replace paeet_lb        = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
        replace paeet_ub        = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
        replace paeet_p         = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))         

        ***********************
        ** Estimate q factor **
        ***********************

        nlcom   2*_pi*(_b[25._at]-_b[26._at])           ///
                /                                       ///
                (                                       ///
                _b[1._at]+_b[2._at]+_b[3._at]       +   ///
                _b[4._at]+_b[5._at]+_b[6._at]       +   ///
                _b[7._at]+_b[8._at]+_b[9._at]       +   ///
                _b[10._at]+_b[11._at]+_b[12._at]    +   ///
                _b[13._at]+_b[14._at]+_b[15._at]    +   ///
                _b[16._at]+_b[17._at]+_b[18._at]    +   ///
                _b[19._at]+_b[20._at]+_b[21._at]    +   ///
                _b[22._at]+_b[23._at]+_b[24._at]        ///
                )

        replace qfactor_hat   = r(b)[1,1]                                           
        replace qfactor_se   = sqrt(r(V)[1,1])                                      
        replace qfactor_lb   = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])         
        replace qfactor_ub   = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])         
        replace qfactor_p   = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))             

        drop _est_m1
        
        restore

        asdf

    }

    else drop if ID == "`curID'"

}

drop hour paee_hour sin24 cos24
duplicates drop

save "C:\Users\tg421\OneDrive - University of Cambridge\Fenland diurnal PA and Met risk\PAEE_estimates.dta", replace


asdf











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