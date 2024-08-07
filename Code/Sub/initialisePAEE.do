version 17.0

********************************************************************************
** Redo Cosinor Modelling
********************************************************************************

capture confirm file "C:\Users\tg421\OneDrive - University of Cambridge\Fenland diurnal PA and Met risk\cosinorEstimates.dta"

if _rc != 0{

    *************************
    ** Initialise postfile **
    *************************

    capture postutil clear
    #delimit ;
    postfile paeepost  
            str9 ID
            double(sin24_hat     sin24_se     sin24_lb     sin24_ub     sin24_p)
            double(cos24_hat     cos24_se     cos24_lb     cos24_ub     cos24_p)
            double(acrophase_hat acrophase_se acrophase_lb acrophase_ub acrophase_p)   
            double(max_hat       max_se       max_lb       max_ub       max_p)                     
            double(min_hat       min_se       min_lb       min_ub       min_p) 
            double(amplitude_hat amplitude_se amplitude_lb amplitude_ub amplitude_p) 
            double(paeet_hat     paeet_se     paeet_lb     paeet_ub     paeet_p)  
            double(qfactor_hat   qfactor_se   qfactor_lb   qfactor_ub   qfactor_p)    
            using
            "C:\Users\tg421\OneDrive - University of Cambridge\Fenland diurnal PA and Met risk\cosinorEstimates.dta"
            ,
            replace
            ;
    #delimit cr

    ******************************************************
    ** Initialise dataset + outcome modelling variables **
    ******************************************************
    
    keep ID paee_hour*
    reshape long paee_hour, i(ID) j(hour)
    drop if paee_hour == .
    qui levelsof ID, local(IDList)

    gen sin24 = sin(hour*2*_pi/24)
    gen cos24 = cos(hour*2*_pi/24)

    ***************************************
    ** Perform IQR outlier analysis loop **
    ***************************************

    gen upperOutlier    = 0
    gen lowerOutlier    = 0

    qui foreach curID of local IDList{

        noisi di "IQR outlier detection: `curID'"
        
        su paee_hour if ID == "`curID'", detail
        replace upperOutlier = cond((paee_hour - r(p75)) / (r(p75) - r(p25)) >=  2 ,1,0) if ID == "`curID'"
        replace lowerOutlier = cond((paee_hour - r(p25)) / (r(p75) - r(p25)) <= -2 ,1,0) if ID == "`curID'"

    }

    qui drop if upperOutlier == 1
    qui drop if lowerOutlier == 1
    qui drop upperOutlier lowerOutlier
    
    ******************************
    ** Begin main analysis loop **
    ******************************

    qui foreach curID of local IDList{

        count if ID == "`curID'"

        if r(N) > 18{

            noisi di "Cosinor modelling: `curID'"

            ***********************************************************
            ** Restrict sample to current ID and initialize postlist **
            ***********************************************************

            preserve
            keep if ID == "`curID'"
            local postlist ("`curID'")
            
            ************************************************************************************************
            ** Apply cosinor GLM model, using log-link function to ensure paee values are always positive **
            ************************************************************************************************

            glm paee_hour sin24 cos24, family(gaussian) link(log)
            estimates store m1

            *******************************************************
            ** Estimate real-values for sin and cos coefficients **
            *******************************************************

            margin, dydx(sin24 cos24) post

            local sin24_hat       = _b[sin24]                                         
            local sin24_se        = _se[sin24]                                        
            local sin24_lb        = _b[sin24] + invnormal(0.025)*_se[sin24]           
            local sin24_ub        = _b[sin24] + invnormal(0.975)*_se[sin24]           
            local sin24_p         = 2*normal(-abs(_b[sin24]/_se[sin24]))

            local postlist `postlist' (`sin24_hat') (`sin24_se') (`sin24_lb') (`sin24_ub') (`sin24_p')              

            local cos24_hat       = _b[cos24]                                         
            local cos24_se        = _se[cos24]                                        
            local cos24_lb        = _b[cos24] + invnormal(0.025)*_se[cos24]           
            local cos24_ub        = _b[cos24] + invnormal(0.975)*_se[cos24]           
            local cos24_p         = 2*normal(-abs(_b[cos24]/_se[cos24]))              

            local postlist `postlist' (`cos24_hat') (`cos24_se') (`cos24_lb') (`cos24_ub') (`cos24_p') 

            ************************
            ** Estimate acrophase **
            ************************

            nlcom cond(_b[sin24]<0,24,0)+atan2(_b[sin24],_b[cos24])*24/(2*_pi) 
            
            local acrophase_hat   = r(b)[1,1]                                         
            local acrophase_se    = sqrt(r(V)[1,1])                                   
            local acrophase_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
            local acrophase_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
            local acrophase_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

            local postlist `postlist' (`acrophase_hat') (`acrophase_se') (`acrophase_lb') (`acrophase_ub') (`acrophase_p')         

            *************************************************************************************************************
            ** Get predicted values at max and min hour (at 25 and 26), as well as hours 1 through 24 (at 1 through 24) **
            *************************************************************************************************************

            estimates restore m1

            #delimit ;
            margin, at(sin24=(`=sin(1*2*_pi/24)')   cos24=(`=cos(1*2*_pi/24)'))           
                    at(sin24=(`=sin(2*2*_pi/24)')   cos24=(`=cos(2*2*_pi/24)'))           
                    at(sin24=(`=sin(3*2*_pi/24)')   cos24=(`=cos(3*2*_pi/24)'))         
                    at(sin24=(`=sin(4*2*_pi/24)')   cos24=(`=cos(4*2*_pi/24)'))

                    at(sin24=(`=sin(5*2*_pi/24)')   cos24=(`=cos(5*2*_pi/24)'))            
                    at(sin24=(`=sin(6*2*_pi/24)')   cos24=(`=cos(6*2*_pi/24)'))            
                    at(sin24=(`=sin(7*2*_pi/24)')   cos24=(`=cos(7*2*_pi/24)'))            
                    at(sin24=(`=sin(8*2*_pi/24)')   cos24=(`=cos(8*2*_pi/24)'))

                    at(sin24=(`=sin(9*2*_pi/24)')   cos24=(`=cos(9*2*_pi/24)'))             
                    at(sin24=(`=sin(10*2*_pi/24)')  cos24=(`=cos(10*2*_pi/24)'))          
                    at(sin24=(`=sin(11*2*_pi/24)')  cos24=(`=cos(11*2*_pi/24)'))            
                    at(sin24=(`=sin(12*2*_pi/24)')  cos24=(`=cos(12*2*_pi/24)'))

                    at(sin24=(`=sin(13*2*_pi/24)')  cos24=(`=cos(13*2*_pi/24)'))          
                    at(sin24=(`=sin(14*2*_pi/24)')  cos24=(`=cos(14*2*_pi/24)'))            
                    at(sin24=(`=sin(15*2*_pi/24)')  cos24=(`=cos(15*2*_pi/24)'))            
                    at(sin24=(`=sin(16*2*_pi/24)')  cos24=(`=cos(16*2*_pi/24)'))

                    at(sin24=(`=sin(17*2*_pi/24)')  cos24=(`=cos(17*2*_pi/24)'))            
                    at(sin24=(`=sin(18*2*_pi/24)')  cos24=(`=cos(18*2*_pi/24)'))           
                    at(sin24=(`=sin(19*2*_pi/24)')  cos24=(`=cos(19*2*_pi/24)'))           
                    at(sin24=(`=sin(20*2*_pi/24)')  cos24=(`=cos(20*2*_pi/24)'))
                                
                    at(sin24=(`=sin(21*2*_pi/24)')  cos24=(`=cos(21*2*_pi/24)'))           
                    at(sin24=(`=sin(22*2*_pi/24)')  cos24=(`=cos(22*2*_pi/24)'))           
                    at(sin24=(`=sin(23*2*_pi/24)')  cos24=(`=cos(23*2*_pi/24)'))            
                    at(sin24=(`=sin(24*2*_pi/24)')  cos24=(`=cos(24*2*_pi/24)'))  

                    at(sin24=(`=sin(`acrophase_hat'*2*_pi/24)')     cos24=(`=cos(`acrophase_hat'*2*_pi/24)'))    
                    at(sin24=(`=sin(`acrophase_hat'*2*_pi/24-_pi)') cos24=(`=cos(`acrophase_hat'*2*_pi/24-_pi)')) 
                    post
                    ;
            #delimit cr

            ************************
            ** Estimate max value **
            ************************

            nlcom  _b[25._at]

            local max_hat         = r(b)[1,1]                                         
            local max_se          = sqrt(r(V)[1,1])                                   
            local max_lb          = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
            local max_ub          = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
            local max_p           = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1]))) 

            local postlist `postlist' (`max_hat') (`max_se') (`max_lb') (`max_ub') (`max_p')        

            ************************
            ** Estimate min value **
            ************************

            nlcom  _b[26._at]

            local min_hat         = r(b)[1,1]                                         
            local min_se          = sqrt(r(V)[1,1])                                   
            local min_lb          = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
            local min_ub          = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
            local min_p           = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1]))) 

            local postlist `postlist' (`min_hat') (`min_se') (`min_lb') (`min_ub') (`min_p')        

            ************************
            ** Estimate amplitude **
            ************************

            nlcom (_b[25._at]-_b[26._at])/2

            local amplitude_hat   = r(b)[1,1]                                         
            local amplitude_se    = sqrt(r(V)[1,1])                                   
            local amplitude_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
            local amplitude_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
            local amplitude_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

            local postlist `postlist' (`amplitude_hat') (`amplitude_se') (`amplitude_lb') (`amplitude_ub') (`amplitude_p')       

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

            local paeet_hat       = r(b)[1,1]                                         
            local paeet_se        = sqrt(r(V)[1,1])                                   
            local paeet_lb        = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
            local paeet_ub        = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
            local paeet_p         = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

            local postlist `postlist' (`paeet_hat') (`paeet_se') (`paeet_lb') (`paeet_ub') (`paeet_p')        

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

            local qfactor_hat   = r(b)[1,1]                                           
            local qfactor_se    = sqrt(r(V)[1,1])                                      
            local qfactor_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])         
            local qfactor_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])         
            local qfactor_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

            local postlist `postlist' (`qfactor_hat') (`qfactor_se') (`qfactor_lb') (`qfactor_ub') (`qfactor_p')

            ******************
            ** Post results **
            ******************

            post paeepost `postlist'

            estimates clear
            restore
        }

        else noisi di "Insufficient data: `curID'"
    }

    postclose paeepost

}

else joinby ID using "C:\Users\tg421\OneDrive - University of Cambridge\Fenland diurnal PA and Met risk\cosinorEstimates.dta"
drop _merge




/*

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

*/

