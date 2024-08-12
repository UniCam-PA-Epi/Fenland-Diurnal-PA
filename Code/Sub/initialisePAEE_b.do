version 17.0

********************************************************************************
** Redo Cosinor Modelling
********************************************************************************

args rootPath cosinorEstimatesFile

//capture erase "C:\Users\tg421\OneDrive - University of Cambridge\Fenland diurnal PA and Met risk\cosinorEstimates.dta"
capture confirm file "`rootPath'/`cosinorEstimatesFile'"

if _rc == 0{

    frame copy dataset tempset
    frame change tempset

    *************************
    ** Initialise postfile **
    *************************

    capture postutil clear
    #delimit ;
    postfile cosinorpost  
            str9 ID
            double(sin24_hat     sin24_se     sin24_lb     sin24_ub     sin24_p     )
            double(cos24_hat     cos24_se     cos24_lb     cos24_ub     cos24_p     )
            double(acrophase_hat acrophase_se acrophase_lb acrophase_ub acrophase_p )   
            double(max_hat       max_se       max_lb       max_ub       max_p       )                     
            double(min_hat       min_se       min_lb       min_ub       min_p       ) 
            double(amplitude_hat amplitude_se amplitude_lb amplitude_ub amplitude_p ) 
            double(paeet_hat     paeet_se     paeet_lb     paeet_ub     paeet_p     )  
            double(qfactor_hat   qfactor_se   qfactor_lb   qfactor_ub   qfactor_p   )    
            using
            "`rootPath'/temp.dta" , replace
            //"`rootPath'/`cosinorEstimatesFile'"
            ;
    #delimit cr

    ******************************************************
    ** Initialise dataset + outcome modelling variables **
    ******************************************************
    
    keep ID paee_hour*
    reshape long paee_hour, i(ID) j(hour)
    drop if paee_hour == .

    gen sin24 = sin(hour*2*_pi/24)
    gen cos24 = cos(hour*2*_pi/24)

    gen sin12 = sin(hour*2*_pi/12)
    gen cos12 = cos(hour*2*_pi/12)

    gen sin8 = sin(hour*2*_pi/8)
    gen cos8 = cos(hour*2*_pi/8)

    ***************************************
    ** Perform IQR outlier analysis loop **
    ***************************************

    frame copy tempset tempsub
    frame change tempsub

    collapse (p25) p25=paee_hour (p75) p75=paee_hour (iqr) iqr=paee_hour, by(ID)

    frame change tempset
    frlink m:1 ID, frame(tempsub)
    frget *, from(tempsub)
    frame drop tempsub

    gen upperOutlier = cond((paee_hour - p75) / (iqr) >=  2 ,1,0)
    gen lowerOutlier = cond((paee_hour - p25) / (iqr) <= -2 ,1,0)
    
    drop if upperOutlier == 1 | lowerOutlier == 1
    drop upperOutlier lowerOutlier tempsub p25 p75 iqr

    frame copy tempset tempsub
    frame change tempsub

    collapse (count) notout=paee_hour, by(ID)

    frame change tempset
    frlink m:1 ID, frame(tempsub)
    frget *, from(tempsub)
    frame drop tempsub

    drop if notout <20
    drop notout tempsub

    ******************************
    ** Begin main analysis loop **
    ******************************
    
    qui levelsof ID, local(IDList)
    qui foreach curID of local IDList{
        local curID =  "R16600025"
        noisi di "Cosinor modelling: `curID'"

        ***********************************************************
        ** Restrict sample to current ID and initialize postlist **
        ***********************************************************

        frame put * if ID == "`curID'", into(tempsub)
        frame change tempsub

        local postlist ("`curID'")
        
        ************************************************************************************************
        ** Apply cosinor GLM model, using log-link function to ensure paee values are always positive **
        ************************************************************************************************

        noisi glm paee_hour sin24 cos24 sin12 cos12 sin8 cos8, family(gamma) link(power 0.5)
        estimates store m1

        asdf

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
                _b[22._at]+_b[23._at]+_b[24._at]    -   ///
                (24*_b[26._at])                         ///
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

        post cosinorpost `postlist'

        frame change tempset
        frame drop tempsub

    }

    postclose cosinorpost
    
    frame change dataset
    frame drop tempset
}

joinby ID using "`rootPath'/`cosinorEstimatesFile'"