version 17.0

***********************
** Cosinor Modelling **
***********************

**  Performs cosinor analysis to model rhythmic patterns in physical 
**  activity energy expenditure (PAEE) data. It fits a gamma GLM with a log link function 
**  to estimate the parameters of the cosinor model, including the mesor, amplitudes, and 
**  acrophases for 24-hour, 12-hour, and 8-hour cycles. The script also estimates the hour 
**  of the day where the maximum PAEE occurs and the total PAEE over 24 hours.
**  Performs IQR-based outlier removal before fitting the cosinor model. 
**  Uses a multi-start optimization approach to estimate the hour of maximum PAEE.
**  Total PAEE is estimated as the average of marginal values computed at each hour within 24-hours.


****************************************************************
** Define cosinor model specification for max hour estimation **
****************************************************************

cap mata: mata drop cosinorModel()
mata
    void cosinorModel(real scalar todo      , 
                        real vector b         , 
                        real vector sinCoef   , 
                        real vector cosCoef   , 
                        real scalar mesorCoef , 
                        val                   , 
                        grad                  , 
                        hess)
    {   

        // Designation for model.  Note the "exp", which corresponds to "log link function" in glm used

        val =   exp(
                sinCoef[1]:*sin(b:*2:*pi():/24) :+ cosCoef[1]:*cos(b:*2:*pi():/24)  :+ 
                sinCoef[2]:*sin(b:*2:*pi():/12) :+ cosCoef[2]:*cos(b:*2:*pi():/12)  :+
                sinCoef[3]:*sin(b:*2:*pi():/8)  :+ cosCoef[3]:*cos(b:*2:*pi():/8)   :+
                mesorCoef
                )
    }
end

****************************************************************************
** Define function to find hour of global maximum of fitted cosinor model **
****************************************************************************

cap mata: mata drop estimateMaxHour()
mata
    void estimateMaxHour(real vector sinCoef    , 
                            real vector cosCoef     , 
                            real scalar mesorCoef)
    {

        transmorphic S
        vector bh, maxInd, maxCount, maxHour_hat, maxHour_se, maxValue_hat
        scalar i

        maxHour_hat  = J(1,0,.)
        maxHour_se   = J(1,0,.) 
        maxValue_hat = J(1,0,.)

        // Run optimisation loop, using different initial guess for hour (i) each time

        for (i=0; i<=23; i++){

            S  = optimize_init()
            optimize_init_argument(S, 1, sinCoef)
            optimize_init_argument(S, 2, cosCoef)
            optimize_init_argument(S, 3, mesorCoef)
            optimize_init_evaluator(S, &cosinorModel())
            optimize_init_params(S, J(1,1,i))
            bh = optimize(S)

            // Only keep results that are within 0 to 24 hours
            
            if (optimize_result_params(S):>=0 :& optimize_result_params(S):<=24){
                maxHour_hat = maxHour_hat, optimize_result_params(S)
                maxHour_se  = maxHour_se, sqrt(diagonal(optimize_result_V(S)))'
                maxValue_hat = maxValue_hat, optimize_result_value(S)
            }
        }

        // Output max hour hat and se

        maxindex(maxValue_hat,2,maxInd,maxCount)
        st_local("maxHour_hat" , strofreal(maxHour_hat[maxInd[1]]))
        st_local("maxHour_se"  , strofreal(maxHour_se[maxInd[1]]))

    }
end





capture confirm file "`rootPath'/cosinorEstimates.dta"

if _rc != 0{

    frame copy dataset tempset
    frame change tempset

    *************************
    ** Initialise postfile **
    *************************

    local postList mesor_hat mesor_se mesor_lb mesor_ub mesor_p

    // Start with cosinor submodel vars

    foreach curVar in sin cos acrophase amplitude{
        foreach p in 24 12 8{
            local postList `postList' `curVar'`p'_hat `curVar'`p'_se `curVar'`p'_lb `curVar'`p'_ub `curVar'`p'_p
        }
    }

    // Next with cosinor whole-model vars

    foreach curVar in maxHour maxValue totalPAEE{
        local postList `postList' `curVar'_hat `curVar'_se `curVar'_lb `curVar'_ub `curVar'_p
    }

    // Set postfile

    capture postutil clear
    postfile cosinorpost str9 ID double(`postList') using "`rootPath'/cosinorEstimates.dta" , replace

    qui{

        ******************************************************
        ** Initialise dataset + outcome modelling variables **
        ******************************************************
        
        keep ID paee_hour* pwear_hour*
        reshape long paee_hour pwear_hour, i(ID) j(hour)
        drop if paee_hour == .

        // Generate sin and cos axis variables for 24, 12, and 8 hour periods

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

        // Collapse data set on percentiles, then exclude values that are + or - 2 IQR

        collapse (p25) p25=paee_hour (p75) p75=paee_hour (iqr) iqr=paee_hour, by(ID)

        frame change tempset
        frlink m:1 ID, frame(tempsub)
        frget *, from(tempsub)
        frame drop tempsub

        gen upperOutlier = cond((paee_hour - p75) / (iqr) >=  2 ,1,0)
        gen lowerOutlier = cond((paee_hour - p25) / (iqr) <= -2 ,1,0)
        
        drop if upperOutlier == 1 | lowerOutlier == 1
        drop upperOutlier lowerOutlier tempsub p25 p75 iqr

        // Drop those with more than 4 data points excluded by IQR outlier analysis

        frame copy tempset tempsub
        frame change tempsub

        collapse (count) notout=paee_hour, by(ID)

        frame change tempset
        frlink m:1 ID, frame(tempsub)
        frget *, from(tempsub)
        frame drop tempsub

        drop if notout <20
        drop notout tempsub

    }

    ******************************
    ** Begin main analysis loop **
    ******************************
    
    qui levelsof ID, local(IDList)
    qui foreach curID of local IDList{

        noisi di "Cosinor modelling: `curID'"

        ***********************************
        ** Restrict sample to current ID **
        ***********************************

        frame put * if ID == "`curID'", into(tempsub)
        frame change tempsub

        *****************************************************************
        ** Apply cosinor GLM model, using gamma with log-link function **
        *****************************************************************

        glm paee_hour sin24 cos24 sin12 cos12 sin8 cos8 [aw=pwear_hour], family(gamma) link(log)
        
        // If model converges, continue.

        if e(converged)==1{

            estimates store m1
            
            local postlist ("`curID'") 

            **********************************
            ** Get sin and cos coefficients **
            **********************************

            // Get "midline-estimating statistic of rhythm" (i.e. mesor) value, which is just the constant of the model

            local mesor_hat = _b[_cons]                                         
            local mesor_se  = _se[_cons]                                        
            local mesor_lb  = _b[_cons] + invnormal(0.025)*_se[_cons]           
            local mesor_ub  = _b[_cons] + invnormal(0.975)*_se[_cons]           
            local mesor_p   = 2*normal(-abs(_b[_cons]/_se[_cons]))

            local postlist `postlist' (`mesor_hat') (`mesor_se') (`mesor_lb') (`mesor_ub') (`mesor_p')

            foreach curAxis in sin cos{

                foreach p in 24 12 8{
                    
                    local `curAxis'`p'_hat  = .                                       
                    local `curAxis'`p'_se   = .                                    
                    local `curAxis'`p'_lb   = .          
                    local `curAxis'`p'_ub   = .          
                    local `curAxis'`p'_p    = .

                    capture di _b[`curAxis'`p']

                    // Get coefficients for each cosinor submodel

                    if _rc == 0{

                        local `curAxis'`p'_hat  = _b[`curAxis'`p']                                         
                        local `curAxis'`p'_se   = _se[`curAxis'`p']                                        
                        local `curAxis'`p'_lb   = _b[`curAxis'`p'] + invnormal(0.025)*_se[`curAxis'`p']           
                        local `curAxis'`p'_ub   = _b[`curAxis'`p'] + invnormal(0.975)*_se[`curAxis'`p']           
                        local `curAxis'`p'_p    = 2*normal(-abs(_b[`curAxis'`p']/_se[`curAxis'`p']))
                    }

                    local postlist `postlist' (``curAxis'`p'_hat') (``curAxis'`p'_se') (``curAxis'`p'_lb') (``curAxis'`p'_ub') (``curAxis'`p'_p')
                }
            }

            *************************
            ** Estimate acrophases **
            *************************

            foreach p in 24 12 8{
                
                local acrophase`p'_hat   = .                                       
                local acrophase`p'_se    = .                            
                local acrophase`p'_lb    = .     
                local acrophase`p'_ub    = . 
                local acrophase`p'_p     = .

                // For acrophase computation, reference quadrant is changed depending on + or - for sin beta

                capture nlcom cond(_b[sin`p']<0,`p',0) + atan2(_b[sin`p'],_b[cos`p'])*`p'/(2*_pi) 

                if _rc == 0{

                    local acrophase`p'_hat   = r(b)[1,1]                                         
                    local acrophase`p'_se    = sqrt(r(V)[1,1])                                   
                    local acrophase`p'_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
                    local acrophase`p'_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
                    local acrophase`p'_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

                }

                else local acrophase`p'_hat = cond(_b[sin`p']<0,`p',0) + atan2(_b[sin`p'],_b[cos`p'])*`p'/(2*_pi) 

                local postlist `postlist' (`acrophase`p'_hat') (`acrophase`p'_se') (`acrophase`p'_lb') (`acrophase`p'_ub') (`acrophase`p'_p') 
                
            }

            *************************
            ** Estimate amplitudes **
            *************************

            foreach p in 24 12 8{
                
                local amplitude`p'_hat   = .                                       
                local amplitude`p'_se    = .
                local amplitude`p'_lb    = .     
                local amplitude`p'_ub    = .     
                local amplitude`p'_p     = .  

                capture nlcom sqrt(_b[sin`p']^2+_b[cos`p']^2)

                if _rc == 0{

                    local amplitude`p'_hat   = r(b)[1,1]                                       
                    local amplitude`p'_se    = sqrt(r(V)[1,1])                                   
                    local amplitude`p'_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
                    local amplitude`p'_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
                    local amplitude`p'_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

                }

                else local amplitude`p'_hat = sqrt(_b[sin`p']^2+_b[cos`p']^2)
                    
                local postlist `postlist' (`amplitude`p'_hat') (`amplitude`p'_se') (`amplitude`p'_lb') (`amplitude`p'_ub') (`amplitude`p'_p')          
            }
            
            ***************************************************
            ** Estimate max attained value, the hour it occurrs
            ***************************************************

            local maxHour_hat   = .                                      
            local maxHour_se    = .                                   
            local maxHour_lb    = .    
            local maxHour_ub    = .     
            local maxHour_p     = .

            local maxValue_hat   = .                                     
            local maxValue_se    = .                                  
            local maxValue_lb    = .     
            local maxValue_ub    = .   
            local maxValue_p     = .

            // Populate matrices with model coeffients to send into mata

            mat sinCoef     = `sin24_hat' , `sin12_hat' , `sin8_hat' 
            mat cosCoef     = `cos24_hat' , `cos12_hat' , `cos8_hat'
            mat mesorCoef   = `mesor_hat'

            // Run mata code to do multi-start optimization of hour where max value occurs in cosinor model

            qui mata: estimateMaxHour(st_matrix("sinCoef"), st_matrix("cosCoef"), st_matrix("mesorCoef"))
            
            if `maxHour_hat' != .{

                local maxHour_hat   = `maxHour_hat'                                        
                local maxHour_se    = `maxHour_se'                                    
                local maxHour_lb    = `maxHour_hat' + invnormal(0.025)*sqrt(`maxHour_se' )      
                local maxHour_ub    = `maxHour_hat' + invnormal(0.975)*sqrt(`maxHour_se' )      
                local maxHour_p     = 2*normal(-abs(`maxHour_hat'/sqrt(`maxHour_se')))

                // Initialize "at statement" for marginal hour value where global max occurs

                foreach p in 24 12 8{
                    
                    local per`p' sin`p'=(`=sin(`maxHour_hat'*2*_pi/`p')')cos`p'=(`=cos(`maxHour_hat'*2*_pi/`p')')

                }

                // Compute marginal value at hour where global max occurs

                capture margin, at(`per24'`per12'`per8') post

                if _rc == 0{

                    local maxValue_hat   = r(b)[1,1]                                         
                    local maxValue_se    = sqrt(r(V)[1,1])                                   
                    local maxValue_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
                    local maxValue_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
                    local maxValue_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

                }
  
            }

            local postlist `postlist' (`maxHour_hat')  (`maxHour_se')  (`maxHour_lb')  (`maxHour_ub')  (`maxHour_p')  
            local postlist `postlist' (`maxValue_hat') (`maxValue_se') (`maxValue_lb') (`maxValue_ub') (`maxValue_p')  

            *************************
            ** Estimate total PAEE **
            *************************

            local totalPAEE_hat   = .                                      
            local totalPAEE_se    = .                                   
            local totalPAEE_lb    = .    
            local totalPAEE_ub    = .     
            local totalPAEE_p     = .

            // Initilize at margin list as each hour of the day, across each cosinor submodel

            local marginList
            forvalues h = 0(1)23{
                foreach p in 24 12 8{
                    local per`p' sin`p'=(`=sin(`h'*2*_pi/`p')')cos`p'=(`=cos(`h'*2*_pi/`p')')
                }
                local marginList `marginList' at(`per24'`per12'`per8')
            }

            // Get marginal values

            estimates restore m1
            capture margin, `marginList' post

            if _rc == 0{
                
                // Compute sum of marginal values at each hour to compute total PAEE

                capture nlcom   _b[1._at]  + _b[2._at]  + _b[3._at]  + _b[4._at]  +     ///
                                _b[5._at]  + _b[6._at]  + _b[7._at]  + _b[8._at]  +     ///
                                _b[9._at]  + _b[10._at] + _b[11._at] + _b[12._at] +     ///
                                _b[13._at] + _b[14._at] + _b[15._at] + _b[16._at] +     ///
                                _b[17._at] + _b[18._at] + _b[19._at] + _b[20._at] +     ///
                                _b[21._at] + _b[22._at] + _b[23._at] + _b[24._at]
                
                if _rc == 0{

                    local totalPAEE_hat   = r(b)[1,1]                                       
                    local totalPAEE_se    = sqrt(r(V)[1,1])                                   
                    local totalPAEE_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
                    local totalPAEE_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
                    local totalPAEE_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))
                }

                else{
                    
                    predict x
                    su x, detail
                    local totalPAEE_hat = r(sum)

                }

            }

            local postlist `postlist' (`totalPAEE_hat') (`totalPAEE_se') (`totalPAEE_lb') (`totalPAEE_ub') (`totalPAEE_p')  
                                       
            *********************************
            ** Output results to post file **
            *********************************

            post cosinorpost `postlist'
            
        }

        frame change tempset
        frame drop tempsub

    }
    
    postclose cosinorpost
    
    frame change dataset
    frame drop tempset
}

joinby ID using "`rootPath'/cosinorEstimates.dta"
