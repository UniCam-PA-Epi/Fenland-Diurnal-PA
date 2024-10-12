version 17.0

graph drop _all

frame copy dataset tempset
frame change tempset

*************************
** Initilaise analysis **
*************************


cap mata: mata drop shiftModel()
mata
    void shiftModel(real scalar todo            , 
                    real vector b               , 
                    real scalar amp24initial    , 
                    real scalar amp12initial    , 
                    real scalar amp8initial     , 
                    real scalar acr8initial     , 
                    real scalar mesorInitial    ,
                    real scalar acr24target     ,
                    real scalar acr12target     ,
                    real scalar initialPAEE     ,
                    val                         , 
                    grad                        , 
                    hess)
    {   
        
        real scalar targetPAEE, i

        targetPAEE = 0
        for (i=1; i<=128; i++){
            targetPAEE =    targetPAEE +
                            exp(
                            amp24initial:*cos((2:*pi():/24):*((24*i/128)-acr24target))   :+
                            amp12initial:*cos((2:*pi():/12):*((24*i/128)-acr12target))   :+
                            amp8initial:*cos((2:*pi():/8):*((24*i/128)-acr8initial))     :+
                            (mesorInitial:+b[1])
                            )
        }

        val = (targetPAEE:-initialPAEE):^2  

    }
end


cap mata: mata drop estimateShift()
mata
    void estimateShift(real scalar amp24initial     , 
                       real scalar amp12initial     , 
                       real scalar amp8initial      , 
                       real scalar acr24initial     , 
                       real scalar acr12initial     , 
                       real scalar acr8initial      , 
                       real scalar mesorInitial     ,
                       real scalar acr24target      ,
                       real scalar acr12target      
                       )
    {

        
        transmorphic S
        vector bh
        real scalar initialPAEE, i

        initialPAEE = 0
        for (i=1; i<=128; i++){
            initialPAEE =   initialPAEE +
                            exp(
                            amp24initial:*cos((2:*pi():/24):*((24*i/128)-acr24initial))  :+
                            amp12initial:*cos((2:*pi():/12):*((24*i/128)-acr12initial))  :+
                            amp8initial:*cos((2:*pi():/8):*((24*i/128)-acr8initial))     :+
                            mesorInitial
                            )
        }
        
        S  = optimize_init()
        optimize_init_argument(S, 1, amp24initial) 
        optimize_init_argument(S, 2, amp12initial)
        optimize_init_argument(S, 3, amp8initial)
        optimize_init_argument(S, 4, acr8initial)
        optimize_init_argument(S, 5, mesorInitial)
        optimize_init_argument(S, 6, acr24target)
        optimize_init_argument(S, 7, acr12target)
        optimize_init_argument(S, 8, initialPAEE)
        optimize_init_evaluator(S, &shiftModel())
        optimize_init_params(S, J(1,1,0))
        optimize_init_which(S,"min")
        bh = optimize(S)

        st_local("mesorOffset" , strofreal(optimize_result_params(S)[1]))
    

    }
end








// Exclude participants for whom the cosinor model did not fit, those with insufficient wear time, and those without fat mass measurement.

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05

drop if mesor_hat == .
drop if amplitude24_hat == . | amplitude12_hat == . | amplitude8_hat == .
drop if acrophase24_hat == . | acrophase12_hat == . | acrophase8_hat  == .

drop if P_Pwear_consolidated<72
drop if fatMass == .
drop if fatFreeMass == .

drop *_se *_lb *_ub *_p
rename *_hat *

// Initialise analysis vars by deconstructing time-based features into their sin and cos axes

gen maxHour24_sin   = sin(maxHour*2*_pi/24)
gen maxHour24_cos   = cos(maxHour*2*_pi/24)

gen acrophase24_sin = sin(acrophase24*2*_pi/24)
gen acrophase24_cos = cos(acrophase24*2*_pi/24)

gen acrophase12_sin = sin(acrophase12*2*_pi/12)
gen acrophase12_cos = cos(acrophase12*2*_pi/12)

gen acrophase8_sin  = sin(acrophase8*2*_pi/8)
gen acrophase8_cos  = cos(acrophase8*2*_pi/8)

// Convert J/min/kg to kJ/hour/kg

replace totalPAEE   = totalPAEE  * 60/1000

****************************************************************************************
** Outcome variables, continuous control variables, and categorical control variables **
****************************************************************************************

#delimit ;
local outcomeVars   //glucose120
                    insulin
                    //leptin
                    //nefa
                    //adiponectin
                    //ldl
                    //hdl
                    //mbpsys
                    //mbpdia
                    //crp
                    ;

local contCovVars   c.age
                    c.age#c.age
                    c.rhr
                    c.rhr#c.rhr
                    c.testTime_sin
                    c.testTime_cos
                    c.testDay_sin
                    c.testDay_cos
                    ;

local catCovVars    i.sex
                    i.diet 
                    i.alcohol
                    i.ethnic
                    i.education
                    i.income
                    i.smoke
                    i.work_status
                    i.marital_status
                    i.cardiometabol_med
                    i.testsite
                    ;

#delimit cr

// Initialise two model levels, with the second level controlling for adiposity  (i.e.fat mass)

local baseModel c.mesor c.mesor#c.mesor c.mesor#i.sex 
local modelLevel1 `contCovVars' `catCovVars'
local modelLevel2 `contCovVars' `catCovVars' c.fatMass c.fatMass#c.fatMass c.fatFreeMass c.fatFreeMass#c.fatFreeMass c.fatMass#c.fatFreeMass


forvalues i = 1/2{

    capture mkdir Results
    putexcel set "Results/6_jointAssociations.xlsx", sheet("lrTest_m`i'") modify
    putexcel A1 = ("outcomeVar")
    putexcel B1 = ("Count: Women")
    putexcel C1 = ("Count: Men")
    putexcel D1 = ("Pvalue: base totalPAEE model")
    putexcel E1 = ("AIC: base totalPAEE model")
    putexcel F1 = ("Pvalue: cosinor model main effects")
    putexcel G1 = ("AIC: cosinor model main effects")
    putexcel H1 = ("Pvalue: cosinor model totalPAEE interactions")
    putexcel I1 = ("AIC: cosinor model totalPAEE interactions")
    putexcel J1 = ("Pvalue: cosinor model sex interactions")
    putexcel K1 = ("AIC: cosinor model sex interactions")
    putexcel L1 = ("Pvalue: full cosinor model vs base totalPAEE model")

}


egen sin24_std = std(sin24)
egen cos24_std = std(cos24)
egen sin12_std = std(sin12)
egen cos12_std = std(cos12)
egen sin8_std  = std(sin8)
egen cos8_std  = std(cos8)
egen mesor_std  = std(mesor)

cluster kmeans sin24_std cos24_std sin12_std cos12_std sin8_std cos8_std mesor_std, k(7) gen(kGroup) measure(L2)

drop *_std

********************************************************************************
** Begin analysis loop, looping through each outcome variable and model level **
********************************************************************************



local curRow = 2
qui foreach curOutcomeVar of local outcomeVars{

    forvalues i = 1/2{
        
        **********************************************************
        ** Apply nested GLM gaussian linear model with log link **
        **********************************************************

        // Note that when "crp" is the outcome var, we apply an inverse-Gaussian GLM due to the extreme right tail.
        // This choice helps resolve convergence issues that may arise when using a Gaussian GLM for "crp".
        // This choice was verified using BIC for all outcome vars.

        #delimit ;
        
        noisi
        glm `curOutcomeVar' 
            `modelLevel`i'' 
            `baseModel' 
            if
            `curOutcomeVar' != .         
            ,
            family(`=cond("`curOutcomeVar'"=="crp","igaussian","gaussian")')
            link(log)
            ;
        estimates store m1  ;

        noisi
        nestreg, lr: glm   `curOutcomeVar' 
                            (
                            `modelLevel`i''
                            )
                            (
                            `baseModel'           
                            )             
                            ( 
                            c.amplitude24
                            c.amplitude12
                            c.amplitude8

                            c.acrophase24_sin c.acrophase24_cos
                            c.acrophase12_sin c.acrophase12_cos
                            c.acrophase8_sin  c.acrophase8_cos
                            )
                            (
                            i.sex#c.amplitude24
                            i.sex#c.amplitude12
                            i.sex#c.amplitude8

                            i.sex#c.acrophase24_sin  i.sex#c.acrophase24_cos
                            i.sex#c.acrophase12_sin  i.sex#c.acrophase12_cos
                            i.sex#c.acrophase8_sin   i.sex#c.acrophase8_cos
                            )
                            (
                            c.mesor#c.amplitude24
                            c.mesor#c.amplitude12
                            c.mesor#c.amplitude8

                            c.acrophase24_sin#c.mesor  c.acrophase24_cos#c.mesor
                            c.acrophase12_sin#c.mesor  c.acrophase12_cos#c.mesor
                            c.acrophase8_sin#c.mesor   c.acrophase8_cos#c.mesor
                            )
                            
                            (
                            c.amplitude24#c.amplitude12
                            c.amplitude24#c.amplitude8
                            c.amplitude12#c.amplitude8
                            )

                            
                            (
                            c.acrophase24_sin#c.acrophase24_cos
                            c.acrophase24_sin#c.acrophase12_cos
                            c.acrophase24_sin#c.acrophase8_cos

                            c.acrophase12_sin#c.acrophase24_cos
                            c.acrophase12_sin#c.acrophase12_cos
                            c.acrophase12_sin#c.acrophase8_cos

                            c.acrophase8_sin#c.acrophase24_cos
                            c.acrophase8_sin#c.acrophase12_cos
                            c.acrophase8_sin#c.acrophase8_cos
                            )
                            
                            (
                            i.sex#c.amplitude24#c.amplitude12
                            i.sex#c.amplitude24#c.amplitude8
                            i.sex#c.amplitude12#c.amplitude8
                            )
                            (
                            i.sex#c.acrophase24_sin#c.acrophase24_cos
                            i.sex#c.acrophase24_sin#c.acrophase12_cos
                            i.sex#c.acrophase24_sin#c.acrophase8_cos

                            i.sex#c.acrophase12_sin#c.acrophase24_cos
                            i.sex#c.acrophase12_sin#c.acrophase12_cos
                            i.sex#c.acrophase12_sin#c.acrophase8_cos

                            i.sex#c.acrophase8_sin#c.acrophase24_cos
                            i.sex#c.acrophase8_sin#c.acrophase12_cos
                            i.sex#c.acrophase8_sin#c.acrophase8_cos
                            )
                            
                            if
                            `curOutcomeVar' != .         
                            ,
                            family(`=cond("`curOutcomeVar'"=="crp","igaussian","gaussian")')
                            link(`=cond("`curOutcomeVar'"=="cmrs","identity","log")')
                            ;
        #delimit cr
            
        capture mkdir Models
        estimates save Models/`curOutcomeVar'_m`i' , replace
        estimates store m2

        mat lrMat = r(lr)

        // Store results of nested lr tests for the contribution of each block to the model
        
        putexcel set "Results/6_jointAssociations.xlsx", sheet("lrTest_m`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")

        count if e(sample) == 1 & sex == 0
        putexcel B`curRow' = (r(N))

        count if e(sample) == 1 & sex == 1
        putexcel C`curRow' = (r(N))

        forvalues j = 0/4{
            local curCol = char(68+`=`j'*2')  
            putexcel `curCol'`curRow' = (lrMat[`=`j'+2',4])
            local curCol = char(69+`=`j'*2')
            putexcel `curCol'`curRow' = (lrMat[`=`j'+2',5])
        }

        noisi lrtest m1
        putexcel L`curRow' = (r(p))
        putexcel clear
        
        if r(p) < 0.001 local curP = "p<0.001"
        else            local curP = "p=`=trim("`: display %10.3f r(p)'")'"

        local pooledDeltaAIC = "{&Delta}AIC=`=trim("`: display %10.0f lrMat[4,5]-lrMat[2,5]'")'"
        local stratDeltaAIC  = "{&Delta}AIC=`=trim("`: display %10.0f lrMat[6,5]-lrMat[2,5]'")'"
        



        local marginsList
        forvalues k = 1/7{

            su mesor if kGroup == `k'
            local mesorInitial = r(mean)
            
            noisi foreach p in 24 12 8{

                su sin`p' if  kGroup == `k'
                local sin`p'_mean = r(mean)

                su cos`p' if  kGroup == `k'
                local cos`p'_mean = r(mean)

                local amp`p'initial = sqrt((`sin`p'_mean')^2+(`cos`p'_mean')^2)
                local acr`p'initial = cond((`sin`p'_mean')<0,`p',0) + atan2((`sin`p'_mean'),(`cos`p'_mean'))*`p'/(2*_pi)
                    
                noisi di `amp`p'initial'
                noisi di `acr`p'initial'
            }


            local marginsList `marginsList' at(
            local marginsList `marginsList' amplitude24=`amp24initial'
            local marginsList `marginsList' amplitude12=`amp12initial'
            local marginsList `marginsList' amplitude8=`amp8initial'
            local marginsList `marginsList' mesor=`mesorInitial'
            local marginsList `marginsList' acrophase24_sin=(`=sin(`acr24initial'*2*_pi/24)')
            local marginsList `marginsList' acrophase24_cos=(`=cos(`acr24initial'*2*_pi/24)')
            local marginsList `marginsList' acrophase12_sin=(`=sin(`acr12initial'*2*_pi/12)')
            local marginsList `marginsList' acrophase12_cos=(`=cos(`acr12initial'*2*_pi/12)')
            local marginsList `marginsList' acrophase8_sin=(`=sin(`acr8initial'*2*_pi/8)')
            local marginsList `marginsList' acrophase8_cos=(`=cos(`acr8initial'*2*_pi/8)')
            local marginsList `marginsList' ) 

        }

        noisi margins, over(sex) `marginsList' asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")')

        marginsplot, plotdimension(sex)

        asdf
        










        forvalues curSex = 0/1{

            su mesor if sex == `curSex', detail
            local mesorInitial = r(p25)

            noisi foreach curPeriod in 24 12 8{

                su sin`curPeriod' if  mesor <= `mesorInitial' & sex == `curSex'
                local sin`curPeriod'_mean = r(mean)

                su cos`curPeriod' if  mesor <= `mesorInitial' & sex == `curSex' 
                local cos`curPeriod'_mean = r(mean)

                local amp`curPeriod'initial = sqrt((`sin`curPeriod'_mean')^2+(`cos`curPeriod'_mean')^2)
                local acr`curPeriod'initial = cond((`sin`curPeriod'_mean')<0,`curPeriod',0) + atan2((`sin`curPeriod'_mean'),(`cos`curPeriod'_mean'))*`curPeriod'/(2*_pi)
                
                noisi di `amp`curPeriod'initial'
                noisi di `acr`curPeriod'initial'
            }

            local marginsList
            forvalues j = -8(1)8{

                local acr24target     = `acr24initial'+0.25*`j'
                local acr12target     = `acr12initial'-0.125*`j'
                local mesorOffset = 0

                if `j' != 0{

                    #delimit ;
                    mata: estimateShift(    `amp24initial'   , 
                                            `amp12initial'   , 
                                            `amp8initial'    ,
                                            `acr24initial'   , 
                                            `acr12initial'   , 
                                            `acr8initial'    ,
                                            `mesorInitial'   , 
                                            `acr24target'    , 
                                            `acr12target'    
                                            )
                                            ;
                    #delimit cr
                }
                
                local marginsList `marginsList' at(
                local marginsList `marginsList' amplitude24=`amp24initial'
                local marginsList `marginsList' amplitude12=`amp12initial'
                local marginsList `marginsList' amplitude8=`amp8initial'
                local marginsList `marginsList' mesor=`=`mesorInitial'+`mesorOffset''
                local marginsList `marginsList' acrophase24_sin=(`=sin(`acr24target'*2*_pi/24)')
                local marginsList `marginsList' acrophase24_cos=(`=cos(`acr24target'*2*_pi/24)')
                local marginsList `marginsList' acrophase12_sin=(`=sin(`acr12target'*2*_pi/12)')
                local marginsList `marginsList' acrophase12_cos=(`=cos(`acr12target'*2*_pi/12)')
                local marginsList `marginsList' acrophase8_sin=(`=sin(`acr8initial'*2*_pi/8)')
                local marginsList `marginsList' acrophase8_cos=(`=cos(`acr8initial'*2*_pi/8)')
                local marginsList `marginsList' ) 

            }

            noisi margins if sex==`curSex', `marginsList' asobserved

            #delimit ;
            marginsplot,    

                            title("")
                            xlab(   , labsize(2.5) labcolor(black) angle(0) nogrid)
                            ylab(#4 , labsize(2.5) labcolor(black) angle(0) nogrid)

                            graphregion(color(white%0))
                            plotregion(color(white%0))
                            legend(off)
                            name(g`curSex', replace)
                            nolabels  
                            ;

            #delimit cr
        }

        graph combine g0 g1, imargins(0 0 0 0) ycommon xcommon

        asdf

        asdf











    
        #delimit ;
        noisi margins,
                    over(sex)

                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.616
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2-2)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2-2)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7+1)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7+1)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )
                    
                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.697
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2-1.5)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2-1.5)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7+0.75)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7+0.75)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )
                    
                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.761
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2-1)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2-1)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7+0.5)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7+0.5)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )
                    
                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.794
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2-0.5)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2-0.5)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7+0.25)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7+0.25)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )
                    
                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.789
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2-0)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2-0)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7-0)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7-0)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )
                    
                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.747
                    //totalPAEE=41.92  

                    acrophase24_sin=(`=sin((15.2+0.5)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2+0.5)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7-0.25)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7-0.25)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )
                    
                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.677
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2+1)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2+1)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7-0.50)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7-0.50)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )

                    at(
                    //amplitude24=1.89 
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.594
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2+1.5)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2+1.5)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7-0.75)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7-0.75)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )

                    at(
                    //amplitude24=1.89
                    //amplitude12=0.93 
                    //amplitude8=0.31
                    //(p80) mesor amplitude24 amplitude12 amplitude8
                    //mesor=2.511
                    //totalPAEE=41.92 

                    acrophase24_sin=(`=sin((15.2+2)*2*_pi/24)')
                    acrophase24_cos=(`=cos((15.2+2)*2*_pi/24)')
                    acrophase12_sin=(`=sin((9.7-1)*2*_pi/12)')
                    acrophase12_cos=(`=cos((9.7-1)*2*_pi/12)')
                    acrophase8_sin=(`=sin(0.18*2*_pi/8)')
                    acrophase8_cos=(`=cos(0.18*2*_pi/8)')
                    )
                    
                    asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")')
                    ;
        #delimit cr


        #delimit ;
        marginsplot,    plotdimension(sex) 

                        title("")
                        xlab(   , labsize(2.5) labcolor(black) angle(0) nogrid)
                        ylab(#4 , labsize(2.5) labcolor(black) angle(0) nogrid)

                        graphregion(color(white))
                        legend(off)

                        nolabels  
                        ;

        #delimit cr


        asdf

        ***********************************************
        ** PA volume stratified time-response curves **
        ***********************************************

        foreach curPercentile in p25 p50 p75{

            local plotResolution = 64
            local marginsList
            forvalues j = 0/`plotResolution'{
                    
                local marginsList `marginsList' at(
                local marginsList `marginsList' acrophase24_sin=(`=sin(1*`j'*2*_pi/`plotResolution')')
                local marginsList `marginsList' acrophase24_cos=(`=cos(1*`j'*2*_pi/`plotResolution')')
                local marginsList `marginsList' acrophase12_sin=(`=sin(2*`j'*2*_pi/`plotResolution')')
                local marginsList `marginsList' acrophase12_cos=(`=cos(2*`j'*2*_pi/`plotResolution')')
                local marginsList `marginsList' acrophase8_sin=(`=sin(3*`j'*2*_pi/`plotResolution')')
                local marginsList `marginsList' acrophase8_cos=(`=cos(3*`j'*2*_pi/`plotResolution')')
                local marginsList `marginsList' (`curPercentile') totalPAEE
                local marginsList `marginsList' )
        
            }

            margins, over(sex) at((`curPercentile') totalPAEE) asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")') 
            local curRef_s0 = r(b)[1,1]
            local curRef_s1 = r(b)[1,2]

            margins, over(sex) `marginsList' asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")') 


            #delimit ;
            marginsplot,    plotdimension(sex) 

                            recast(line)
                            plot1opts(lcolor(navy))
                            plot2opts(lcolor(maroon))


                            recastci(rarea) 
                            ci1opts(fcolor(navy%30) lcolor(navy%0)) 
                            ci2opts(fcolor(maroon%30) lcolor(maroon%0)) 
                            
                            title("At `=regexr("`curPercentile'","p","")'th percentile" "of total PAEE", size(2.5) color(black) nospan)

                            xtitle("Hour of day", size(2.5) color(black))
                            ytitle("Predicted `curOutcomeVar'", size(2.5) color(black))

                            xlab(1 "0" 17 "6" 33 "12" 49 "18" 65 "24"   , labsize(2.5) labcolor(black) angle(0) nogrid)
                            ylab(#4                                     , labsize(2.5) labcolor(black) angle(0) nogrid)

                            graphregion(color(white))
                            legend(off)
                            yline(`curRef_s0', lcolor(navy%40))
                            yline(`curRef_s1', lcolor(maroon%40))

                            name(`curPercentile', replace)

                            nolabels  
                            ;

            #delimit cr

        }

        #delimit ;
        graph   combine p25 p50 p75
                , 
                rows(1) 
                cols(3) 
                ycommon 
                xcommon 
                graphregion(color(white) 
                margin(l=17 r=17 t=32 b=32)) 
                name("Strat", replace)
                note("`stratDeltaAIC'", ring(0) position(11) size(1.7))
                ;
        #delimit cr
        
        capture mkdir Plots
        graph save "Strat" Plots/`curOutcomeVar'_strat_m`i'.gph , replace
        graph close _all

    }
    

    local curRow = `curRow'+1   
}

frame change dataset
frame drop tempset

