version 17.0

frame copy dataset tempset
frame change tempset

*************************
** Initilaise analysis **
*************************

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

//replace mesor       = exp(mesor) * 60/1000
replace totalPAEE   = totalPAEE  * 60/1000
//replace maxValue    = maxValue   * 60/1000

****************************************************************************************
** Outcome variables, continuous control variables, and categorical control variables **
****************************************************************************************

#delimit ;
local outcomeVars   glucose120
                    insulin
                    leptin
                    nefa
                    adiponectin
                    ldl
                    hdl
                    mbpsys
                    mbpdia
                    crp
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

local baseModel c.totalPAEE c.totalPAEE#c.totalPAEE c.totalPAEE#i.sex
local modelLevel1 `contCovVars' `catCovVars'
local modelLevel2 `contCovVars' `catCovVars' c.fatMass c.fatMass#c.fatMass c.fatFreeMass#c.fatFreeMass c.fatMass#c.fatFreeMass


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
    putexcel L1 = ("Pvalue: cosinor model quadratic")
    putexcel M1 = ("AIC: cosinor model quadratic")
    putexcel N1 = ("Pvalue: full cosinor model vs base totalPAEE model")


}

********************************************************************************
** Begin analysis loop, looping through each outcome variable and model level **
********************************************************************************

local curRow = 2
qui foreach curOutcomeVar of local outcomeVars{

    forvalues i = 1/2{
        
        **********************************************************
        ** Apply nested GLM gaussian linear model with log link **
        *****************************************************s*****

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
                            c.mesor
                            c.amplitude24
                            c.amplitude12
                            c.amplitude8
                            c.acrophase24_sin c.acrophase24_cos
                            c.acrophase12_sin c.acrophase12_cos
                            c.acrophase8_sin  c.acrophase8_cos
                            )
                            (
                            c.mesor#c.totalPAEE
                            c.amplitude24#c.totalPAEE
                            c.amplitude12#c.totalPAEE
                            c.amplitude8#c.totalPAEE
                            c.acrophase24_sin#c.totalPAEE  c.acrophase24_cos#c.totalPAEE
                            c.acrophase12_sin#c.totalPAEE  c.acrophase12_cos#c.totalPAEE
                            c.acrophase8_sin#c.totalPAEE   c.acrophase8_cos#c.totalPAEE
                            )  
                            (
                            c.mesor#i.sex
                            c.amplitude24#i.sex
                            c.amplitude12#i.sex
                            c.amplitude8#i.sex
                            c.acrophase24_sin#i.sex  c.acrophase24_cos#i.sex
                            c.acrophase12_sin#i.sex  c.acrophase12_cos#i.sex 
                            c.acrophase8_sin#i.sex   c.acrophase8_cos#i.sex
                            )
                            (
                            c.mesor#c.mesor
                            c.amplitude24#c.amplitude24
                            c.amplitude12#c.amplitude12
                            c.amplitude8#c.amplitude8

                            c.acrophase24_sin#c.acrophase24_sin 
                            c.acrophase12_sin#c.acrophase12_sin
                            c.acrophase8_sin#c.acrophase8_sin

                            c.acrophase24_sin#c.acrophase12_sin
                            c.acrophase24_sin#c.acrophase8_sin
                            c.acrophase12_sin#c.acrophase8_sin

                            c.acrophase24_cos#c.acrophase12_cos
                            c.acrophase24_cos#c.acrophase8_cos
                            c.acrophase12_cos#c.acrophase8_cos
                            )
                            if
                            `curOutcomeVar' != .         
                            ,
                            family(`=cond("`curOutcomeVar'"=="crp","igaussian","gaussian")')
                            link(log)
                            ;
        #delimit cr
        
        capture mkdir Models
        estimates save Models/`curOutcomeVar'_m`i' , replace 
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
        putexcel N`curRow' = (r(p))
        putexcel clear

        if r(p) < 0.001 local curP = "p<0.001"
        else            local curP = "p=`=trim("`: display %10.3f r(p)'")'"

        local curDeltaAIC = "{&Delta}AIC=`=trim("`: display %10.1f lrMat[6,5]-lrMat[2,5]'")'"

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

            margins, over(sex) at((`curPercentile') totalPAEE)
            local curRef_s0 = r(b)[1,1]
            local curRef_s1 = r(b)[1,2]

            margins, over(sex) `marginsList' asobserved

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
                name(`curOutcomeVar'_m`i', replace)
                note("`curP'" "`curDeltaAIC'", ring(0) position(11) size(2))
                ;
        #delimit cr
        
        capture mkdir Plots
        graph save `curOutcomeVar'_m`i' Plots/`curOutcomeVar'_m`i'.gph , replace
        graph close _all

    }

    local curRow = `curRow'+1   
}

frame change dataset
frame drop tempset



       */

        /*

        *********************************************************
        ** Get point estimates for the "linear" curve features **
        *********************************************************
        
        foreach curExposure in mesor amplitude24 amplitude12 amplitude8{

            putexcel set "Results/6_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
            putexcel A`curRow' = ("`curOutcomeVar'")
            
            estimates restore fullModel
            margins, over(sex) eydx(`curExposure') asobserved post

            forvalues curSex = 0/1{
                
                nlcom _b[`curExposure':`curSex'.sex]

                local curEstimate   = r(b)[1,1]
                local curLB         = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  
                local curUB         = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])
                local curPvalue     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1]))) 

                local curEstimateCI =   `"`=trim("`: display %10.3f `curEstimate''")'"' +   ///
                                        " ("                                            +   ///
                                        `"`=trim("`: display %10.3f `curLB''")'"'       +   ///
                                        ", "                                            +   ///
                                        `"`=trim("`: display %10.3f `curUB''")'"'       +   ///
                                        ")"

                local curSigSymbol = cond(`curPvalue'<0.01,"**",cond(`curPvalue'<0.05,"*",""))
    
                local curCol = char(66+`curSex')
                putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")
            }       
        }


        ********************************************************************************
        ** Get amplitude and "ideal normalised acrophase values" for each subharmonic **
        ********************************************************************************

        foreach curExposure in acrophase24 acrophase12 acrophase8{

            putexcel set "Results/6_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
            putexcel A`curRow' = ("`curOutcomeVar'")

            estimates restore fullModel
            margins, over(sex) eydx(`curExposure'_sin `curExposure'_cos) asobserved post

            forvalues curSex = 0/1{

                // Compute the amplitude: sqrt((sinCoeff)^2 + cosCoeff)^2)

                nlcom sqrt((_b[`curExposure'_sin:`curSex'.sex])^2+(_b[`curExposure'_cos:`curSex'.sex])^2)

                local curEstimate   = r(b)[1,1]
                local curLB         = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  
                local curUB         = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])
                local curPvalue     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1]))) 

                local curEstimateCI =   `"`=trim("`: display %10.3f `curEstimate''")'"' +   ///
                                        " ("                                            +   ///
                                        `"`=trim("`: display %10.3f `curLB''")'"'       +   ///
                                        ", "                                            +   ///
                                        `"`=trim("`: display %10.3f `curUB''")'"'       +   ///
                                        ")"

                local curSigSymbol = cond(`curPvalue'<0.01,"**",cond(`curPvalue'<0.05,"*",""))

                local curCol = char(66+2*`curSex')
                putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")

                // Compute the "ideal normalised acrophase value"
                // This is a value from 0 to 1 indicating the relative time point on the current time period that maximises its effect on the current outcome
                // For example, a value of "0.5" for "acrophase24" would indicate that a "acrophase24" of "12 hours" would maximise its effect on the current outcome.
                // This should be interpretted based on the intended direction for the current outcome.  For example, maximising "insulin" would not be good in this context!
                // In that example, finding the minimum would be desirable (which would just be the normalised valued + 0.5)

                nlcom (cond(_b[`curExposure'_sin:`curSex'.sex]<0,2*_pi,0)+atan2(_b[`curExposure'_sin:`curSex'.sex], _b[`curExposure'_cos:`curSex'.sex]))*1/(2*_pi)

                local curEstimate   = r(b)[1,1]
                local curLB         = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  
                local curUB         = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])
                local curPvalue     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1]))) 

                local curEstimateCI =   `"`=trim("`: display %10.3f `curEstimate''")'"' +   ///
                                        " ("                                            +   ///
                                        `"`=trim("`: display %10.3f `curLB''")'"'       +   ///
                                        ", "                                            +   ///
                                        `"`=trim("`: display %10.3f `curUB''")'"'       +   ///
                                        ")"

                local curSigSymbol = cond(`curPvalue'<0.01,"**",cond(`curPvalue'<0.05,"*",""))

                local curCol = char(67+2*`curSex')
                putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")

            } 
        }

        */
