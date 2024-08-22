version 17.0

frame copy dataset tempset
frame change tempset

*************************
** Initilaise analysis **
*************************

// Exclude participants for whom the cosinor model did not fit, those with insufficient wear time, and those without fat mass measurement.

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05
drop if P_Pwear_consolidated<72
drop if fatMass == .

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
//replace totalPAEE   = totalPAEE  * 60/1000
//replace maxValue    = maxValue   * 60/1000

****************************************************************************************
** Outcome variables, continuous control variables, and categorical control variables **
****************************************************************************************

#delimit ;
local outcomeVars   //glucose120
                    //insulin
                    //nefa
                    //leptin
                    //adiponectin
                    //ldl
                    //hdl
                    //fatFreeMass
                    //mbpsys
                    //mbpdia
                    crp
                    ;

local contCovVars   c.age
                    c.rhr
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
                    i.season
                    i.cardiometabol_med
                    i.testsite
                    ;

#delimit cr

// Initialise two model levels, with the second level controlling for adiposity  (i.e.fat mass)

local modelLevel1 `contCovVars' `catCovVars'
local modelLevel2 `contCovVars' `catCovVars' c.fatMass
local modelLevel3 `contCovVars' `catCovVars' c.fatMass c.totalPAEE
local modelLevel4 `contCovVars' `catCovVars' c.fatMass c.totalPAEE c.maxValue
local modelLevel5 `contCovVars' `catCovVars' c.fatMass c.totalPAEE c.maxValue c.maxHour24_sin c.maxHour24_cos

**********************************
** Initialise excel spreadsheet **
**********************************

capture mkdir Results
//capture erase "Results/6_jointAssociations.xlsx"

qui forvalues i = 1/5{

    putexcel set "Results/6_jointAssociations.xlsx", sheet("waldBlockTest_m`i'") modify
    putexcel A1 = ("outcomeVar")
    putexcel B1 = ("Pvalue: mesor")
    putexcel C1 = ("Pvalue: mesor#sex")
    putexcel D1 = ("Pvalue: amplitude24")
    putexcel E1 = ("Pvalue: amplitude24#sex")
    putexcel F1 = ("Pvalue: amplitude12")
    putexcel G1 = ("Pvalue: amplitude12#sex")
    putexcel H1 = ("Pvalue: amplitude8")
    putexcel I1 = ("Pvalue: amplitude8#sex")
    putexcel J1 = ("Pvalue: acrophase24")
    putexcel K1 = ("Pvalue: acrophase24#sex")
    putexcel L1 = ("Pvalue: acrophase12")
    putexcel M1 = ("Pvalue: acrophase12#sex")
    putexcel N1 = ("Pvalue: acrophase8")
    putexcel O1 = ("Pvalue: acrophase8#sex")
    putexcel P1 = ("Count: Women")
    putexcel Q1 = ("Count: Men")

    foreach curExposure in mesor amplitude24 amplitude12 amplitude8{

        putexcel set "Results/6_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
        putexcel A1 = ("outcomeVar")
        putexcel B1 = ("Estimate: Women")
        putexcel C1 = ("Estimate: Men")
    }

    foreach curExposure in acrophase24 acrophase12 acrophase8{

        putexcel set "Results/6_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
        putexcel A1 = ("outcomeVar")
        putexcel B1 = ("Amplitude: Women")
        putexcel C1 = ("Normalised acrophase: Women")
        putexcel D1 = ("Amplitude: Men")
        putexcel E1 = ("Normalised acrophase: Men")

    }
}


********************************************************************************
** Begin analysis loop, looping through each outcome variable and model level **
********************************************************************************

local curRow = 1
foreach curOutcomeVar of local outcomeVars{

    local curRow = `curRow'+1

    forvalues i = 4/5{
        
        **********************************************************
        ** Apply nested GLM gaussian linear model with log link **
        *****************************************************s*****

        // Note that when "crp" is the outcome var, we apply an inverse-Gaussian GLM due to the extreme right tail.
        // This choice helps resolve convergence issues that may arise when using a Gaussian GLM for "crp".
        // This choice was verified using BIC for all outcome vars.

        #delimit ;
        
        nestreg : glm   `curOutcomeVar' 
                        (`modelLevel`i'')

                        (c.mesor)
                        (c.mesor#i.sex)
                        
                        (c.amplitude24)
                        (c.amplitude24#i.sex)

                        (c.amplitude12)
                        (c.amplitude12#i.sex)

                        (c.amplitude8)
                        (c.amplitude8#i.sex)

                        (c.acrophase24_sin        c.acrophase24_cos)
                        (c.acrophase24_sin#i.sex  c.acrophase24_cos#i.sex)

                        (c.acrophase12_sin        c.acrophase12_cos)
                        (c.acrophase12_sin#i.sex  c.acrophase12_cos#i.sex)

                        (c.acrophase8_sin          c.acrophase8_cos)
                        (c.acrophase8_sin#i.sex    c.acrophase8_cos#i.sex)                      
                        ,
                        family(`=cond("`curOutcomeVar'"=="crp","igaussian","gaussian")')
                        link(log)
                        ;
        #delimit cr
        
        capture mkdir Models
        estimates save Models/`curOutcomeVar'_m`i' , replace
        estimates store fullModel
        
        // Store results of nested Wald tests for the contribution of each block to the model
        
        putexcel set "Results/6_jointAssociations.xlsx", sheet("waldBlockTest_m`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")

        forvalues j = 1/18{
            local curCol = char(65+`j')
            putexcel `curCol'`curRow' = (r(wald)[`=`j'+1',3])

        }
        count if e(sample) == 1 & sex == 0
        putexcel P`curRow' = ("`r(N)'")

        count if e(sample) == 1 & sex == 1
        putexcel Q`curRow' = ("`r(N)'")

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
    }   
}

frame change dataset
frame drop tempset
