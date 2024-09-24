version 17.0

graph drop _all

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

replace totalPAEE   = totalPAEE  * 60/1000

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
                            c.mesor
                            c.amplitude24
                            c.amplitude12
                            c.amplitude8
                            c.acrophase24_sin c.acrophase24_cos
                            c.acrophase12_sin c.acrophase12_cos
                            c.acrophase8_sin  c.acrophase8_cos
                            )
                            (
                            c.amplitude24#c.amplitude12
                            c.amplitude24#c.amplitude8
                            c.amplitude12#c.amplitude8

                            c.mesor#c.amplitude24
                            c.mesor#c.amplitude12
                            c.mesor#c.amplitude8
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

                            if
                            `curOutcomeVar' != .         
                            ,
                            family(`=cond("`curOutcomeVar'"=="crp","igaussian","gaussian")')
                            link(`=cond("`curOutcomeVar'"=="cmrs","identity","log")')
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
        putexcel L`curRow' = (r(p))
        putexcel clear
        
        if r(p) < 0.001 local curP = "p<0.001"
        else            local curP = "p=`=trim("`: display %10.3f r(p)'")'"

        local pooledDeltaAIC = "{&Delta}AIC=`=trim("`: display %10.0f lrMat[4,5]-lrMat[2,5]'")'"
        local stratDeltaAIC  = "{&Delta}AIC=`=trim("`: display %10.0f lrMat[6,5]-lrMat[2,5]'")'"
        






        ********************************
        ** Pooled time-response curve **
        ********************************

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
            local marginsList `marginsList' )
    
        }

        margins, asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")')
        local curRef = r(b)[1,1] 

        margins, `marginsList' asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")') 

        #delimit ;
        marginsplot,    recast(line)
                        plot1opts(lcolor(navy))

                        recastci(rarea) 
                        ci1opts(fcolor(navy%30) lcolor(navy%0)) 
                        
                        title("Pooled", size(2.5) color(black) nospan)

                        xtitle("Hour of day", size(2.5) color(black))
                        ytitle("Predicted `curOutcomeVar'", size(2.5) color(black))

                        xlab(1 "0" 17 "6" 33 "12" 49 "18" 65 "24"   , labsize(2.5) labcolor(black) angle(0) nogrid)
                        ylab(#4                                     , labsize(2.5) labcolor(black) angle(0) nogrid)

                        graphregion(color(white) margin(l=45 r=45 t=32 b=32))
                        legend(off)

                        name("Pooled", replace)
                        note("`pooledDeltaAIC'", ring(0) position(11) size(2))
                        yline(`curRef', lcolor(gs8%40))
                        nolabels  
                        ;

        #delimit cr

        capture mkdir Plots
        graph save "Pooled" Plots/`curOutcomeVar'_pooled_m`i'.gph , replace
        graph close _all



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

