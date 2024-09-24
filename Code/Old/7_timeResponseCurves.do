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



#delimit ;
local outcomeVars   //glucose120
                    insulin
                    //nefa
                    //leptin
                    //adiponectin
                    //ldl
                    //hdl
                    //fatFreeMass
                    //mbpsys
                    //mbpdia
                    //crp
                    ;
#delimit cr

foreach curOutcomeVar of local outcomeVars{

    forvalues i = 4/5{

        ********************************************************************************
        ** Construct acrophase time-response curves for each phase value (24, 12, 18) **
        ********************************************************************************

        local dothis = 1

        if `dothis' == 1{

        foreach curExposure in acrophase24 acrophase12 acrophase8{
            
            // Get p-values for the statistical significance of acrophase time-response curves for each sex
            // This is determined from the signficance of the derived amplitude feature of the cosinor model
            // If the amplitude is not significant, the curve is not significant.

            estimates use Models/`curOutcomeVar'_m`i'
            estimates esample:
            margins, over(sex) eydx(`curExposure'_sin `curExposure'_cos) asobserved post

            forvalues curSex = 0/1{

                nlcom sqrt((_b[`curExposure'_sin:`curSex'.sex])^2+(_b[`curExposure'_cos:`curSex'.sex])^2)
                local Pvalue_s`curSex' = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

                if `Pvalue_s`curSex'' < 0.001   local Pvalue_s`curSex' = "p<0.001"
                else                            local Pvalue_s`curSex' = "p=`=trim("`: display %10.3f `Pvalue_s`curSex'''")'"

            }

            // Initialise marginal value points. 32 points resolution produces good looking plots, so going with that.
            
            local plotResolution = 32
            local marginsList
            forvalues j = 0/`plotResolution'{
                
                local marginsList `marginsList' at(`curExposure'_sin=(`=sin(`j'*2*_pi/`plotResolution')')`curExposure'_cos=(`=cos(`j'*2*_pi/`plotResolution')'))

            }

            forvalues curSex = 0/1{
                
                // Compute marginal values at the selected points: point estimates (hat) and upper and lower bounds  (lb, ub)
                // Also compute the mean overall values for the current sex, which is to compute relative differences.

                gen hat = .
                gen lb = .
                gen ub = .

                estimates use Models/`curOutcomeVar'_m`i'
                estimates esample:
                margins if sex==`curSex', at() `marginsList' asobserved post

                // Compute relative difference between each point estimate (i.e. currentValue) and the overall mean (i.e. referenceValue):
                // (currentValue - referenceValue)/referenceValue

                forvalues j = 1/`=1+`plotResolution''{

                    nlcom (_b[`=`j'+1'._at]-_b[1._at])/_b[1._at]
                    replace hat = r(b)[1,1]                                     in `j'
                    replace lb  = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  in `j'
                    replace ub  = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])  in `j'

                }

                // Setup time-axis labeling

                if regexm("`curExposure'","24") == 1 local curTime = 24 
                if regexm("`curExposure'","12") == 1 local curTime = 12 
                if regexm("`curExposure'","8")  == 1 local curTime = 8
                            
                gen timePLot = `curTime'*(_n-1)/`plotResolution' if hat !=.
                local timeLab 0 `=1*`curTime'/4' `=2*`curTime'/4' `=3*`curTime'/4' `curTime'

                // Plot current acrophase-response curve based on all values computed above

                set graphics on

                #delimit ;
                twoway  (
                        rarea lb ub timePLot if hat !=.
                        ,
                        fcolor(navy%40)
                        lcolor(navy%0)
                        )
                        (
                        line hat timePLot  if hat !=.
                        ,
                        lcolor(navy)
                        )
                        (
                        func y=0
                        ,
                        lcolor(gs8)
                        lpattern(dash)
                        lwidth(0.2)
                        range(0 `curTime')
                        )
                        (
                        ,
                        xlab(`timeLab'   , nogrid angle(0))
                        ylab(-.5(0.25).5 , nogrid angle(0))
                        graphregion(color(white))
                        name("`curOutcomeVar'_`curExposure'_s`curSex'_m`i'", replace)

                        note("`Pvalue_s`curSex''", ring(0) position(1))
                        legend(off)                       
                        )
                        ;                      
                #delimit cr

                capture mkdir Plots
                graph save "`curOutcomeVar'_`curExposure'_s`curSex'_m`i'" "Plots/`curOutcomeVar'_`curExposure'_s`curSex'_m`i'.gph", replace
                graph close "`curOutcomeVar'_`curExposure'_s`curSex'_m`i'"

                drop hat lb ub timePLot
                
            }      
        }

        }

        *****************************************************************************************************
        ** Construct acrophase time-response curve illustrating the effect of transitioning from the least **  
        ** ideal to the most ideal combination of acrophases for the 24-hour, 12-hour, and 8-hour cycles   **
        *****************************************************************************************************

        // In essence, the "least ideal" combination of acrophases represents the scenario where the timing of 
        // peak PA for each of the 24-hour, 12-hour, and 8-hour cycles is least favorable for metabolic health, 
        // based on the model's coefficients.


        // Calculate the "ideal" acrophase (in radians) for each component cycle, separately for each sex

        foreach curExposure in acrophase24 acrophase12 acrophase8{
            
            estimates use Models/`curOutcomeVar'_m`i'
            estimates esample:
            margins, over(sex) eydx(`curExposure'_sin `curExposure'_cos) asobserved post

            forvalues curSex = 0/1{

                nlcom (cond(_b[`curExposure'_sin:`curSex'.sex]<0,2*_pi,0)+atan2(_b[`curExposure'_sin:`curSex'.sex], _b[`curExposure'_cos:`curSex'.sex]))
                local `curExposure'_s`curSex' = r(b)[1,1]

            }
        }

        // For each sex, generate the combined acrophase time-response curve

        forvalues curSex = 0/1{

            // Define points along the curve, with resolution for smooth plotting.
            // Origin is at the "least ideal" combination (estimated acrophases).

            local plotResolution = 32
            local marginsList
            forvalues j = 0/`plotResolution'{
                
                // Calculate sin and cos components for each phase incrementally shifting through one cycle
                // Here we are holding the 12 hour cycle at its worst location.

                local p24sin acrophase24_sin=(`=sin(`acrophase24_s`curSex''+(`j'*2*_pi/`plotResolution'))')
                local p24cos acrophase24_cos=(`=cos(`acrophase24_s`curSex''+(`j'*2*_pi/`plotResolution'))')

                //local p24sin acrophase24_sin=(`=sin(`acrophase24_s`curSex'')')
                //local p24cos acrophase24_cos=(`=sin(`acrophase24_s`curSex'')')
                //local p24all `p24sin'`p24cos'

                local p12sin acrophase12_sin=(`=sin(`acrophase12_s`curSex''+(2*`j'*2*_pi/`plotResolution'))')
                local p12cos acrophase12_cos=(`=cos(`acrophase12_s`curSex''+(2*`j'*2*_pi/`plotResolution'))')
                local p12all `p12sin'`p12cos'

                local p8sin acrophase8_sin=(`=sin(`acrophase8_s`curSex''+(3*`j'*2*_pi/`plotResolution'))')
                local p8cos acrophase8_cos=(`=cos(`acrophase8_s`curSex''+(3*`j'*2*_pi/`plotResolution'))')
                local p8all `p8sin'`p8cos'

                //local marginsList `marginsList' at(`p24all'`p12all'`p8all')
                local marginsList `marginsList' at(`p12all'`p8all')
            }

            gen hat = .
            gen lb = .
            gen ub = .

            // Calculate predicted outcomes & CIs at each point on the curve

            estimates use Models/`curOutcomeVar'_m`i'
            estimates esample:
            margins if sex==`curSex', `marginsList' asobserved post
            
                // Calculate the relative change in outcome compared to the "least ideal" point
            forvalues j = 1/`=`plotResolution'+1'{

                nlcom (_b[`j'._at]-_b[1._at])/_b[1._at]
                replace hat = r(b)[1,1]                                     in `j'
                replace lb  = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  in `j'
                replace ub  = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])  in `j'

            }
            
            // Generate the time axis for plotting (0 to 1 represents the full transition)

            gen timePLot = 24*(_n-1)/`plotResolution' if hat !=.

            set graphics on

            #delimit ;
            twoway  (
                    rarea lb ub timePLot if hat !=.
                    ,
                    fcolor(navy%40)
                    lcolor(navy%0)
                    )
                    (
                    line hat timePLot  if hat !=.
                    ,
                    lcolor(navy)
                    )
                    (
                    func y=0
                    ,
                    lcolor(gs8)
                    lpattern(dash)
                    lwidth(0.2)
                    range(0 24)
                    )
                    (
                    ,
                    xlab(0(6)24    , nogrid angle(0))
                    ylab(-.5(0.25).5 , nogrid angle(0))
                    graphregion(color(white))
                    name("`curOutcomeVar'_acrophaseAll_s`curSex'_m`i'", replace)

                    //note("`Pvalue_s`curSex''", ring(0) position(1))
                    legend(off)                       
                    )
                    ;                      
            #delimit cr

            capture mkdir Plots
            graph save "`curOutcomeVar'_acrophaseAll_s`curSex'_m`i'" "Plots/`curOutcomeVar'_acrophaseAll_s`curSex'_m`i'.gph", replace
            //graph close "`curOutcomeVar'_acrophaseAll_s`curSex'_m`i'"

            drop hat lb ub timePLot

        }
    }
}

frame change dataset
frame drop tempset