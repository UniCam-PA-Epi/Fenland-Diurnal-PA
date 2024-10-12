version 17.0

frame copy dataset tempset
frame change tempset


replace totalPAEE_hat = totalPAEE_hat*60/1000

#delimit ;

/* List of outcome variables to be analyzed. */
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

/* List of continuous covariates, including quadratic terms and cyclical time variables. */ 
local contCovVars   c.age#i.sex
                    c.age#c.age#i.sex

                    c.rhr#i.sex
                    c.rhr#c.rhr#i.sex

                    c.testTime_sin
                    c.testTime_cos
                    c.testTime_sin#c.testTime_cos
                    
                    c.testDay_sin
                    c.testDay_cos
                    c.testDay_sin#c.testDay_cos
                    ;

/* List of categorical covariates. */ 
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


/* Define two sets of covariates for the nested models. */
/* modelLevel1: Includes baseline covariates (continuous and categorical). */
local modelLevel1   `contCovVars' 
                    `catCovVars'
                    ;

/* modelLevel2: Extends modelLevel1 by adding body composition variables (fat mass and fat-free mass), */
/* including their interactions with sex and quadratic terms. */
local modelLevel2   `contCovVars' 
                    `catCovVars' 
                    c.fatMass#i.sex 
                    c.fatMass#c.fatMass#i.sex 
                    c.fatFreeMass#i.sex 
                    c.fatFreeMass#c.fatFreeMass#i.sex
                    ;

/* baseModel: Includes the mesor (average value) of the outcome, its quadratic term, and interaction with sex. */
local baseModel c.totalPAEE_hat
                c.totalPAEE_hat#c.totalPAEE_hat 
                c.totalPAEE_hat#i.sex
                c.totalPAEE_hat#c.totalPAEE_hat#i.sex
                ; 

#delimit cr

/* Initialize Excel sheet for storing likelihood ratio test results. */

capture mkdir Tables
capture erase "Tables/9_totalPAEEAnalysis.xlsx"

qui forvalues i = 1/2{

    putexcel set "Tables/9_totalPAEEAnalysis.xlsx", sheet("lrTest_m`i'") modify
    putexcel A1 = ("outcomeVar")
    putexcel B1 = ("Count: Women")
    putexcel C1 = ("Count: Men")
    putexcel D1 = ("Pvalue: full model vs base model")
    putexcel E1 = ("deltaAIC: full model vs base model")

}

********************************************************************************
** Begin analysis loop, looping through each outcome variable and model level **
********************************************************************************

local curRow = 2
qui foreach curOutcomeVar of local outcomeVars{
    
    qui forvalues i = 1/2{
        
        *************************************************************************************
        ** Apply nested GLM (gaussian/inverse-gaussian) linear model with log/identity link **
        **************************************************************************************
        
        /* Note that when "crp" is the outcome var, we apply an inverse-Gaussian GLM due to its extreme right tail. */
        /* This choice helps resolve convergence issues that arise when using a Gaussian GLM for "crp". */
        /* For "cmrs" we use the identity link function. */
        /* For all other outcome vars, we use the log link function. */

        local curFamily    = cond("`curOutcomeVar'"=="crp","igaussian","gaussian") 
        local curLink      = cond("`curOutcomeVar'"=="cmrs","identity","log")

        noisi nestreg, lr: glm `curOutcomeVar' (`modelLevel`i'')(`baseModel') if `curOutcomeVar' != . , family(`curFamily') link(`curLink')

        capture mkdir Models
        capture mkdir Models/totalPAEEAnalysis
        estimates save Models/totalPAEEAnalysis/`curOutcomeVar'_m`i' , replace

        mat lrMat = r(lr)

        /* Store results of nested lr tests for the contribution of each block to the model */

        putexcel set "Tables/9_totalPAEEAnalysis.xlsx", sheet("lrTest_m`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")

        count if e(sample) == 1 & sex == 0
        putexcel B`curRow' = (r(N))

        count if e(sample) == 1 & sex == 1
        putexcel C`curRow' = (r(N))

        putexcel D`curRow' = (lrMat[2,4])
        if r(p) < 0.001 local curP = "p<0.001"
        else            local curP = "p=`=trim("`: display %10.3f r(p)'")'"

        putexcel E`curRow' = (lrMat[2,5]-lrMat[1,5])
        local deltaAIC = "{&Delta}AIC=`=trim("`: display %10.0f lrMat[2,5]-lrMat[1,5]'")'"

        putexcel clear

        ********************************************************************************
        ** Generate and Plot Predicted Outcome Curves for different total PAEE values **
        ********************************************************************************

        local marginsList
        forvalues j = 20(2.5)100{
            local marginsList `marginsList' at(totalPAEE_hat=`j')
        }

        ** Extract predicted outcome values for women (ref0) and men (ref1) at median total PAEE value
        ** These values will be used as reference lines in the plot.
        noisi margins, over(sex) at((p50) totalPAEE_hat)  asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")')
        local ref0 = r(b)[1,1]
        local ref1 = r(b)[1,2]

         ** Calculate predicted outcomes using the margins command with the specified total PAEE values.
        noisi margins, over(sex) `marginsList' asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")')
        
         ** Construct the predicted outcome curve and save to file.
        local plotName "Pooled_`curOutcomeVar'_m`i'"
        #delimit ;
        marginsplot,    plotdimension(sex) 

                        recast(line)
                        plot1opts(lcolor("214 40  40 %100"))
                        plot2opts(lcolor("0   47  73 %100"))

                        recastci(rarea) 
                        ci1opts(fcolor("214 40  40 %30") lcolor("214 40  40 %0")) 
                        ci2opts(fcolor("0   47  73 %30") lcolor("0   47  73 %0")) 
                        
                        title("", size(2.5) color(black) nospan)

                        xtitle("Total PAEE", size(2.5) color(black))
                        ytitle("Predicted `curOutcomeVar'", size(2.5) color(black))

                        xlab(1 "20" 9 "40" 17 "60" 25 "80" 33 "100" , labsize(2.5) labcolor(black) angle(0) nogrid)
                        ylab(#4 , labsize(2.5) labcolor(black) angle(0) nogrid)

                        plotregion(color(white))
                        graphregion(color(white))
                        legend(off)

                        yline(`ref0', lcolor("214 40  40 %50") lpattern(dash) noextend)
                        yline(`ref1', lcolor("0   47  73 %50") lpattern(dash) noextend)

                        name("`plotName'", replace)

                        nolabels  
                        ;
        #delimit cr

        capture mkdir Plots
        capture mkdir Plots/totalPAEEAnalysis
        graph save "`plotName'" Plots/totalPAEEAnalysis/`plotName'.gph , replace
        graph close "`plotName'"
        graph drop "`plotName'"
        
    }
    ** Increment Excel spreadsheet to next row for the the next outcome variable.
    local curRow = `curRow'+1  
}

frame change dataset
frame drop tempset


