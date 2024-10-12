version 17.0

/**************************************************************************************************************
*  Cardiometabolic Risk Factor Associations with Multi-component Cosinor Model Features:
*  A Nested GLM Approach with Predictions at Shifted PAEE Profiles 
*
*  Author: Tomas Isaac Gonzales
*
*  This script investigates the associations between various cardiometabolic risk factors and 
*  features of a multi-component (24h, 12h, 8h) cosinor model. The analysis involves:
*
*    1.  Decomposing acrophases into sine and cosine components for linear modeling.
*
*    2.  Constructing nested generalized linear models (GLMs) with different sets of covariates,
*        including baseline characteristics, body composition measures, and interactions.
*
*    3.  Fitting the nested models and conducting likelihood ratio tests to evaluate the 
*        contribution of different predictor blocks (main effects, interactions, etc.).
*
*    4.  Generating predicted outcomes for time-shifted Physical Activity Energy Expenditure (PAEE) 
*        profiles by systematically adjusting the acrophases of the cosinor model components.
*
*    5.  Visualizing the impact of shifting the PAEE profile on predicted outcomes.
*
*  The analysis uses both Gaussian and inverse-Gaussian GLMs with log or identity link functions, 
*  selected based on the characteristics of each outcome variable.
*
****************************************************************************************************************/


frame copy dataset tempset
frame change tempset

*************************
** Initilaise analysis **
*************************

** Decompose the circular acrophase variables (acrophase24, acrophase12, acrophase8) 
** into their corresponding sine and cosine components. This transformation is necessary 
** to incorporate the cyclical nature of the acrophase into GLM models.

gen acrophase24_sin = sin(acrophase24_hat*2*_pi/24)
gen acrophase24_cos = cos(acrophase24_hat*2*_pi/24)

gen acrophase12_sin = sin(acrophase12_hat*2*_pi/12)
gen acrophase12_cos = cos(acrophase12_hat*2*_pi/12)

gen acrophase8_sin  = sin(acrophase8_hat*2*_pi/8)
gen acrophase8_cos  = cos(acrophase8_hat*2*_pi/8)


/* Define Outcome variables, continuous control variables, and categorical control variables */

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

/* Define blocks of variables to be added sequentially to the nested models. */
/* baseModel: Includes the mesor (average value) of the outcome, its quadratic term, and interaction with sex. */
local baseModel c.mesor_hat 
                c.mesor_hat#c.mesor_hat 
                c.mesor_hat#i.sex
                c.mesor_hat#c.mesor_hat #i.sex
                ;

/* block1: Main effects of amplitude and acrophase for each cosinor component (24h, 12h, 8h). */
local block1    c.amplitude24_hat
                c.amplitude12_hat
                c.amplitude8_hat

                c.acrophase24_sin 
                c.acrophase24_cos

                c.acrophase12_sin 
                c.acrophase12_cos

                c.acrophase8_sin  
                c.acrophase8_cos
                ;

/* block2: Interactions between sex and the amplitude/acrophase terms. */
local block2    i.sex#c.amplitude24_hat
                i.sex#c.amplitude12_hat
                i.sex#c.amplitude8_hat

                i.sex#c.acrophase24_sin  
                i.sex#c.acrophase24_cos

                i.sex#c.acrophase12_sin  
                i.sex#c.acrophase12_cos

                i.sex#c.acrophase8_sin   
                i.sex#c.acrophase8_cos
                ;

/* block3: Interactions between mesor and amplitude/acrophase terms. */               
local block3    c.mesor_hat#c.amplitude24_hat
                c.mesor_hat#c.amplitude12_hat
                c.mesor_hat#c.amplitude8_hat

                c.acrophase24_sin#c.mesor_hat  
                c.acrophase24_cos#c.mesor_hat

                c.acrophase12_sin#c.mesor_hat  
                c.acrophase12_cos#c.mesor_hat

                c.acrophase8_sin#c.mesor_hat   
                c.acrophase8_cos#c.mesor_hat
                ;

/* block4: Interactions between the acrophase components of different periods (24h, 12h, 8h). */
local block4    c.acrophase24_sin#c.acrophase24_cos
                c.acrophase24_sin#c.acrophase12_cos
                c.acrophase24_sin#c.acrophase8_cos

                c.acrophase12_sin#c.acrophase24_cos
                c.acrophase12_sin#c.acrophase12_cos
                c.acrophase12_sin#c.acrophase8_cos

                c.acrophase8_sin#c.acrophase24_cos
                c.acrophase8_sin#c.acrophase12_cos
                c.acrophase8_sin#c.acrophase8_cos
                ;

/* block5: Interactions between the amplitude components of different periods. */
local block5    c.amplitude24_hat#c.amplitude12_hat
                c.amplitude24_hat#c.amplitude8_hat
                c.amplitude12_hat#c.amplitude8_hat
                ;

/* block6: Acrophase component standard error values */
local block6    c.acrophase24_se
                c.acrophase12_se
                c.acrophase8_se
                ;

#delimit cr

/* Initialize Excel sheet for storing likelihood ratio test results. */

capture mkdir Tables
capture erase "Tables/8_cosinorFeatureAnalysis.xlsx"

qui forvalues i = 1/2{

    putexcel set "Tables/8_cosinorFeatureAnalysis.xlsx", sheet("lrTest_m`i'") modify
    putexcel A1 = ("outcomeVar")
    putexcel B1 = ("Count: Women")
    putexcel C1 = ("Count: Men")
    putexcel D1 = ("Pvalue: base mesor model")
    putexcel E1 = ("AIC: base mesor model")
    putexcel F1 = ("Pvalue: main effects")
    putexcel G1 = ("AIC: main effects")
    putexcel H1 = ("Pvalue: sex interactions")
    putexcel I1 = ("AIC: sex interactions")
    putexcel J1 = ("Pvalue: mesor interactions")
    putexcel K1 = ("AIC: mesor interactions")
    putexcel L1 = ("Pvalue: acrophase interactions")
    putexcel M1 = ("AIC: acrophase interactions")
    putexcel N1 = ("Pvalue: amplitude interactions")
    putexcel O1 = ("AIC: amplitude interactions")
    putexcel P1 = ("Pvalue: acrophase se main effects")
    putexcel Q1 = ("AIC: acrophase se main effects")
    putexcel R1 = ("Pvalue: full model vs base model")
    putexcel S1 = ("deltaAIC: full model vs base model")

}

********************************************************************************
** Begin analysis loop, looping through each outcome variable and model level **
********************************************************************************

local curRow = 2
qui foreach curOutcomeVar of local outcomeVars{
    
    qui forvalues i = 1/2{
        
        **************************************************************************************
        ** Apply nested GLM (gaussian/inverse-gaussian) linear model with log/identity link **
        **************************************************************************************
        
        /* Note that when "crp" is the outcome var, we apply an inverse-Gaussian GLM due to its extreme right tail. */
        /* This choice helps resolve convergence issues that arise when using a Gaussian GLM for "crp". */
        /* For "cmrs" we use the identity link function. */
        /* For all other outcome vars, we use the log link function. */

        local curFamily    = cond("`curOutcomeVar'"=="crp","igaussian","gaussian") 
        local curLink      = cond("`curOutcomeVar'"=="cmrs","identity","log")

        /* Fit the base model. */
        noisi glm `curOutcomeVar' `modelLevel`i'' `baseModel'    if `curOutcomeVar' != . , family(`curFamily') link(`curLink')
        estimates store baseModel

        /* Fit the full model by sequentially adding blocks of predictors. */

        local nestedModel (`modelLevel`i'')(`baseModel')(`block1')(`block2')(`block3')(`block4')(`block5')(`block6')

        noisi nestreg, lr: glm `curOutcomeVar' `nestedModel' if `curOutcomeVar' != . , family(`curFamily') link(`curLink')
        estimates store fullModel


        ****************************************
        ** Store GLM model diagnostic outputs **
        ****************************************

        capture mkdir Models
        capture mkdir Models/cosinorFeatureAnalysis
        estimates save Models/cosinorFeatureAnalysis/`curOutcomeVar'_m`i' , replace

        mat lrMat = r(lr)

        /* Store results of nested lr tests for the contribution of each block to the model */

        putexcel set "Tables/8_cosinorFeatureAnalysis.xlsx", sheet("lrTest_m`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")

        count if e(sample) == 1 & sex == 0
        putexcel B`curRow' = (r(N))

        count if e(sample) == 1 & sex == 1
        putexcel C`curRow' = (r(N))

        /* Store p values and AIC values */

        forvalues j = 0/7{
            local curCol = char(68+`=`j'*2')  
            putexcel `curCol'`curRow' = (lrMat[`=`j'+2',4])
            local curCol = char(69+`=`j'*2')
            putexcel `curCol'`curRow' = (lrMat[`=`j'+2',5])
        }

        /* Perform LR test between base model and full model */
         /* Store strings for p and delta AIC values for plotting */
         
        noisi lrtest baseModel fullModel

        putexcel R`curRow' = (r(p))
        if r(p) < 0.001 local curP = "p<0.001"
        else            local curP = "p=`=trim("`: display %10.3f r(p)'")'"

        putexcel S`curRow' = (lrMat[8,5]-lrMat[2,5])
        local deltaAIC = "{&Delta}AIC=`=trim("`: display %10.0f lrMat[8,5]-lrMat[2,5]'")'" 
    
        putexcel clear

        
        
        **************************************************************************
        ** Generate and Plot Predicted Outcome Curves for Shifted PAEE Profiles **
        **************************************************************************

        /* Loop through each time profile extracted from k-means cluster analysis. */
        /* "0" represents the pooled sample (all clusters combined). */
        /* Note that "5_clusterAnalysis.do" must have been run previously before the code block below can be run. */

        levelsof kGroup, local(kVals)
        local kVals 0 `kVals'

        foreach k of local kVals{
            
            /* Calculate the mean mesor for the current cluster (or pooled sample). */ 
            /* The same logic is applied to conditional statements below. */
            su mesor_hat if kGroup `=cond(`k'==0,">","==")' `k'
            local mesor_`k' = r(mean)
            
            /* Compute amplitude and acrophase for each cosinor component (24h, 12h, 8h). */
            foreach p in 24 12 8{

                /* Calculate the mean sine and cosine values for the current period. */
                su sin`p'_hat if  kGroup `=cond(`k'==0,">","==")' `k'
                local sin`p'_mean = r(mean)

                su cos`p'_hat if  kGroup `=cond(`k'==0,">","==")' `k'
                local cos`p'_mean = r(mean)

                /* Calculate amplitude and acrophase using the mean sine and cosine values. */
                local amp`p'_`k' = sqrt((`sin`p'_mean')^2+(`cos`p'_mean')^2)
                local acr`p'_`k' = cond((`sin`p'_mean')<0,`p',0) + atan2((`sin`p'_mean'),(`cos`p'_mean'))*`p'/(2*_pi)
                    
            }

            ** Generate predicted outcome values for shifted PAEE profiles.
            ** This involves shifting the acrophases of the 24h, 12h, and 8h cosinor components 
            ** by a range of hourly values and calculating predicted outcomes at each shift.
            local marginsList
            forvalues shift = -1.5(0.125)1.5{
                
                ** Construct the list of values for the margins command, including shifted acrophases.
                local marginsList `marginsList' at(
                local marginsList `marginsList' amplitude24_hat=`amp24_`k''
                local marginsList `marginsList' amplitude12_hat=`amp12_`k''
                local marginsList `marginsList' amplitude8_hat=`amp8_`k''
                local marginsList `marginsList' mesor_hat=`mesor_`k''
                local marginsList `marginsList' acrophase24_sin=(`=sin((`acr24_`k''+`shift')*2*_pi/24)')
                local marginsList `marginsList' acrophase24_cos=(`=cos((`acr24_`k''+`shift')*2*_pi/24)')
                local marginsList `marginsList' acrophase12_sin=(`=sin((`acr12_`k''+`shift')*2*_pi/12)')
                local marginsList `marginsList' acrophase12_cos=(`=cos((`acr12_`k''+`shift')*2*_pi/12)')
                local marginsList `marginsList' acrophase8_sin=(`=sin((`acr8_`k''+`shift')*2*_pi/8)')
                local marginsList `marginsList' acrophase8_cos=(`=cos((`acr8_`k''+`shift')*2*_pi/8)')
                local marginsList `marginsList' ) 

            }
            
            ** Calculate predicted outcomes using the margins command with the specified shifts.
            margins, over(sex) `marginsList' asobserved predict(`=cond("`curOutcomeVar'"=="crp","xb","")')

            ** Extract predicted outcome values for women (ref0) and men (ref1) at the zero-shift point (i.e., no shift applied to the PAEE profile). 
            ** These values will be used as reference lines in the plot.
            local ref0 = r(b)[1,25]
            local ref1 = r(b)[1,26]

            ** Construct the predicted outcome curve and save to file.
            local plotName `=cond(`k'==0,"Pooled_`curOutcomeVar'_m`i'","k`k'_`curOutcomeVar'_m`i'")' 
            #delimit ;
            marginsplot,    plotdimension(sex) 

                            recast(line)
                            plot1opts(lcolor("214 40  40 %100"))
                            plot2opts(lcolor("0   47  73 %100"))

                            recastci(rarea) 
                            ci1opts(fcolor("214 40  40 %30") lcolor("214 40  40 %0")) 
                            ci2opts(fcolor("0   47  73 %30") lcolor("0   47  73 %0")) 
                            
                            title("", size(2.5) color(black) nospan)

                            xtitle("Shift in PAEE profile (hour)", size(2.5) color(black))
                            ytitle("Predicted `curOutcomeVar'", size(2.5) color(black))

                            xlab(1 "-1.5" 5 "-1" 9 "-0.5" 13 "0" 17 "0.5" 21 "1" 25 "1.5" , labsize(2.5) labcolor(black) angle(0) nogrid)
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
            capture mkdir Plots/cosinorFeatureAnalysis
            graph save "`plotName'" Plots/cosinorFeatureAnalysis/`plotName'.gph , replace
            graph close "`plotName'"
            graph drop "`plotName'"
        }
        
    }

    ** Increment Excel spreadsheet to next row for the the next outcome variable.
    local curRow = `curRow'+1  
}

frame change dataset
frame drop tempset
