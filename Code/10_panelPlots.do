version 17.0

frame copy dataset tempset
frame change tempset

#delimit ;

/* List of outcome variables that were analyzed. */
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

#delimit cr

// Loop through each analysis permutation, fetching plots from each, and construct panel figures.

foreach curGroup in Pooled k1 k2 k3 k4 k5 k6{

    foreach curAnalysis in totalPAEEAnalysis cosinorFeatureAnalysis{

        foreach curModel in m1 m2{

            local plotList
            foreach curVar of local outcomeVars{

                capture graph use Plots/`curAnalysis'/`curGroup'_`curVar'_`curModel'
                local plotList `plotList' `curGroup'_`curVar'_`curModel'
            }

            #delimit ;

            capture graph combine `plotList'
                    , 
                    row(2) 
                    col(5) 
                    imargin(1 1 1 1)  
                    graphregion(color(white) margin(l=15 r=15 t=30 b=30))
                    name(`curGroup'_`curAnalysis'_`curModel', replace)
                    ;

            #delimit cr
            
            if _rc == 0{

                capture mkdir Figures
                graph export "Figures/`curGroup'_`curAnalysis'_`curModel'.png" , height(2000) width(2750) replace
                graph close `plotList' `curGroup'_`curAnalysis'_`curModel'
                graph drop  `plotList' `curGroup'_`curAnalysis'_`curModel'

            }
        }
    }
}

frame change dataset
frame drop tempset















