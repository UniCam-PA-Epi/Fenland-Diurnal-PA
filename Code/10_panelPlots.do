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

foreach curAnalysis in totalPAEEAnalysis cosinorFeatureAnalysis{

    foreach curModel in m1 m2{

        local plotList
        foreach curVar of local outcomeVars{

            graph use Plots/`curAnalysis'/Pooled_`curVar'_`curModel'
            local plotList `plotList' Pooled_`curVar'_`curModel'
        }

        #delimit ;

        graph combine `plotList'
                , 
                row(2) 
                col(5) 
                imargin(1 1 1 1)  
                graphregion(color(white) margin(l=5 r=5 t=20 b=20))
                name(`curAnalysis'_`curModel', replace)
                ;

        #delimit cr

        capture mkdir Figures
        graph export "Figures/`curAnalysis'_`curModel'.png" , height(2000) width(2750) replace
        graph close `plotList' `curAnalysis'_`curModel'
        graph drop  `plotList' `curAnalysis'_`curModel'

    }
}

asdf



frame change dataset
frame drop tempset















