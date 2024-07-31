version 17.0

clear
set more off
set seed 1234

args filePath


********************************************************************************
** Descriptives Tables (Table 1)
********************************************************************************

use "`filePath'"

qui do Code/Sub/initialiseCovariates.do
qui do Code/Sub/initialiseOutcomes.do
qui do Code/Sub/initialisePAEE.do

// Define continuous and categorical variables

#delimit ;
local contVars  age
                smoking
                diet
                alcohol
                paeeTt
                bodyfat
                insulin
                leptin
                nefa
                ldl
                adiponectin
                crp
                mbpdia
                mbpsys
                ;

local catVars   ethnic
                education
                income
                work_s
                marital_s
                season
                cardiometabol_med
                testsite
                glucose_cat
                ;

#delimit cr



levelsof sex, local(sexLevels)
foreach curSex of local sexLevels{
    
    local curSexLabel : label (sex) `curSex'
    putexcel set "Results/1_descriptiveTables", sheet("`curSexLabel'") modify

    preserve
    keep if sex == `curSex'
    local curRow = 0

    // Get total count for current sex

    count
    local curTotalCount = r(N)

    local curRow = `curRow' + 1
    putexcel A`curRow' = ("N")
    putexcel B`curRow' = `curTotalCount'

    // Loop through continuous vars

    foreach curVar of local contVars{

        su `curVar'

		local curMeanSD = 	`"`=trim("`: display %10.1f r(mean)'")'"' 	+   ///
							`" `=ustrunescape("\u00B1")' "' 			+   ///
							`"`=trim("`: display %10.1f r(sd)'")'"'		

        local curRow = `curRow' + 1
        putexcel A`curRow' = ("`curVar'")
        putexcel B`curRow' = ("`curMeanSD'")

    }

    // Loop through categorical vars

    foreach curVar of local catVars{

        local curRow = `curRow' + 1
        levelsof `curVar', local(curVarLevels)
        local curVarLevels `curVarLevels' . // Add "Missing" category

        foreach curLevel of local curVarLevels{
        
            count if `curVar' == `curLevel'
            local curCount = r(N)
            local curPercentage = `"`=trim("`: display %10.1f `=(100 * `r(N)'/`curTotalCount')''")'"'
            
            if `curLevel' != .  local curLabel : label (`curVar') `curLevel'
            else                local curLabel = "Missing"
                           
            local curRow = `curRow' + 1
            putexcel A`curRow' = ("`curVar': `curLabel'")
            putexcel B`curRow' = ("`curPercentage'% (`curCount')")
        
        }

    }

    restore
    putexcel close

}