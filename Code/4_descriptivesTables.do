version 17.0

frame copy dataset tempset
frame change tempset

********************************************************************************
** Descriptives Tables (Table 1)
********************************************************************************

// Define continuous and categorical variables

#delimit ;
local contVars  age
                fatMass
                fatFreeMass
                insulin
                leptin
                nefa
                ldl
                adiponectin
                crp
                mbpdia
                mbpsys
                glucose0
                glucose120
                ;

local catVars   alcohol
                ethnic
                education
                income
                smoke
                work_status
                marital_status
                season
                cardiometabol_med
                testsite
                ;

#delimit cr

capture erase "Results/2_descriptiveTables.xlsx"

levelsof sex, local(sexLevels)
foreach curSex of local sexLevels{
    
    local curSexLabel : label (sex) `curSex'
    putexcel set "Results/2_descriptiveTables.xlsx", sheet("`curSexLabel'") modify

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

        su `curVar', detail

        local curMedIQR =   `"`=trim("`: display %10.1f r(p50)'")'"'    +   ///
                            " ("                                        +   ///
                            `"`=trim("`: display %10.1f r(p25)'")'"'    +   ///
                            "-"                                         +   ///
                            `"`=trim("`: display %10.1f r(p75)'")'"'    +   ///
                            ")"

        local curRow = `curRow' + 1
        putexcel A`curRow' = ("`curVar'")
        putexcel B`curRow' = ("`curMedIQR'")

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

frame change dataset
frame drop tempset
