version 17.0

frame copy dataset tempset
frame change tempset

********************************************************************************
** Descriptives Tables (Table 1)
********************************************************************************

// Define continuous and categorical variables

// Convert J/min/kg to kJ/hour/kg
replace totalPAEE_hat   = totalPAEE_hat   * 60/1000

#delimit ;
local contVars  age
                bmi
                acrophase24_hat
                acrophase12_hat
                acrophase8_hat
                amplitude24_hat
                amplitude12_hat
                amplitude8_hat
                mesor_hat
                totalPAEE_hat
                glucose120
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

local catVars   sex
                agecat
                bmicat
                ethnic
                alcohol
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

capture mkdir Tables
capture erase "Tables/6_descriptiveTables.xlsx"

gen pooled = 0

foreach curGroup in pooled sex kGroup{

    levelsof `curGroup', local(groupLevels)

    foreach curGroupLevel of local groupLevels{
    
        local curLabel `curGroup'`curGroupLevel'
        
        putexcel set "Tables/6_descriptiveTables.xlsx", sheet("`curLabel'") modify

        preserve
        keep if `curGroup' == `curGroupLevel'
        local curRow = 0

        // Get total count for current group level

        count
        local curTotalCount = r(N)

        local curRow = `curRow' + 1
        putexcel A`curRow' = ("N")
        putexcel B`curRow' = `curTotalCount'

        // Loop through continuous vars

        foreach curVar of local contVars{

            if regexm("`curVar'","acrophase")==1{

                local p = regexr(regexr("`curVar'","acrophase",""),"_hat","")

                gen acrophase`p'_sin  = sin(acrophase`p'_hat*2*_pi/`p')
                gen acrophase`p'_cos  = cos(acrophase`p'_hat*2*_pi/`p')

                su acrophase`p'_sin, detail
                local sin_p25 = r(p25)
                local sin_p50 = r(p50)
                local sin_p75 = r(p75)

                su acrophase`p'_cos, detail
                local cos_p25 = r(p25)
                local cos_p50 = r(p50)
                local cos_p75 = r(p75)

                local p50   =  cond((`sin_p50')<0,`p',0) + atan2((`sin_p50'),(`cos_p50'))*`p'/(2*_pi)

                local pVal1 = cond((`sin_p75')<0,`p',0) + atan2((`sin_p75'),(`cos_p25'))*`p'/(2*_pi)  
                local pVal2 = cond((`sin_p25')<0,`p',0) + atan2((`sin_p25'),(`cos_p75'))*`p'/(2*_pi)

                local p25 = cond(`pVal1'<`p50',`pVal1',`pVal2')
                local p75 = cond(`pVal1'>`p50',`pVal1',`pVal2')

                drop acrophase`p'_sin acrophase`p'_cos

            }

            else{
                
                su `curVar', detail
                local p25 = r(p25)
                local p50 = r(p50)
                local p75 = r(p75)

            }

            if regexm("`curVar'","acrophase|amplitude|mesor")==1 local sig = 2
            else local sig = 1

            local curMedIQR =   `"`=trim("`: display %10.`sig'f `p50''")'"'     +   ///
                                " ("                                            +   ///
                                `"`=trim("`: display %10.`sig'f `p25''")'"'     +   ///
                                "-"                                             +   ///
                                `"`=trim("`: display %10.`sig'f `p75''")'"'     +   ///
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
}

frame change dataset
frame drop tempset
