version 17.0

frame copy dataset tempset
frame change tempset

#delimit ;
local outcomeVars   fatFreeMass
                    insulin
                    mbpsys
                    mbpdia
                    nefa
                    leptin
                    adiponectin
                    crp
                    ldl
                    hdl
                    glucose0
                    ;

local contCovVars   c.log_rel_amplitude
                    c.sin_acro
                    c.cos_acro
                    c.age
                    c.diet 
                    c.alcohol
                    ;

local catCovVars    i.ethnic
                    i.education
                    i.income
                    i.smoke
                    i.work_s
                    i.marital_s
                    i.season
                    i.cardiometabol_med
                    i.testsite
                    ;

#delimit cr

local modelLevel1 `contCovVars' `catCovVars'
local modelLevel2 `contCovVars' `catCovVars' c.fatMass

capture erase "Results/2_jointAssociations.xlsx"
local curRow = 1

forvalues i = 1/2{

    putexcel set "Results/2_jointAssociations.xlsx", sheet("modelLevel`i'") modify
    putexcel A`curRow' = ("outcomeVar")
    putexcel B`curRow' = ("Count")
    putexcel C`curRow' = ("Women")
    putexcel D`curRow' = ("Men")
    putexcel E`curRow' = ("Diff")
    
}


foreach curOutcomeVar of local outcomeVars{

    local curRow = `curRow'+1

    forvalues i = 1/2{

        putexcel set "Results/2_jointAssociations.xlsx", sheet("modelLevel`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")

        #delimit ;
        glm `curOutcomeVar' 
            c.paeeTt##i.sex
            `modelLevel`i''
            if
            fatMass != .
            ,
            family(gaussian) 
            link(log)
            vce(robust)
            ;
        #delimit cr

        count if e(sample) == 1 & sex == 0
        local womenCount = r(N)
        count if e(sample) == 1 & sex == 1
        local menCount = r(N)
        putexcel B`curRow' = ("`womenCount', `menCount'")
        
        margins, over(sex) eyex(paeeTt) atmeans post vce(unconditional)

        local coefTest1 _b[paeeTt:0.sex]
        local coefTest2 _b[paeeTt:1.sex]
        local coefTest3 _b[paeeTt:0.sex]-_b[paeeTt:1.sex]

        forvalues j = 1/3{
            
            lincom `coefTest`j''

            local curCoefCI =   `"`=trim("`: display %10.3f r(estimate)'")'"'   +   ///
                                " ("                                            +   ///
                                `"`=trim("`: display %10.3f r(lb)'")'"'         +   ///
                                ", "                                            +   ///
                                `"`=trim("`: display %10.3f r(ub)'")'"'         +   ///
                                ")"

            local curSig    = cond(r(p)<0.01,"**",cond(r(p)<0.05,"*",""))
            local curCol = char(66+`j')

            putexcel `curCol'`curRow' = ("`curCoefCI'`curSig'")

        }      
    }
}

frame change dataset
frame drop tempset
