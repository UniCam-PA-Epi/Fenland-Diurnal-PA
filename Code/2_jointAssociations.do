version 17.0

clear
set more off
set seed 1234

args filePath

use "`filePath'"

qui do Code/Sub/initialiseCovariates.do
qui do Code/Sub/initialiseOutcomes.do
qui do Code/Sub/initialisePAEE.do

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
                    c.smoking 
                    c.diet 
                    c.alcohol
                    ;

local catCovVars    i.ethnic
                    i.testsite 
                    i.season 
                    i.education 
                    i.income 
                    i.work_s
                    i.cardiometabol_med
                    ;
#delimit cr

local modelLevel1 `contCovVars' `catCovVars'
local modelLevel2 `contCovVars' `catCovVars' fatMass

foreach curOutcomeVar of local outcomeVars{
    forvalues i = 1/2{

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

        margins, over(sex) eyex(paeeTt) vce(unconditional) grand post
        lincom _b[paeeTt:0.sex] - _b[paeeTt:1.sex]    
        asdf

    }
}

