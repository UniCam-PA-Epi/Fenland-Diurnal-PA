version 17.0

frame copy dataset tempset
frame change tempset

drop if P_Pwear_consolidated<72
drop if sin24_p>0.05 & cos24_p>0.05
drop if fatMass == .

drop *_se *_lb *_ub //*_p
rename *_hat *

gen sinCosinorAcro = sin(acrophase*2*_pi/24)
gen cosCosinorAcro = cos(acrophase*2*_pi/24)
gen maxCosinorPaee = max*60/1000
gen minCosinorPaee = min*60/1000 

forvalues curSex = 0/1{

    egen fatFreeMass_std_s`curSex' = std((fatFreeMass)) if sex == `curSex'
    egen insulin_std_s`curSex'     = std((insulin)) if sex == `curSex'
    egen nefa_std_s`curSex'        = std((nefa)) if sex == `curSex'
    egen adiponectin_std_s`curSex'= std((adiponectin)) if sex == `curSex'
    egen leptin_std_s`curSex'      = std((leptin)) if sex == `curSex'
    egen crp_std_s`curSex'         = std((crp)) if sex == `curSex'
    egen mbpsys_std_s`curSex'      = std((mbpsys)) if sex == `curSex'
    egen glucose120_std_s`curSex' = std((glucose120)) if sex == `curSex'

}

egen fatFreeMass_std   = rowtotal(fatFreeMass_std_*)    , missing
egen insulin_std       = rowtotal(insulin_std_*)    , missing
egen nefa_std          = rowtotal(nefa_std_*)    , missing
egen adiponectin_std   = rowtotal(adiponectin_std_*)    , missing
egen leptin_std        = rowtotal(leptin_std_*)    , missing
egen crp_std           = rowtotal(crp_std_*)    , missing
egen mbpsys_std        = rowtotal(mbpsys_std_*)    , missing
egen glucose120_std    = rowtotal(glucose120_std_*)    , missing

replace insulin_std     = -insulin_std
replace nefa_std        = -nefa_std 
replace leptin_std      = -leptin_std 
replace crp_std         = -crp_std 
replace mbpsys_std      = -mbpsys_std 
replace glucose120_std  = -glucose120_std 

egen notmiss = rowmiss( insulin_std nefa_std mbpsys_std leptin_std crp_std glucose120_std)
egen CCMR = rowmean( insulin_std nefa_std mbpsys_std leptin_std crp_std glucose120_std) if notmiss ==0

drop notmiss *_std*

#delimit ;
local outcomeVars   //CCMR
                    fatFreeMass
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
                    glucose120
                    ;

local contCovVars   c.age
                    c.diet 
                    c.alcohol
                    ;

local catCovVars    i.sex
                    i.ethnic
                    i.education
                    i.income
                    i.smoke
                    i.work_status
                    i.marital_status
                    i.season
                    i.cardiometabol_med
                    i.testsite
                    ;

#delimit cr

local modelLevel1 `contCovVars' `catCovVars'
local modelLevel2 `contCovVars' `catCovVars' c.fatMass

capture erase "Results/3_jointAssociations.xlsx"

forvalues i = 1/2{

    putexcel set "Results/3_jointAssociations.xlsx", sheet("waldBlockTest_m`i'") modify
    putexcel A1 = ("outcomeVar")
    putexcel B1 = ("Pvalue: CosinorPaee block")
    putexcel C1 = ("Pvalue: CosinorPaee#sex block")
    putexcel D1 = ("Pvalue: CosinorAcro block")
    putexcel E1 = ("Pvalue: CosinorAcro#sex block")
    putexcel F1 = ("Count: Women")
    putexcel G1 = ("Count: Men")

    foreach curExposure in maxCosinorPaee minCosinorPaee{

        putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
        putexcel A1 = ("outcomeVar")
        putexcel B1 = ("Estimate: Women")
        putexcel C1 = ("Estimate: Men")
    }
}

forvalues i = 1/2{

    putexcel set "Results/3_jointAssociations.xlsx", sheet("CosinorAcro_m`i'") modify
    putexcel A1 = ("outcomeVar")
    putexcel B1 = ("Estimate: Women")
    putexcel C1 = ("Hour at max: Women")
    putexcel D1 = ("Hour at min: Women")
    putexcel E1 = ("Midpoint value: Women")
    putexcel F1 = ("Estimate: Men")
    putexcel G1 = ("Hour at max: Men")
    putexcel H1 = ("Hour at min: Men")
    putexcel I1 = ("Midpoint value: Men")
}
   
    
local curRow = 1
foreach curOutcomeVar of local outcomeVars{

    local curRow = `curRow'+1

    forvalues i = 1/2{

        #delimit ;
        nestreg: glm `curOutcomeVar' 
                     (`modelLevel`i'') 
                     (c.maxCosinorPaee       c.minCosinorPaee)
                     (c.maxCosinorPaee#i.sex c.minCosinorPaee#i.sex)
                     (c.sinCosinorAcro       c.cosCosinorAcro)
                     (c.sinCosinorAcro#i.sex c.cosCosinorAcro#i.sex)
                     ,
                     family(gaussian)
                     link(`=cond("`curOutcomeVar'"=="CCMR","identity","log")')
                     ;
        #delimit cr
        
        estimates store fullModel

        putexcel set "Results/3_jointAssociations.xlsx", sheet("waldBlockTest_m`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")
        putexcel B`curRow' = (r(wald)[2,3])
        putexcel C`curRow' = (r(wald)[3,3])
        putexcel D`curRow' = (r(wald)[4,3])
        putexcel E`curRow' = (r(wald)[5,3])

        count if e(sample) == 1 & sex == 0
        putexcel F`curRow' = ("`r(N)'")

        count if e(sample) == 1 & sex == 1
        putexcel G`curRow' = ("`r(N)'")
        
        foreach curExposure in maxCosinorPaee minCosinorPaee{

            putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
            putexcel A`curRow' = ("`curOutcomeVar'")
            
            estimates restore fullModel
            if "`curOutcomeVar'" == "CCMR"  margins, over(sex) dydx(`curExposure') asobserved post
            else                            margins, over(sex) eydx(`curExposure') asobserved post

            forvalues curSex = 0/1{
                
                nlcom _b[`curExposure':`curSex'.sex]

                local curEstimate   = r(b)[1,1]
                local curLB         = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  
                local curUB         = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])
                local curPvalue     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1]))) 

                local curEstimateCI =   `"`=trim("`: display %10.3f `curEstimate''")'"' +   ///
                                        " ("                                            +   ///
                                        `"`=trim("`: display %10.3f `curLB''")'"'       +   ///
                                        ", "                                            +   ///
                                        `"`=trim("`: display %10.3f `curUB''")'"'       +   ///
                                        ")"

                local curSigSymbol = cond(`curPvalue'<0.01,"**",cond(`curPvalue'<0.05,"*",""))
    
                local curCol = char(66+`curSex')
                putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")
            }       
        }
        
        putexcel set "Results/3_jointAssociations.xlsx", sheet("CosinorAcro_m`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")

        forvalues curSex = 0/1{

            estimates restore fullModel

            if "`curOutcomeVar'" == "CCMR"  margins if sex == `curSex', dydx(sinCosinorAcro cosCosinorAcro) asobserved post
            else                            margins if sex == `curSex', eydx(sinCosinorAcro cosCosinorAcro) asobserved post

            nlcom sqrt((_b[sinCosinorAcro])^2+(_b[cosCosinorAcro])^2)

            local curEstimate   = r(b)[1,1]
            local curLB         = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  
            local curUB         = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])
            local curPvalue     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1]))) 

            local curEstimateCI =   `"`=trim("`: display %10.3f `curEstimate''")'"' +   ///
                                    " ("                                            +   ///
                                    `"`=trim("`: display %10.3f `curLB''")'"'       +   ///
                                    ", "                                            +   ///
                                    `"`=trim("`: display %10.3f `curUB''")'"'       +   ///
                                    ")"

            local curSigSymbol = cond(`curPvalue'<0.01,"**",cond(`curPvalue'<0.05,"*",""))

            local curCol = char(66+4*`curSex')
            putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")

            nlcom cond(_b[sinCosinorAcro]<0,24,0)+atan2(_b[sinCosinorAcro], _b[cosCosinorAcro])*24/(2*_pi)
            local acroMax = r(b)[1,1]
            local curCol = char(67+4*`curSex')
            putexcel `curCol'`curRow' = (`"`=trim("`: display %10.2f `acroMax''")'"')

            nlcom cond(_b[sinCosinorAcro]<0,12,12)+atan2(_b[sinCosinorAcro], _b[cosCosinorAcro])*24/(2*_pi)
            local acroMin = r(b)[1,1]
            local curCol = char(68+4*`curSex')
            putexcel `curCol'`curRow' = (`"`=trim("`: display %10.2f `acroMin''")'"')

            nlcom cond(_b[sinCosinorAcro]<0,18,6)+atan2(_b[sinCosinorAcro], _b[cosCosinorAcro])*24/(2*_pi)
            local acroMid = r(b)[1,1]
            
            estimates restore fullModel
            margins if sex == `curSex', at(sinCosinorAcro=`=sin(`acroMid'*2*_pi/24)' cosCosinorAcro=(`=cos(`acroMid'*2*_pi/24)')) asobserved post

            nlcom _b[_cons]
            local midPoint = r(b)[1,1]
            local curCol = char(69+4*`curSex')
            putexcel `curCol'`curRow' = (`"`=trim("`: display %10.2f `midPoint''")'"')

            
        }
    }
}

frame change dataset
frame drop tempset



/*

      asdf

        forvalues curSex = 0/1{
            
            su acrophase if sex == `curSex'
            local curSinMax = sin(r(mean)*2*_pi/24)
            local curCosMax = cos(r(mean)*2*_pi/24)

            local curSinMin = sin(r(mean)*2*_pi/24-_pi)
            local curCosMin = cos(r(mean)*2*_pi/24-_pi)

            estimates restore fullModel
            margins if sex == `curSex', eydx(sinCosinorAcro cosCosinorAcro) asobserved post

            local acroMax = 24+atan2(_b[sinCosinorAcro],_b[cosCosinorAcro])*24/(2*_pi)
            local acroMid = 18+atan2(_b[sinCosinorAcro],_b[cosCosinorAcro])*24/(2*_pi)
            local acroMin = 12+atan2(_b[sinCosinorAcro],_b[cosCosinorAcro])*24/(2*_pi)

            estimates restore fullModel
            #delimit ;
            margins if sex == `curSex'
                    ,  
                    at(sinCosinorAcro=(`=sin(0*2*_pi/24)')   cosCosinorAcro=(`=cos(0*2*_pi/24)')) 
                    at(sinCosinorAcro=(`=sin(1*2*_pi/24)')   cosCosinorAcro=(`=cos(1*2*_pi/24)'))           
                    at(sinCosinorAcro=(`=sin(2*2*_pi/24)')   cosCosinorAcro=(`=cos(2*2*_pi/24)'))           
                    at(sinCosinorAcro=(`=sin(3*2*_pi/24)')   cosCosinorAcro=(`=cos(3*2*_pi/24)'))         
                        
                    at(sinCosinorAcro=(`=sin(4*2*_pi/24)')   cosCosinorAcro=(`=cos(4*2*_pi/24)'))
                    at(sinCosinorAcro=(`=sin(5*2*_pi/24)')   cosCosinorAcro=(`=cos(5*2*_pi/24)'))            
                    at(sinCosinorAcro=(`=sin(6*2*_pi/24)')   cosCosinorAcro=(`=cos(6*2*_pi/24)'))            
                    at(sinCosinorAcro=(`=sin(7*2*_pi/24)')   cosCosinorAcro=(`=cos(7*2*_pi/24)'))            
                        
                    at(sinCosinorAcro=(`=sin(8*2*_pi/24)')   cosCosinorAcro=(`=cos(8*2*_pi/24)'))
                    at(sinCosinorAcro=(`=sin(9*2*_pi/24)')   cosCosinorAcro=(`=cos(9*2*_pi/24)'))             
                    at(sinCosinorAcro=(`=sin(10*2*_pi/24)')  cosCosinorAcro=(`=cos(10*2*_pi/24)'))          
                    at(sinCosinorAcro=(`=sin(11*2*_pi/24)')  cosCosinorAcro=(`=cos(11*2*_pi/24)'))            
                       
                    at(sinCosinorAcro=(`=sin(12*2*_pi/24)')  cosCosinorAcro=(`=cos(12*2*_pi/24)'))
                    at(sinCosinorAcro=(`=sin(13*2*_pi/24)')  cosCosinorAcro=(`=cos(13*2*_pi/24)'))          
                    at(sinCosinorAcro=(`=sin(14*2*_pi/24)')  cosCosinorAcro=(`=cos(14*2*_pi/24)'))            
                    at(sinCosinorAcro=(`=sin(15*2*_pi/24)')  cosCosinorAcro=(`=cos(15*2*_pi/24)'))            
                        
                    at(sinCosinorAcro=(`=sin(16*2*_pi/24)')  cosCosinorAcro=(`=cos(16*2*_pi/24)'))
                    at(sinCosinorAcro=(`=sin(17*2*_pi/24)')  cosCosinorAcro=(`=cos(17*2*_pi/24)'))            
                    at(sinCosinorAcro=(`=sin(18*2*_pi/24)')  cosCosinorAcro=(`=cos(18*2*_pi/24)'))           
                    at(sinCosinorAcro=(`=sin(19*2*_pi/24)')  cosCosinorAcro=(`=cos(19*2*_pi/24)'))           
                        
                    at(sinCosinorAcro=(`=sin(20*2*_pi/24)')  cosCosinorAcro=(`=cos(20*2*_pi/24)'))            
                    at(sinCosinorAcro=(`=sin(21*2*_pi/24)')  cosCosinorAcro=(`=cos(21*2*_pi/24)'))           
                    at(sinCosinorAcro=(`=sin(22*2*_pi/24)')  cosCosinorAcro=(`=cos(22*2*_pi/24)'))           
                    at(sinCosinorAcro=(`=sin(23*2*_pi/24)')  cosCosinorAcro=(`=cos(23*2*_pi/24)'))  

                    at(sinCosinorAcro=`=sin(`acroMax'*2*_pi/24)' cosCosinorAcro=(`=cos(`acroMax'*2*_pi/24)'))
                    at(sinCosinorAcro=`=sin(`acroMid'*2*_pi/24)' cosCosinorAcro=(`=cos(`acroMid'*2*_pi/24)'))
                    at(sinCosinorAcro=`=sin(`acroMin'*2*_pi/24)' cosCosinorAcro=(`=cos(`acroMin'*2*_pi/24)'))
                    //at(sinCosinorAcro=(`=sin(11*2*_pi/24)')  cosCosinorAcro=(`=cos(11*2*_pi/24)')) 
                    //at(sinCosinorAcro=(`curSinMean')         cosCosinorAcro=(`curCosMean'))
                    asobserved
                    post
                    ;
            #delimit cr

            nlcom (_b[25._at]- _b[27._at])/2
            nlcom (abs(_b[25._at]- _b[26._at])+abs(_b[27._at]- _b[26._at]))/2
            asdf
            forvalues curHour = 1/24{
               
               nlcom _b[`curHour'._at]- _b[25._at] 
               
                //nlcom  100*(exp((_b[`curHour'._at] - _b[25._at])/_b[25._at])-1)

                //if `curHour' == 24  nlcom 100*(exp((_b[1._at]-_b[24._at])/_b[24._at])-1)
                //else                nlcom 100*(exp((_b[`=`curHour'+1'._at]-_b[`curHour'._at])/_b[`curHour'._at])-1)
            }

            
            asdf
        }      

        */