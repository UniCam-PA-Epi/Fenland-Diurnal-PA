version 17.0

frame copy dataset tempset
frame change tempset

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05
drop if P_Pwear_consolidated<72
drop if fatMass == .

drop *_se *_lb *_ub *_p
rename *_hat *

gen maxHour24_sin = sin(maxHour*2*_pi/24)
gen maxHour24_cos = cos(maxHour*2*_pi/24)

gen acrophase24_sin = sin(acrophase24*2*_pi/24)
gen acrophase24_cos = cos(acrophase24*2*_pi/24)

gen acrophase12_sin = sin(acrophase12*2*_pi/12)
gen acrophase12_cos = cos(acrophase12*2*_pi/12)

gen acrophase8_sin = sin(acrophase8*2*_pi/8)
gen acrophase8_cos = cos(acrophase8*2*_pi/8)

replace mesor = exp(mesor) * 60/1000
replace totalPAEE = totalPAEE * 60/1000
replace maxValue = maxValue * 60/1000


#delimit ;
local outcomeVars   insulin
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
                    fatFreeMass
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

qui forvalues i = 1/2{

    putexcel set "Results/3_jointAssociations.xlsx", sheet("waldBlockTest_m`i'") modify
    putexcel A1 = ("outcomeVar")
    putexcel B1 = ("Pvalue: mesor")
    putexcel C1 = ("Pvalue: mesor#sex")
    putexcel D1 = ("Pvalue: maxValue")
    putexcel E1 = ("Pvalue: maxValue#sex")
    putexcel F1 = ("Pvalue: amplitude24")
    putexcel G1 = ("Pvalue: amplitude24#sex")
    putexcel H1 = ("Pvalue: amplitude12")
    putexcel I1 = ("Pvalue: amplitude12#sex")
    putexcel J1 = ("Pvalue: amplitude8")
    putexcel K1 = ("Pvalue: amplitude8#sex")
    putexcel L1 = ("Pvalue: maxHour24")
    putexcel M1 = ("Pvalue: maxHour24#sex")
    putexcel N1 = ("Pvalue: acrophase24")
    putexcel O1 = ("Pvalue: acrophase24#sex")
    putexcel P1 = ("Pvalue: acrophase12")
    putexcel Q1 = ("Pvalue: acrophase12#sex")
    putexcel R1 = ("Pvalue: acrophase8")
    putexcel S1 = ("Pvalue: acrophase8#sex")
    putexcel T1 = ("Count: Women")
    putexcel U1 = ("Count: Men")

    foreach curExposure in mesor maxValue amplitude24 amplitude12 amplitude8{

        putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
        putexcel A1 = ("outcomeVar")
        putexcel B1 = ("Estimate: Women")
        putexcel C1 = ("Estimate: Men")
    }

    foreach curExposure in maxHour24 acrophase24 acrophase12 acrophase8{

        putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
        putexcel A1 = ("outcomeVar")
        putexcel B1 = ("Amplitude: Women")
        putexcel C1 = ("Normalised acrophase: Women")
        putexcel D1 = ("Amplitude: Men")
        putexcel E1 = ("Normalised acrophase: Men")

    }
}
    
local curRow = 1
foreach curOutcomeVar of local outcomeVars{

    local curRow = `curRow'+1

    forvalues i = 1/2{

        #delimit ;
        nestreg : glm    `curOutcomeVar' 
                            (`modelLevel`i'')

                            (c.mesor)
                            (c.mesor#i.sex)

                            (c.maxValue)
                            (c.maxValue#i.sex)
                            
                            (c.amplitude24)
                            (c.amplitude24#i.sex)

                            (c.amplitude12)
                            (c.amplitude12#i.sex)

                            (c.amplitude8)
                            (c.amplitude8#i.sex)

                            (c.maxHour24_sin        c.maxHour24_cos)
                            (c.maxHour24_sin#i.sex  c.maxHour24_cos#i.sex)

                            (c.acrophase24_sin        c.acrophase24_cos)
                            (c.acrophase24_sin#i.sex  c.acrophase24_cos#i.sex)

                            (c.acrophase12_sin        c.acrophase12_cos)
                            (c.acrophase12_sin#i.sex  c.acrophase12_cos#i.sex)

                            (c.acrophase8_sin          c.acrophase8_cos)
                            (c.acrophase8_sin#i.sex    c.acrophase8_cos#i.sex)                      
                            ,
                            family(gaussian)
                            link(log)
                            ;
        #delimit cr
        
        estimates store fullModel
        
        putexcel set "Results/3_jointAssociations.xlsx", sheet("waldBlockTest_m`i'") modify
        putexcel A`curRow' = ("`curOutcomeVar'")

        forvalues j = 1/18{
            local curCol = char(65+`j')
            putexcel `curCol'`curRow' = (r(wald)[`=`j'+1',3])

        }
        count if e(sample) == 1 & sex == 0
        putexcel T`curRow' = ("`r(N)'")

        count if e(sample) == 1 & sex == 1
        putexcel U`curRow' = ("`r(N)'")

        
        foreach curExposure in mesor maxValue amplitude24 amplitude12 amplitude8{

            putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
            putexcel A`curRow' = ("`curOutcomeVar'")
            
            estimates restore fullModel
            margins, over(sex) eydx(`curExposure') asobserved post

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

        foreach curExposure in maxHour24 acrophase24 acrophase12 acrophase8{

            putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
            putexcel A`curRow' = ("`curOutcomeVar'")

            estimates restore fullModel
            margins, over(sex) eydx(`curExposure'_sin `curExposure'_cos) asobserved post

            forvalues curSex = 0/1{

                nlcom sqrt((_b[`curExposure'_sin:`curSex'.sex])^2+(_b[`curExposure'_cos:`curSex'.sex])^2)

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

                local curCol = char(66+2*`curSex')
                putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")


                nlcom cond(_b[`curExposure'_sin:`curSex'.sex]<0,1,0)+atan2(_b[`curExposure'_sin:`curSex'.sex], _b[`curExposure'_cos:`curSex'.sex])*1/(2*_pi)

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

                local curCol = char(67+2*`curSex')
                putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")


            }
        }

        gen obs = _n/64
        gen hat = .
        gen lb = .
        gen ub = .

        local marginsList
        foreach curExposure in maxHour24 acrophase24 acrophase12 acrophase8{

            forvalues j = 1/64{
                
                local marginsList `marginsList' at(`curExposure'_sin=(`=sin(`j'*2*_pi/64)')`curExposure'_cos=(`=cos(`j'*2*_pi/64)'))

            }

            forvalues curSex = 0/1{

                estimates restore fullModel
                margins if sex==`curSex', at()  `marginsList' asobserved post

                forvalues j = 1/64{

                    nlcom (_b[`=`j'+1'._at]-_b[1._at])/_b[1._at]
                    replace hat = r(b)[1,1]                                     in `j'
                    replace lb  = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])  in `j'
                    replace ub  = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])  in `j'

                }
                
                #delimit ;
                twoway  (
                        rarea lb ub obs if hat !=.
                        ,
                        fcolor(navy%40)
                        lcolor(navy%0)
                        )
                        (
					    line hat obs  if hat !=.
					    ,
					    lcolor(navy)
					    )
                        (
                        func y=0
                        ,
					    lcolor(gs8)
					    lpattern(dash)
					    lwidth(0.2)
					    range(0 1)
					    )
                        (
                        ,
 					    graphregion(color(white))
					    legend(off)                       
                        )
                        ;                      
                #delimit cr
         
                asdf
            }
            
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