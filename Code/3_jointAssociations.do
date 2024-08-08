version 17.0

frame copy dataset tempset
frame change tempset

drop if sin24_p>0.05 & cos24_p>0.05

drop *_se *_lb *_ub //*_p
rename *_hat *

gen sinCosinorAcro = sin(acrophase*2*_pi/24)
gen cosCosinorAcro = cos(acrophase*2*_pi/24)
gen maxCosinorPaee = max*60/1000
gen minCosinorPaee = min*60/1000 


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
local curRow = 1
foreach curExposure in maxCosinorPaee minCosinorPaee{
    forvalues i = 1/2{


        putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
        putexcel A`curRow' = ("outcomeVar")
        putexcel B`curRow' = ("Estimate: Women")
        putexcel C`curRow' = ("Count: Women")
        putexcel D`curRow' = ("Estimate: Men")
        putexcel E`curRow' = ("Count: Men")
        
    }
}

foreach curOutcomeVar of local outcomeVars{

    local curRow = `curRow'+1

    forvalues i = 1/2{

        #delimit ;
        glm `curOutcomeVar'
            c.maxCosinorPaee#i.sex
            c.minCosinorPaee#i.sex
            c.sinCosinorAcro#i.sex
            c.cosCosinorAcro#i.sex
            `modelLevel`i''
            if
            fatMass != .
            ,
            family(gaussian) 
            link(log)
            vce(robust)
            ;
        #delimit cr
        estimates store m1

        foreach curExposure in maxCosinorPaee minCosinorPaee{

            putexcel set "Results/3_jointAssociations.xlsx", sheet("`curExposure'_m`i'") modify
            putexcel A`curRow' = ("`curOutcomeVar'")
            
            margins, over(sex) eydx(`curExposure') atmeans vce(unconditional) post

            forvalues curSex = 0/1{

                nlcom 100*(exp(_b[`curExposure':`curSex'.sex])-1)

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
    
                local curCol = char(66+`curSex'*2)
                putexcel `curCol'`curRow' = ("`curEstimateCI'`curSigSymbol'")

                count if e(sample) == 1 & sex == `curSex'
                local curCol = char(67+`curSex'*2)
                putexcel `curCol'`curRow' = ("`r(N)'")

            }
            estimates restore m1
        }

        forvalues curSex = 0/1{
            
            su sinCosinorAcro if sex == `curSex'
            local curSinMean = r(mean)
            su cosCosinorAcro if sex == `curSex'
            local curCosMean = r(mean)

            #delimit ;
            margins,    at(sinCosinorAcro=(`=sin(0*2*_pi/24)')   cosCosinorAcro=(`=cos(0*2*_pi/24)')     sex==`curSex') 
                        at(sinCosinorAcro=(`=sin(1*2*_pi/24)')   cosCosinorAcro=(`=cos(1*2*_pi/24)')     sex==`curSex')           
                        at(sinCosinorAcro=(`=sin(2*2*_pi/24)')   cosCosinorAcro=(`=cos(2*2*_pi/24)')     sex==`curSex')           
                        at(sinCosinorAcro=(`=sin(3*2*_pi/24)')   cosCosinorAcro=(`=cos(3*2*_pi/24)')     sex==`curSex')         
                        
                        at(sinCosinorAcro=(`=sin(4*2*_pi/24)')   cosCosinorAcro=(`=cos(4*2*_pi/24)')     sex==`curSex')
                        at(sinCosinorAcro=(`=sin(5*2*_pi/24)')   cosCosinorAcro=(`=cos(5*2*_pi/24)')     sex==`curSex')            
                        at(sinCosinorAcro=(`=sin(6*2*_pi/24)')   cosCosinorAcro=(`=cos(6*2*_pi/24)')     sex==`curSex')            
                        at(sinCosinorAcro=(`=sin(7*2*_pi/24)')   cosCosinorAcro=(`=cos(7*2*_pi/24)')     sex==`curSex')            
                        
                        at(sinCosinorAcro=(`=sin(8*2*_pi/24)')   cosCosinorAcro=(`=cos(8*2*_pi/24)')     sex==`curSex')
                        at(sinCosinorAcro=(`=sin(9*2*_pi/24)')   cosCosinorAcro=(`=cos(9*2*_pi/24)')     sex==`curSex')             
                        at(sinCosinorAcro=(`=sin(10*2*_pi/24)')  cosCosinorAcro=(`=cos(10*2*_pi/24)')    sex==`curSex')          
                        at(sinCosinorAcro=(`=sin(11*2*_pi/24)')  cosCosinorAcro=(`=cos(11*2*_pi/24)')    sex==`curSex')            
                       
                        at(sinCosinorAcro=(`=sin(12*2*_pi/24)')  cosCosinorAcro=(`=cos(12*2*_pi/24)')    sex==`curSex')
                        at(sinCosinorAcro=(`=sin(13*2*_pi/24)')  cosCosinorAcro=(`=cos(13*2*_pi/24)')    sex==`curSex')          
                        at(sinCosinorAcro=(`=sin(14*2*_pi/24)')  cosCosinorAcro=(`=cos(14*2*_pi/24)')    sex==`curSex')            
                        at(sinCosinorAcro=(`=sin(15*2*_pi/24)')  cosCosinorAcro=(`=cos(15*2*_pi/24)')    sex==`curSex')            
                        
                        at(sinCosinorAcro=(`=sin(16*2*_pi/24)')  cosCosinorAcro=(`=cos(16*2*_pi/24)')    sex==`curSex')
                        at(sinCosinorAcro=(`=sin(17*2*_pi/24)')  cosCosinorAcro=(`=cos(17*2*_pi/24)')    sex==`curSex')            
                        at(sinCosinorAcro=(`=sin(18*2*_pi/24)')  cosCosinorAcro=(`=cos(18*2*_pi/24)')    sex==`curSex')           
                        at(sinCosinorAcro=(`=sin(19*2*_pi/24)')  cosCosinorAcro=(`=cos(19*2*_pi/24)')    sex==`curSex')           
                        
                        at(sinCosinorAcro=(`=sin(20*2*_pi/24)')  cosCosinorAcro=(`=cos(20*2*_pi/24)')    sex==`curSex')            
                        at(sinCosinorAcro=(`=sin(21*2*_pi/24)')  cosCosinorAcro=(`=cos(21*2*_pi/24)')    sex==`curSex')           
                        at(sinCosinorAcro=(`=sin(22*2*_pi/24)')  cosCosinorAcro=(`=cos(22*2*_pi/24)')    sex==`curSex')           
                        at(sinCosinorAcro=(`=sin(23*2*_pi/24)')  cosCosinorAcro=(`=cos(23*2*_pi/24)')    sex==`curSex')            

                        at(sinCosinorAcro=(`curSinMean')         cosCosinorAcro=(`curCosMean')           sex==`curSex')
                        atmeans
                        post
                        ;
            #delimit cr

            forvalues curHour = 1/24{
               
               nlcom _b[`curHour'._at]- _b[25._at] 
               
                //nlcom  100*(exp((_b[`curHour'._at] - _b[25._at])/_b[25._at])-1)

                //if `curHour' == 24  nlcom 100*(exp((_b[1._at]-_b[24._at])/_b[24._at])-1)
                //else                nlcom 100*(exp((_b[`=`curHour'+1'._at]-_b[`curHour'._at])/_b[`curHour'._at])-1)
            }

            

            asdf
        }        
    }
}

frame change dataset
frame drop tempset