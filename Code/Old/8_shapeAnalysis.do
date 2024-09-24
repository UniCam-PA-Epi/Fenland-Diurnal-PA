version 17.0

mat drop _all




cap mata: mata drop conditionalShape()
mata
    void conditionalShape(real scalar applyShift, real scalar numRandomDraws)
    {   

        if (applyShift :== 0) shift = "unshifted"
        if (applyShift :== 1) shift = "shifted"
        
        Mu_out      = st_matrix("age")                      ,
                      st_matrix("bmi")                      ,
                      st_matrix("rhr")                      ,
                      st_matrix("totalPAEE")                ,
                      st_matrix("acrophase24_"+shift):*24   ,
                      st_matrix("acrophase12_"+shift):*12   ,
                      st_matrix("acrophase8_"+shift):*8     ,
                      st_matrix("amplitude24")              ,
                      st_matrix("amplitude12")              ,
                      st_matrix("amplitude8")               ,
                      st_matrix("mesor")
        
        Sigma_out   = st_matrix("covMatrix")

        X_1                 = J(1,0,.)
        conditioningIndices = J(1,0,.)
        dependentIndices    = J(1,0,.)

        for (i=1; i<=cols(Mu_out); i++){
            if (Mu_out[i] :!= .){
                X_1 = X_1 , Mu_out[i]
                conditioningIndices = conditioningIndices, i
            }
            else{
                dependentIndices = dependentIndices, i
            }
        }
   
        
        if (cols(dependentIndices):!=0){

            meanVector  = quadcross(st_matrix("catVector")',st_matrix("meanMatrix")) 
        
            Mu_1 = meanVector[conditioningIndices]
            Mu_2 = meanVector[dependentIndices]

            Sigma_11 = Sigma_out[conditioningIndices,conditioningIndices]
            Sigma_12 = Sigma_out[conditioningIndices,dependentIndices]
            Sigma_22 = Sigma_out[dependentIndices,dependentIndices]

            solved      = cholsolve(Sigma_11,Sigma_12)

            Mu_hat      = Mu_2 :+ quadcross((X_1:-Mu_1)',solved)
            Sigma_hat   = Sigma_22 - quadcross(Sigma_12,solved)

            Mu_out[dependentIndices] = Mu_hat
            //Sigma_out[dependentIndices,dependentIndices] = Sigma_hat
        }
        
        if (numRandomDraws:>0){

            standardNormalRandoms       = rnormal(cols(Mu_out),numRandomDraws,0,0.5)
            correlatedNormalRandoms     = quadcross(standardNormalRandoms,cholesky(Sigma_out))
            randomDraws                 = correlatedNormalRandoms + J(numRandomDraws, 1, Mu_out)

            randomDraws[.,5] = mod(randomDraws[.,5],24)
            randomDraws[.,6] = mod(randomDraws[.,6],12)
            randomDraws[.,7] = mod(randomDraws[.,7],8)
            randomDraws = randomDraws[.,5::11]
        }

        st_matrix("Mu_out", Mu_out[5::11])
        st_matrix("randomDraws", randomDraws)

    }
end






*************************
** Initilaise analysis **
*************************

// Exclude participants for whom the cosinor model did not fit, those with insufficient wear time, and those without fat mass measurement.

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05
drop if P_Pwear_consolidated<72
drop if fatMass == .

drop *_se *_lb *_ub *_p
rename *_hat *

// Initialise analysis vars by deconstructing time-based features into their sin and cos axes

gen maxHour24_sin   = sin(maxHour*2*_pi/24)
gen maxHour24_cos   = cos(maxHour*2*_pi/24)

gen acrophase24_sin = sin(acrophase24*2*_pi/24)
gen acrophase24_cos = cos(acrophase24*2*_pi/24)

gen acrophase12_sin = sin(acrophase12*2*_pi/12)
gen acrophase12_cos = cos(acrophase12*2*_pi/12)

gen acrophase8_sin  = sin(acrophase8*2*_pi/8)
gen acrophase8_cos  = cos(acrophase8*2*_pi/8)


#delimit ;
local featureVec age
                 bmi
                 rhr
                 totalPAEE
                 amplitude24
                 amplitude12
                 amplitude8
                 mesor
                 ;
   
local acroVec   acrophase24
                acrophase12
                acrophase8
                ;
#delimit cr







forvalues curSex = 0/1{

    #delimit ;
    local contVars  age
                    bmi
                    rhr
                    totalPAEE
                    acrophase24
                    acrophase12
                    acrophase8
                    amplitude24
                    amplitude12 
                    amplitude8
                    mesor
                    ;

    local catVars   i.work_status
                    i.marital_status
                    i.education
                    i.income
                    ;
    #delimit cr






    qui mi set mlong
    qui mi register imputed `contVars'
    mi impute mvn `contVars' = `catVars' if sex == `curSex', emonly(iterate(1000)) rseed(1234)
    
    mat covMatrix  = r(Sigma_em)
    mat meanMatrix = r(Beta_em)
    mi unregister `contVars'






    capture mat drop catVector
    foreach curCat of local catVars{
        
        local firstLevel = 1
        local curCat = regexr("`curCat'","i.","")
        levelsof `curCat', local(levelsList)

        foreach curLevel of local levelsList{

            if `firstLevel' == 1{

                capture mat catVector = catVector , 0
                if _rc != 0 mat catVector = 0

                count if `curCat' != . & sex == `curSex'
                local countTotal = r(N)
                local firstLevel = 0
            }
            else{
                count if `curCat' == `curLevel' & sex == `curSex'
                mat catVector = catVector, r(N)/`countTotal'
            }
        }
    }
    mat catVector = catVector, 1




    foreach curVar of local contVars{
        mat `curVar' = .
    }

    foreach curVar of local featureVec{
        qui su `curVar' if sex == `curSex'
        mat `curVar' = r(mean)
    }



    mat acrophase24_unshifted = .
    mat acrophase12_unshifted = .
    mat acrophase8_unshifted  = .

    mat acrophase24_shifted = .
    mat acrophase12_shifted = .
    mat acrophase8_shifted  = .

    qui foreach curAcro of local acroVec{
    
        estimates use Models/insulin_m5
        estimates esample:
        margins if sex==`curSex', eydx(`curAcro'_sin `curAcro'_cos) asobserved post

        noisi nlcom (cond(_b[`curAcro'_sin]<0,2*_pi,0)+atan2(_b[`curAcro'_sin], _b[`curAcro'_cos]))*1/(2*_pi)


        mat `curAcro'_unshifted = r(b)[1,1]

        if regexm("`curAcro'","24") == 1 local shiftVal = 5/24      
        if regexm("`curAcro'","12") == 1 local shiftVal =  5/12   
        if regexm("`curAcro'","8")  == 1 local shiftVal =  5/8       

        mat `curAcro'_shifted   = mod(r(b)[1,1] + `shiftVal',1)

        /*

        nlcom sqrt((_b[`curAcro'_sin])^2+(_b[`curAcro'_cos])^2)
        local curP = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

        if `curP' < 0.05{

            nlcom (cond(_b[`curAcro'_sin]<0,2*_pi,0)+atan2(_b[`curAcro'_sin], _b[`curAcro'_cos]))*1/(2*_pi)

            //mat `curAcro'_unshifted = r(b)[1,1]
            //mat `curAcro'_shifted   = mod(r(b)[1,1] + cond(r(b)[1,1]<0.5, 0.25, -0.25),1)

            mat `curAcro'_unshifted = mod(r(b)[1,1] + cond(r(b)[1,1]<0.5, 0.075, -0.075),1)
            mat `curAcro'_shifted   = mod(r(b)[1,1] + cond(r(b)[1,1]<0.5, 0.325, -0.325),1)

        }

        */
    }


    forvalues i = 0/1{
        
        mata: conditionalShape(`i',128)

        forvalues j = 0/128{
            
            if `j' == 0{
                local f24 Mu_out[1,4]*cos((2*_pi/24)*(x-Mu_out[1,1]))
                local f12 Mu_out[1,5]*cos((2*_pi/12)*(x-Mu_out[1,2]))
                local f8  Mu_out[1,6]*cos((2*_pi/8)*(x-Mu_out[1,3]))
                local mesor Mu_out[1,7]

                local plotCommand (func y=exp(`f24'+`f12'+`f8'+`mesor'), range(0 24) lcolor(maroon) lwidth(0.5)) 
            }

            else{
                local f24 randomDraws[`j',4]*cos((2*_pi/24)*(x-randomDraws[`j',1]))
                local f12 randomDraws[`j',5]*cos((2*_pi/12)*(x-randomDraws[`j',2]))
                local f8  randomDraws[`j',6]*cos((2*_pi/8)*(x-randomDraws[`j',3]))
                local mesor randomDraws[`j',7]

                local plotCommand (func y=exp(`f24'+`f12'+`f8'+`mesor'), range(0 24) lcolor(navy%15) lwidth(0.25)) `plotCommand' 
            }

        }

        set graphics on

        twoway `plotCommand' (, legend(off) name(s`curSex'_i`i', replace))

        //twoway (func y=exp(`f24'+`f12'+`f8'+`mesor'),range(0 24)) //name(s`curSex'_i`i', replace))
        
    }
    
}

set graphics on

graph combine s0_i0 s0_i1 s1_i0 s1_i1, xcommon ycommon

