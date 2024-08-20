version 17.0


cap mata: mata drop conditionalShape()
mata
    void conditionalShape(  real matrix covMatrix   , 
                            real matrix meanMatrix  , 
                            real vector catVector   ,
                            real scalar age         ,
                            real scalar bmi         ,
                            real scalar rhr         ,
                            real scalar totalPAEE   ,
                            real scalar acrophase24 ,
                            real scalar acrophase12 ,
                            real scalar acrophase8                   
                        )
    {   
        
        meanVector = quadcross(catVector',meanMatrix)

       // X_1 = age , bmi, rhr, totalPAEE, 24*(acrophase24) , 12*(acrophase12), 8*(acrophase8)
        
        X_1 = age , bmi, rhr, totalPAEE, 12*(acrophase12), 8*(acrophase8), 2.8

        conditioningIndices = 1,2,3,4,6,7,11
        dependentIndices    = 5,8,9,10
    
        Mu_1 = meanVector[conditioningIndices]
        Mu_2 = meanVector[dependentIndices]

        Sigma_11 = covMatrix[conditioningIndices,conditioningIndices]
        Sigma_12 = covMatrix[conditioningIndices,dependentIndices]
        Sigma_22 = covMatrix[dependentIndices,dependentIndices]

        solved      = cholsolve(Sigma_11,Sigma_12)
        Mu_hat      = Mu_2 :+ quadcross((X_1:-Mu_1)',solved)
        Sigma_hat   = Sigma_22 - quadcross(Sigma_12,solved)
        
        st_matrix("Mu_hat", Mu_hat)
        Mu_hat
        standardNormalRandoms       = rnormal(cols(Mu_hat),5,0,1)
        correlatedNormalRandoms     = quadcross(standardNormalRandoms,cholesky(Sigma_hat))
        multivariateNormalRandoms   = correlatedNormalRandoms + J(5, 1, Mu_hat) 

        st_matrix("randomDraws", multivariateNormalRandoms)

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
                    ;
    #delimit cr

    mi set mlong
    mi register imputed `contVars'
    mi impute mvn `contVars' = `catVars' if sex == `curSex', emonly(iterate(1000)) rseed(1234)
    
    mat covMatrix_s`curSex'  = r(Sigma_em)
    mat meanVector_s`curSex'   = r(Beta_em)

    mi unregister `contVars'

    mat workDummy_s`curSex'   = 0,1,0,0,0,0
    mat catVector_s`curSex'   = workDummy , 1


    foreach curVar in age bmi rhr totalPAEE{

        su `curVar' if sex == `curSex'
        local `curVar'_s`curSex' = r(mean) //+ cond("`curVar'"=="totalPAEE",-2*r(sd),2*r(sd))

    }
}


foreach curExposure in acrophase24 acrophase12 acrophase8{
    
    graph use Plots/insulin_acrophaseAll_s0_m3
    graph use Plots/insulin_acrophaseAll_s1_m3

    estimates use Models/insulin_m3
    estimates esample:
    margins, over(sex) eydx(`curExposure'_sin `curExposure'_cos) asobserved post

    forvalues curSex = 0/1{

        nlcom sqrt((_b[`curExposure'_sin:`curSex'.sex])^2+(_b[`curExposure'_cos:`curSex'.sex])^2)
        nlcom (cond(_b[`curExposure'_sin:`curSex'.sex]<0,2*_pi,0)+atan2(_b[`curExposure'_sin:`curSex'.sex], _b[`curExposure'_cos:`curSex'.sex]))*1/(2*_pi)
        

        if "`curExposure'" == "acrophase24" local `curExposure'_s`curSex' = 15.23/24
        else                                local `curExposure'_s`curSex' = mod(r(b)[1,1],1)

        //if "`curExposure'" == "acrophase24" local `curExposure'_s`curSex' = r(b)[1,1]
        //else                                local `curExposure'_s`curSex' = mod(r(b)[1,1]-0.5,1)

    }
}



forvalues curSex = 0/1{

    #delimit ;          
    mata: conditionalShape( st_matrix("covMatrix_s`curSex'")    , 
                            st_matrix("meanVector_s`curSex'")   , 
                            st_matrix("catVector_s`curSex'")    ,
                            `age_s`curSex''                     ,
                            `bmi_s`curSex''                     ,
                            `rhr_s`curSex''                     ,
                            `totalPAEE_s`curSex''               ,
                            `acrophase24_s`curSex''             ,
                            `acrophase12_s`curSex''             ,
                            `acrophase8_s`curSex''              
                       )
                       ;
    #delimit cr

    local f24 Mu_hat[1,2]*cos((2*_pi/24)*(x-Mu_hat[1,1]))
    local f12 Mu_hat[1,3]*cos((2*_pi/12)*(x-12*`acrophase12_s`curSex''))
    local f8  Mu_hat[1,4]*cos((2*_pi/8)*(x-8*`acrophase8_s`curSex''))
    local mesor 2.8 //Mu_hat[1,5]

    //twoway (func y=exp(`f24'+`mesor'),range(0 24) name(s`curSex', replace))
    twoway (func y=exp(`f24'+`f12'+`f8'+`mesor'),range(0 24) name(s`curSex'_b, replace))

}

