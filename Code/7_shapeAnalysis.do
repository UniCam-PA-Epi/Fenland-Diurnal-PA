version 17.0

// Exclude participants for whom the cosinor model did not fit, those with insufficient wear time, and those without fat mass measurement.

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05
drop if P_Pwear_consolidated<72
drop if fatMass == .

drop *_se *_lb *_ub *_p
rename *_hat *

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
mi impute mvn `contVars' = `catVars' if sex == 1   , emonly(iterate(1000)) rseed(1234)

mat covMatrix  = r(Sigma_em)
mat meanVector = r(Beta_em)

mat workDummy   = 0,1,0,0,0,0
mat catVector   = workDummy , 1

mata
    
    covMatrix   = st_matrix("covMatrix")
    meanVector  = st_matrix("meanVector")
    catVector   = st_matrix("catVector")

    meanVector = cross(catVector',meanVector)

    X_1 = 48.56946 , 27.31037, 61.40449, 994.194, 24*.3387948, 12*.9222235, 8*.5243812

    conditioningIndices = 1,2,3,4,5,6,7
    dependentIndices = 8, 9, 10, 11
 
    Mu_1 = meanVector[conditioningIndices]
    Mu_2 = meanVector[dependentIndices]

    Sigma_11 = covMatrix[conditioningIndices,conditioningIndices]
    Sigma_12 = covMatrix[conditioningIndices,dependentIndices]
    Sigma_22 = covMatrix[dependentIndices,dependentIndices]

    solved = cholsolve(Sigma_11,Sigma_12)
    Mu_hat      = Mu_2 :+ quadcross((X_1:-Mu_1)',solved)
    Sigma_hat   = Sigma_22 - quadcross(Sigma_12,solved)
    
    Mu_hat
    Sigma_hat

    st_matrix("Mu_hat", Mu_hat)

    standardNormalRandoms = rnormal(cols(Mu_hat),5,0,1)
    correlatedNormalRandoms = quadcross(standardNormalRandoms,cholesky(Sigma_hat))
    multivariateNormalRandoms = correlatedNormalRandoms + J(5, 1, Mu_hat) 
    multivariateNormalRandoms

    st_matrix("randomDraws", multivariateNormalRandoms)

end

local f24 Mu_hat[1,1]*cos((2*_pi/24)*(x-24*.3387948))
local f12 Mu_hat[1,2]*cos((2*_pi/12)*(x-12*.9222235))
local f8  Mu_hat[1,3]*cos((2*_pi/8)*(x-8*.5243812))
local mesor Mu_hat[1,4]

twoway (func y=exp(`f24'+`f12'+`f8'+`mesor'),range(0 24))




