version 17.0

********************************************************************************
** Redo Cosinor Modelling
********************************************************************************

args rootPath cosinorEstimatesFile


cap mata: mata drop getEstimates()
mata:
    void getEstimates(real matrix marginMat)
    {   
        real vector i
        real vector w

        i = J(0,1,.) 
        w = J(0,1,.) 
        maxindex(marginMat,1, i, w)
        st_local("maxValue", strofreal(marginMat[i]))
        st_local("maxHour", strofreal(0.1*(i-1)))
        st_local("totalPAEE",strofreal(colsum(marginMat)*0.1))

    }
end

capture confirm file "`rootPath'/`cosinorEstimatesFile'"

if _rc != 0{

    frame copy dataset tempset
    frame change tempset

    *************************
    ** Initialise postfile **
    *************************

    local postList mesor_hat mesor_se mesor_lb mesor_ub mesor_p

    foreach curVar in sin cos acrophase amplitude{
        foreach i in 24 12 8{

            local postList `postList' `curVar'`i'_hat `curVar'`i'_se `curVar'`i'_lb `curVar'`i'_ub `curVar'`i'_p

        }
    }

    local postList `postList' maxValue maxHour totalPAEE

    capture postutil clear
    postfile cosinorpost str9 ID double(`postList') using "`rootPath'/`cosinorEstimatesFile'" , replace

    ******************************************************
    ** Initialise dataset + outcome modelling variables **
    ******************************************************
    
    keep ID paee_hour*
    reshape long paee_hour, i(ID) j(hour)
    drop if paee_hour == .

    gen sin24 = sin(hour*2*_pi/24)
    gen cos24 = cos(hour*2*_pi/24)

    gen sin12 = sin(hour*2*_pi/12)
    gen cos12 = cos(hour*2*_pi/12)

    gen sin8 = sin(hour*2*_pi/8)
    gen cos8 = cos(hour*2*_pi/8)

    ***************************************
    ** Perform IQR outlier analysis loop **
    ***************************************

    frame copy tempset tempsub
    frame change tempsub

    collapse (p25) p25=paee_hour (p75) p75=paee_hour (iqr) iqr=paee_hour, by(ID)

    frame change tempset
    frlink m:1 ID, frame(tempsub)
    frget *, from(tempsub)
    frame drop tempsub

    gen upperOutlier = cond((paee_hour - p75) / (iqr) >=  2 ,1,0)
    gen lowerOutlier = cond((paee_hour - p25) / (iqr) <= -2 ,1,0)
    
    drop if upperOutlier == 1 | lowerOutlier == 1
    drop upperOutlier lowerOutlier tempsub p25 p75 iqr

    frame copy tempset tempsub
    frame change tempsub

    collapse (count) notout=paee_hour, by(ID)

    frame change tempset
    frlink m:1 ID, frame(tempsub)
    frget *, from(tempsub)
    frame drop tempsub

    drop if notout <20
    drop notout tempsub

    ******************************
    ** Begin main analysis loop **
    ******************************
    
    qui levelsof ID, local(IDList)
    qui foreach curID of local IDList{

        noisi di "Cosinor modelling: `curID'"

        ***********************************************************
        ** Restrict sample to current ID and initialize postlist **
        ***********************************************************

        frame put * if ID == "`curID'", into(tempsub)
        frame change tempsub

        ***********************************************************************************************
        ** Apply cosinor GLM model, using log-link function to ensure paee values are always positive **
        ************************************************************************************************

        local modelLevel1 sin24 cos24 sin12 cos12 sin8 cos8
        local modelLevel2 sin24 cos24 sin12 cos12
        local modelLevel3 sin24 cos24

        forvalues m = 1/3{

            glm paee_hour `modelLevel`m'', family(gamma) link(power 0.5)
            if e(converged)==1 continue, break
        }

        if e(converged)==1{
            
            local postlist ("`curID'") 

            **********************************
            ** Get sin and cos coefficients **
            **********************************

            local mesor_hat = _b[_cons]                                         
            local mesor_se  = _se[_cons]                                        
            local mesor_lb  = _b[_cons] + invnormal(0.025)*_se[_cons]           
            local mesor_ub  = _b[_cons] + invnormal(0.975)*_se[_cons]           
            local mesor_p   = 2*normal(-abs(_b[_cons]/_se[_cons]))

            local postlist `postlist' (`mesor_hat') (`mesor_se') (`mesor_lb') (`mesor_ub') (`mesor_p')

            foreach curAxis in sin cos{

                foreach i in 24 12 8{
                    
                    capture di _b[`curAxis'`i']

                    if _rc == 0{

                        local `curAxis'`i'_hat  = _b[`curAxis'`i']                                         
                        local `curAxis'`i'_se   = _se[`curAxis'`i']                                        
                        local `curAxis'`i'_lb   = _b[`curAxis'`i'] + invnormal(0.025)*_se[`curAxis'`i']           
                        local `curAxis'`i'_ub   = _b[`curAxis'`i'] + invnormal(0.975)*_se[`curAxis'`i']           
                        local `curAxis'`i'_p    = 2*normal(-abs(_b[`curAxis'`i']/_se[`curAxis'`i']))
                    }

                    else{

                        local `curAxis'`i'_hat  = .                                       
                        local `curAxis'`i'_se   = .                                    
                        local `curAxis'`i'_lb   = .          
                        local `curAxis'`i'_ub   = .          
                        local `curAxis'`i'_p    = .
                    }

                    local postlist `postlist' (``curAxis'`i'_hat') (``curAxis'`i'_se') (``curAxis'`i'_lb') (``curAxis'`i'_ub') (``curAxis'`i'_p')
                }
            }

            *************************
            ** Estimate acrophases **
            *************************

            foreach i in 24 12 8{

                capture nlcom cond(_b[sin`i']<0,`i',0) + atan2(_b[sin`i'],_b[cos`i'])*`i'/(2*_pi) 

                if _rc == 0{

                    local acrophase`i'_hat   = r(b)[1,1]                                         
                    local acrophase`i'_se    = sqrt(r(V)[1,1])                                   
                    local acrophase`i'_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
                    local acrophase`i'_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
                    local acrophase`i'_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

                }

                else{
                    
                    local acrophase`i'_hat   = .                                       
                    local acrophase`i'_se    = .                            
                    local acrophase`i'_lb    = .     
                    local acrophase`i'_ub    = . 
                    local acrophase`i'_p     = .
                }

                local postlist `postlist' (`acrophase`i'_hat') (`acrophase`i'_se') (`acrophase`i'_lb') (`acrophase`i'_ub') (`acrophase`i'_p') 
                
            }

            *************************
            ** Estimate amplitudes **
            *************************

            foreach i in 24 12 8{

                capture nlcom sqrt(_b[sin`i']^2+_b[cos`i']^2)

                if _rc == 0{

                    local amplitude`i'_hat   = cond(_rc == 0, r(b)[1,1],.)                                        
                    local amplitude`i'_se    = sqrt(r(V)[1,1])                                   
                    local amplitude`i'_lb    = r(b)[1,1] + invnormal(0.025)*sqrt(r(V)[1,1])      
                    local amplitude`i'_ub    = r(b)[1,1] + invnormal(0.975)*sqrt(r(V)[1,1])      
                    local amplitude`i'_p     = 2*normal(-abs(r(b)[1,1]/sqrt(r(V)[1,1])))

                }

                else{

                    local amplitude`i'_hat   = .                                       
                    local amplitude`i'_se    = .
                    local amplitude`i'_lb    = .     
                    local amplitude`i'_ub    = .     
                    local amplitude`i'_p     = .        
                }

                local postlist `postlist' (`amplitude`i'_hat') (`amplitude`i'_se') (`amplitude`i'_lb') (`amplitude`i'_ub') (`amplitude`i'_p')          
            }
            
            **********************************************************************
            ** Estimate max attained value, the hour it occurrs, and total paee **
            **********************************************************************

            local marginList
            forvalues i = 0(0.1)24{
                foreach j in 24 12 8{
                    capture di _b[sin`j'] _b[cos`j']
                    if _rc == 0 local per`j' sin`j'=(`=sin(`i'*2*_pi/`j')')cos`j'=(`=cos(`i'*2*_pi/`j')')
                    else        local per`j'

                }
                local marginList `marginList' at(`per24'`per12'`per8')
            }

            margin, `marginList' post
            mata: getEstimates(st_matrix("r(table)")[1,.]')
            local postlist `postlist' (`maxValue') (`maxHour') (`totalPAEE')
            
            post cosinorpost `postlist'
        
        }

        frame change tempset
        frame drop tempsub

    }
    
    postclose cosinorpost
    
    frame change dataset
    frame drop tempset
}

asdf

joinby ID using "`rootPath'/`cosinorEstimatesFile'"



            /*
            noisi nlcom _b[135._at]

            local marginList
            forvalues i = 1/`e(k_at)'{

                local marginList = "`marginList'"+"_b[`i'._at]"+cond(`i'!=`e(k_at)',"+","")

            }
            noisi nlcom (`marginList')*0.1
            */