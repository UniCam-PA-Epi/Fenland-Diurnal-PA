version 17.0

frame copy dataset tempset
frame change tempset

set seed 1234

*************************
** Initilaise analysis **
*************************

// Exclude participants for whom the cosinor model did not fit, those with insufficient wear time, and those without fat mass measurement.

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05

drop if mesor_hat == .
drop if amplitude24_hat == . | amplitude12_hat == . | amplitude8_hat == .
drop if acrophase24_hat == . | acrophase12_hat == . | acrophase8_hat  == .

drop if P_Pwear_consolidated<72
drop if fatMass == .
drop if fatFreeMass == .

drop *_se *_lb *_ub *_p
rename *_hat *

#delimit ;
local modelSpec c.sin24
                c.cos24 
                c.sin12 
                c.cos12 
                c.sin8 
                c.cos8 
                //c.mesor 

                i.sex
                c.totalPAEE

                i.sex#c.totalPAEE

                c.sin24#i.sex
                c.cos24#i.sex
                c.sin12#i.sex
                c.cos12#i.sex 
                c.sin8#i.sex 
                c.cos8#i.sex 
                //c.mesor#i.sex

                c.sin24#c.totalPAEE
                c.cos24#c.totalPAEE
                c.sin12#c.totalPAEE
                c.cos12#c.totalPAEE 
                c.sin8#c.totalPAEE 
                c.cos8#c.totalPAEE 
                //c.mesor#c.totalPAEE
                ;
#delimit cr

local fullModel
foreach curVar in sin24 cos24 sin12 cos12 sin8 cos8 mesor{

    local subModel
    foreach curPredictor of local modelSpec{

        if regexm("`curPredictor'", "`curVar'")==1 continue
        else local subModel `subModel' `curPredictor'

    }
    local fullModel `fullModel' (`curVar' `subModel')
}

sureg `fullModel'
estimates store m1

foreach curVar in sin24 cos24 sin12 cos12 sin8 cos8 mesor{

    margins, over(sex) at((p20) totalPAEE) at((p40) totalPAEE) at((p60) totalPAEE) at((p80) totalPAEE) predict(xb equation(`curVar')) asobserved post
    
    forvalues curSex = 0/1{
        forvalues curQuart = 1/4{
        
            local `curVar'_`curQuart'_`curSex' = _b[`curQuart'._at#`curSex'.sex]
        }
    }

    estimates restore m1
}

forvalues curSex = 0/1{
    forvalues curQuart = 1/4{

        local per24 `sin24_`curQuart'_`curSex''*sin(x*2*_pi/24)+`cos24_`curQuart'_`curSex''*cos(x*2*_pi/24)
        local per12 `sin12_`curQuart'_`curSex''*sin(x*2*_pi/12)+`cos12_`curQuart'_`curSex''*cos(x*2*_pi/12)
        local per8  `sin8_`curQuart'_`curSex''*sin(x*2*_pi/8)+`cos8_`curQuart'_`curSex''*cos(x*2*_pi/8)
        local mesor 2.5 //`mesor_`curQuart'_`curSex''

        local cosinor_`curQuart'_`curSex' y=(60/1000)*exp(`per24'+`per12'+`per8'+`mesor')
    }

    #delimit ;
    twoway  (func `cosinor_4_`curSex'', range(0 24) recast(area) lcolor("0   47  73 %0") fcolor("0   47  73 %100")) 
            (func `cosinor_3_`curSex'', range(0 24) recast(area) lcolor("214 40  40 %0") fcolor("214 40  40 %100")) 
            (func `cosinor_2_`curSex'', range(0 24) recast(area) lcolor("247 128 0  %0") fcolor("247 128 0  %100")) 
            (func `cosinor_1_`curSex'', range(0 24) recast(area) lcolor("252 192 73 %0") fcolor("252 192 73 %100")) 
            (
            ,
            title(`=cond(`curSex'==0,"Women","Men")',size(3) color(black))
            yscale(`=cond(`curSex'==0,"","alt")')
            xtitle("Clocktime",size(3) color(black))
            ytitle("PAEE",size(3) color(black))
            xlab(0 8 16 24, nogrid angle(0))
            ylab(#3       , nogrid angle(0))
            graphregion(color(white))
            name(g`curSex', replace)
            legend(off)   
            )
            ;
    #delimit cr
    
}

graph combine g0 g1, row(1) col(2) ycommon graphregion(color(white) margin(l=32 r=32 t=30 b=30)) imargin(2 2 0 0)






asdf
















#delimit ;
reg sin24   c.cos24 
            c.sin12 
            c.cos12 
            c.sin8 
            c.cos8 
            c.mesor 

            i.sex
            c.totalPAEE

            i.sex#c.totalPAEE

            c.cos24#i.sex
            c.sin12#i.sex
            c.cos12#i.sex 
            c.sin8#i.sex 
            c.cos8#i.sex 
            c.mesor#i.sex

            c.cos24#c.totalPAEE
            c.sin12#c.totalPAEE
            c.cos12#c.totalPAEE 
            c.sin8#c.totalPAEE 
            c.cos8#c.totalPAEE 
            c.mesor#c.totalPAEE
            ;
#delimit cr









asdf


keep ID paee_hour* pwear_hour* sex totalPAEE_hat
reshape long paee_hour pwear_hour, i(ID) j(hour)
drop if paee_hour == .

// Generate sin and cos axis variables for 24, 12, and 8 hour periods

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

// Collapse data set on percentiles, then exclude values that are + or - 2 IQR

collapse (p25) p25=paee_hour (p75) p75=paee_hour (iqr) iqr=paee_hour, by(ID)

frame change tempset
frlink m:1 ID, frame(tempsub)
frget *, from(tempsub)
frame drop tempsub

gen upperOutlier = cond((paee_hour - p75) / (iqr) >=  2 ,1,0)
gen lowerOutlier = cond((paee_hour - p25) / (iqr) <= -2 ,1,0)

drop if upperOutlier == 1 | lowerOutlier == 1
drop upperOutlier lowerOutlier tempsub p25 p75 iqr

// Drop those with more than 4 data points excluded by IQR outlier analysis

frame copy tempset tempsub
frame change tempsub

collapse (count) notout=paee_hour, by(ID)

frame change tempset
frlink m:1 ID, frame(tempsub)
frget *, from(tempsub)
frame drop tempsub

drop if notout <20
drop notout tempsub


rename *_hat *
replace totalPAEE = totalPAEE *60/1000

egen totalPAEE_std = std(totalPAEE)

#delimit ;
glm paee_hour   c.sin24 
                c.cos24 
                c.sin12 
                c.cos12 
                c.sin8 
                c.cos8

                i.sex
                
                c.sin24#i.sex
                c.cos24#i.sex 
                c.sin12#i.sex 
                c.cos12#i.sex 
                c.sin8#i.sex 
                c.cos8#i.sex
                        

                [aw=pwear_hour]
                , 
                family(gamma) link(log) vce(cluster ID)
                ;
#delimit cr

margins, over(sex) predict(xb)                              at((p75) sin24 (p75) cos24) asobserved 
margins, over(sex) eydx(sin24 cos24 sin12 cos12 sin8 cos8)  at((p75) sin24 (p75) cos24) asobserved

asdf

forvalues curSex = 0/2{

    su totalPAEE if sex == `curSex', detail
    
    local lb1 = r(p1)
    local ub1 = r(p25)
    
    local lb2 = r(p25)
    local ub2 = r(p50)

    local lb3 = r(p50)
    local ub3 = r(p75)

    local lb4 = r(p75)
    local ub4 = r(p99)

     forvalues i = 1/4{

        glm paee_hour sin24 cos24 sin12 cos12 sin8 cos8 [aw=pwear_hour] if sex == `curSex' & totalPAEE >= `lb`i'' & totalPAEE < `ub`i'', family(gamma) link(log) vce(cluster ID)

        foreach curVar in sin24 cos24 sin12 cos12 sin8 cos8{
            local `curVar'_`i' = _b[`curVar']
        }
        local mesor_`i' = _b[_cons]

        local per24_`i' `sin24_`i''*sin(x*2*_pi/24)+`cos24_`i''*cos(x*2*_pi/24)
        local per12_`i' `sin12_`i''*sin(x*2*_pi/12)+`cos12_`i''*cos(x*2*_pi/12)
        local per8_`i'  `sin8_`i''*sin(x*2*_pi/8)+`cos8_`i''*cos(x*2*_pi/8)

        local funcQuart_`i' y=(60/1000)*exp(`per24_`i''+`per12_`i''+`per8_`i''+`mesor_`i'')

     }

    #delimit ;
    twoway  (func `funcQuart_4', range(0 24) recast(area) lcolor("0   47  73 %0") fcolor("0   47  73 %100")) 
            (func `funcQuart_3', range(0 24) recast(area) lcolor("214 40  40 %0") fcolor("214 40  40 %100")) 
            (func `funcQuart_2', range(0 24) recast(area) lcolor("247 128 0  %0") fcolor("247 128 0  %100")) 
            (func `funcQuart_1', range(0 24) recast(area) lcolor("252 192 73 %0") fcolor("252 192 73 %100")) 
            (
            ,
            title(`=cond(`curSex'==0,"Women","Men")',size(3) color(black))
            yscale(`=cond(`curSex'==0,"","alt")')
            xtitle("Clocktime",size(3) color(black))
            ytitle("PAEE",size(3) color(black))
            xlab(0 8 16 24, nogrid angle(0))
            ylab(#3       , nogrid angle(0))
            graphregion(color(white))
            name(g`curSex', replace)
            legend(off)   
            )
            ;
    #delimit cr
}

graph combine g0 g1, row(1) col(2) ycommon graphregion(color(white) margin(l=32 r=32 t=30 b=30)) imargin(2 2 0 0)














asdf


asdf




// Exclude participants for whom the cosinor model did not fit, those with insufficient wear time, and those without fat mass measurement.

drop if sin24_p>0.05 & cos24_p>0.05 & sin12_p>0.05 & cos12_p>0.05 & sin8_p>0.05 & cos8_p>0.05

drop if mesor_hat == .
drop if amplitude24_hat == . | amplitude12_hat == . | amplitude8_hat == .
drop if acrophase24_hat == . | acrophase12_hat == . | acrophase8_hat  == .

drop if P_Pwear_consolidated<72
drop if fatMass == .
drop if fatFreeMass == .

drop *_se *_lb *_ub *_p
rename *_hat *

replace totalPAEE = totalPAEE *60/1000

egen sin24_std = std(sin24)
egen cos24_std = std(cos24)
egen sin12_std = std(sin12)
egen cos12_std = std(cos12)
egen sin8_std  = std(sin8)
egen cos8_std  = std(cos8)
egen mesor_std = std(mesor)

forvalues curSex = 0/1{
    
    su totalPAEE if sex == `curSex', detail
    
    local lb1 = r(p1)
    local ub1 = r(p25)
    
    local lb2 = r(p25)
    local ub2 = r(p50)

    local lb3 = r(p50)
    local ub3 = r(p75)

    local lb4 = r(p75)
    local ub4 = r(p99)

    forvalues i = 1/4{

        asdf


        
        #delimit ;
        sureg   (sin24 c.cos24_std c.sin12_std c.cos12_std c.sin8_std  c.cos8_std c.mesor_std)
                (cos24 c.sin24_std c.sin12_std c.cos12_std c.sin8_std  c.cos8_std c.mesor_std)
                (sin12 c.sin24_std c.cos24_std c.cos12_std c.sin8_std  c.cos8_std c.mesor_std)
                (cos12 c.sin24_std c.cos24_std c.sin12_std c.sin8_std  c.cos8_std c.mesor_std)
                (sin8  c.sin24_std c.cos24_std c.sin12_std c.cos12_std c.cos8_std c.mesor_std)
                (cos8  c.sin24_std c.cos24_std c.sin12_std c.cos12_std c.sin8_std c.mesor_std)
                (mesor c.sin24_std c.cos24_std c.sin12_std c.cos12_std c.sin8_std c.cos8_std )
                if
                sex == `curSex' & totalPAEE >= `lb`i'' & totalPAEE < `ub`i''
                ;
        #delimit cr

        foreach curVar in sin24 cos24 sin12 cos12 sin8 cos8 mesor{
            local `curVar'_`i' = _b[`curVar':_cons]
        }

        local func24_`i' `sin24_`i''*sin(x*2*_pi/24)+`cos24_`i''*cos(x*2*_pi/24)
        local func12_`i' `sin12_`i''*sin(x*2*_pi/12)+`cos12_`i''*cos(x*2*_pi/12)
        local func8_`i'  `sin8_`i''*sin(x*2*_pi/8)+`cos8_`i''*cos(x*2*_pi/8)

        local funcAll_`i' y=(60/1000)*exp(`func24_`i''+`func12_`i''+`func8_`i''+`mesor_`i'')
    }

    #delimit ;
    twoway  (func `funcAll_4', range(0 24) recast(area) lcolor("0   47  73 %0") fcolor("0   47  73 %100")) 
            (func `funcAll_3', range(0 24) recast(area) lcolor("214 40  40 %0") fcolor("214 40  40 %100")) 
            (func `funcAll_2', range(0 24) recast(area) lcolor("247 128 0  %0") fcolor("247 128 0  %100")) 
            (func `funcAll_1', range(0 24) recast(area) lcolor("252 192 73 %0") fcolor("252 192 73 %100")) 
            (
            ,
            title(`=cond(`curSex'==0,"Women","Men")',size(3) color(black))
            yscale(`=cond(`curSex'==0,"","alt")')
            xtitle("Clocktime",size(3) color(black))
            ytitle("PAEE",size(3) color(black))
            xlab(0 8 16 24, nogrid angle(0))
            ylab(#3       , nogrid angle(0))
            graphregion(color(white))
            name(g`curSex', replace)
            legend(off)   
            )
            ;
    #delimit cr
}

graph combine g0 g1, row(1) col(2) ycommon graphregion(color(white) margin(l=32 r=32 t=30 b=30)) imargin(2 2 0 0)


asdf

#delimit ;
local contVars  sin24
                cos24
                sin12
                cos12
                sin8 
                cos8
                mesor
                //age
                //fatFreeMass
                //fatMass
                //rhr
                ;

local catVars   //b1.sex
                //b1.work_status
                //b2.marital_status
                //b1.education
                //b2.income
                ;
#delimit cr


qui mi set mlong
qui mi register imputed `contVars'
mi impute mvn `contVars' = `catVars', emonly(iterate(1000)) rseed(1234)

asdf

mat covMatrix  = r(Sigma_em)
mat meanMatrix = r(Beta_em)
mi unregister `contVars'


foreach curVar of local contVars{
    mat `curVar' = .
}



mata: conditionalShape(`i',128)

asdf










asdf

egen acrophase24_std = std(acrophase24)
egen amplitude24_std = std(amplitude24)

egen acrophase12_std = std(acrophase12)
egen amplitude12_std = std(amplitude12)

egen acrophase8_std = std(acrophase8)
egen amplitude8_std = std(amplitude8)

egen mesor_std = std(mesor)
egen totalPAEE_std = std(totalPAEE)


** Lets cluster! **



// ktest


/*
qui forvalues i = 1/20{

    cluster kmedians acrophase24_std amplitude24_std acrophase12_std amplitude12_std acrophase8_std amplitude8_std mesor_std totalPAEE_std, k(`i') keepcenters gen(k`i') measure(L1)

    forvalues j = 1/`i'{

        cluster measures acrophase24_std amplitude24_std acrophase12_std amplitude12_std acrophase8_std amplitude8_std mesor_std totalPAEE_std  if k`i'==`j', compare(`=12064+`j'') gen(k`i'_j`j')

    }

    egen k`i'_dt = rowtotal(k`i'_j*)
    drop in 12065/`=12064+`i''
    drop k`i'_j*

    su k`i'_dt, detail
    noisi di "`r(sum)'"

}
*/

// after test, k = 8

cluster kmedians acrophase24_std amplitude24_std acrophase12_std amplitude12_std acrophase8_std amplitude8_std mesor_std totalPAEE_std, k(8) gen(k8) measure(L1)

forvalues i = 1/8{

    su acrophase8 if k8 == `i'

}