version 17.0

set seed 1234

** Control whether to run the optimal k selection process (1 = yes, 0 = no) **
local getBestK = 0
local plotClusters = 0


** Standardize variables for k-means clustering (improves L2 distance computation) **
egen sin24_std = std(sin24_hat)
egen cos24_std = std(cos24_hat)
egen sin12_std = std(sin12_hat)
egen cos12_std = std(cos12_hat)
egen sin8_std  = std(sin8_hat)
egen cos8_std  = std(cos8_hat)
egen mesor_std  = std(mesor_hat)

if `getBestK' == 1{

    ** Calculate the within-cluster sum of squares (WCSS) for k = 1 to 20 **
    ** This will be used to estimate the optimal k using the elbow method **

    gen kValue = .
    gen WCSS = .

    qui forvalues i = 1/20{

        ** Perform k-means clustering with the current k **
        cluster kmeans sin24_std cos24_std sin12_std cos12_std sin8_std cos8_std mesor_std , k(`i') gen(k`i') measure(L2) start(krandom(1234))

        ** Calculate WCSS for the current k **
        local errorSum = 0
        forvalues j = 1/`i'{
    
            foreach curVar in sin24_std cos24_std sin12_std cos12_std sin8_std cos8_std mesor_std{

                su `curVar'                             if k`i' == `j' 
                gen diff = (`curVar'-r(mean))^2         if k`i' == `j' 
                su diff                                 if k`i' == `j' , detail
                local errorSum = `errorSum' + r(sum) 
                drop diff
            } 
        }

        drop k`i'      
        replace kValue = `i'        in `i'
        replace WCSS = `errorSum'   in `i'
        
    }
    
    ** Use linear splines regression to find the "elbow" in the WCSS curve **
    ** The optimal k is estimated as the midpoint of the best spline knots **

    gen knot1 = .
    gen ll = .

    local curIndex = 1
    forvalues i = 1(0.1)20{

        ** Generate linear splines with two knots **
        mkspline spline1 `i' spline2 = kValue  

        ** Fit a regression model with an inverse square link function **
        glm WCSS spline*, family(gaussian) link(power -2) 
        drop spline*
        
        replace knot1 = `i'     in `curIndex'
        replace ll    =  e(ll)  in `curIndex'

        local curIndex = `curIndex' + 1
    
    }

    ** Select k based on the spline model with the lowest log-likelihood **
    preserve
    gsort -ll
    su knot1 in 1/20
    local kSelected = round(r(mean))
    restore

    noisi di "Optimal K : `kSelected'"
    
    ** Clean up temporary variables **
    cluster drop _all
    drop kValue WCSS knot1 knot2 ll

}

** Use pre-determined k if getBestK is off. The above process yields k = 6 **
else local kSelected = 6

** Perform k-means clustering with the selected k **

cluster kmeans sin24_std cos24_std sin12_std cos12_std sin8_std cos8_std mesor_std, k(`kSelected') gen(kGroup) measure(L2) start(krandom(1234)) name(kMeans)
drop sin24_std cos24_std sin12_std cos12_std sin8_std cos8_std mesor_std



if `plotClusters' == 1 {

    forvalues k = 1/`kSelected'{

        foreach p in 24 12 8{

            su sin`p'_hat if kGroup == `k'
            local curSin`p' = r(mean)
            su cos`p'_hat if kGroup == `k'
            local curCos`p' = r(mean)

            local f`p' `curSin`p''*sin((2*_pi/`p')*x)+`curCos`p''*cos((2*_pi/`p')*x)

        }

        su mesor_hat if kGroup == `k'
        local curMesor = r(mean)

        #delimit ;
        twoway (func    y=exp(`f24'+`f12'+`f8'+`curMesor')
                        ,
                        recast(area) 
                        lcolor("214 40  40 %0") 
                        fcolor("214 40  40 %100") 
                        
                        title("k`k'",color(black) size(2.5)) 
                        
                        xtitle("Clocktime (hour)", size(2.5)) 
                        range(0 24) 

                        yscale(off) 

                        ylab(,nogrid) 
                        xlab(0 12 24, labsize(2.5)) 
                        graphregion(color(white)) 
                        name(k`k', replace)
                )
                ;
        #delimit cr

    }

    graph combine k1 k2 k3 k4 k5 k6, row(2) col(3) ycommon xcommon imargin(1 1 1 1) graphregion(color(white) margin(t=26 b=26 l=42 r=42) ) name(kClusterPAEEProfiles, replace)
    
    capture mkdir Figures
    graph export "Figures/kClusterPAEEProfiles.png" , height(2000) width(2750) replace
    graph close k1 k2 k3 k4 k5 k6 kClusterPAEEProfiles
    graph drop  k1 k2 k3 k4 k5 k6 kClusterPAEEProfiles

}

foreach p in 24 12 8 {

    forvalues k=1/6{

        ellip amplitude`p'_hat acrophase`p'_hat if kGroup == `k', g(y`k' x`k') nograph c(chi2) level(95) 

    }

    #delimit ;

    twoway
    (scatter amplitude`p'_hat acrophase`p'_hat if kGroup == 1, msize(0.4) mcolor("168 199 253  %10")   mlwidth(0) )
    (scatter amplitude`p'_hat acrophase`p'_hat if kGroup == 2, msize(0.4) mcolor("252 63  246  %10")   mlwidth(0) )
    (scatter amplitude`p'_hat acrophase`p'_hat if kGroup == 3, msize(0.4) mcolor("254 186 122  %10")   mlwidth(0) )  
    (scatter amplitude`p'_hat acrophase`p'_hat if kGroup == 4, msize(0.4) mcolor("178 214 0    %10")   mlwidth(0) ) 
    (scatter amplitude`p'_hat acrophase`p'_hat if kGroup == 5, msize(0.4) mcolor("37  175 22   %10")   mlwidth(0) ) 
    (scatter amplitude`p'_hat acrophase`p'_hat if kGroup == 6, msize(0.4) mcolor("50  223 176  %10")   mlwidth(0) )
    (line    y1 x1, lcolor("168 199 253  %100"))
    (line    y2 x2, lcolor("252 63  246  %100"))
    (line    y3 x3, lcolor("254 186 122  %100"))
    (line    y4 x4, lcolor("178 214 0    %100"))
    (line    y5 x5, lcolor("37  175 22   %100"))
    (line    y6 x6, lcolor("50  223 176  %100"))
    (,
    xscale(off)
    yscale(off)
    xlab(,nogrid)
    ylab(,nogrid)
    plotregion(color(black))
    graphregion(color(black))
    legend(off)
    name(g`p', replace)
    )
    ;

    #delimit cr

    drop y1 y2 y3 y4 y5 y6 x1 x2 x3 x4 x5 x6

}

graph combine g24 g12 g8, row(1) col(3) graphregion(color(black) margin(l=15 r=15 t=20 b=20))



asdf

/*
 || 


scatter sin24_hat cos24_hat if kGroup == 2, msize(0.2) mcolor(blue%20)      || 
scatter sin24_hat cos24_hat if kGroup == 3, msize(0.2) mcolor(orange%20)    || scatter sin24_hat cos24_hat if kGroup == 4, msize(0.2) mcolor(purple%20) || scatter sin24_hat cos24_hat if kGroup == 5 , msize(0.2) mcolor(black%20) || scatter sin24_hat cos24_hat if kGroup == 6 , msize(0.2) mcolor(yellow%20)

*/

    (func y= sqrt(1^2-x^2)  , range(-1 1) lcolor(gs16%30))
    (func y=-sqrt(1^2-x^2)  , range(-1 1) lcolor(gs16%30))
    (func y= sqrt(2^2-x^2)  , range(-2 2) lcolor(gs16%30))
    (func y=-sqrt(2^2-x^2)  , range(-2 2) lcolor(gs16%30))
    (func y= sqrt(3^2-x^2)  , range(-3 3) lcolor(gs16%30))
    (func y=-sqrt(3^2-x^2)  , range(-3 3) lcolor(gs16%30))
    (func y= sqrt(4^2-x^2)  , range(-4 4) lcolor(gs16%30))
    (func y=-sqrt(4^2-x^2)  , range(-4 4) lcolor(gs16%30))
    (func y= x              , range(-3.02 3.02) lcolor(gs16%30))
    (func y=-x              , range(-3.02 3.02) lcolor(gs16%30))
    (, yline(0, lcolor(gs16%30)))
    (, xline(0, lcolor(gs16%30)))