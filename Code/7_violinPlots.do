version 17.0

frame copy dataset tempset
frame change tempset

if "`getPackages'" == "" local getPackages = 0


if `getPackages' == 1{

    capture ssc install violinplot
    capture ssc install dstat
    capture ssc install moremata
    capture ssc install palettes
    capture ssc install colrspace
    capture net install grc1leg.pkg
}


#delimit ;
local   plotVars    glucose120
                    insulin
                    leptin
                    nefa
                    adiponectin
                    ldl
                    hdl
                    mbpsys
                    mbpdia
                    crp
                    ;

local   plotXLabs   ln(glucose120)
                    ln(insulin)
                    ln(leptin)
                    ln(nefa)
                    ln(adiponectin)
                    ln(ldl)
                    ln(hdl)
                    ln(mbpsys)
                    ln(mbpdia)
                    ln(crp)
                    ;


#delimit cr


foreach curVar of local plotVars{

    gen ln_`curVar' = ln(`curVar')

}




noisi di `plotVarCount'

xtile paeeTx_0 = totalPAEE_hat if sex==0, nq(3)
xtile paeeTx_1 = totalPAEE_hat if sex==1, nq(3)

egen paeeTx = rowtotal(paeeTx_*)
drop paeeTx_*
label define quarLab 1 "T1" 2 "T2" 3 "T3"
label values paeeTx quarLab


local plotVarCount: word count `plotVars'

forvalues i = 1/`plotVarCount'{

    local curPlotVar:   word `i' of `plotVars'
    local curPlotXLab:  word `i' of `plotXLabs'

    #delimit ;
    violinplot  ln_`curPlotVar'
                , 
                over(paeeTx)
                split(sex) 
                nobox

                ytitle("`curPlotXLab'")
                ylab(#3,nogrid angle(0))

                lcolors("214 40 40" "0 47 73")
                fcolors("214 40 40" "0 47 73")
                medcolors("214 40 40" "0 47 73")

                graphregion(color(white))

                name(`curPlotVar', replace)
                ;
    #delimit cr
    
    gr_edit .xaxis1.style.editstyle majorstyle(tickstyle(textstyle(color(black)))) editcopy
    gr_edit .xaxis1.style.editstyle majorstyle(tickstyle(show_ticks(yes))) editcopy

    capture mkdir Plots
    capture mkdir Plots/violinPlots
    graph save "`curPlotVar'" Plots/violinPlots/`curPlotVar'.gph , replace

}


#delimit ;

grc1leg `plotVars'
        , 
        row(2) 
        col(5) 
        imargin(1 1 1 1)  
        graphregion(color(white) margin(l=5 r=5 t=20 b=20)) 
        legendfrom("glucose120") 
        position(3)
        ;

set graphics on ;
gr_edit .legend.Edit 
        , 
        style(  rows(2) 
                cols(0) 
                key_xsize(3) 
                labelstyle(size(vsmall)) 
                row_gap(minuscule) 
                col_gap(minuscule) 
                boxstyle(linestyle(color(white))) 
                margin(zero))
        ;

#delimit cr


capture mkdir Figures
graph export "Figures/violinPlots.png" , height(2000) width(2750) replace
graph close _all
graph drop _all

frame change dataset
frame drop tempset
