version 17.0

frame copy dataset tempset
frame change tempset

args getPackages
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
local   plotVars    fatMass
                    fatFreeMass
                    mbpsys
                    mbpdia
                    glucose0
                    insulin
                    adiponectin
                    crp
                    ldl
                    hdl
                    leptin
                    nefa
                    ;

local   plotXLabs   ln(FM)
                    ln(FFM)
                    ln(SBP)
                    ln(DBP)
                    ln(FPG)
                    ln(FI)
                    ln(Adadiponectin)
                    ln(CRP)
                    ln(LDL)
                    ln(HDL)
                    ln(Leptin)
                    ln(NEFA)
                    ;

local   plotYLabs   1(2)5
                    3(1)5
                    4(.75)5.5
                    3.5(.75)5
                    1(1)3
                    0(4)8
                    -1(2.5)4
                    -4(5)6
                    -1(2)3
                    -2(2)2
                    -4(5)6
                    3(3)9
                    ;

#delimit cr



noisi di `plotVarCount'

xtile paeeTx_0 = paeeTt if sex==0, nq(4)
xtile paeeTx_1 = paeeTt if sex==1, nq(4)

egen paeeTx = rowtotal(paeeTx_*)
drop paeeTx_*
label define quarLab 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4"
label values paeeTx quarLab


set graphics off
local plotVarCount: word count `plotVars'

forvalues i = 1/`plotVarCount'{

    local curPlotVar:   word `i' of `plotVars'
    local curPlotXLab:  word `i' of `plotXLabs'
    local curPlotYLab:  word `i' of `plotYLabs'

    #delimit ;
    violinplot  log_`curPlotVar'
                , 
                over(paeeTx)
                split(sex) 
                nobox

                ytitle("`curPlotXLab'")
                ylab(`curPlotYLab',nogrid angle(0))

                lcolors(navy maroon) 
                graphregion(color(white))

                name(`curPlotVar', replace)
                ;
    #delimit cr

    gr_edit .xaxis1.style.editstyle majorstyle(tickstyle(textstyle(color(white)))) editcopy
    gr_edit .xaxis1.style.editstyle majorstyle(tickstyle(show_ticks(yes))) editcopy
    
}

#delimit ;

grc1leg `plotVars'
        , 
        row(3) 
        col(4) 
        imargin(1 1 1 1)  
        graphregion(color(white) margin(l=15 r=15 t=12 b=12)) 
        legendfrom("fatMass") 
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

graph export "Results/3_violinPlots.png" , height(2000) width(2750) replace
graph close _all

frame change dataset
frame drop tempset
