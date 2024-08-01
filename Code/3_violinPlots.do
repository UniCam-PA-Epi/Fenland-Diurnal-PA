version 17.0

clear
set more off
set seed 1234

args filePath

use "`filePath'"

qui do Code/Sub/initialiseCovariates.do
qui do Code/Sub/initialiseOutcomes.do
qui do Code/Sub/initialisePAEE.do

//capture ssc install violinplot
//capture ssc install dstat
//capture ssc install moremata
//capture ssc install palettes
//capture ssc install colrspace

xtile paeeTx_0 = paeeTt if sex==0, nq(4)
xtile paeeTx_1 = paeeTt if sex==1, nq(4)

egen paeeTx = rowtotal(paeeTx_*)
drop paeeTx_*
label define quarLab 1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4"
label values paeeTx quarLab

set graphics off

violinplot log_leptin, over(paeeTx) split(sex) ylab(,nogrid angle(0)) lcolors(navy maroon) graphregion(color(white)) legend(off) name(g1, replace)
violinplot log_glucose0, over(paeeTx) split(sex) ylab(,nogrid angle(0)) lcolors(navy maroon) graphregion(color(white)) legend(off) name(g2, replace)
violinplot log_fatMass, over(paeeTx) split(sex) ylab(,nogrid angle(0)) lcolors(navy maroon) graphregion(color(white)) legend(off) name(g3, replace)
violinplot log_fatFreeMass, over(paeeTx) split(sex) ylab(,nogrid angle(0)) lcolors(navy maroon) graphregion(color(white)) legend(off) name(g4, replace)

set graphics on

graph combine g1 g2 g3 g4, row(2) col(2)

asdd

graph box fatMass, over(sex) over(paeeTx) asy 
