version 17.0

clear
set more off
set seed 1234

args filePath

use "`filePath'"

qui do Code/Sub/initialiseCovariates.do
qui do Code/Sub/initialiseOutcomes.do
qui do Code/Sub/initialisePAEE.do


capture ssc install violinplot
capture ssc install dstat
capture ssc install moremata
capture ssc install palettes
capture ssc install colrspace

xtile paeeTx_0 = paeeTt if sex==0, nq(4)
xtile paeeTx_1 = paeeTt if sex==1, nq(4)

egen paeeTx = rowtotal(paeeTx_*)
drop paeeTx_*

violinplot fatMass, over(paeeTx) split(sex)

asdd

graph box fatMass, over(sex) over(paeeTx) asy 
