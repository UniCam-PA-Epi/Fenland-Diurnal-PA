version 17.0

clear
estimates clear
set more off
set seed 1234

args rootPath fenlandReleaseFile cosinorEstimatesFile

do Code/1_initialiseDataset.do "`rootPath'" "`fenlandReleaseFile'" "`cosinorEstimatesFile'"
//do Code/2_descriptivesTables.do 
do Code/3_jointAssociations.do  
//do Code/4_violinPlots.do        