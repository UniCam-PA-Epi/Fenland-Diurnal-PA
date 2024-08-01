version 17.0

clear
set more off
set seed 1234

args filePath

do Code/1_descriptivesTables.do "`filePath'"
do Code/2_jointAssociations.do  "`filePath'"
//do Code/3_violinPlot.do "`filePath'"