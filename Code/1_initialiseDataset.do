version 17.0

frames reset
frame rename default dataset

args rootPath               ///
     fenlandReleaseFile     ///
     cosinorEstimatesFile

use "`rootPath'/`fenlandReleaseFile'"

do Code/Sub/initialiseCovariates.do
do Code/Sub/initialiseOutcomes.do
do Code/Sub/initialisePAEE.do "`rootPath'" "`cosinorEstimatesFile'"  

