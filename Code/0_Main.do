version 17.0

clear
estimates clear
set more off
set seed 1234

args rootPath

frames reset
frame rename default dataset

use "`rootPath'/fenlandRaw.dta"

include 1_initialiseOutcomes.do
include 2_initialiseCovariates.do
include 3_applyCosinorModel.do
include 4_descriptivesTables.do
include 5_violinPlots.do
include 6_jointAssociations.do
