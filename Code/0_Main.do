version 17.0

clear
set more off
set seed 1234

args rootPath getPackages

frames reset
frame rename default dataset

use "`rootPath'/fenlandRaw.dta"

include 1_initialiseCovariates.do
include 2_initialiseOutcomes.do
include 3_applyCosinorModel.do
include 4_applyExclusions.do
include 5_clusterAnalysis.do

include 6_descriptivesTables.do
//include 7_violinPlots.do
//include 8_cosinorFeatureAnalysis.do
//include 9_totalPAEEAnalysis.do
//include 10_panelPlots.do

