clear all
macro drop _all

// Instal required packages

cap net install psidtools, from("http://fmwww.bc.edu/RePEc/bocode/p") // PSID tools (make sure you're using a recent version to use the 2021 files)
cap net install st0110_1, from("http://www.stata-journal.com/software/sj15-3")  // Freduse (Stata<14)
cap net install blindschemes, from("http://fmwww.bc.edu/RePEc/bocode/b") // plotplainblind
cap net install winsor2.pkg, from("http://fmwww.bc.edu/RePEc/bocode/w") // Winsor2
cap net install regsave, from("https://raw.githubusercontent.com/reifjulian/regsave/master") replace // Regsave
cap net install grc1leg.pkg, from("http://www.stata.com/users/vwiggins") // grc1leg

set varabbrev off
set more off
set scheme plotplainblind, permanently
global basepath "/mcr/res-m1eeb00/Projects/delete/HomeownershipBankMomDad.jl"               // Path to root directory of the git repo.
global PSIDpath "$basepath/data/PSID"

global realvars "valhouse wealthwoequit wealth income valstocks homeequity transf2val houseexpend rent expend_healthcare expend_food expend_transport" // List of all variables that should be deflated
global winvars = "wealth income valhouse wealth_prnt" // variables to winsorize
global baseyear "2015" // Which year to use as reference year
global baseunit "1000" // Divide by thousand, to get thousands of 2015 dollars"

// switches/values
global mingift   = .5 // Value for the minimum size of gifts that are given
global firstyear = 1999
global agestep 	 = 6
global firstage  = 25
global lastage   = 85
global retage    = 67
global minagepar = 55
global period    = "benchmark"
global nv        = 3
global nvp       = 3
global minyear = 1999
global minmatch = 2
global Nxk = 4
global Nxp = 4
global wealth_threshold = 0 // Left-censoring for initial wealth (e.g., can't start with negative values since households can't start with debt)
global maxagediff = 40
global minagediff = 15
global polyorder = 4
global biennialize = 1 // Do you want to biennialize
global biennializerate = 1. 	// Increase transfer rate by 
global biennializesize = 1.     // Increase transfer size by
global Nboot=100

cd $basepath
// Setup data
do src/Stata/FRED.do // Download FRED data (used for deflators)
do src/Stata/PSID.do // Unzip, verify, and merge PSID waves:
 * First run: Unzip all PSID files (slow, say 3-5 minutes)
 * Later runs: Much faster
 * Note: As described in the readme, the RT13 files have to be manually unzipped and the dofiles modified
do src/Stata/PSID_panel.do // Sample selection and cleaning, merging with parent info etc

// Create tables, figures, and moments
do src/Stata/descr_table.do
do src/Stata/emp_eventstudy.do 
do src/Stata/regr_hypotheses.do
do src/Stata/PSID_calibration.do 
do src/Stata/PSID_moments.do

// Just to plot plots in introduction
do src/Stata/AHS.do
do src/Stata/SHED.do
do src/Stata/SHED_AHS_plots.do





