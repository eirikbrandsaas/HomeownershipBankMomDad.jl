*** Code to create the datafile from the SCF, waves 1999-2016
** Only requires the summary files

// First download SHED data

clear
cd "$basepath/data/SHED/" //  have to change dir since Stata doesn't allow unzipping to destination"

forvalues y = 2013/2017 {
	cap confirm file "public`y'.dta" // Check if file exists
	if _rc { // if not download it
		local file = "SHED_public_use_data_`y'_(STATA).zip"
		if `y' ==  2013 {
			local file = "SHED_data_`y'_(STATA).zip"
		}
		copy "https://www.federalreserve.gov/consumerscommunities/files/`file'" .
		unzipfile "`file'", replace
	}
}
cd "$basepath" // change directory back to root

****************
** Open and clean the data
****************

// Variables that you want to rename

local commonvars "ppmsacat ppreg9 ppstaten ppethm ppgender pphhsize ppfs0596"

local oldnames2013 "R1_a 	R1_b 	R1_c 	R1_d 	R1_e 	R1_f 	R1_g 	R1_h"
local newnames2013 "rcheap	rconv	rmove	rqual	rdown	rpref	rlook	rothe"
local oldnames2014 "R1_a 	R1_b 	R1_c 	R1_d 	R1_e 	R1_f 	R1_g 	R1_h"
local newnames2014 "rcheap	rconv	rmove	rqual	rdown	rpref	rlook	rothe"
local oldnames2015 " "
local newnames2015 " "
local oldnames2016 "R1_a 	R1_c 	R1_d 	R1_e 	R1_f 	R1_g 	R1_h	R1_i	R1_b"
local newnames2016 "rcheap	rconv	rmove	rqual	rdown	rpref	rlook	rothe	rrisk"
local oldnames2017 "	"
local newnames2017 "	"

local oldfinan2013 " "
local newfinan2013 " "
local oldfinan2014 "H7_a	H7_b	H7_c	H7_d	H7_e	H7_f"
local newfinan2014 "hprev	hsave	hgift	hscnd	hgovt	hnodo"
local oldfinan2015 "H7_a	H7_b	H7_c	H7_d	H7_e	H7_f	H7_g 	H6"
local newfinan2015 "hprev	hsave	hgift	hscnd	hgovt	hnodo	hothe	prevown"
local oldfinan2016 "H7_a	H7_b	H7_c	H7_d	H7_e	H7_f	H7_g	H6"
local newfinan2016 "hprev	hsave	hgift	hscnd	hgovt	hnodo	hothe	prevown"
local oldfinan2017 " "
local newfinan2017 " "


local oldvars	"ppage	ppincimp ppfs0596	SL1	GH2	FS10	FS20_a 		FS20_b	FS20_c 	FS20_d 	FS20_e	FS30_a 	FS30_b 	FS30_c 	FS30_d 	FS30_e 	FS40"
local newvars	"age	totinc	 totsavings	studln	h_year	g	g_rentmort	g_educ	g_car	g_bills g_gnrl	gs_prnt gs_chld	gs_rltv	gs_frnd	gs_oth	gp"		 	 	 

local captolow2013 "PPAGE PPINCIMP PPMSACAT PPREG9 PPSTATEN PPETHM PPGENDER PPHHSIZE"

local capgenvars 	"ppfs0596 SL1" // Variables that you want to generate for each year, even though they may be missing

local renameold2013	"H0 	weight_orig"
local renamenew2013	"GH2	wgt"

local renameold2014	"D6	weight3	"
local renamenew2014	"GH2	wgt"

local renameold2015	"weight3"
local renamenew2015	"wgt"

local renameold2016	"weight3b"
local renamenew2016	"wgt"

local renameold2017	"weight5b"
local renamenew2017	"wgt"

global years "2013 2014 2015 2016  2017"
qui  foreach year of global years {

	use "data/SHED/public`year'.dta", clear
	noi disp `year'
	
	if `year' == 2013 {
		rename `captolow2013', lower
	}
	
	foreach var of local capgenvars {
		  cap gen `var' = .
	}
	
	local oldnames "oldnames`year'"
	local newnames "newnames`year'"
	local oldfinan "oldfinan`year'"
	local newfinan "newfinan`year'"
	local oldrename "renameold`year'"
	local newrename "renamenew`year'"
	cap rename (``oldnames'') (``newnames'')
	rename (``oldrename'') (``newrename'')
	cap rename (``oldfinan'') (``newfinan'')
	

	local Nvars: word count `oldvars'
	
	forv i=1/`Nvars' {
		local oldvar `: word `i' of `oldvars''
		local newvar `:word `i' of `newvars''
		gen `newvar' = .
		cap replace `newvar' = `oldvar'
	}
	
	keep wgt ``newnames'' ``newfinan'' `commonvars' `newvars'  //  `newfinan`year''
	des ``newnames'' 
	
	gen year = `year'
	
	save "data/SHED/SHED`year'", replace
}


clear all
set obs 0
foreach year of global years {
	append using "data/SHED/SHED`year'.dta"
}


// Cleaning/making variables comparable across years
replace totinc = 19 if totinc>=19 & !missing(totinc)

// Label variables
label var year		"Year"

label var rcheap 	"Cheaper to rent"
label var rconv  	"More convenient to rent"
label var rmove		"Plan to move"
label var rqual		"Doesn't qualify for mortgage"
label var rdown		"Can't afford down payment"
label var rpref		"Prefer to rent"
label var rlook		"Looking to buy"
label var rrisk 	"Safer to rent"
label var rothe		"Other"

label var hprev 	"Sale of previous home"
label var hsave  	"Personal savings"
label var hgift		"Loan or gift from family"
label var hscnd		"Second mortgage"
label var hgovt		"Goverment assistance"
label var hnodo		"No down payment"
label var hothe		"Other"

label var prevown 	"Previously owned"
label var h_year	"Year bought/moved into current home"

label var age		"Age"
label var totinc	"Total income"
label var studln 	"Has student loan"
label var totsavings	"Total value of savings and investments"

label var g		"Recieves regular financial support"
label var g_rentmort	"Rent/Mortgage"
label var g_educ	"Educ exp./stud. loan"
label var g_car		"Car payments"
label var g_bills	"Other bills"
label var g_gnrl	"General exp./other"

label var gs_prnt	"Parent src"
label var gs_chld	"Adult chlild src"
label var gs_rltv	"Other relative src"
label var gs_frnd	"Friends"
label var gs_oth	"Other"

label var gp		"Regular provides support to others?"



save "data/SHED/SHED_stacked.dta", replace
