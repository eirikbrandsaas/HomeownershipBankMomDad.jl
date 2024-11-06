local datapath $basepath/data

*********************************
** Create a stacked file with all id's and their wealth (just to match with givers in the next block)
*********************************

local cons = "expend_food  + expend_transport + expend_healthcare +  rent*12 + 0.06*valhouse"
foreach id in hd wf {
	use "$basepath/data/PSID/PSID17_HouseholdPanel.dta", clear
	gen cons = 0 // `cons'
	noi disp "Here you used to find consumption (PSID_panel.do)"
	keep id_`id' wealth income age year cons
	rename id_`id' id
	save temp_`id',replace
}

use temp_hd, clear
append using temp_wf
rename wealth Wealth
rename income Income
rename cons Cons
rename age Age
drop if id == .
save "$basepath/data/PSID/wealthyearind.dta", replace

rm temp_hd.dta
rm temp_wf.dta

*********************************
** Merge parental observations
*********************************
use "$basepath/data/PSID/PSID17_HouseholdPanel.dta", clear
gen cons = 0 // `cons'
noi disp "Here you used to find consumption (PSID_panel.do)"
local gen = "dad"
local sp = "hd"
sort id_hd year

// Merge in income and wealth of all four parents
foreach sp in hd wf{
foreach gen in dad mom {
	rename id`gen'_`sp' id
	merge m:1 id year using "$basepath/data/PSID/wealthyearind.dta", nogen keep(1 3) keepusing(Wealth Income Age Cons)
	rename (Wealth Income Age Cons) (wealth_`gen'_`sp' income_`gen'_`sp' age_`gen'_`sp' cons_`gen'_`sp' )
	rename id id`gen'_`sp'
	
}
}

// Use wealth of moms, not dads, unless mom's wealth is missing
foreach sp in hd wf {
	foreach var in wealth income age cons {
		replace `var'_mom_`sp' = `var'_dad_`sp' if missing(`var'_mom_`sp' )
		gen `var'_prnt_`sp' = `var'_mom_`sp'
	}
	
}

foreach var in wealth income age cons {
	gen `var'_prnt = `var'_prnt_hd
	replace `var'_prnt = `var'_prnt_wf if `var'_prnt_hd == .
	drop `var'_prnt_hd `var'_prnt_wf
}

drop wealth_???_* income_???_* age_???_* cons_???_* 
sort id_hd year

egen matchedobs = total(!missing(wealth_prnt)) , by(id_hd)
egen maxage = max(age), by(id_hd)
egen minage = min(age), by(id_hd)
gen agediff = age_prnt - age


// Create variables that you use later
gen rent2inc = rent*12/income
replace rent2inc = . if owner == 1
replace rent2inc = . if rent2inc < 0 | rent2inc > 1
gen hval2wealth = valhouse/wealth
replace hval2wealth = . if owner == 0
replace hval2wealth = . if hval2wealth < 0  | hval2wealth > 15
replace hval2wealth = . if wealth <= 0 

save "$basepath/data/PSID/panel_temp", replace

*********************************
** Merge the individual level transfer data
*********************************
local datapath $basepath/data
local basepath $basepath
global keepvars "totgiven tother thouse tschool rectype id_giver"
use "$basepath/data/PSID/PSID17_HouseholdPanel.dta", clear


sort id_hd t
keep if year == 2013

keep id_hd id_wf year age
rename id_hd id_rcver 
ds id_rcver, not
des `r(varlist)', varlist
global vars "`r(varlist)'" // Find all variables in pre-merge dataset
merge 1:m id_rcver using "$basepath/data/PSID/transfers_clean.dta", keep(1 3) keepusing($keepvars)
rename _merge merge1
drop if rectype == 2 // Only want child records (parent giving to kids)

rename id_giver id // Merge in the givers wealth
merge m:1 id year using "$basepath/data/PSID/wealthyearind.dta", nogen keep(1 3) keepusing(Wealth)
rename id id_giver

collapse (first) $vars (sum) totgiven (first) merge1 (mean) Wealth, by(id_rcver)
// replace totgiven = . if merge1 != 3 // Stata thinks sums of reals and missings are 0...
// replace totrcved = . if merge1 != 3
rename tot* tot*_hd 
rename id_rcver id_hd

replace id_wf = -_n if id_wf == . // Negative unique numbers for households without wifes (So Stata doesnt complain about non-unique )
rename id_wf id_rcver
ds id_rcver, not
des `r(varlist)', varlist
global vars "`r(varlist)'" // Find all variables in pre-merge dataset
merge 1:m id_rcver using "$basepath/data/PSID/transfers_clean.dta", keep(1 3) keepusing($keepvars)
rename _merge merge2
// replace totgiven = . if merge2 != 3 // Stata thinks sums of reals and missings are 0...
// replace totrcved = . if merge2 != 3
drop if rectype == 2

rename id_giver id
merge m:1 id year using "$basepath/data/PSID/wealthyearind.dta", nogen keep(1 3) keepusing(Wealth)
rename id id_giver

collapse (first) $vars (sum) totgiven tother thouse tschool (first) merge2, by(id_rcver rectype) 
rename id_rcver id_wf
replace id_wf =. if id_wf < 0
rename (totgiven) (totgiven_wf)

// drop if totrcved_hd == . | totrcved_wf == .
keep if merge1>1 | merge2>1 // Only keep households who have parents in the ParChd file
drop if merge1 == 3 & merge2 == 3 // Drop if any households has two PSID genes

gen totgiven_fromparent = totgiven_hd + totgiven_wf
label var totgiven_fromparent "Transfer From Parents"


rename Wealth wealth_giver
drop tot?????_hd tot?????_wf

gen rcvr = (totgiven_fromparent>0.0) 
rename *given* *rcved* 

foreach var in rcved {
	gen ctot`var'_fromparent = tot`var'_fromparent
	replace ctot`var'_fromparent = . if tot`var'_fromparent == 0
}

save "$basepath/data/PSID/indtransfer_17", replace

********************************************************
** Import that transfer data into the main file:
********************************************************
keep if age >=$firstage & age<=$lastage
local basepath $basepath
use "$basepath/data/PSID/panel_temp", replace
merge 1:1 id_hd year using "$basepath/data/PSID/indtransfer_17", keepusing(wealth_giver rcvr tot* ctot* tother tschool thouse) nogen keep(1 3)

order id_hd id_wf wealth wealth_*
save "$basepath/data/PSID/panel_ind", replace

********************************************************
** Winsorize the main file
********************************************************
foreach var of global winvars {
	winsor2 `var', cuts(1 99) by(age) replace
}

save "$basepath/data/PSID/panel_winsor", replace


*********************************
** Merge the family level transfer data
*********************************

use "$basepath/data/PSID/PSID17_HouseholdPanel.dta", clear

rename fuid fuid_store
gen fuid = -_n 				// give a negative number so that you merge 1:1
replace fuid = fuid_store if year == 2013 // replace the negative numbers with the actual interview number in 2013 (>0 in the using file)
merge 1:1 fuid using $basepath/data/PSID/transfers_family_clean.dta, nogen keep(1 3) // Drop households only in the family transfer file)

// Some cleaning and treatment of the gift variables
replace gve_chld = . if n_chld_records == 0
replace gve_chld = . if gve_chld < 5.0
gen haskids = (n_chld_records != 0)
gen giver = (gve_chld != 0 & gve_chld !=.)*100 
gen rcver = (rcv_chld != 0 & rcv_chld !=.)*100 
gen giver_house = (gve_home_amt!= 0 & gve_home_amt !=.)*100
gen cgift = gve_chld
gen crcvd = rcv_chld
keep if haskids == 1
replace	cgift = . if gve_chld == 0 | gve_chld <= $mingift  
replace	crcvd = . if rcv_chld == 0 | rcv_chld <= $mingift  
gen gift2wealth = gve_chld/wealth
replace gift2wealth = . if wealth<0 | gift2wealth>1
replace gift2wealth = . if cgift == .

foreach ind in hd wf {
	replace wtr_`ind'_rcv_home = . if  wtr_`ind'_rcv_home == 8 | wtr_`ind'_rcv_home == 9
	replace wtr_`ind'_rcv_home = 0 if wtr_`ind'_rcv_home == 5
}

gen wtr_rcv_home = (wtr_hd_rcv_home == 1 | wtr_wf_rcv_home == 1)


foreach var in gve_chld haskids giver giver_house cgift gift2wealth {
	replace `var' = . if year!=2013
}

save "$basepath/data/PSID/panel_family", replace


