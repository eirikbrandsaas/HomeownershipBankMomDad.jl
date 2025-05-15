local minage = 50
local maxage = 80
local step = 1


global tabvars "totrcved_fromparent wealth wealth_prnt income coll white  owner Lowner age" // Variables for statistisc

use  "$basepath/data/PSID/panel_ind", clear

sort id_hd t
gen Lowner = l.owner
gen Lwealth = l.wealth
gen Lwealth_prnt = l.wealth_prnt
keep if year == 2013 // Only year with the transfer supplement data
local agestep = 7 // Pick agestep 
keep if age <= 44
keep if  age >= $firstage
noi sum totrcved_fromparent if rcvr == 1 [aw = famwgt]
noi sum rcvr [aw = famwgt]
replace totrcved_fromparent = 0 if totrcved_fromparent < $mingift
replace rcvr = 0 if totrcved_fromparent < $mingift
// Cross-tab over age and recieving

** Labeling
label var totrcved_fromparent "Transfer"
label var wealth "Wealth"
label var Lwealth "Wealth t-2"
label var wealth_prnt "Wealth Parent"
label var Lwealth_prnt "Wealth Parent t-2"
label var income "Income"
label var age "Age"
label var coll "College"
label var hs "High School"
label var white "White"
label var owner "Owner"
label var Lowner "Owner t-2"

local j = 0
forv age = 25(`agestep')39 {
	preserve
	qui {
	local j =  `j' + 1
	keep if age >= `age' & age < `age' + `agestep'
	sum age
	noi display "`r(min)' to `r(max)'"
	estpost tabstat $tabvars if rcvr==0  [aw = famwgt], statistics(mean) columns(statistics)    
	est store Tab_0_age`age'_mean
	estpost tabstat $tabvars if rcvr==1  [aw = famwgt], statistics(mean) columns(statistics)    
	est store Tab_1_age`j'_mean
	}
	restore
}

// Cross-tab over wealth quintiles and recieving
cap drop qwealth
xtile qwealth = wealth, nq(3)
local wealthstep = 1
forv nq = 1/3 {
	preserve
	qui {
	keep if qwealth == `nq'
	estpost tabstat $tabvars if rcvr==0   [aw = famwgt] ,statistics(mean) columns(statistics)
	est store Tab_0_w`nq'_mean
	estpost tabstat $tabvars if rcvr==1   [aw = famwgt], statistics(mean) columns(statistics)    
	est store Tab_1_w`nq'_mean
	}
	restore
}
// Cross-tab over house tenure and recieving
qui {	
forv owner = 0/1 {
	preserve
	keep if owner == `owner'
	estpost tabstat $tabvars if rcvr==0   [aw = famwgt] ,statistics(mean) columns(statistics)
	est store Tab_0_owner`owner'_mean
	estpost tabstat $tabvars if rcvr==1   [aw = famwgt], statistics(mean) columns(statistics)    
	est store Tab_1_owner`owner'_mean
	restore
}

estpost tabstat $tabvars if rcvr==0   [aw = famwgt] ,statistics(mean) columns(statistics)
est store Tab_0_all_mean
estpost tabstat $tabvars if rcvr==1   [aw = famwgt] ,statistics(mean) columns(statistics)
est store Tab_1_all_mean
}

** Store tables to tex	
esttab Tab*_age*_mean using "tabfig/descr/descr_trsnf_age_mean.tex", cells("mean(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap
esttab Tab*_w*_mean   using "tabfig/descr/descr_trsnf_wlt_mean.tex", cells("mean(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap
esttab Tab_*all_mean Tab*_owner*_mean using "tabfig/descr/descr_trsnf_disc_mean.tex", cells("mean(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap
