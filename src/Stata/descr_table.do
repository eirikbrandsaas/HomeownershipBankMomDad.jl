local minage = 50
local maxage = 80
local step = 1


global tabvars "totrcved_fromparent wealth wealth_prnt income coll white  owner Lowner age" // Variables for statistisc

use  "$basepath/data/PSID/panel_ind", clear

sort id_hd t
gen Lowner = l.owner
gen Lwealth = l.wealth
gen Lwealth_prnt = l.wealth_prnt
keep if year == 2013
keep if age < 45 
keep if  age >= $firstage
noi sum totrcved_fromparent if rcvr == 1 [aw = famwgt]
noi sum rcvr [aw = famwgt]
replace totrcved_fromparent = 0 if totrcved_fromparent < $mingift
replace rcvr = 0 if totrcved_fromparent < $mingift
noi sum totrcved_fromparent if rcvr == 1 [aw = famwgt]
noi sum rcvr [aw = famwgt]

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


// Cross-tab over age and recieving
keep if age <= 55 
foreach stat in mean p50 {
	local agestep = 4 // Pick agestep and start age to generate 4 categories
	forv age = 25(`agestep')41 {
		preserve
		 {
		noi display `age'
		keep if age >= `age' & age < `age' + `agestep'
		estpost tabstat $tabvars if rcvr==0 [aw = famwgt], statistics(`stat') columns(statistics)    
		est store Tab_0_age`age'_`stat'
		estpost tabstat $tabvars if rcvr==1 [aw = famwgt], statistics(`stat') columns(statistics)    
		est store Tab_1_age`age'_`stat'
		}
		restore
	}

	// Cross-tab over wealth quintiles and recieving
	cap drop qwealth
	xtile qwealth = wealth, nq(5)
	local wealthstep = 1
	forv nq = 1/5 {
		preserve
		qui {
		keep if qwealth == `nq'
		estpost tabstat $tabvars if rcvr==0 [aw = famwgt] ,statistics(`stat') columns(statistics)
		est store Tab_0_w`nq'_`stat'
		estpost tabstat $tabvars if rcvr==1 [aw = famwgt], statistics(`stat') columns(statistics)    
		est store Tab_1_w`nq'_`stat'
		}
		restore
	}


	// Cross-tab over maritial status and house tenure and recieving
	qui {	
	forv marr = 0/1 {
		preserve
		keep if married == `marr'
		estpost tabstat $tabvars if rcvr==0 [aw = famwgt] ,statistics(`stat') columns(statistics)
		est store Tab_0_marr`marr'_`stat'
		estpost tabstat $tabvars if rcvr==1 [aw = famwgt], statistics(`stat') columns(statistics)    
		est store Tab_1_marr`marr'_`stat'
		restore
	}
	forv owner = 0/1 {
		preserve
		keep if owner == `owner'
		estpost tabstat $tabvars if rcvr==0 [aw = famwgt] ,statistics(`stat') columns(statistics)
		est store Tab_0_owner`owner'_`stat'
		estpost tabstat $tabvars if rcvr==1 [aw = famwgt], statistics(`stat') columns(statistics)    
		est store Tab_1_owner`owner'_`stat'
		restore
	}

	estpost tabstat $tabvars if rcvr==0 [aw = famwgt] ,statistics(`stat') columns(statistics)
	est store Tab_0_all_`stat'
	estpost tabstat $tabvars if rcvr==1 [aw = famwgt] ,statistics(`stat') columns(statistics)
	est store Tab_1_all_`stat'
	}

** Store tables to tex	
esttab Tab*_age*_`stat' using "tabfig/descr/descr_trsnf_age_`stat'.tex", cells("`stat'(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap
esttab Tab*_w*_`stat'   using "tabfig/descr/descr_trsnf_wlt_`stat'.tex", cells("`stat'(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap
esttab Tab_*all_`stat' Tab*_marr*_`stat' Tab*_owner*_`stat' using "tabfig/descr/descr_trsnf_disc_`stat'.tex", cells("`stat'(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap

// esttab Tab*_age25_`stat' Tab*_age30_`stat' Tab*_age35_`stat' Tab*_age40_`stat' using "tabfig/descr/descr_trsnf_age_`stat'_narrow.tex", cells("`stat'(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap
esttab Tab_*all_`stat' Tab*_owner*_`stat' using "tabfig/descr/descr_trsnf_disc_`stat'_narrow.tex", cells("`stat'(fmt(%5.2f))") label booktabs nonum f collabels(none) plain replace nomtitles nogap
}
